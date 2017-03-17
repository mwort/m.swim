import os
import sys
import numpy as np
import grass.script as grass
import grass.script.core as gcore
import mygrass as mg
import datetime as dt

grun = grass.run_command
gm = grass.message


def routing(units, accumulation, outname, searchradius=1.5):
    grun('r.mask',raster=units,overwrite=True)

    # make outlet raster points, border subbasins with neg accumul. not correct
    grass.message('Searching outlets...')
    grun('r.statistics', base=units,cover=accumulation, method='max',
         output='maxaccum__')
    exp = "outlets__=if('%s' == @maxaccum__,%s,null())" %(accumulation,units)
    grass.mapcalc(exp)

    grass.message('Growing outlets...')
    # join outlets (set all to 1) and then grow outlets
    grass.mapcalc("outlets1__=if(isnull(outlets__),null(),1)",overwrite=True)
    # grow to all 8 neighbours (1.5 times cell size)
    grun('r.grow', input='outlets1__', output='outlets__grown', radius=searchradius)
    grun('r.clump', input='outlets__grown', output='outlets__grown__clumped')

    # make inlets
    grass.message('Searching inlets...')
    grun('r.statistics', base='outlets__grown__clumped', cover=accumulation,
         method='max', output='maxoutlet__clumps')
    grass.mapcalc("inlets__=if(%s == @maxoutlet__clumps,%s,null())" %(accumulation,units))

    # transfer inlet units to clumps
    grun('r.stats.zonal', base='outlets__grown__clumped', cover='inlets__',
         method='max', output='clumps__unitID')
    # transfer to original units again (if n outlets >1, take largest inlet ID)
    grun('r.stats.zonal', base='outlets__', cover='clumps__unitID',
        method='max', output='next__unitID__outlets')
    grun('r.stats.zonal', base=units, cover='next__unitID__outlets',
         method='max', output='nextID__units__float')
    grass.mapcalc(outname+'=int(nextID__units__float)')

    return

def valley(slope, output, thresh=15, filter=300):
    '''Make a clean valley=1 and slope/rockface=0 raster'''
    # make valley 1/0 raster
    grass.mapcalc(exp='valleys__=%s<=%s'%(slope,thresh))
    # clean
    print 'Cleaning...'
    grun('r.resamp.filter', input='valleys__', output='valleys__smooth',
         filter='box', radius=filter)
    grass.mapcalc(exp='%s=int(round(valleys__smooth))'%output)

    return output


def cleanRast(input,output,areathresh,fill=True,lenthresh=None,quiet=False):

    # RASTER BASED ONLY
    grun('r.clump', input=input,output=output+'__clumped',quiet=quiet)
    grass.mapcalc('help__1=1')
    grun('r.stats.zonal',base=output+'__clumped',cover='help__1',method='sum',
         output='n__cells',quiet=quiet)
    # output by thresholding cells
    mincells=int(areathresh/float(grass.region()['nsres'])**2)
    print 'Removing all units with less than %s cells'%mincells
    interout = output+'__filtered'
    grass.mapcalc(exp='%s=if(n__cells>=%s,%s,null())'%(interout,mincells,input))

    if fill:
        # close/fill gaps
        grun('r.grow', input=interout, output=output+'__float', radius=mincells,
             quiet=quiet)
        interout = output+'__float'
    # make int
    grass.mapcalc(output+'=int(%s)'%interout)

#    # VECTOR BASED, v.clean
#    grun('r.to.vect', input=input, type='area', output=input+'__',#flags='t',
#         quiet=quiet,overwrite=True)
#    if lenthresh:
#        grun('v.clean', input=input+'__', output=input+'__clean',overwrite=True,quiet=quiet,
#             tool='rmarea,prune', thresh='%s,%s'%(areathresh,lenthresh))
#    else:
#        grun('v.clean', input=input+'__', output=input+'__clean',overwrite=True,quiet=quiet,
#             tool='rmarea', thresh=areathresh)
#    grun('v.to.rast', input=input+'__clean', output=output,
#         use='attr', attrcolumn='value',overwrite=True,quiet=quiet)
    return

def cleanGunits(gunits,subbasins,outname,smallarea):

    # get area and subbasins
    print 'Reading gunits for all subbasins...'
    gu = mg.zonalStats(subbasins,gunits)
    # largest gunit
    maxu = gu['non_null_cells'].max()

    print 'Largest unit: %s cells' %maxu

    # subbasins
    sbs = np.unique(gu['mean'])
    print len(sbs),'subbasins found'
    cunits = []
    maxn = 0
    print 'Looping over each subbasin...'
    for s in sbs:
        print 'Setting region and mask...'
        grun('v.extract',input=subbasins,cats=s,output='sb__region',quiet=True)
        grun('g.region',vect='sb__region')
        grun('r.mask',raster=subbasins,maskcats=s)

        outn = 'clean__gu__%04i' %s

        cleanRast(gunits,'first__iteration',smallarea,quiet=True)
        # second iteration
        cleanRast('first__iteration','second__iteration',smallarea,quiet=True)

        # add maxid to avoid neighbouring merge when patching
        grass.mapcalc('{0}=second__iteration+{1}'.format(outn,maxn))
        # save for patching
        cunits+=[outn]
        maxn = max(map(int,grass.read_command('r.stats',input=outn,flags='n').split()))+1
        grun('r.mask',flags='r')

        grun('g.region',flags='d')
        print 'Subbasin %s done!' %s

    # patch together
    if len(cunits)<500:
        grun('r.patch',input=','.join(cunits),output=outname)
    else:
    # patch in batches due to open files limit of r.patch
        patched=[]
        batches=range(0,len(cunits),500)+[len(cunits)]
        for i in range(len(batches)-1):
            n = outname+'__%s'%i
            grun('r.patch',input=','.join(cunits[batches[i]:batches[i+1]]),output=n)
            patched+=[n]
        grun('r.patch',input=','.join(patched),output=outname)
    grass.mapcalc(exp='%s=int(%s)'%(outname,outname))

    return

def checkRouting(nextID,gunits,outname):
    import pandas as pa
    tbl=grass.read_command('r.univar',map=nextID,zones=gunits,flags='t').split('\n')
    cols = tbl[0].split('|')
    tbl=pa.DataFrame([tuple(l.split('|')) for l in tbl[1:-1]],columns=cols)
    tbl.set_index(tbl['zone'].astype(int),inplace=True)
    tbl=tbl[['mean','non_null_cells']].astype(int)
    # group by 'inlet ids' and accum cells
    incells = tbl.groupby('mean').agg('sum')['non_null_cells']
    receivingcells= tbl.ix[incells.index]['non_null_cells']
    # ratio
    routratio = incells/receivingcells
    print routratio.describe()
    mg.areclass(gunits,np.array(zip(routratio.index,routratio)),outname)

    return 0


def aspect(elevation,output,nclasses):
    '''Make nice aspect classes: 360 should be devisible by nclasses'''
    ### ASPECT
    grun('r.slope.aspect', elevation=elevation, aspect='aspect__', precision='CELL')
    # resulting map: degrees from E ccw, ie. E=0, N=270

    # make aspect classes
    step=360/nclasses
    start=270+(step/2) # centered around north

    # from start to 0 (-1 for python)
    breaks=range(start,-1,-step)
    exp=['if(aspect__ <= %s & aspect__ > %s,%s,0)' %(breaks[i],breaks[i]-step,i+1) for i in range(len(breaks))][:-1]

    # from start to 360
    breaks2=range(start,360,step)[::-1]
    # if not breaking at 0deg/E, add class going across 0
    if breaks[-1]!=0:
        exp+=['if(aspect__ <= %s | aspect__ > %s,%s,0)' %(breaks[i],360-breaks[i],i+1)]
        breaks2 = breaks2[1:]

    exp+=['if(aspect__ <= %s & aspect__ > %s,%s,0)' %(breaks2[i]+step,breaks2[i],len(exp)+i+1) for i in range(len(breaks2))]

    # add all rules together
    exp = output+'='+('+'.join(exp))
    gm(exp)
    grass.mapcalc(exp=exp,overwrite=True)

    return output

def hydAddressRast(hydrotopes,subbasins,output):
    '''Create a hydrotope address map, i.e. the n'th hydrotope for each subbasin
    as they appear in the structure file.
    '''
    # get subbasin for each hydrotope
    tbl=grass.read_command('r.univar', map=subbasins,zones=hydrotopes,
                           flags='gt').split('\n')[:-1] #:-1 as last line hast line break]
    tbl = [tuple(l.split('|')) for l in tbl]
    tbl = np.array(tbl[1:],dtype=zip(tbl[0],['S250']*len(tbl[0])))
    sub = np.array(tbl['mean'],dtype=int)
    subIDs = np.unique(sub)
    # get counts
    freq,bins = np.histogram(sub,len(subIDs))
    # cumulative and shift by one subbasin
    cumfreq = np.array([subIDs,[0]+freq.cumsum().tolist()[:-1]]).T
    # reclass subbasin map to this count
    tmpf = grass.tempfile()
    np.savetxt(tmpf,cumfreq,delimiter='=',fmt='%i')
    grass.run_command('r.reclass',input=subbasins,output='subbasins__cumfreq',rules=tmpf)
    # subtract from hydrotopes
    grass.mapcalc('%s=%s - subbasins__cumfreq' %(output,hydrotopes))
    # remove help stuff
    grass.run_command('g.remove',type='rast',name='subbasins__cumfreq',flags='f')
    grass.message('Created hydrotope address raster: %s' %output)
    return

def classedContours(elevation,valleys, output,valinter=40, slopeinter=400):
    '''Make nice aspect classes: 360 should be devisible by nclasses'''
    val = valleys

    # get stats of DEM for nice breaks
    stats=grass.parse_command('r.univar',map=elevation,flags='g')
    # function to return the r.mapcalc condition with nice elevation breaks
    def niceBreaks(interval,last):
        # make expression
        exp = 'int({elev}/{interval})+{last}'.format(elev=elevation,last=last,
                                                     interval=interval)
        # calculate next last/highest value
        last = int(float(stats['max'])/interval)+last
        return exp, last

#    # make classed condition
#    exp = []
#    last=1
#    for n,c in enumerate(contours):
#        con,last = niceBreaks(c,last)
#        condi = 'if({sl} >= {lo} & {sl} < {up}, {contours})'
#        exp += [condi.format(sl='slope__',lo=classes[n],up=classes[n+1],contours=con)]
#
    # valleys
    vcon,last = niceBreaks(valinter,0)
    scon,last = niceBreaks(slopeinter,last)
    exp = output+'=if(%s,%s,%s)'%(val,vcon,scon)
#    # add all rules together
#    exp = output+'='+('+'.join(exp))
    print(exp)
    grass.mapcalc(exp=exp,overwrite=True)

    return output

def writeStr(outname, unitmap, subbasins, hydAddress, nextID, *othercols):
    # get subbasin, hydAddress, nextID
    subb   = meanUnit(unitmap,subbasins,int)
    subbasins = subbasins.split('@')[0]
    hydAdd = meanUnit(unitmap,hydAddress,int)
    nxt    = meanUnit(unitmap,nextID,int)

    # make sure none are flowing to unit 0
    nxt[nextID][nxt[nextID]==0] = nxt['id'][nxt[nextID]==0]

    # get other columns
    ocols = ()
    for c in othercols: ocols += (meanUnit(unitmap, c, float)[c.split('@')[0]],)

    # get all columns
    data = zip(subb['id'],subb[subbasins],hydAdd[hydAddress],
               nxt[nextID],*ocols)
    columnnames = ['gunitID','subbasinID','hydrotopeID','nextID']+[c.split('@')[0] for c in othercols]
    data = np.array(data, dtype=zip(columnnames,4*[int]+len(othercols)*[float]))
#    # dont take those that are outflows
#    data = data[nxt['id']!=nxt[nextID]]
    # save to txt file
    f = file(outname,'w')
    f.write((4*'%-10s '+len(othercols)*'%-14s '+'\n')%tuple(columnnames))
    np.savetxt(f,data,fmt=4*'%10i '+len(othercols)*'%16.6f ')
    f.close()
    print 'Wrote %s' %outname
    return data

def meanUnit(units,rast,outtype=str):
    '''Get mean values over units and return as array with two column
    (units int, rast with outtype type)
    '''
    tbl=grass.read_command('r.univar', map=rast,zones=units,flags='gt')
    tbl   = [tuple(l.split('|')) for l in tbl.split('\n')[:-1]]#:-1 as last line hast line break]
    array = np.array(tbl[1:],dtype=zip(tbl[0],['S250']*len(tbl[0])))
    out   =  np.array(zip(array['zone'],array['mean'],
                          array['non_null_cells'],array['null_cells']),
                      dtype=[('id',int),(rast.split('@')[0],outtype),
                             ('non_null_cells',int),('null_cells',int)])
    print 'Got %r x %s' %(out.dtype.names,out.shape[0])
    return out


def loadGWE(file, prefix, reclassrast):
    # load gla file
    gla=pa.read_table(file,index_col=[0,1],names=['id','value'],delim_whitespace=True)
    styear=1979
    # colormap
    max = 2100
    colormap='''-%s   magenta
                 -200 red
                   -1 orange
                    0 yellow
                    1 green
                  200 blue
                 %s   cyan''' %(max,max)
    # loop over indeces (year,index) and reclass units raster
    dom = [31,28,31,30,31,30,31,31,30,31,30,31]
    outmaps = []
    for l1,l2 in np.unique(gla.index.values):
        outname = prefix+'_%04i-%02i'%(l1,l2)
        grass.message('Processing %s '%(outname))
        print gla.ix[l1,l2].as_matrix().shape
        mg.areclass(reclassrast,gla.ix[l1,l2].as_matrix(),outname)
        # nice colour map
        grass.write_command('r.colors',map=outname,rules='-',stdin=colormap)
        # timestamp
        grass.run_command('r.timestamp',map=outname,
                          # get last date of month by month+1-1day
                          date=dt.date(styear-1+l1,l2,dom[l2-1]).strftime('%d %b %Y'))
        outmaps += [outname]
    # make nice timeseries out of it
    grass.run_command('t.create', output=prefix, title=prefix, description=prefix)
    grass.run_command('t.register', input=prefix, maps=','.join(outmaps))
    return


if __name__ == '__main__':
    st = dt.datetime.now()
    # get options and flags
    o, f = parser()
    fmt = lambda d: '\n'.join(['%s: %s' % (k, v) for k, v in d.items()])+'\n'
    grass.message('GIS Environment:\n'+fmt(grass.gisenv()))
    grass.message('Parameters:\n'+fmt(o)+fmt(f))


    # Parameters
    glacierarea = 'glacierarea'  # existing raster
    minarea = 100000  # m2
    # all hydrotope inputs are needed including elevation and contours
    # could maybe be parsed through grass.script.core._parse_opts(lines)

    valleythreshold = 19  # dg
    valleysmoothing = 200  # some sort of radius
    countoursvalley = 30  # m
    countoursslopes = 200  # m (maybe make dependent on mean valley/slope slopes?)
    naspectclasses = 4  # number of aspect classes
    hemisphere = 'N'  # N | S for summer/winter

    ## Output
    d = False  # skip recreating slope_dg, aspect and sunhours
    gunits = 'gunits'  # glacier units raster
    slope_dg = 'slope_dg'  # slope in degrees
    aspect = 'aspect'  # aspect CELL precision
    sunhoursprefix = 'sun_hr'  # prefix for r.sun Output
    valleys = 'valleys'  # valley/slope raster
    glaciercontours = 'glaciercontours'  # valley/slope classed contours
    slopeaspect = 'slope_aspect'  # slope aspects in naspectclasses
    routingratio = 'routingratio'  # optional: number of cells upstream/downstream
    gunitssoils = 'gunitssoils'  # soils mapped onto glacier units
    gunitslanduse = 'gunitslanduse'  # landuse mapped onto glacier units

    ############################################################

    # mask entire catchment
    r.mask subbasins@subbasins

    # clean small bits from garea
    cleanRast glacierarea garea $minarea False

    r.mask garea --o
    # reduce region to garea

    ### some DEM derived maps
    # make slope and aspect raster
    r.slope.aspect elevation=aster_iclean@DEM slope=slope_dg aspect=aspect precision=CELL
    # calculate sun hours on June 21 and Dec 21
    r.sun elevation=aster_iclean@DEM day=172 insol_time=sun_hr_doy172
    r.sun elevation=aster_iclean@DEM day=355 insol_time=sun_hr_doy355

    # distinguish valleys and slopes (with smoothing)
    valley() slope_dg valleys_smooth 19 200
    # clean
    cleanRast valleys_smooth valleys $minarea

    # make contours (also writes slope and valleys)
    classedContours aster_iclean@DEM valleys contours_classed 30 200
    # clean
    cleanRast contours_classed contours_clean $minarea

    # make aspect classes (only on slopes)
    aspect aster_iclean@DEM aspect_classed 4
    # only take those not in valleys and rest null
    r.mapcalc exp='aspect_slopes=int(if(!valleys,aspect_classed,0))'
    # filter out small slopes and add them to neighbouring slopes
    cleanRast aspect_slopes aspect_clean $minarea

    # make raw glacier units
    r.cross input=subbasins@subbasins,contours_clean,aspect_clean out=gunits_raw
    # spatially explicit
    r.clump input=gunits_raw output=gunits_clumped

    ### clean gunits and spatially explicit
    # then clean each subbasin individually to fill each subbasin
    cleanGunits gunits_clumped subbasins@subbasins gunits_clean $minarea

    # remove small overspilled units at glacierarea fringes and reduce glacierarea accordingingly
    r.mask garea
    cleanRast gunits_clean gunits_cleaner $minarea 0
    r.mapcalc exp='garea=if(!isnull(gunits_cleaner),1,null())' --o
    r.mask garea --o
    # clump again to renumber (enhancing performance, would also work without)
    r.clump -d input=gunits_cleaner output=gunits

    # routing
    routing gunits accumulation@subbasins nextID
    # double check
    checkRouting nextID gunits routing_ratio

    # make glacier units part of hydrotope map and include contours outside of glacier area
    r.mask raster=subbasins@subbasins --o
    r.mapcalc exp='gcontours=if(isnull(gunits),int(aster_iclean@DEM/200),gunits)'

    # Prepare hydrotope maps
    # get valley/slope glacier units
    r.statistics base=gunits cover=valleys method=mode output=gunits_valleys
    r.mapcalc exp='gunits_valleys=int(@gunits_valleys)' --o

    # map soils to gunits with additional slope soil unit
    r.statistics base=gunits@glaciers cover=HWSD@soil method=mode output=gunits_soils
    r.mapcalc exp='gunits_soils=int(@gunits_soils)' --o
    r.mapcalc exp='soils_HWSD=if(isnull(gunits_valleys@glaciers),HWSD@soil,if(gunits_valleys@glaciers,gunits_soils,null()))'

    # map landuse to gunits
    r.statistics base=gunits@glaciers cover=GLC_EUv2@landuse method=mode output=gunits_landuse
    r.mapcalc exp='GLC_gunits=if(!isnull(garea@glaciers),int(@gunits_landuse),GLC_EUv2@landuse)'

    # make hydrotopes and write str file
    m.swim.hydrotopes subbasins=subbasins@subbasins landuse=landuse@landuse \
                        soil=soils_SWIM@soil -c contourrast=gcontours elevation=aster_iclean@DEM \
                        hydrotopes=hydrotopes strfilepath=proSWIM/Input/rhone.str

    # hydrotope addresses
    hydAddressRast hydrotopes@hydrotopes subbasins@subbasins hydAddress

    # write glacier str file (with 0 glacier height, 0 debris cover (cant have the same null0 name))
    r.mask --o garea
    writeStr ~/wortmann/Alps/proSWIM/Input/glaciers.str gunits \
                        subbasins@subbasins hydAddress nextID gla0 slope_dg \
                        sun_hr_doy172 sun_hr_doy355 debris0

    ############################################################


    # report time it took
    delta = dt.datetime.now()-st
    grass.message('Execution took %s hh:mm:ss' %delta)
