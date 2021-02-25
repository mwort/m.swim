#!/usr/bin/env python
#
############################################################################
#
# MODULE:      m.swim.glaciers v1.5
# AUTHOR(S):   Michel Wortmann, wortmann@pik-potsdam.de
# PURPOSE:     Preprocessing suit for the Soil and Water Integrated Model
#              Glacier dynamics (SWIM-G)
# COPYRIGHT:   (C) 2012-2021 by Wortmann/PIK
#
#              This program is free software under the GNU General Public
#              License (>=v2). Read the file COPYING that comes with GRASS
#              for details.
#
#############################################################################

#%Module
#% description: Preprocessor for the Soil and Water Integrated Model - Glacier dynamics (SWIM-G), requires pandas.
#%End
#%Option
#% guisection: Required
#% key: elevation
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% description: Elevation raster
#% gisprompt: old,cell,raster
#%end
#%Option
#% guisection: Required
#% key: accumulation
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% description: Accumulation raster, e.g. from m.swim.subbasins
#% gisprompt: old,cell,raster
#%end
#%Option
#% guisection: Required
#% key: glacierarea
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% description: Glacier modelling domain raster, e.g. created using a min. elevation
#% gisprompt: old,cell,raster
#%end
#%Option
#% guisection: Required
#% key: subbasins
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% description: Subbasin raster map
#% gisprompt: old,cell,raster
#%end
#%Option
#% guisection: Required
#% key: landuse
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% description: Land use raster map in SWIM land use classification
#% gisprompt: old,cell,raster
#%end
#%Option
#% guisection: Required
#% key: soil
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% description: Soil raster map classified to SWIM soil files
#% gisprompt: old,cell,raster
#%end
#%Option
#% guisection: Glacier units
#% key: minarea
#% type: double
#% answer: 100000
#% multiple: no
#% key_desc: sq. m
#% description: Minimum unit area in sq. m
#%end
#%Option
#% guisection: Glacier units
#% key: valleythreshold
#% type: double
#% answer: 19
#% multiple: no
#% key_desc: dg
#% description: Valley slope threshold in dg
#%end
#%Option
#% guisection: Glacier units
#% key: valleysmoothing
#% type: double
#% answer: 200
#% multiple: no
#% key_desc: m
#% description: Valley smoothing radius in m
#%end
#%Option
#% guisection: Glacier units
#% key: contoursvalleys
#% type: double
#% answer: 30
#% multiple: no
#% key_desc: m
#% description: Elevation contour intervals in valleys in m
#%end
#%Option
#% guisection: Glacier units
#% key: contoursslopes
#% type: double
#% answer: 200
#% multiple: no
#% key_desc: m
#% description: Elevation contour intervals on slopes in m
#%end
#%Option
#% guisection: Glacier units
#% key: naspectclasses
#% type: integer
#% answer: 4
#% multiple: no
#% key_desc: n classes
#% description: Number of aspect classes
#%end
#%Option
#% guisection: Glacier units
#% key: avalanchethreshold
#% type: double
#% answer: 45
#% multiple: no
#% key_desc: degrees
#% description: Slope threshold to delineate avalanche area in dg. If not given, avalanchearea raster must exist.
#%end
#%Option
#% guisection: Glacier units
#% key: hemisphere
#% type: string
#% answer: N
#% multiple: no
#% key_desc: N|S
#% description: Summer/winter differentiation
#%end
#%Flag
#% guisection: Output
#% key: d
#% label: Skip recreating slope_dg, aspect and sunhours (must all exist)
#%end
#%Flag
#% guisection: Output
#% key: s
#% label: Only rewrite structure file, all maps must exist
#%end
#%Option
#% guisection: Output
#% key: gunits
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% answer: gunits
#% gisprompt: new,cell,raster
#% description: Name of glacier units raster
#%end
#%Option
#% guisection: Output
#% key: downstreamgunits
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% answer: downstreamgunits
#% gisprompt: new,cell,raster
#% description: Name of next downstream glacier unit ID raster
#%end
#%Option
#% guisection: Output
#% key: slope_dg
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% answer: slope_dg
#% description: Name of slope raster
#%end
#%Option
#% guisection: Output
#% key: aspect
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% answer: aspect
#% description: Name of aspect raster
#%end
#%Option
#% guisection: Output
#% key: avalanchearea
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% answer: avalanchearea
#% description: Name of raster with avalanche area. Must exist if avalanchethreshold not given.
#%end
#%Option
#% guisection: Output
#% key: gunitsglacierarea
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% answer: gunitsglacierarea
#% description: Cleaned glacier area raster fitting to gunits
#%end
#%Option
#% guisection: Output
#% key: sunhoursprefix
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% answer: sun_hr
#% description: Prefix for sun hours rasters *_summer and *_winter
#%end
#%Option
#% guisection: Output
#% key: valleys
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% answer: valleys
#% gisprompt: new,cell,raster
#% description: Name of valleys raster
#%end
#%Option
#% guisection: Output
#% key: glaciercontours
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% answer: glaciercontours
#% gisprompt: new,cell,raster
#% description: Name of glacier valley and slope classed contours raster
#%end
#%Option
#% guisection: Output
#% key: slopeaspect
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% answer: slopeaspect
#% gisprompt: new,cell,raster
#% description: Name of slope aspect raster
#%end
#%Option
#% guisection: Output
#% key: gunitssoil
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% answer: gunitssoil
#% gisprompt: new,cell,raster
#% description: Name of soils raster mapped to glacier units
#%end
#%Option
#% guisection: Output
#% key: gunitslanduse
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% answer: gunitslanduse
#% gisprompt: new,cell,raster
#% description: Name of landuse raster mapped to glacier units
#%end

#%Option
#% guisection: Hydrotopes
#% key: hydrotopes
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% answer: hydrotopes
#% gisprompt: new,cell,raster
#% description: Name of hydrotope raster to be created
#%end
#%Option
#% guisection: Hydrotopes
#% key: contours
#% type: integer
#% required: yes
#% answer: 100
#% multiple: no
#% key_desc: int m
#% description: Elevation contour intervals outside of glacier area
#%end
#%Option
#% guisection: Hydrotopes
#% key: contourrast
#% type: string
#% required: no
#% multiple: no
#% key_desc: string
#% answer: contours
#% description: Elevation contours to be created with glacier units and contours outside of glacier area
#%end
#%Option
#% guisection: Hydrotopes
#% key: management
#% type: string
#% required: no
#% multiple: no
#% key_desc: name
#% description: Management raster name (if not given, column is filled with default)
#% gisprompt: old,cell,raster
#%end
#%Option
#% guisection: Hydrotopes
#% key: wetland
#% type: string
#% required: no
#% multiple: no
#% key_desc: name
#% description: Wetlands raster name (if not given, column is filled with default)
#% gisprompt: old,cell,raster
#%end

#%Option
#% guisection: Optional
#% key: strfilepath
#% type: string
#% multiple: no
#% key_desc: path
#% description: Path where structure file will be written if given, glacier input files will be created in the same directory as glaciers.str.
#% gisprompt: new,file,file
#%end
#%Flag
#% guisection: Optional
#% key: a
#% label: Add new soil ID for slope soils (recommended)
#%end
#%Flag
#% guisection: Optional
#% key: k
#% label: Keep intermediat files (include __ in names)
#%end
#%Option
#% guisection: Optional
#% key: initialdebris
#% type: string
#% multiple: no
#% key_desc: name
#% gisprompt: new,cell,raster
#% description: Initial debris concentration (fraction of total volume), if not given set to 0
#%end
#%Option
#% guisection: Optional
#% key: routingratio
#% type: string
#% multiple: no
#% key_desc: name
#% gisprompt: new,cell,raster
#% description: Create raster with number of cells upstream/downstream raster
#%end
#%Option
#% guisection: Optional
#% key: mswimhydrotopespath
#% type: string
#% multiple: no
#% key_desc: name
#% gisprompt: old,file,file
#% description: Path to custom m.swim.hydrotopes script or if not installed (e.g. in testing)
#%end

from __future__ import print_function
import os
import datetime as dt
import numpy as np


import grass.script as grass
import grass.script.core as gcore
gm = grass.message


def grun(*args, **kwargs):
    if 'quiet' not in kwargs:
        kwargs['quiet'] = True
    return grass.run_command(*args, **kwargs)


def routing(units, accumulation, outname, searchradius=1.5):
    grun('r.mask', raster=units, overwrite=True)

    # make outlet raster points, border subbasins with neg accumul. not correct
    grass.message('Searching outlets...')
    grun('r.statistics', base=units, cover=accumulation, method='max',
         output='maxaccum__')
    exp = "outlets__=if('%s' == @maxaccum__,%s,null())" % (accumulation, units)
    grass.mapcalc(exp)

    grass.message('Growing outlets...')
    # join outlets (set all to 1) and then grow outlets
    grass.mapcalc("outlets1__=if(isnull(outlets__),null(),1)", overwrite=True)
    # grow to all 8 neighbours (1.5 times cell size)
    grun('r.grow', input='outlets1__',
         output='outlets__grown', radius=searchradius)
    grun('r.clump', input='outlets__grown', output='outlets__grown__clumped')

    # make inlets
    grass.message('Searching inlets...')
    grun('r.statistics', base='outlets__grown__clumped', cover=accumulation,
         method='max', output='maxoutlet__clumps')
    grass.mapcalc("inlets__=if(%s == @maxoutlet__clumps,%s,null())" %
                  (accumulation, units))

    # transfer inlet units to clumps
    grun('r.stats.zonal', base='outlets__grown__clumped', cover='inlets__',
         method='max', output='clumps__unitID')
    # transfer to original units again (if n outlets >1, take largest inlet ID)
    grun('r.stats.zonal', base='outlets__', cover='clumps__unitID',
         method='max', output='next__unitID__outlets')
    grun('r.stats.zonal', base=units, cover='next__unitID__outlets',
         method='max', output='nextID__units__float')
    grass.mapcalc(outname + '=int(nextID__units__float)')

    return


def valley(slope, output, thresh=15, filter=300):
    '''Make a clean valley=1 and slope/rockface=0 raster'''
    # make valley 1/0 raster
    grass.mapcalc(exp=output+'__=%s<=%s' % (slope, thresh))
    # clean
    print('Cleaning...')
    grun('r.resamp.filter', input=output+'__', output=output+'__smooth',
         filter='box', radius=filter)
    grass.mapcalc(exp=output+'=int(round(%s__smooth))' % output)

    return output


def cleanRast(input, output, areathresh, fill=True, lenthresh=None,
              quiet=True):
    # RASTER BASED ONLY
    grun('r.clump', input=input, output=output + '__clumped')
    grass.mapcalc('help__1=1')
    grun('r.stats.zonal', base=output + '__clumped', cover='help__1',
         method='sum', output='n__cells')
    # output by thresholding cells
    mincells = int(areathresh / float(grass.region()['nsres'])**2)
    if not quiet:
        print('Removing all units with less than %s cells' % mincells)
    interout = output + '__filtered'
    grass.mapcalc(exp='%s=if(n__cells>=%s,%s,null())' %
                  (interout, mincells, input))

    if fill:
        # close/fill gaps
        grun('r.grow', input=interout, output=output + '__float',
             radius=mincells)
        interout = output + '__float'
    # make int
    grass.mapcalc(output + '=int(%s)' % interout)

    tmprast = ['help__1', output + '__clumped', 'n__cells', interout]
    grun('g.remove', type='rast', name=tmprast, flags='f')

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


def zonalStats(rast, zones, **runivarargs):
    '''Copied from mwpy.grass !!!
    Return statistics of the raster for the given zones as a pd.Dataframe.
    '''
    # make flags
    if 'flags' in runivarargs:
        runivarargs['flags'] += 't'
    else:
        runivarargs['flags'] = 't'
    # call r.univar
    lines = grass.read_command(
        'r.univar', map=rast, zones=zones, **runivarargs).split('\n')
    # convert to str array
    stats = np.array([tuple(s.split('|')) for s in lines[1:-1]], dtype=str)
    colnames = lines[0].split('|')
    # parse cols as either int,float or str
    cols = []
    for c in stats.T:
        try:
            cols += [list(map(int, c))]
        except ValueError:
            try:
                cols += [list(map(float, c))]
            except ValueError:
                cols += [c]
    # make dataframe and set index to zone
    df = pd.DataFrame(zip(*cols), columns=colnames)
    df.set_index('zone', inplace=True)
    # convert cells to area
    reg = grass.region()
    df['non_null_area'] = df['non_null_cells'] * reg['ewres'] * reg['nsres']
    df['null_area'] = df['null_cells'] * reg['ewres'] * reg['nsres']
    return df


def cleanGunits(gunits, subbasins, outname, smallarea):

    # get area and subbasins
    gm('Reading gunits for all subbasins...')
    gu = zonalStats(subbasins, gunits)
    # largest gunit
    maxu = gu['non_null_cells'].max()

    gm('Largest unit: %s cells' % maxu)

    # subbasins
    sbs = np.unique(gu['mean'])
    nsb = len(sbs)
    cunits = []
    maxn = 0
    grun('r.mask', flags='r')
    gm('Looping over %s subbasin...' % nsb)
    for i, s in enumerate(sbs):
        grun('r.mask', raster=subbasins, maskcats=s)
        grun('g.region', zoom='MASK')

        outn = 'clean__gu__%04i' % s

        cleanRast(gunits, outn+'__first__iteration', smallarea)
        # second iteration
        cleanRast(outn+'__first__iteration', outn+'__second__iteration',
                  smallarea)
        # add maxid to avoid neighbouring merge when patching
        grass.mapcalc('{0}={0}__second__iteration+{1}'.format(outn, maxn))
        # save for patching
        cunits += [outn]
        maxn = int(grass.raster_info(outn)['max']) + 1
        # tidy
        grun('r.mask', flags='r')
        grun('g.region', raster=subbasins)
        grun('g.remove', type='raster', flags='f',
             name=[outn+'__first__iteration', outn+'__second__iteration'])
        print('%s/%s' % (i+1, nsb))

    # patch together
    if len(cunits) < 500:
        grun('r.patch', input=','.join(cunits), output=outname+'__patched')
    else:
        # patch in batches due to open files limit of r.patch
        patched = []
        batches = range(0, len(cunits), 500) + [len(cunits)]
        for i in range(len(batches) - 1):
            n = outname + '__%s' % i
            grun('r.patch', input=','.join(
                cunits[batches[i]:batches[i + 1]]), output=n)
            patched += [n]
        grun('r.patch', input=','.join(patched), output=outname+'__patched')
    grass.mapcalc(exp='%s=int(%s)' % (outname, outname+'__patched'))

    # tidy
    grun('g.remove', type='raster', pattern='clean__gu__*', flags='f')
    grun('g.remove', type='raster', name=outname+'__patched', flags='f')
    return


def checkRouting(nextID, gunits, outname):
    import mygrass as mg

    tbl = grass.read_command('r.univar', map=nextID,
                             zones=gunits, flags='t').split('\n')
    cols = tbl[0].split('|')
    tbl = pd.DataFrame([tuple(l.split('|')) for l in tbl[1:-1]], columns=cols)
    tbl.set_index(tbl['zone'].astype(int), inplace=True)
    tbl = tbl[['mean', 'non_null_cells']].astype(int)
    # group by 'inlet ids' and accum cells
    incells = tbl.groupby('mean').agg('sum')['non_null_cells']
    receivingcells = tbl.ix[incells.index]['non_null_cells']
    # ratio
    routratio = incells / receivingcells
    print(routratio.describe())
    mg.areclass(gunits, np.array(zip(routratio.index, routratio)), outname)

    return 0


def aspect(elevation, output, nclasses):
    '''Make nice aspect classes: 360 should be devisible by nclasses'''
    # ASPECT
    grun('r.slope.aspect', elevation=elevation,
         aspect='aspect__', precision='CELL')
    # resulting map: degrees from E ccw, ie. E=0, N=270

    # make aspect classes
    step = int(360 / nclasses)
    start = 270 + int(step / 2)  # centered around north

    # from start to 0 (-1 for python)
    breaks = range(start, -1, -step)
    nbr = len(breaks)
    exp = ['if(aspect__ <= %s & aspect__ > %s,%s,0)' % (
        breaks[i], breaks[i] - step, i + 1) for i in range(nbr)][:-1]
    # from start to 360
    breaks2 = range(start, 360, step)[::-1]
    # if not breaking at 0deg/E, add class going across 0
    if breaks[-1] != 0:
        exp += ['if(aspect__ <= %s | aspect__ > %s,%s,0)' %
                (breaks[-1], 360 - breaks[-1], nbr)]
        breaks2 = breaks2[1:]

    exp += ['if(aspect__ <= %s & aspect__ > %s,%s,0)' % (
            breaks2[-1] + step, breaks2[-1], len(exp) + nbr)
            for i in range(len(breaks2))]

    # add all rules together
    exp = output + '=' + ('+'.join(exp))
    gm(exp)
    grass.mapcalc(exp=exp, overwrite=True)

    return output


def hydAddressRast(hydrotopes, subbasins, output):
    '''Create a hydrotope address map, i.e. the n'th hydrotope for each subbasin
    as they appear in the structure file.
    '''
    # get subbasin for each hydrotope
    tbl = grass.read_command('r.univar', map=subbasins, zones=hydrotopes,
                             flags='gt').split('\n')[:-1]  # :-1 as last line hast line break]
    tbl = [tuple(l.split('|')) for l in tbl]
    tbl = np.array(tbl[1:], dtype=list(zip(tbl[0], ['S250'] * len(tbl[0]))))
    sub = np.array(tbl['mean'], dtype=int)
    subIDs = np.unique(sub)
    # get counts
    freq, bins = np.histogram(sub, len(subIDs))
    # cumulative and shift by one subbasin
    cumfreq = np.array([subIDs, [0] + freq.cumsum().tolist()[:-1]]).T
    # reclass subbasin map to this count
    tmpf = grass.tempfile()
    np.savetxt(tmpf, cumfreq, delimiter='=', fmt='%i')
    grass.run_command('r.reclass', input=subbasins,
                      output='subbasins__cumfreq', rules=tmpf)
    # subtract from hydrotopes
    grass.mapcalc('%s=%s - subbasins__cumfreq' % (output, hydrotopes))
    # remove help stuff
    grass.run_command('g.remove', type='rast',
                      name='subbasins__cumfreq', flags='f')
    grass.message('Created hydrotope address raster: %s' % output)
    return


def classedContours(elevation, valleys, output, valinter=40, slopeinter=400):
    '''Make nice aspect classes: 360 should be devisible by nclasses'''
    val = valleys

    # get stats of DEM for nice breaks
    stats = grass.parse_command('r.univar', map=elevation, flags='g')
    # function to return the r.mapcalc condition with nice elevation breaks

    def niceBreaks(interval, last):
        # make expression
        exp = 'int({elev}/{interval})+{last}'.format(elev=elevation, last=last,
                                                     interval=interval)
        # calculate next last/highest value
        last = int(float(stats['max']) / interval) + last
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
    vcon, last = niceBreaks(valinter, 0)
    scon, last = niceBreaks(slopeinter, last)
    exp = output + '=if(%s,%s,%s)' % (val, vcon, scon)
#    # add all rules together
#    exp = output+'='+('+'.join(exp))
    gm(exp)
    grass.mapcalc(exp=exp, overwrite=True)

    return output


def writeStr(outname, unitmap, subbasins, hydAddress, nextID, *othercols):
    # get subbasin, hydAddress, nextID
    subb = meanUnit(unitmap, subbasins, int)
    subbasins = subbasins.split('@')[0]
    hydAdd = meanUnit(unitmap, hydAddress, int)
    nxt = meanUnit(unitmap, nextID, int)

    # make sure none are flowing to unit 0
    nxt[nextID][nxt[nextID] == 0] = nxt['id'][nxt[nextID] == 0]

    # get other columns
    ocols = ()
    for c in othercols:
        if type(c) in [int, float]:
            ocols += (np.ones((len(subb),)) * c,)
        else:
            ocols += (meanUnit(unitmap, c, float)[c.split('@')[0]],)

    # get all columns
    data = list(zip(subb['id'], subb[subbasins], hydAdd[hydAddress],
                    nxt[nextID], *ocols))
    columnnames = ['gunitID', 'subbasinID', 'hydrotopeID', 'nextID']
    columnnames += [c.split('@')[0] if type(c) == str else 'raster_%s' % i
                    for i, c in enumerate(othercols)]
    data = np.array(data, dtype=list(zip(columnnames, 4 *
                                         [int] + len(othercols) * [float])))
#    # dont take those that are outflows
#    data = data[nxt['id']!=nxt[nextID]]
    # save to txt file
    f = open(outname, 'w')
    f.write((4 * '%-10s ' + len(othercols) * '%-14s ' + '\n') %
            tuple(columnnames))
    np.savetxt(f, data, fmt=4 * '%10i ' + len(othercols) * '%16.6f ')
    f.close()
    gm('Wrote %s' % outname)
    return data


def meanUnit(units, rast, outtype=str):
    '''Get mean values over units and return as array with two column
    (units int, rast with outtype type)
    '''
    tbl = grass.read_command('r.univar', map=rast, zones=units, flags='gt')
    # :-1 as last line hast line break]
    tbl = [tuple(l.split('|')) for l in tbl.split('\n')[:-1]]
    array = np.array(tbl[1:], dtype=list(zip(tbl[0], ['S250'] * len(tbl[0]))))
    out = np.array(list(zip(array['zone'], array['mean'],
                            array['non_null_cells'], array['null_cells'])),
                   dtype=[('id', int), (rast.split('@')[0], outtype),
                          ('non_null_cells', int), ('null_cells', int)])
    gm('Got %r x %s' % (out.dtype.names, out.shape[0]))
    return out


def loadGWE(file, prefix, reclassrast):
    # load gla file
    gla = pd.read_table(file, index_col=[0, 1], names=[
                        'id', 'value'], delim_whitespace=True)
    styear = 1979
    # colormap
    max = 2100
    colormap = '''-%s   magenta
                 -200 red
                   -1 orange
                    0 yellow
                    1 green
                  200 blue
                 %s   cyan''' % (max, max)
    # loop over indeces (year,index) and reclass units raster
    dom = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    outmaps = []
    for l1, l2 in np.unique(gla.index.values):
        outname = prefix + '_%04i-%02i' % (l1, l2)
        grass.message('Processing %s ' % (outname))
        print(gla.ix[l1, l2].as_matrix().shape)
        mg.areclass(reclassrast, gla.ix[l1, l2].as_matrix(), outname)
        # nice colour map
        grass.write_command('r.colors', map=outname, rules='-', stdin=colormap)
        # timestamp
        grass.run_command('r.timestamp', map=outname,
                          # get last date of month by month+1-1day
                          date=dt.date(styear - 1 + l1, l2, dom[l2 - 1]).strftime('%d %b %Y'))
        outmaps += [outname]
    # make nice timeseries out of it
    grass.run_command('t.create', output=prefix,
                      title=prefix, description=prefix)
    grass.run_command('t.register', input=prefix, maps=','.join(outmaps))
    return


class Main:

    def __init__(self, **optionsandflags):
        '''Process all arguments and prepare processing'''
        # add all options and flags as attributes (only nonempty ones)
        self.options = {}
        for o in optionsandflags:
            if optionsandflags[o] != '':
                self.options[o] = optionsandflags[o]
        self.__dict__.update(self.options)

        # input checking
        assert self.hemisphere in ['N', 'S'], 'hemisphere must be N or S'

        if not hasattr(self, 'initialdebris'):
            self.initialdebris = 0

        convertnumber = ['minarea', 'naspectclasses', 'valleythreshold',
                         'valleysmoothing', 'contoursvalleys',
                         'contoursslopes', 'avalanchethreshold']
        convertnumber = [n for n in convertnumber if hasattr(self, n)]
        for n in convertnumber:
            v = self.__dict__[n]
            try:
                self.__dict__[n] = int(v)
            except ValueError:
                try:
                    self.__dict__[n] = float(v)
                except ValueError:
                    raise 'Cant convert %s to a number. %s' % (n, v)

        # defaults
        self.hydrotope_address = self.hydrotopes + '__address'
        self.sunhours_summer = self.sunhoursprefix + '_summer'
        self.sunhours_winter = self.sunhoursprefix + '_winter'
        self.avalanchearea_frac = self.avalanchearea + '__frac'
        self.strcolumns = [self.gunits,
                           self.subbasins,
                           self.hydrotope_address,
                           self.downstreamgunits,
                           0,  # initial glacier thickness must be 0
                           self.slope_dg,
                           self.avalanchearea_frac,
                           self.sunhours_summer,
                           self.sunhours_winter,
                           self.initialdebris,
                           0]  # initial snow
        if hasattr(self, 'strfilepath'):
            self.glacier_structure_file = os.path.join(
                            os.path.dirname(self.strfilepath), 'glaciers.str')
        return

    def _raster_exists(self, name):
        if len(name.split('@')) == 1:
            mapset = grass.gisenv()['MAPSET']
        else:
            name, mapset = name.split('@')

        rasters = grass.read_command('g.list', type='raster', mapset=mapset)
        return name in rasters.split()

    def process(self):
        gm('Cleaning and zooming to glacier area...')
        gm('''Resolution of subbasins raster is used for all output with
           minimal non-null extent of glacierarea raster.''')

        if not self.s:
            self.create_garea()

        if not self.d and not self.s:
            gm('Deriving statistics from DEM...')
            self.process_dem()

        if not self.s:
            gm('Creating glacier units...')
            self.create_units()

        if hasattr(self, 'routingratio'):
            gm('Calculating routing ratio...')
            self.check_routing()

        if not self.s:
            gm('Creating hydrotopes...')
            self.create_hydrotopes()

        if hasattr(self, 'strfilepath'):
            gm('Wrting glacier structure file...')
            self.glacier_structure()

        # clean
        grun('r.mask', flags='r')
        if not main.k:
            grass.run_command('g.remove', type='raster,vector', pattern='*__*',
                              flags='fb', quiet=True)

        return

    def mask(self, raster):
        # zoom to area and mask
        grun('g.region', raster=raster)
        grun('r.mask', raster=raster, overwrite=True)
        return

    def create_garea(self):
        # mask entire catchment
        self.mask(self.subbasins)
        # clean small bits from garea
        garea_clean = self.gunitsglacierarea + '__clean'
        cleanRast(self.glacierarea, garea_clean, self.minarea, fill=False)
        grun('g.region', zoom=garea_clean)
        self.mask(garea_clean)
        self.glacierarea_clean = garea_clean + '__zoomed'
        grass.mapcalc(self.glacierarea_clean + '=' + garea_clean)
        return

    def process_dem(self):
        """Create some DEM derived maps."""
        # make slope and aspect raster
        grun('r.slope.aspect', elevation=self.elevation, precision='CELL',
             slope=self.slope_dg, aspect=self.aspect, quiet=False)

        doysum, doywin = (172, 355) if self.hemisphere == 'N' else (355, 172)

        # calculate sun hours on June 21 and Dec 21
        grun('r.sun', elevation=self.elevation, day=doysum,
             insol_time=self.sunhours_summer, quiet=False)
        grun('r.sun', elevation=self.elevation, day=doywin,
             insol_time=self.sunhours_winter, quiet=False)
        return

    def create_units(self):
        """
        Main processing function to create the glacier units. Will be run
        every time.
        """
        # distinguish valleys and slopes (with smoothing)
        valleys_smooth = self.valleys + '__smooth'
        valley(self.slope_dg, valleys_smooth,
               self.valleythreshold, self.valleysmoothing)
        # clean
        cleanRast(valleys_smooth, self.valleys, self.minarea)

        # make contours (also writes slope and valleys)
        contours_classed = self.glaciercontours + '__classed'
        classedContours(self.elevation, self.valleys, contours_classed,
                        self.contoursvalleys, self.contoursslopes)
        # clean
        cleanRast(contours_classed, self.glaciercontours, self.minarea)

        # make aspect classes (only on slopes)
        aspect_classed = self.slopeaspect + '__classed'
        aspect(self.elevation, aspect_classed, self.naspectclasses)
        # only take those not in valleys and rest null
        aspect_slopes = aspect_classed + '__slopes'
        grass.mapcalc(aspect_slopes + '=int(if(!%s,%s,0))' % (self.valleys,
                                                              aspect_classed))
        # filter out small slopes and add them to neighbouring slopes
        cleanRast(aspect_slopes, self.slopeaspect, self.minarea)

        # make raw glacier units
        gunits_maps = [self.subbasins, self.glaciercontours, self.slopeaspect]
        gunits_raw = self.gunits + '__raw'
        grun('r.cross', input=','.join(gunits_maps), output=gunits_raw)
        # spatially explicit
        gunits_clumped = self.gunits + '__clumped'
        grun('r.clump', input=gunits_raw, output=gunits_clumped)

        # clean gunits and make spatially explicit
        # then clean each subbasin individually to fill each subbasin
        gunits_clean = self.gunits + '__clean'
        cleanGunits(gunits_clumped, self.subbasins, gunits_clean, self.minarea)
        # mask garea again (cleanGunits loses mask)
        self.mask(self.glacierarea_clean)
        # remove small overspilled units at glacierarea fringes and reduce
        # glacierarea accordingingly
        gunits_cleaner = self.gunits + '__cleaner'
        cleanRast(gunits_clean, gunits_cleaner, self.minarea,
                  quiet=False, fill=False)
        exp = 'if(!isnull(%s),1,null())' % gunits_cleaner
        grass.mapcalc(self.gunitsglacierarea + '=' + exp)
        grun('r.mask', raster=self.gunitsglacierarea, overwrite=True)
        # clump again to renumber (enhancing performance, would also work w/o)
        grun('r.clump', flags='d', input=gunits_cleaner, output=self.gunits)

        # routing
        routing(self.gunits, self.accumulation, self.downstreamgunits)

        return

    def check_routing(self):
        # double check
        checkRouting(self.downstreamgunits, self.gunits, self.routingratio)
        return

    def map_gunits(self, inrast, outrast):
        grun('r.statistics', base=self.gunits, cover=inrast,
             method='mode', output=outrast + '__labeled')
        grass.mapcalc(outrast + '=int(@%s)' % (outrast + '__labeled'))
        return

    def create_hydrotopes(self):
        """
        Create hydrotopes that include the glacier units and get addresses.
        """
        # make glacier units part of hydrotope map and include contours outside
        # of glacier area
        self.mask(self.subbasins)

        # create contours with glacier units
        exp = 'if(isnull($gunits),int($elevation/$interval),$gunits)'
        grass.mapcalc(self.contourrast + '=' + exp, gunits=self.gunits,
                      elevation=self.elevation, interval=self.contours)

        # Prepare hydrotope maps
        # map soils to gunits
        gunits_soil = self.gunits + '__soil'
        self.map_gunits(self.soil, gunits_soil)

        if self.a:  # add slope soilID
            # get valley/slope glacier units
            gunits_soilslopes = gunits_soil + '__slopes'
            gunits_valleys = self.gunits + '__valleys'
            self.map_gunits(self.valleys, gunits_valleys)
            # get next higher soil for slopes
            ssid = int(grass.raster_info(self.soil)['max'])+1
            gm('Assigning soilID=%s to slopes inside glacier area.' % ssid)
            grass.mapcalc('$output=if($gvalleys,$gsoil,$slopesoil)',
                          output=gunits_soilslopes, gvalleys=gunits_valleys,
                          gsoil=gunits_soil, slopesoil=ssid)
            gunits_soil = gunits_soilslopes

        ex = '$output=if(isnull($garea),$soil,$gsoil)'
        grass.mapcalc(ex, output=self.gunitssoil, garea=self.gunitsglacierarea,
                      soil=self.soil, gsoil=gunits_soil)

        # map landuse to gunits
        gunits_landuse = self.gunits + '__landuse'
        self.map_gunits(self.landuse, gunits_landuse)
        exp = 'if(!isnull($garea),int($glanduse),$landuse)'
        grass.mapcalc(self.gunitslanduse + '=' + exp,
                      garea=self.gunitsglacierarea,
                      glanduse=gunits_landuse, landuse=self.landuse)

        # make hydrotopes and write str file
        # all these are options of this module that are parsed on to
        # m.swim.hydrotopes if given
        optionkeys = ['subbasins', 'elevation', 'hydrotopes', 'strfilepath',
                      'wetland', 'managment']
        options = {i: self.__dict__[i]
                   for i in optionkeys
                   if hasattr(self, i)}
        options.update(dict(
            flags='k',
            contours=self.contourrast,  # contains gunits
            landuse=self.gunitslanduse,
            soil=self.gunitssoil,
        ))
        gm('Running m.swim.hydrotopes...')
        if hasattr(self, "mswimhydrotopespath"):
            mshp = self.mswimhydrotopespath
            gm("...using custom m.swim.hydrotopes %s" % mshp)
        else:
            mshp = 'm.swim.hydrotopes'
        grun(mshp, quiet=False, **options)

        return

    def create_avalanche(self):
        """
        Create the avalanche area if necessary and the avalanche fraction.
        """
        if hasattr(self, 'avalanchethreshold'):
            exp = '$output=if($slope >= $threshold, 1,null())'
            grass.mapcalc(exp=exp, slope=self.slope_dg,
                          output=self.avalanchearea,
                          threshold=self.avalanchethreshold)
        elif not self._raster_exists(self.avalanchearea):
            grass.fatal('''avalanchearea raster must exist if
                        avalanchethreshold is not given.''')
        # make fraction
        grun('r.stats.zonal', base=self.gunits, cover=self.avalanchearea,
             method='sum', output='n__avalanche__cells')
        grass.mapcalc(exp='help__1=1', overwrite=True)
        grun('r.stats.zonal', base=self.gunits, cover='help__1', method='sum',
             output='n__cells', overwrite=True)
        grass.mapcalc(exp="$output=n__avalanche__cells/n__cells",
                      output=self.avalanchearea_frac)
        return

    def glacier_structure(self):
        """Prepare all input and write the glacier.str file."""

        # hydAddress needs to be for all subbasins
        self.mask(self.subbasins)
        hydAddressRast(self.hydrotopes, self.subbasins, self.hydrotope_address)

        self.mask(self.gunitsglacierarea)
        self.create_avalanche()

        # check if all input exists
        maps = [i for i in self.strcolumns if type(i) == str]
        for m in maps:
            if not self._raster_exists(m):
                grass.fatal('%s does not exist to write structure file.' % m)

        writeStr(self.glacier_structure_file, *self.strcolumns)
        return


if __name__ == '__main__':
    st = dt.datetime.now()
    # get options and flags
    o, f = grass.parser()

    def fmt(d): return '\n'.join(['%s: %s' % (k, v)
                                  for k, v in d.items()]) + '\n'
    grass.message('GIS Environment:\n' + fmt(grass.gisenv()))
    grass.message('Parameters:\n' + fmt(o) + fmt(f))

    try:
        import pandas as pd
    except ModuleNotFoundError:
        raise ImportError('Cant import pandas. Is it installed?')

    # Parameters
    # glacierarea = 'glacierarea'  # existing raster
    # accumulation = 'accumulation'  # input accumulation map
    # minarea = 100000  # m2
    # valleythreshold = 19  # dg
    # valleysmoothing = 200  # some sort of radius
    # contoursvalleys = 30  # m
    # contoursslopes = 200  # m (maybe make dependent on mean valley/slope slopes?)
    # naspectclasses = 4  # number of aspect classes
    # hemisphere = 'N'  # N | S for summer/winter
    # initialdebris = 0  # optional raster

    # ## Output
    # d = False  # skip recreating slope_dg, aspect and sunhours (must all exist)
    # gunits = 'gunits'  # glacier units raster
    # downstreamgunits = 'downstreamgunits'  # raster with downstream gunit id
    # slope_dg = 'slope_dg'  # slope in degrees
    # aspect = 'aspect'  # aspect CELL precision
    # sunhoursprefix = 'sun_hr'  # prefix for r.sun Output
    # gunitsglacierarea = 'gunitsglacierarea'  # glacier area raster fitting to gunits
    # valleys = 'valleys'  # valley/slope raster
    # glaciercontours = 'glaciercontours'  # valley/slope classed contours
    # slopeaspect = 'slopeaspect'  # slope aspects in naspectclasses
    # routingratio = 'routingratio'  # optional: number of cells upstream/downstream
    # gunitssoil = 'gunitssoil'  # soils mapped onto glacier units
    # gunitslanduse = 'gunitslanduse'  # landuse mapped onto glacier units

    ############################################################

    keywords = o
    keywords.update(f)
    main = Main(**keywords)

    main.process()

    # report time it took
    delta = dt.datetime.now() - st
    grass.message('Execution took %s hh:mm:ss' % delta)
