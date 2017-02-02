#!/usr/bin/env python
#
############################################################################
#
# MODULE:      m.swim.climate v1.0
# AUTHOR(S):   Michel Wortmann, wortmann@pik-potsdam.de
# PURPOSE:     Preprocessing suit for the Soil and Water Integrated Model (SWIM)
# COPYRIGHT:   (C) 2012-2016 by Wortmann/PIK
#
#              This program is free software under the GNU General Public
#              License (>=v2). Read the file COPYING that comes with GRASS
#              for details.
#
#############################################################################

#%Module
#% description: Soil and Water Integrated Model (SWIM) climate input preprocessor
#% keywords: hydrological modelling, SWIM, climate interpolation
#%End

#%Option
#% guisection: Subbasins
#% key: subbasins
#% type: string
#% multiple: no
#% key_desc: name
#% label: Subbasins vector or raster
#% description: Subbasins for which climate data will be interpolated
#% gisprompt: old,raster,raster
#%end

#%Option
#% guisection: Subbasins
#% key: elevation
#% type: string
#% multiple: no
#% key_desc: name
#% description: Name of column with subbasin elevation
#%end

#%Option
#% guisection: Climate data
#% key: climstations
#% type: string
#% multiple: no
#% key_desc: name
#% label: Climate stations vector
#% description: Climate stations from which to interpolate
#% gisprompt: old,vector,vector
#%end

#%Option
#% guisection: Climate data
#% key: stationelevation
#% type: string
#% multiple: no
#% key_desc: name
#% description: Name of column with elevations of climate stations
#%end

#%Option
#% guisection: Climate data
#% key: fnames
#% type: string
#% multiple: no
#% key_desc: name
#% description: Name of column with filenames/ids (filenames: fnames.ext)
#%end

#%Option
#% guisection: Climate data
#% key: ext
#% type: string
#% multiple: no
#% key_desc: name
#% description: Extension of climate input files
#%end

#%Option
#% guisection: Climate data
#% key: datadir
#% type: string
#% multiple: no
#% key_desc: path
#% description: Folder with climate data files
#% gisprompt: old,dir,dir
#% answer: .
#%end

#%Option
#% guisection: Climate data
#% key: nodata
#% type: string
#% multiple: no
#% key_desc: value
#% description: No data value
#%end

#%flag
#% key: e
#% label: Elevation correction with daily lapse rates
#% guisection: Interpolation
#%end

#%flag
#% key: p
#% label: Precipitation correction
#% guisection: Interpolation
#%end

#%Option
#% guisection: Interpolation
#% key: method
#% type: string
#% multiple: no
#% key_desc: method
#% description: Interpolation method
#% options: IDW, Kriging
#%end

#%Option
#% guisection: Interpolation
#% key: start
#% type: string
#% multiple: no
#% key_desc: DD.MM.YYYY
#% description: Start date
#%end
#%Option
#% guisection: Interpolation
#% key: end
#% type: string
#% multiple: no
#% key_desc: DD.MM.YYYY
#% description: End date
#%end

#%Option
#% guisection: Interpolation
#% key: minnb
#% type: integer
#% multiple: no
#% key_desc: int
#% description: Minimum neighbours
#%end

#%Option
#% guisection: Interpolation
#% key: maxnb
#% type: integer
#% multiple: no
#% key_desc: int
#% description: Maximum neighbours
#%end

#%Option
#% guisection: Interpolation
#% key: maxdist
#% type: double
#% multiple: no
#% key_desc: float
#% description: Maximum distance of search radius [km]
#%end

#%Option
#% guisection: Output
#% key: outdir
#% type: string
#% multiple: no
#% key_desc: path
#% label: Path where all output will be written
#% description: Output files written: stationfile.dat, virtualstationfile.dat
#% gisprompt: old,dir,dir
#% answer: .
#%end

#%Option
#% guisection: Grid
#% key: grid
#% type: string
#% multiple: no
#% required: no
#% key_desc: name
#% label: Name of new or old grid vector
#% description: Resulting climate grid area vector with centroids and table or if -d existing grid vector with lon & lat columns
#%end

#%Option
#% guisection: Grid
#% key: res
#% type: double
#% multiple: no
#% key_desc: float
#% answer: 0.5
#% description: grid size in degrees
#%end


#%Option
#% guisection: Grid
#% key: ncinfopath
#% type: string
#% required: no
#% multiple: no
#% key_desc: path/myproj.dat
#% description: path/name of the nc info file that will be written
#% gisprompt: new,file,file
#%end

#%Flag
#% guisection: Optional
#% key: d
#% label: Predefined grid, will be 'grown' onto subbasin area if smaller
#%end

#%Flag
#% guisection: Optional
#% key: k
#% label: Keep intermediate files (those named *__*)
#%end

import grass.script as grass
import numpy as np
import sys, os
import datetime as dt
grun = grass.run_command
gread= grass.read_command
gm   = grass.message

class Interpolation:
    def __init__(self,**optionsandflags):
        '''Process all arguments and prepare processing'''
        # add all options and flags as attributes (only nonempty ones)
        self.options = {}
        for o in optionsandflags:
            if optionsandflags[o]!='':
                try: self.options[o] = int(optionsandflags[o])             # int
                except ValueError:
                    try: self.options[o] = float(optionsandflags[o])       # float
                    except ValueError: self.options[o] = optionsandflags[o]# str
        self.__dict__.update(self.options)
        # save region for convenience
        self.region = grass.region()

        # INPUT CHECK

        # input dir
        if not os.path.exists(self.datadir): grass.fatal('%s doesnt exisist!' %self.datadir)

        # climate stations columns
        cols=grass.vector_columns(self.climstations)
        for c in [self.fnames, self.stationelevation]:
            if c not in cols: grass.fatal('Cant find %s in table %s' %(c,self.climatestations))
        # subbasins
        if self.elevation not in grass.vector_columns(self.subbasins):
            grass.fatal('Cant find %s in table %s' %(self.elevation,self.subbasins))

        # no extension
        if 'ext' not in self.options: self.ext=''

        # check if INTERPOL.PAR can be written
        self.par = True
        parneeded = ['method','start','end','minnb','maxnb','maxdist','nodata']
        if any([e not in self.options for e in parneeded]):
            grass.warning('''Won't write INTERPOL.PAR as any of these arguments
            is not set %s''' %(parneeded,))
            self.par = False
        return

    def writeClimStationfile(self):
        # get coordinates and table
        coor, a = vectCoordsTbl(self.climstations)
        # add filepath to name/2nd column
        fnames = []
        for i,n in enumerate(a[self.fnames]):
            fnames += [os.path.join(self.datadir,str(n)+self.ext)]
            if not os.path.exists(fnames[-1]): grass.warning('%s doesnt exisist!' %fnames[-1])
        # get needed cols
        cols = np.array(zip(coor['cat'],
                            fnames,coor['x'],coor['y'],
                            a[self.stationelevation]),
                        dtype=[('cat',int),('fnames','S1000'),('x',float),('y',float),('z',float)])
        # write file
        sfpath = os.path.join(self.outdir,'stationfile.dat')
        self.stfile = sfpath
        f = file(sfpath,'w')
        f.write('ID  FILE                       POINT_X     POINT_Y            ELEV\n')
        np.savetxt(f, cols, fmt='%6i  %-'+str(len(self.datadir)+32)+'s%12.1f%12.1f%10.1f')
        gm( 'Saved climate stations to %s' %sfpath)

        return

    def writeVirtualStationfile(self):
        # get coordinates and table
        coor, a = vectCoordsTbl(self.subbasins)
        coor['z'] = a[self.elevation]

        # save to file
        sfpath = os.path.join(self.outdir,'virtualstationfile.dat')
        f = file(sfpath,'w')
        f.write('id    x     y      z\n')
        np.savetxt(f, coor, fmt='%6i %12.1f%12.1f%10.1f')
        gm( 'Saved virtual stations to %s' %sfpath)
        self.vstfile=sfpath
        return

    def stationDistance(self):
        '''Calculate distances from all virtualstations to all climate stations
        within the max distance threshold'''
        # create table with distances
        dist = gread('v.distance', flags='pa', from_=self.subbasins, from_type='centroid',
               to=self.climstations, dmax=self.maxdist*1e3, upload='dist')
        dist = [tuple(l.split('|')) for l in dist.split()]
        dist = np.array(dist[1:],dtype=zip(['from_cat','to_cat','distance'],[int,int,float]))
        n = len(np.unique(dist['from_cat']))
        if n==0: grass.fatal('No subbasins have corresponding climstation within the maxdist.')
        grass.warning('%s subbasins have %3.1f climstations on average within the maxdist.' \
                      %(n,len(dist)/float(n)))
        # m to km
        dist['distance'] = dist['distance']*1e-3
        # check n station and mean distance
        aggregated = []
        for i in np.unique(dist['from_cat']):
            subb = dist[dist['from_cat']==i]
            aggregated += [(i,len(subb),subb['distance'].min(),
                            subb['distance'].mean(),subb['distance'].max())]
        # make rec array
        aggregated = np.array(aggregated,dtype=[('subbasinID',int),('n',int),
                                ('min',float),('mean',float),('max',float)])
        # save
        statspath = os.path.join(self.outdir,'stationstats.dat')
        np.savetxt(statspath,aggregated,fmt='%7i %7i %9.1f %9.1f %9.1f',
                   header='subbasin n      min       mean      max      ')
        # report
        gm('Wrote table with number of climstations, min, mean, max distance to %s'%statspath)
        gm('Summary statistics:')
        gm('Stats  min    mean    max')
        for i,u in zip(['n','min','mean','max'],['stations']+3*['km']):
            gm(('%-4s'+3*' %7.1f'+' %s')%(i,aggregated[i].min(),aggregated[i].mean(),
                                    aggregated[i].max(),u))

        return aggregated

    def writePARfile(self):

        f = file(os.path.join(self.outdir,'INTERPOL.PAR'),'w')
        f.write('ClimateStationFile %r \n' %self.stfile)
        f.write('VirtualStations %r \n' %self.vstfile)
        f.write('OutputPath %r \n' %(self.outdir+'/'))
        f.write('PrecCorr %s \n' %(str(bool(self.p)).upper()))
        f.write('NoDataValue %r \n' %str(self.nodata))
        f.write('BeginDate %r \n' %self.start)
        f.write('EndDate %r \n' %self.end)
        # Method
        method = {'IDW':1, 'Kriging': 3}[self.method]
        if self.e: method += 1
        f.write('InterpolationMethod %s \n' %(' '.join(6*[str(method)])))
        f.write('NumNeighboursMax %s \n' %self.maxnb)
        f.write('NumNeighboursMin %s \n' %self.minnb)
        f.write('Maxdist %s \n' %(self.maxdist*1e3))

        f.write('''distributed  FALSE \nSubbasinMap ""  \nDemMap ""''')


def vectCoordsTbl(vect):
    '''Return the v.report table with coordinates as an array with
    column names in the entries for a point vector'''
    # get table, returned columns:
    # table columns, x, y, z
    t=grass.read_command('v.report', map=vect, option='coor').split('\n')[:-1]
    t=[tuple(l.split('|')) for l in t]
    colnames = t[0]
    t=np.array(t[1:])
    # set empty cell to nan
    t[t=='']=np.nan
    columns = []
    for c in t.T: # loop over columns
        try: columns += [np.array(c,int)]
        except ValueError:
            try: columns += [np.array(c,float)]
            except ValueError: columns += [c]
    # now fuse to make recArray of table and of coordinates
    tblcols  = columns[:-3]
    tbl   = np.array(zip(*tblcols), dtype=zip(colnames,[c.dtype for c in tblcols]))
    coorcols = [columns[0]]+columns[-3:]
    coor  = np.array(zip(*coorcols), dtype=zip(['cat','x','y','z'],[c.dtype for c in coorcols]))

    return coor,tbl


class Grid:
    def __init__(self,**optionsandflags):
        '''Process all arguments and prepare processing'''
        # add all options and flags as attributes (only nonempty ones)
        self.options = {}
        for o in optionsandflags:
            if optionsandflags[o]!='': self.options[o] = optionsandflags[o]
        self.__dict__.update(self.options)

        # get some location infos
        self.env=grass.gisenv()
        self.region=grass.region()
        self.proj=grass.parse_command('g.proj',flags='g')

        # convert res
        self.res = float(self.res)
        return

    def mkLonLatGrid(self):
        # make temporary roi
        roi='roi__'
        grass.run_command('v.in.region',output=roi,quiet=True)
        # create temporary lonlat location
        tmpdir=grass.tempdir()
        tmploc='lonlat'
        grass.core.create_location(tmpdir,tmploc,epsg=4326)
        grass.run_command('g.mapset',mapset='PERMANENT',location=tmploc,
                          dbase=tmpdir,quiet=True)
        # reproj roi, smax in meters = 200km per degree
        grass.run_command('v.proj',input=roi,mapset=self.env['MAPSET'],
                          location=self.env['LOCATION_NAME'],dbase=self.env['GISDBASE'],
                          quiet=True)
        grass.run_command('g.region',vector=roi,quiet=True)
        llregion = grass.region()

        # bounds to extend to next resolution break
        extent = {c:int(float(llregion[c])/self.res)*self.res for c in ('s','w')}
        extent.update({c:int((float(llregion[c])+self.res)/self.res)*self.res for c in ('n','e')})
        # set region
        grass.run_command('g.region',res=self.res,**extent)
        grass.message(('Lon/Lat extent of region:',extent))

        # make grid
        grass.run_command('v.mkgrid',map=self.grid, type='area')
        grass.run_command('v.db.addcolumn',map=self.grid,columns='lon double,lat double')
        grass.run_command('v.to.db', map=self.grid, type='centroid',
                          option='coor', columns='lon,lat',quiet=True)

        # back to origional location and reproj
        grass.run_command('g.mapset',mapset=self.env['MAPSET'],location=self.env['LOCATION_NAME'],
                          dbase=self.env['GISDBASE'],quiet=True)
        grass.run_command('v.proj',input=self.grid,mapset='PERMANENT',
                          location=tmploc,dbase=tmpdir,quiet=True,
                          smax=float(self.region['nsres'])+float(self.region['ewres']))
        return 0


    def getSubbasinGridProportions(self):
        # make raster
        grass.run_command('v.to.rast',input=self.grid,output=self.grid+'__',
                          use='cat',type='area')

        # grow grid if predefined
        if self.d:
            grass.run_command('r.grow.distance',input=self.grid+'__',
                              value = self.grid+'__grown')
            # make int and overwrite grid__
            grass.mapcalc(exp='{0}=int({1})'.format(self.grid+'__',self.grid+'__grown'))

        # cross subbasin raster and grid raster and get areas
        crossrast = self.subbasins+','+self.grid+'__'
        cross = grass.read_command('r.stats',flags='an',input=crossrast,
                                   separator='pipe').split()
        # cross = ( subbasinID, gridID, area )
        dtype = [('subbasinID',int),('gridID',int),('area',float)]
        tbl = np.array([tuple(l.split('|')) for l in cross],dtype=dtype)
        # sort by subbasin and area
        tbl.sort(order=('subbasinID','area'))

        return tbl

    def writeNCInfo(self,tbl):

        # get lon,lats (ditionary with cell ids as keys and [lon,lat] as entry
        lonlatmap = grass.vector_db_select(self.grid,columns='lon,lat')['values']

        # pickout lonlats for gridids
        lons = np.array([lonlatmap[i][0] for i in tbl['gridID']],float)
        lats = np.array([lonlatmap[i][1] for i in tbl['gridID']],float)

        # get proportions in each subbasin
        props = np.zeros(len(tbl))
        for sb in np.unique(tbl['subbasinID']):
            isb = tbl['subbasinID']==sb
            props[isb] = tbl['area'][isb]/tbl['area'][isb].sum()

        # make out array
        outtbl = np.column_stack((tbl['subbasinID'],lons,lats,props))
        # write out
        fmt = '%12i %12.3f %12.3f %12.6f'
        head= '%10s '%'subbasinID' + '%12s %12s %12s'%('lon','lat','weight')
        np.savetxt(self.ncinfopath,outtbl,fmt=fmt,header=head)
        grass.message('Wrote %s lines and %s columns to %s'%(outtbl.shape+(self.ncinfopath,)))
        return 0


if __name__=='__main__':
    st = dt.datetime.now()
    # get options and flags
    o, f = grass.parser()
    fmt = lambda d: '\n'.join(['%s: %s' % (k, v) for k, v in d.items()])+'\n'
    grass.message('GIS Environment:\n'+fmt(grass.gisenv()))
    grass.message('Parameters:\n'+fmt(o)+fmt(f))
    keywords = o; keywords.update(f)

    # decide what climate module
    if 'grid'!='':
        main=Grid(**keywords)

        # make grid
        if not main.d:
            main.mkLonLatGrid()
        # calculate proportional coverage with subbasins
        tbl = main.getSubbasinGridProportions()

        if 'ncinfopath' in main.options:
            main.writeNCInfo(tbl)

        # clean
        if not main.k:
            grass.run_command('g.remove',type='raster,vector', pattern='*__*',
                              flags='fb',quiet=True)
    else:
        main=Interpolation(**keywords)
        # process
        main.writeClimStationfile()
        main.writeVirtualStationfile()
        if main.par:
            main.writePARfile()
        if 'maxdist' in main.options:
            stats = main.stationDistance()

    # report time it took
    delta = dt.datetime.now()-st
    grass.message('Execution took %s hh:mm:ss' %delta)
