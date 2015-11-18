#!/usr/bin/env python
#
############################################################################
#
# MODULE:      m.swim.climate.grid
# AUTHOR(S):   Michel Wortmann, wortmann@pik-potsdam.de
# PURPOSE:     Preprocessing suit for the Soil and Water Integrated Model (SWIM)
# COPYRIGHT:   (C) 2012 by Wortmann/PIK
#
#              This program is free software under the GNU General Public
#              License (>=v2). Read the file COPYING that comes with GRASS
#              for details.
#
#############################################################################
 
#%Module
#% description: Soil and Water Integrated Model (SWIM) climate input preprocessor for netCDF gridded climate input
#% keywords: hydrological modelling, SWIM, climate, grid, netCDF
#%End
 
#%Option
#% key: subbasins
#% type: string
#% multiple: no
#% required: yes
#% key_desc: name
#% label: Subbasins raster
#% description: Subbasins of grass project
#% gisprompt: old,raster,raster
#%end

#%Option
#% key: grid
#% type: string
#% multiple: no
#% required: yes
#% key_desc: name
#% label: Name of grid vector
#% description: Resulting climate grid area vector with centroids and table
#% gisprompt: new,vector,vector
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
#% guisection: Output
#% key: maxcells
#% type: integer
#% multiple: no
#% key_desc: int
#% answer: 4
#% description: Maximum cells to write into output
#%end

#%Option
#% guisection: Output
#% key: outpath
#% type: string
#% required: no
#% multiple: no
#% key_desc: path/myproj.dat
#% description: path/name of the nc info file that will be written
#% gisprompt: new,file,file
#%end

#%Flag
#% guisection: Optional
#% key: k
#% label: Keep intermediate files (those named *__*)
#%end

import os, sys
import datetime as dt
import numpy as np
import grass.script as grass




class main:
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
        
        # convert res and maxcells
        self.res = float(self.res)
        self.maxcells = int(self.maxcells)
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
        # get unique subbasins
        subb = np.unique(tbl['subbasinID'])
        
        areas  = np.zeros((len(subb),self.maxcells))
        gridids= np.zeros((len(subb),self.maxcells),int)
        many = {}
        # format
        for i,sbid in enumerate(subb):
            sb = tbl[['gridID','area']][tbl['subbasinID']==sbid][::-1]
            n = len(sb)
            # save those areas that are more than self.maxcells
            if n>self.maxcells: many[sbid]=sb
            # fit into data arrays
            nn = min(n,self.maxcells)
            areas[i,:nn]   = sb['area'][:nn]
            gridids[i,:nn] = sb['gridID'][:nn]
        
        # report on cell counts
        counts = (gridids>0).sum(axis=1)
        grass.message('Cells per subbasin [subbasin count]:')
        for i in range(1,self.maxcells+1):
            grass.message('%s: %s'%(i,len(counts[counts==i])))
        
        # report on areas not taken into account
        if len(many)>0:
            grass.message('There are %s subbasins with more than %s climate cells:'%(len(many),self.maxcells))
            grass.message('SubbasinID _ n grid cells _ proportion of subbasin (%)')
            for s in sorted(many.keys()):
                # proportion of areas more than self.maxcells in
                prop = many[s]['area'][self.maxcells:].sum() * 100/many[s]['area'].sum()
                print('%11i %14i %12.1f'%(s,len(many[s]),prop))
    
        # make ratios out of areas (must be transposed to make division work along columns)
        ratios = (areas.T/areas.sum(axis=1)).T
    
        # get lon,lats (ditionary with cell ids as keys and [lon,lat] as entry
        lonlats = grass.vector_db_select(self.grid,columns='lon,lat')['values']
        lonlats[0] = [0,0] # lonlat for grids with 0 weight
        # pickout lonlats for gridids
        lon = np.array([[lonlats[i][0] for i in l] for l in gridids],float)
        lat = np.array([[lonlats[i][1] for i in l] for l in gridids],float)
        
        # fill data array interweaving subb + lon1, lon1, ratio1 ...
        data = np.empty((len(subb),self.maxcells*3+1), dtype=float)
        data[:,0]    = subb
        data[:,1::3] = lon
        data[:,2::3] = lat
        data[:,3::3] = ratios
        # write out
        fmt = '%12i '+' '.join(['%12.3f %12.3f %12.6f']*self.maxcells)
        head= '%10s '%'subbasinID' + ' '.join(['%12s %12s %12s'%('lon','lat','weight')]*self.maxcells)
        np.savetxt(self.outpath,data,fmt=fmt,header=head)
        grass.message('Wrote %s lines and %s columns to %s'%(data.shape+(self.outpath,)))
        return 0


if __name__=='__main__':
    st = dt.datetime.now()
    # get options and flags
    o, f = grass.parser()
    grass.message(('GIS Environment:',grass.gisenv()))
    grass.message(('Parameters:',o,f))
    
    # send all to main
    keywords = o; keywords.update(f)
    main=main(**keywords)
    
    # make grid
    main.mkLonLatGrid()
    # calculate proportional coverage with subbasins
    tbl = main.getSubbasinGridProportions()
    
    if 'outpath' in main.options:
        main.writeNCInfo(tbl)

    # clean
    if not main.k:
        grass.run_command('g.remove',type='raster,vector', pattern='*__*',
                          flags='fb',quiet=True)

    # report time it took
    delta = dt.datetime.now()-st
    grass.message('Execution took %s hh:mm:ss' %delta)

