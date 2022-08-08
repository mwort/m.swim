#!/usr/bin/env python
#
############################################################################
#
# MODULE:      m.swim.climate v2.0
# AUTHOR(S):   Michel Wortmann, wortmann@pik-potsdam.de
# PURPOSE:     Preprocessing suit for the Soil and Water Integrated Model (SWIM)
# COPYRIGHT:   (C) 2012-2022 by Wortmann/PIK
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
#% guisection: Input
#% key: subbasins
#% type: string
#% multiple: no
#% key_desc: name
#% label: Subbasins raster
#% description: SWIM subbasins
#% required: yes
#% gisprompt: old,raster,raster
#%end

#%Option
#% guisection: Grid
#% key: grid
#% type: string
#% multiple: no
#% required: no
#% key_desc: name
#% required: yes
#% label: (Output/input) Name of new or old grid vector
#% description: Resulting climate grid area vector with centroids and table or if -d existing vector with lon & lat columns
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
#% key: lon_column
#% type: string
#% multiple: no
#% key_desc: column name
#% answer: lon
#% description: If grid is prescribed, use this column as lon value.
#%end

#%Option
#% guisection: Grid
#% key: lat_column
#% type: string
#% multiple: no
#% key_desc: column name
#% answer: lat
#% description: If grid is prescribed, use this column as lat value.
#%end

#%Option
#% guisection: Grid
#% key: lonlat_precision
#% type: integer
#% multiple: no
#% required: no
#% key_desc: name
#% answer: 3
#% label: Floating point precision of longitude and latitude.
#%end

#%Option
#% guisection: Output
#% key: gridfilepath
#% type: string
#% required: no
#% multiple: no
#% key_desc: path
#% description: path/name of the subbasin-grid file that will be written
#% gisprompt: new,file,file
#%end

#%Option
#% guisection: Output
#% key: ncinfopath
#% type: string
#% required: no
#% multiple: no
#% key_desc: path
#% description: Deprecated, use gridfilepath
#% gisprompt: new,file,file
#%end

#%Flag
#% guisection: Optional
#% key: d
#% label: Predefined grid, will be 'grown' onto subbasin area if smaller or the grid are points.
#%end

#%Flag
#% guisection: Optional
#% key: k
#% label: Keep intermediate files (those named *__*)
#%end

#%Flag
#% guisection: Optional
#% key: v
#% label: Show version and change/install date of this module and grass.
#%end


import sys
import grass.script as grass
import numpy as np
import datetime as dt
grun = grass.run_command
gread = grass.read_command
gm = grass.message

# cautious Alpha implementation of the mswim abstraction package
try:
    path = grass.utils.get_lib_path(modname='m.swim', libname='mswim')
    if path:
        sys.path.extend(path.split(':'))
        import mswim
    else:
        grass.warning('Unable to find the mswim python library.')
except Exception as e:
    grass.warning('An error occurred while loading the mswim python library.\n'+str(e))
    mswim = None


class Grid:
    def __init__(self,**optionsandflags):
        '''Process all arguments and prepare processing'''
        # add all options and flags as attributes (only nonempty ones)
        self.options = {}
        for o in optionsandflags:
            if optionsandflags[o]!='': self.options[o] = optionsandflags[o]
        self.__dict__.update(self.options)

        if "ncinfopath" in self.options and "gridfilepath" not in self.options:
            self.gridfilepath = self.ncinfopath
            import warnings
            warnings.warn(
                "ncinfopath is deprecated and will be removed, "
                "use gridfilepath instead.", DeprecationWarning)

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
        grass.run_command('g.region', vector=roi, res=self.res, flags='a')

        # make grid
        grass.run_command('v.mkgrid',map=self.grid, type='area')
        cols = '%s double,%s double' % (self.lon_column, self.lat_column)
        grass.run_command('v.db.addcolumn', map=self.grid, columns=cols)
        grass.run_command('v.to.db', map=self.grid, type='centroid',
                          option='coor', columns='lon,lat', quiet=True, overwrite=True)

        # back to origional location and reproj
        grass.run_command('g.mapset',mapset=self.env['MAPSET'],location=self.env['LOCATION_NAME'],
                          dbase=self.env['GISDBASE'],quiet=True)
        grass.run_command('v.proj',input=self.grid,mapset='PERMANENT',
                          location=tmploc,dbase=tmpdir,quiet=True,
                          smax=float(self.region['nsres'])+float(self.region['ewres']))
        return 0


    def getSubbasinGridProportions(self):
        # make raster
        grass.run_command('v.to.rast', input=self.grid, output=self.grid+'__',
                          use='cat', type="area,point")

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

    def write_grid_file(self,tbl):

        # get lon,lats (ditionary with cell ids as keys and [lon,lat] as entry
        cols = "%s,%s" % (self.lon_column, self.lat_column)
        lonlatmap = grass.vector_db_select(self.grid, columns=cols)['values']

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
        fmts = '%12i %12.{0}f %12.{0}f %12.6f'.format(self.lonlat_precision).split()
        cols = ['subbasin_id', 'lon','lat','weight']
        mswim.inout.write_csv(self.gridfilepath, outtbl,cols, min_width=12, column_formats=fmts)
        grass.message('Wrote %s lines and %s columns to %s'%(outtbl.shape+(self.gridfilepath,)))
        return 0


if __name__=='__main__':
    st = dt.datetime.now()
    # print version/date before doing anything else
    mswim.utils.print_version(__file__) if '-v' in sys.argv else None
    # get options and flags
    o, f = grass.parser()
    fmt = lambda d: '\n'.join(['%s: %s' % (k, v) for k, v in d.items()])+'\n'
    grass.message('GIS Environment:\n'+fmt(grass.gisenv()))
    grass.message('Parameters:\n'+fmt(o)+fmt(f))
    keywords = o; keywords.update(f)

    main=Grid(**keywords)

    # make grid
    if not main.d:
        main.mkLonLatGrid()

    # calculate proportional coverage with subbasins
    tbl = main.getSubbasinGridProportions()

    if 'ncinfopath' in main.options or 'gridfilepath' in main.options:
        main.write_grid_file(tbl)

    # clean
    if not main.k:
        grass.run_command('g.remove',type='raster,vector', pattern='*__*',
                          flags='fb',quiet=True)

    # report time it took
    delta = dt.datetime.now()-st
    grass.message('Execution took %s hh:mm:ss' %delta)
