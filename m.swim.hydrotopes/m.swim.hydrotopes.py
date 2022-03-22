#!/usr/bin/env python
#
############################################################################
#
# MODULE:      m.swim.hydrotopes v1.6
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
#% description: Soil and Water Integrated Model (SWIM) preprocessor: hydrotopes
#%End
#%Option
#% guisection: Required
#% key: subbasin_id
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% description: Subbasin raster map
#% gisprompt: old,cell,raster
#%end
#%Option
#% guisection: Required
#% key: landuse_id
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% description: Land use raster map in SWIM land use classification
#% gisprompt: old,cell,raster
#%end
#%Option
#% guisection: Required
#% key: soil_id
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% description: Soil raster map classified to SWIM soil files
#% gisprompt: old,cell,raster
#%end
#%Option
#% guisection: Required
#% key: output
#% type: string
#% required: yes
#% multiple: no
#% key_desc: path
#% description: Path where structure file will be written
#% gisprompt: new,file,file
#%end
#%Option
#% guisection: Required
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
#% guisection: Contours
#% key: contours
#% type: string
#% required: no
#% multiple: no
#% key_desc: ints or list of ints or raster
#% description: Elevation contours to include in hydrotopes [optional, if set, DEM in Subbasins also needs to be set]
#%end
#%Option
#% guisection: Contours
#% key: contourrast
#% type: string
#% required: no
#% multiple: no
#% answer: contours
#% description: Resultant contours raster if contours option is int/list of ints
#%end
#%Flag
#% guisection: Contours
#% key: c
#% label: deprecated option, left for backwards compatibility
#%end
#%Option
#% guisection: Contours
#% key: elevation
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% description: Elevation raster
#% gisprompt: old,cell,raster
#%end
#%Option
#% guisection: Optional
#% key: crop_management_id
#% type: string
#% required: no
#% multiple: no
#% key_desc: name
#% description: Raster name
#% gisprompt: old,cell,raster
#%end
#%Option
#% guisection: Optional
#% key: wetland
#% type: string
#% required: no
#% multiple: no
#% key_desc: name
#% description: Raster name of 1/0 wetland mask
#% gisprompt: old,cell,raster
#%end
#%Option
#% guisection: Optional
#% key: glacier_thickness
#% type: string
#% required: no
#% multiple: no
#% key_desc: name
#% description: Raster name
#% gisprompt: old,cell,raster
#%end
#%Option
#% guisection: Optional
#% key: more
#% type: string
#% required: no
#% multiple: yes
#% key_desc: name
#% description: More rasters maps included in hydrotopes [optional]
#% gisprompt: old,cell,raster
#%end
#%Flag
#% guisection: Optional
#% key: k
#% label: Keep intermediat files (include __ in names)
#%end
#%Flag
#% guisection: Optional
#% key: v
#% label: Show version and change/install date of this module and grass.
#%end


import os,sys
from collections import OrderedDict
import grass.script as grass
g_run=grass.run_command
gm   =grass.message
import numpy as np
import datetime as dt

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


class main:

    def __init__(self,**optionsandflags):
        '''Process all arguments and prepare processing'''
        # add all options and flags as attributes (only nonempty ones)
        self.options = {}
        for o in optionsandflags:
            if optionsandflags[o]!='': self.options[o] = optionsandflags[o]
        self.__dict__.update(self.options)

        # dictionary of maps where the value will be averaged over all hydrotopes
        self.floatmaps = {}

        # default columns
        self.strcolumns = OrderedDict([(i, getattr(self, i)) for i in ('subbasin_id', 'landuse_id', 'soil_id')])

        # managment and wetland
        for r in ['crop_management_id','wetland']:
            self.__dict__[r] = self._maskOrFalse(r)
            if self.__dict__[r] != False:
                self.strcolumns[r] = self.__dict__[r]

        # check contours and create raster if needed
        if 'contours' in self.options:
            try:
                self.contours = list(map(int,self.contours.split(',')))
                if len(self.contours)==1: self.contours = self.contours[0]
            except ValueError:
                gm(('''Using %s as predefined contours.''' % self.contours))
                self.contourrast = self.contours
            else:
                # create contourrast
                self.mkContours()
            # add to columns
            self.strcolumns['elevation'] = self.contourrast
            self.floatmaps[self.contourrast] = self.elevation

        # glaciers
        self.glacier_thickness = self._maskOrFalse('glacier_thickness')
        if self.glacier_thickness != False:
            self.strcolumns['glacier_thickness']= self.glacier_thickness

        # other raster to include?
        if 'more' in self.options:
            self.more = self.more.split(',')
            self.strcolumns.update({k: k for k in self.more})

        # check if all input maps are int/CELL maps
        gm('Check if all input are integer raster...')
        for m in self.strcolumns.values():
            if grass.raster_info(m)['datatype']!='CELL':
                grass.fatal('%s is not an integer/CELL raster, convert using int() in r.mapcalc' %m)

        # check that all input maps have no NULLs in catchment over subbasins area
        g_run('r.mask', raster=self.subbasin_id, overwrite=True)

        gm('Checking for NULLs in input maps...')
        pow2 = 2**np.arange(len(self.strcolumns))
        ifnull = ' + '.join(['isnull(%s)*%s' %(m,i) for m,i in zip(self.strcolumns.values(), pow2)])
        grass.mapcalc('null__areas=%s' %ifnull, overwrite=True)
        g_run('r.null',map='null__areas',setnull=0,quiet=True)
        null = grass.parse_command('r.stats',input='null__areas',flags='aN',separator='=')

        if len(null)>0:
            gm("In any of the input maps are NULL/no data values over the subbasins area.")
            gm(" See the null__area raster with this legend including their combinations.")
            for m,i in zip(self.strcolumns.values(), pow2): gm("%s = %s" %(i,m))
            gm('null__area   area [km2]')
            for m in sorted(null.keys()): gm("%s        %s" %(m,float(null[m])*10**-6))
            gm("How should they appear in the .str file?")
            gm("Set them with r.null")
            grass.fatal('Exiting!')
        g_run('r.mask', flags='r', quiet=True)

        return

    def _maskOrFalse(self, name):
        '''Check if name option is given, if yes make a mask,
        if not return False.'''
        # if given, then it must be a raster, if not just make blank
        argv = getattr(self, name, None)
        if argv:
            # make mask for DCELL and FCELL
            if not grass.raster_info(argv)['datatype'] == 'CELL':
                outname = '%s__mask' % name
                grass.mapcalc(exp=outname+'=if(isnull(%s), 0, 1)' % argv)
                self.floatmaps[outname] = self.__dict__[name]
            else:
                # CELL is just passed on
                outname = self.__dict__[name]
        else:
            return False

        return outname

    def mkContours(self):
        """Make a contour map of the given DEM, in the given catchment"""

        # check what breaks to take, if int find nice breaks, if list take as breaks
        if type(self.contours) is int:
            interval = self.contours
            # get stats of DEM for nice breaks
            stats=grass.parse_command('r.univar',map=self.elevation,flags='g')
            # make nice breaks
            minelev = int(float(stats['min'])/interval+1)*interval
            maxelev = int(float(stats['max'])/interval)*interval
            breaks  = list(range(minelev, maxelev+1, interval))
        else:
            breaks = self.contours
        if len(breaks)<2:
                grass.fatal('Need at least 2 contour breaks: %s \nRange of elevation: %s - %s' %(breaks,stats['min'],stats['max']))
        grass.message(('Contour breaks:', str(breaks)))
        # form mapcalc expression and execute
        exp = self.contourrast+'= if(%s<%s,1)' %(self.elevation, breaks[0]) # for the first, smaller than
        for b in breaks: exp+='+ if(%s>=%s,1)' %(self.elevation,b) # for all greater eq than
        grass.mapcalc(exp, overwrite=True)

        grass.message(('Created contourmap: %s' %self.contourrast))

        return

    def mkHydrotopes(self):
        """Calculates hydrotope map for the maps in mapslist in the catchment
        and writes a structure file with category values of the maps in the same
        order plus an area and cell column
        """
        g_run('r.mask',raster=self.subbasin_id,overwrite=True)

        # maps list
        ml=','.join(self.strcolumns.values())
        # take cross product of maps
        g_run('r.cross', overwrite=True, flags='z', input=ml,
              output='hydrotopes__rcross')
        # start hydrotope count at 1 instead of 0
        grass.mapcalc('$output=$input+1', output=self.hydrotopes,
                      input='hydrotopes__rcross')
        self.strcolumns['hydrotope_id'] = self.hydrotopes
        self.strcolumns.move_to_end('hydrotope_id', last=False)

        # read basic structure file info from hydrotope map
        struct = readinStr(self.strcolumns.values())
        self.strcolumns['area'] = ''

        # replace all float maps with float values
        for intmap,floatmap in self.floatmaps.items():
            # get mean hydrotope array with columns: hydrotope cats, mean
            catval = hydrotopeQ(floatmap, 'hydrotopes__rcross')
            # exchange column in strct with mean hydrotope value
            struct[intmap][catval['cat']] = catval['mean']

        # write structure file
        self.writeStr(struct)

        # report hydrotope count and check max number of hydrotopes in subbasins
        nmax = np.max(np.bincount(struct[self.subbasin_id].astype(int)))
        grass.message('''%s hydrotopes created, %5.2f per subbasin on average, max.
number of hydrotopes per subbasin %i
                      ''' %(len(struct),len(struct)/len(np.unique(struct[self.subbasin_id])),nmax))
        g_run('r.colors', map=self.hydrotopes, color="random", quiet=True)
        return struct

    def writeStr_no_mswim_lib(self, array):
        """Deprecated in favour of using abstracted mswim.io.write_csv"""
        # formate for structure file
        colwidth = [max(i, 9) for i in map(len, self.strcolumns)]
        heacolfmt = ['%%-%ss' % i for i in colwidth] 
        datcolfmt = [('%%%s.3f' if c in self.floatmaps else '%%%si') % i
                     for i, c in zip(colwidth, self.strcolumns.values())]
        # write out structure file
        with open(self.output, 'w') as strf:
            # header
            strf.write((', '.join(heacolfmt) + '\n') % tuple(self.strcolumns.keys()))
            # data
            np.savetxt(strf, array, fmt=', '.join(datcolfmt))
        grass.message(('Wrote structure file %s' % self.output))
        return

    def writeStr(self, array):
        """Write the array into the structure file path given in strname"""
        floatcols = [k for k, v in self.strcolumns.items() if v in self.floatmaps]
        mswim.io.write_csv(self.output, array, self.strcolumns.keys(),
                           float_columns=floatcols)
        grass.message(('Wrote hydrotope file %s' % self.output))
        return

def readinStr(strcolumns):
    """ Return subbasinID, landuseID, soilID, ..., area, cellcount from
    hydrotope map as an numpy array"""
    # get all values with area and cell, [:-1] as last line hast line break
    tbl = grass.read_command('r.stats', input=','.join(strcolumns), flags='an').split('\n')[:-1]
    tbl = [tuple(l.split()) for l in tbl]
    tbl = np.array(tbl, dtype=list(zip(list(strcolumns)+['area'], [float]*(len(strcolumns) + 1))))

    return tbl

def hydrotopeQ(cover,hydrotopemap):
    """Get mean values of the cover map for the hydrotopes"""
    grass.message(('Get mean hydrotope values for %s' %cover))

    tbl = grass.read_command('r.univar', map=cover, zones=hydrotopemap,
                           flags='gt').split('\n')[:-1] #:-1 as last line hast line break]
    tbl = [tuple(l.split('|')) for l in tbl]
    tbl = np.array(tbl[1:], dtype=list(zip(tbl[0],['S250']*len(tbl[0]))))
    tbl = np.array(list(zip(tbl['zone'],tbl['mean'])), dtype=[('cat',np.int64),('mean',np.float64)])
    return tbl[np.isfinite(tbl['mean'])]


if __name__=='__main__':
    st = dt.datetime.now()
    # print version/date before doing anything else
    mswim.utils.print_version(__file__) if '-v' in sys.argv else None
    # get options and flags
    o, f = grass.parser()
    fmt = lambda d: '\n'.join(['%s: %s' % (k, v) for k, v in d.items()])+'\n'
    grass.message('GIS Environment:\n'+fmt(grass.gisenv()))
    grass.message('Parameters:\n'+fmt(o)+fmt(f))

    # send all to main
    keywords = o; keywords.update(f)
    main=main(**keywords)

    # calculate hydrotopes
    grass.message(('Calculating hydrotopes with: %s...' %main.strcolumns))
    main.mkHydrotopes()

    # clean
    g_run('r.mask', flags='r', quiet=True)
    if not main.k:
        grass.run_command('g.remove',type='raster,vector', pattern='*__*',flags='fb',quiet=True)

    # report time it took
    delta = dt.datetime.now()-st
    grass.message('Execution took %s hh:mm:ss' %delta)
