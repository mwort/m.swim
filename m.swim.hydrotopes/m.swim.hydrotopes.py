#!/usr/bin/env python
#
############################################################################
#
# MODULE:      m.swim.hydrotopes v1.0
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
#% description: Soil and Water Integrated Model (SWIM) preprocessor: hydrotopes
#%End
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
#% guisection: Required
#% key: strfilepath
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
#%Flag
#% guisection: Contours
#% key: c
#% label: include contours, either create if contours given or use contourrast
#%end
#%Option
#% guisection: Contours
#% key: contours
#% type: string
#% required: no
#% multiple: no
#% key_desc: ints or list of ints
#% description: Elevation contours to include in hydrotopes [optional, if set, DEM in Subbasins also needs to be set]
#%end
#%Option
#% guisection: Contours
#% key: contourrast
#% type: string
#% required: no
#% multiple: no
#% key_desc: existing rast of elevation contours or to be created
#% answer: contours
#% description: Elevation contours to include in hydrotopes [optional, if set, DEM in Subbasins also needs to be set]
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
#% key: management
#% type: string
#% required: no
#% multiple: no
#% key_desc: name
#% description: Raster name (if not given, column is filled with default)
#% gisprompt: old,cell,raster
#%end
#%Option
#% guisection: Optional
#% key: wetland
#% type: string
#% required: no
#% multiple: no
#% key_desc: name
#% description: Raster name (if not given, column is filled with default)
#% gisprompt: old,cell,raster
#%end
#%Option
#% guisection: Optional
#% key: glaciers
#% type: string
#% required: no
#% multiple: no
#% key_desc: name
#% description: Raster name (if not given, column is filled with default)
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


import os,sys,collections
import grass.script as grass
g_run=grass.run_command
gm   =grass.message
import numpy as np
import datetime as dt

class main:

    DEFAULTS = {'management': 1,
                'wetland':    0,
                'glaciers':   0,
                'contourrast':   0}

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
        self.strcolumns = [self.subbasins,self.landuse, self.soil]

        # managment and wetland
        for r in ['management','wetland']:
            self.__dict__[r] = self._maskOrBlank(r)
            self.strcolumns += [self.__dict__[r]]


        # check contours and create raster if needed
        if 'contours' in self.options:
            try:
                self.contours = map(int,self.contours.split(','))
                if len(self.contours)==1: self.contours = self.contours[0]
            except:
                grass.fatal(('''Contours should be either an interval [integer]
                or a list of breaks separated by commas.'''))
            # create contourrast
            self.mkContours()
        else:
            self.contourrast = self._maskOrBlank('contourrast')
        # add to columns
        self.strcolumns += [self.contourrast]
        self.floatmaps[self.contourrast] = self.elevation

        # glaciers
        self.glaciers = self._maskOrBlank('glaciers')
        self.strcolumns += [self.glaciers]

        # other raster to include?
        if 'more' in self.options:
            self.more = self.more.split(',')
            self.strcolumns += self.more

        # check if all input maps are int/CELL maps
        gm('Check if all input are integer raster...')
        for m in self.strcolumns:
            if grass.raster_info(m)['datatype']!='CELL':
                grass.fatal('%s is not an integer/CELL raster, convert using int() in r.mapcalc' %m)

        # check that all input maps have no NULLs in catchment over subbasins area
        g_run('r.mask',raster=self.subbasins,overwrite=True)

        gm('Checking for NULLs in input maps...')
        pow2 = 2**np.arange(len(self.strcolumns))
        ifnull = ' + '.join(['isnull(%s)*%s' %(m,i) for m,i in zip(self.strcolumns,pow2)])
        grass.mapcalc('null__areas=%s' %ifnull, overwrite=True)
        g_run('r.null',map='null__areas',setnull=0,quiet=True)
        null = grass.parse_command('r.stats',input='null__areas',flags='aN',separator='=')

        if len(null)>0:
            gm("In any of the input maps are NULL/no data values over the subbasins area.")
            gm(" See the null__area raster with this legend including their combinations.")
            for m,i in zip(self.strcolumns,pow2): gm("%s = %s" %(i,m))
            gm('null__area   area [km2]')
            for m in sorted(null.keys()): gm("%s        %s" %(m,float(null[m])*10**-6))
            gm("How should they appear in the .str file?")
            gm("Set them with r.null")
            grass.fatal('Exiting!')
        g_run('r.mask', flags='r', quiet=True)

        return

    def _maskOrBlank(self, name):
        '''Check if name option is given, if yes make a mask,
        if not make a blank map out of the default value.'''
        # if given, then it must be a raster, if not just make blank
        if (hasattr(self, name) and
            tuple(name.split('@')) in grass.list_pairs('rast')):
            # make mask for DCELL and FCELL
            if not grass.raster_info(self.__dict__[name])['datatype']=='CELL':
                outname = '%s__mask'%name
                grass.mapcalc(exp=outname+'=isnull(%s)'%self.__dict__[name])
                self.floatmaps[outname] = self.__dict__[name]
            else:
                # CELL is just passed on
                outname = self.__dict__[name]
        else:
            # not given, make default blank map
            outname = '%s__default'%name
            grass.mapcalc(exp=outname+'=%s'%self.DEFAULTS[name])

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
            breaks  = range(minelev,maxelev+1,interval)
        else:
            breaks = self.contours
        if len(breaks)<2:
                grass.fatal('Need at least 2 contour breaks: %s \nRange of elevation: %s - %s' %(breaks,stats['min'],stats['max']))
        grass.message(('Contour breaks:',str(breaks)))
        # form mapcalc expression and execute
        exp = self.contourrast+'= if(%s<%s,1)' %(self.elevation, breaks[0]) # for the first, smaller than
        for b in breaks: exp+='+ if(%s>=%s,1)' %(self.elevation,b) # for all greater eq than
        grass.mapcalc(exp, overwrite=True)

        grass.message(('Calculated contourmap: %s' %self.contourrast))

        return

    def mkHydrotopes(self):
        """Calculates hydrotope map for the maps in mapslist in the catchment
        and writes a structure file with category values of the maps in the same
        order plus an area and cell column
        """
        g_run('r.mask',raster=self.subbasins,overwrite=True)

        # maps list
        ml=','.join(self.strcolumns)
        # take cross product of maps
        g_run('r.cross', overwrite=True, flags='z',
              input=ml, output=self.hydrotopes)

        # read basic structure file info from hydrotope map
        struct = readinStr(self.strcolumns)

        # replace all float maps with float values
        for intmap,floatmap in self.floatmaps.items():
            # get mean hydrotope array with columns: hydrotope cats, mean
            catval = hydrotopeQ(floatmap, self.hydrotopes)
            # exchange column in strct with mean hydrotope value
            struct[intmap][catval['cat']] = catval['mean']

        # write structure file
        self.writeStr(struct)

        # report hydrotope count and check max number of hydrotopes in subbasins
        nmax = np.max(np.bincount(struct[self.subbasins].astype(int)))
        grass.message('''%s hydrotopes created, %5.2f per subbasin on average, max.
number of hydrotopes per subbasin %i
                      ''' %(len(struct),len(struct)/len(np.unique(struct[self.subbasins])),nmax))

        return struct

    def writeStr(self, array):
        """Write the array into the structure file path given in strname"""
        # formate for structure file
        ncolumns=len(array[0])
        colwidth=14
        heacolfmt = ('%%-%ss '%colwidth)*ncolumns
        datcolfmt = ('%%%si '%colwidth)*ncolumns
        colnames= tuple([c.split('@')[0] for c in self.strcolumns]) + ('area','cells')
        # write out structure file
        strname = self.strfilepath
        strf=file(strname,'w')
        # header
        strf.write( (heacolfmt + '\n')%colnames )
        # data
        np.savetxt(strf, array, fmt=datcolfmt)
        # 0s as the last line
        strf.write( datcolfmt%((0,)*ncolumns))
        strf.close()
        grass.message(('Wrote structure file %s' %strname))
        return

def readinStr(strcolumns):
    """ Return subbasinID, landuseID, soilID, ..., area, cellcount from
    hydrotope map as an numpy array"""
    # get all values with area and cell, [:-1] as last line hast line break
    tbl = grass.read_command('r.stats',input=','.join(strcolumns),flags='acn').split('\n')[:-1]
    tbl = [tuple(l.split()) for l in tbl]
    tbl = np.array(tbl, dtype=zip(strcolumns+['area','cells'],
                                  [np.int64]*len(strcolumns) + [np.float,np.int64]))

    return tbl

def hydrotopeQ(cover,hydrotopemap):
    """Get mean values of the cover map for the hydrotopes"""
    grass.message(('Get mean hydrotope values for %s' %cover))

    tbl = grass.read_command('r.univar', map=cover, zones=hydrotopemap,
                           flags='gt').split('\n')[:-1] #:-1 as last line hast line break]
    tbl = [tuple(l.split('|')) for l in tbl]
    tbl = np.array(tbl[1:],dtype=zip(tbl[0],['S250']*len(tbl[0])))

    return np.array(zip(tbl['zone'],tbl['mean']),dtype=[('cat',np.int64),('mean',np.float64)])



if __name__=='__main__':
    st = dt.datetime.now()
    # get options and flags
    o, f = grass.parser()
    fmt = lambda d: '\n'.join(['%s: %s' % (k, v) for k, v in d.items()])+'\n'
    grass.message('GIS Environment:\n'+fmt(grass.gisenv()))
    grass.message('Parameters:\n'+fmt(o)+fmt(f))

    # send all to main
    keywords = o; keywords.update(f)
    main=main(**keywords)

    ### EXECUTION OPTIONS
    # include or make contours
    if main.c and 'contours' in main.options:
        main.mkContours()

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
