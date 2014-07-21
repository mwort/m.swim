#!/usr/bin/env python
#
############################################################################
#
# MODULE:      m.swim.hydrotopes
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
#% required: no
#% multiple: no
#% key_desc: name
#% description: Elevation raster
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
 

import os,sys
import grass.script as grass
g_run=grass.run_command
gm   =grass.message
import numpy as np
import datetime as dt

class main:
    def __init__(self,**optionsandflags):
        '''Process all arguments and prepare processing'''
        # add all options and flags as attributes (only nonempty ones)
        self.options = {}
        for o in optionsandflags:
            if optionsandflags[o]!='': self.options[o] = optionsandflags[o]
        self.__dict__.update(self.options)
        
        # default columns
        self.strcolumns = [self.subbasins,self.landuse, self.soil]

        # check what kind of contours
        if self.c:
            if 'elevation' not in self.options:
                grass.fatal(('DEM not specified for contours.'))
            # convert contours to ints
            if 'contours' in self.options:
                try:
                    self.contours = map(int,self.contours.split(','))
                    if len(self.contours)==1: self.contours = self.contours[0]
                except:
                    grass.fatal(('''Contours should be either an interval [integer]
                    or a list of breaks separated by commas.'''))
            # add to columns
            self.strcolumns += [self.contourrast]
        
        # other raster to include?
        if 'more' in self.options:
            self.more = self.more.split(',')
            self.strcolumns += self.more
            
        # check that all input maps have no NULLs in catchment over subbasins area
        g_run('r.mask',raster=self.subbasins,overwrite=True)
        # maps that neednt be null
        maps = [self.landuse,self.soil]
        # optional
        if self.c: maps += [self.elevation]
        if self.c and 'contours' not in self.options: maps += [self.contourrast]
        if 'more' in self.options: maps += self.more
        # check
        gm('Checking for NULLs in input maps...')
        pow2 = 2**np.arange(len(maps))
        ifnull = ' + '.join(['isnull(%s)*%s' %(m,i) for m,i in zip(maps,pow2)])
        grass.mapcalc('null__areas=%s' %ifnull, overwrite=True)
        g_run('r.null',map='null__areas',setnull=0,quiet=True)
        null = grass.parse_command('r.stats',input='null__areas',flags='aN',separator='=')
        
        if len(null)>0:
            gm("In any of the input maps are NULL/no data values over the subbasins area.")
            gm(" See the null__area raster with this legend including their combinations.")
            for m,i in zip(maps,pow2): gm("%s = %s" %(i,m))
            gm('null__area   area [km2]')
            for m in sorted(null.keys()): gm("%s        %s" %(m,float(null[m])*10**-6))
            gm("How should they appear in the .str file?")
            gm("Set them with r.null")
            grass.fatal('Exiting!')
        g_run('r.mask', flags='r', quiet=True)

# LONGER CHECK, possibly wihtout resolution errors
#        for m in maps:
#            areas = grass.parse_command('r.stats',input=','.join([self.subbasins,m]),
#                                        flags='aN',separator='=')
#            # only look for areas that have null * in m and make them floats
#            nullarea = []
#            for sb in areas:
#                na = areas[sb].split('=')
#                if na[0]=='*':
#                    nullarea += [(int(sb),float(na[1])*10**-6)]
#            if len(nullarea)>0:
#                #for i in nullarea: print i
#                grass.fatal('''%s has NULL/no data values over a total area of %s sq km.
#                They are located in the above subbasins. How should they appear in the .str file?
#                Set them with r.null''' %(m,np.array(nullarea)[:,1].sum()))
        
        return
        
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
        g_run('r.cross', overwrite=True, flags='', # also where NULL, MASK takes care of the valid catchment part
              input=ml, output=self.hydrotopes)
        g_run('r.null', map=self.hydrotopes, setnull=0)
        
        # read basic structure file info from hydrotope map
        struct = readinStr(self.hydrotopes,self.strcolumns)
        if self.c:
            # get mean hydrotope elevations array with columns: hydrotope cats, mean elevation
            mhe = hydrotopeQ(self.elevation, self.hydrotopes)
            # exchange contour column in strct with mean hydrotope elevation
            struct[self.contourrast] = mhe[:,1]
        # write structure file
        writeStr(struct,self.strfilepath)
        
        # report hydrotope count and check max number of hydrotopes in subbasins
        nmax = np.max(np.bincount(struct[self.subbasins].astype(int)))
        grass.message('''%s hydrotopes created, %5.2f per subbasin on average, max.
number of hydrotopes per subbasin %i
                      ''' %(len(struct),len(struct)/len(np.unique(struct[self.subbasins])),nmax))

        return struct
        
    def addcovertostrfile(self,cover):
        '''Function to add additional columns to saved .str file'''
        # get old str file    
        strf = np.loadtxt(self.strfilepath)[:-1,:]
        # get cover value for hydrotopes
        carray= hydrotopeQ(cover,self.hydrotopes)
        # assemble new str array
        newstrn=np.column_stack((strf[:,:-2],carray[:,1],strf[:,-2:]))
        # write new strfile
        newname=self.strfilepath+'+'
        writeStr(newstrn,newname)    
        return
        
def readinStr(hydrotopemap,strcolumns):
    """ Return subbasinID, landuseID, soilID, ..., area, cellcount from
    hydrotope map as an numpy array"""
    
    # get mean values of each raster over each hydrotope
    struct = []
    for s in strcolumns:
        tbl=grass.read_command('r.univar', map=s,zones=hydrotopemap,
                               flags='gt').split('\n')[:-1] #:-1 as last line hast line break]
        tbl      = [tuple(l.split('|')) for l in tbl]
        array = np.array(tbl[1:],dtype=zip(tbl[0],['S250']*len(tbl[0])))
        struct += [array['mean']]
        gm('Read hydrotope values for %s' %s)
    # add number of cells and area
    ncells = array['non_null_cells'].astype(int)
    reg = grass.region()
    area   = (ncells*reg['nsres']*reg['ewres']).astype(int)
    struct += [area,ncells]
    # make nice record array
    dtype  = zip(strcolumns+['area','ncells'],[int]*(len(strcolumns)+2))
    struct = np.array(zip(*struct),dtype=dtype)

    return struct

def writeStr(array,strname):
    """Write the array into the structure file path given in strname"""
    # formate for structure file
    ncolumns=len(array[0])
    formats='%14i'*ncolumns
    # write out structure file
    strf=file(strname,'w')
    np.savetxt(strf, array, fmt=formats)
    strf.write(formats %((0,)*ncolumns)) #not sure if SWIM needs that
    strf.close()
    grass.message(('Wrote structure file %s' %strname))
    return strname
    
def hydrotopeQ(cover,hydrotopemap):
    """Get mean values of the cover map for the hydrotopes"""
    grass.message(('Get mean hydrotope values for %s' %cover))
    # convert base map to integers
    coverint=cover.split('@')[0]+'__int'
    grass.mapcalc(coverint+'=int('+cover+')', overwrite=True)
    # run r.statistics
    outname=hydrotopemap.split('@')[0]+'__mean__'+coverint
    g_run('r.statistics', overwrite=True,
          base=hydrotopemap, cover=coverint, method='average', output=outname)
    
    # get category values and labels where statistics are stored
    vals = grass.read_command('r.stats',flags='nl', input=outname).split()
    vals = np.array(np.reshape(vals,(len(vals)/2,2)),dtype=float)
    
    return vals
    

    
if __name__=='__main__':
    st = dt.datetime.now()
    # get options and flags
    o, f = grass.parser()
    grass.message(('GIS Environment:',grass.gisenv()))
    grass.message(('Parameters:',o,f))
    
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
        grass.run_command('g.mremove',rast='*__*',vect='*__*',flags='f',quiet=True)
    
    # report time it took
    delta = dt.datetime.now()-st
    grass.message('Execution took %s hh:mm:ss' %delta)
    

