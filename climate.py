#!/usr/bin/env python
#
############################################################################
#
# MODULE:      m.swim.climate
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
#% description: Soil and Water Integrated Model (SWIM) climate input preprocessor
#% keywords: hydrological modelling, SWIM, climate interpolation
#%End
 
#%Option
#% guisection: Subbasins
#% key: subbasins
#% type: string
#% multiple: no
#% required: yes
#% key_desc: name
#% label: Subbasins vector
#% description: Subbasins for which climate data will be interpolated
#% gisprompt: old,vector,vector
#%end

#%Option
#% guisection: Subbasins
#% key: elevation
#% type: string
#% required: yes
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

import grass.script as grass
import numpy as np
import sys, os
import datetime as dt
grun = grass.run_command
gread= grass.read_command
gm   = grass.message

class main:
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
        ## TODO        
        return
        
    def writeClimStationfile(self):
        # get coordinates and table
        coor, a = vectCoordsTbl(self.climstations)
        # add filepath to name/2nd column
        fnames = []
        for i,n in enumerate(a[self.fnames]):
            fnames += [os.path.join(self.datadir,n+self.ext)]
            if not os.path.exists(fnames[-1]): grass.warning('%s doesnt exisist!' %fnames[-1])
        # get needed cols
        cols = np.array(zip(coor['cat'],
                            fnames,coor['x'],coor['y'],
                            a[self.stationelevation]),
                        dtype=[('cat',int),('fnames','S1000'),('x',float),('y',float),('z',float)])
        # write file
        sfpath = os.path.join(self.outdir,'stationfile.dat')
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
        return
        
def vectCoordsTbl(vect):
    '''Return the v.report table with coordinates as an array with
    column names in the entries for a point vector'''
    # get table, returned columns:
    # table columns, x, y, z
    t=grass.read_command('v.report', map=vect, option='coor').split()
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
    
if __name__=='__main__':
    st = dt.datetime.now()
    # get options and flags
    o, f = grass.parser()
    
    # send all to main
    keywords = o; keywords.update(f)
    main=main(**keywords)

    grass.message(('GIS Environment:',grass.gisenv()))
    grass.message(('Parameters:',o,f))
    
    # process
    main.writeClimStationfile()
    main.writeVirtualStationfile()

    # report time it took
    delta = dt.datetime.now()-st
    grass.message('Execution took %s hh:mm:ss' %delta)
