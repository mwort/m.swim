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
#% key_desc: name
#% label: Subbasins vector
#% description: Subbasins for which climate data will be interpolated
#% gisprompt: old,vector,vector
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
        # no extension
        if 'ext' not in self.options: self.ext=''
## TODO        
        return
        
    def writeClimStationfile(self):
        filepath = self.datadir
        stationvect = self.climstations
        # get columns in correct order into array
        cols = [self.fnames,self.stationelevation]
        a   = grass.vector_db_select(stationvect,columns=','.join(cols))
        a   = np.array(a['values'].values())
        names = map(str,a[:,0])
        z   = map(float,a[:,1])
        # get coordinates
        coor = vectCoords(stationvect)
        # make nice array
        a   = np.array(zip(coor['cat'],names,coor['x'],coor['y'],z),
                       dtype=[('ids',int),('names','S5120'),('x',float), ('y',float), ('z',int)])
        # add filepath to name/2nd column
        for i,n in enumerate(a['names']):
            a['names'][i] = str(os.path.join(self.outdir,n+self.ext))
    
        # write file
        sfpath = os.path.join(self.outdir,'stationfile.dat')
        # put header
        head='ID  FILE                       POINT_X     POINT_Y            ELEV'
        np.savetxt(sfpath, a, fmt='%6i  %-'+str(len(filepath)+32)+'s%12.1f%12.1f%10.1f',
                   header=head)
        gm( 'Saved climate stations to %s' %sfpath)
        
        return
    
    def writeVirtualStationfile(self):
        # get coordinates and cats
        coor = vectCoords(self.subbasins)
        # add z
        z = np.array(grass.vector_db_select(self.subbasins,
                     columns=self.elevation)['values'].values(),dtype=float)
        # table to write out
        tbl = np.column_stack((coor['cat'],coor['x'],coor['y'],z))
        
        # save to file
        sfpath = os.path.join(self.outdir,'virtualstationfile.dat')
        # put header
        head='ID  POINT_X     POINT_Y            ELEV'
        np.savetxt(sfpath, tbl, fmt='%6i %12.1f%12.1f%10.1f',
                   header=head)
        gm( 'Saved virtual stations to %s' %sfpath)
        return
        
def vectCoords(vect):
    '''Return the v.report table with coordinates as an array with
    column names in the entries for a point vector'''
    t=grass.read_command('v.report', map=vect, option='coor').split()
    t=[tuple(l.split('|')) for l in t][1:]
    t=[(l[0],l[-3],l[-2]) for l in t]
    t=np.array(t,dtype=[('cat',int),('x',float),('y',float)])
    
    return t
    
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
