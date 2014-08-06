#!/usr/bin/env python
#
############################################################################
#
# MODULE:      m.swim.subbasins
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
#% description: Soil and Water Integrated Model (SWIM) subbasin preprocessor
#% keywords: hydrological modelling, SWIM, subbasins
#%End

#%Option
#% guisection: Input
#% key: elevation
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% description: Elevation raster (hole-filled)
#% gisprompt: old,cell,raster
#%end

#%Option
#% guisection: Input
#% key: stations
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% label: Station point vector
#% description: Will be snapped to the nearest stream
#% gisprompt: old,vector,vector
#%end

#%Option
#% guisection: Subbasin design
#% key: upthresh
#% type: double
#% required: no
#% multiple: no
#% label: Upper subbasin threshold (in sq km,mostly underestimated)
#% description: ignored if upthreshcolumn is given, mostly underestimated
#% answer: 100
#%end

#%Option
#% guisection: Subbasin design
#% key: lothresh
#% type: double
#% required: no
#% multiple: no
#% key_desc: km2
#% description: Lower threshold of subbasin size in sq km (default 5% of upper)
#%end

#%Option
#% guisection: Subbasin design
#% key: upthreshcolumn
#% type: string
#% required: no
#% multiple: no
#% key_desc: name
#% description: Column with upper subbasin threshold in stations vector
#%end

#%Flag
#% guisection: Topography
#% key: d
#% label: don't process DEM (accumulation, drainage, streams must exist)
#%end
#%Option
#% guisection: Topography
#% key: depression
#% type: string
#% required: no
#% multiple: no
#% key_desc: name
#% description: Raster map of known real depressions in or around catchment (only if not -d)
#% gisprompt: old,cell,raster
#%end
#%Option
#% guisection: Topography
#% key: accumulation
#% type: string
#% required: no
#% multiple: no
#% key_desc: string
#% description: Name of accumlation map to be created (or existing if -d)
#% answer: accumulation
#% gisprompt: new,cell,raster
#%end
#%Option
#% guisection: Topography
#% key: drainage
#% type: string
#% required: no
#% multiple: no
#% key_desc: string
#% description: Name of drainge map to be created (or existing if -d)
#% answer: drainage
#% gisprompt: new,cell,raster
#%end

#%Option
#% guisection: Topography
#% key: streams
#% type: string
#% required: no
#% multiple: no
#% key_desc: string
#% description: Name of streams vector to be created (or existing if -d)
#% answer: streams
#% gisprompt: new,vector,vector
#%end

#%Option
#% guisection: Topography
#% key: streamthresh
#% type: double
#% required: no
#% multiple: no
#% key_desc: name
#% label: Drainage area of smallest stream in km2 (influences station snapping, default: 10% of region)
#% description: Stations will be snapped to these streams, ie. should not be smaller than the smallest catchment
#%end

#%Option
#% guisection: Output
#% key: subbasins
#% type: string
#% required: no
#% multiple: no
#% key_desc: string
#% description: Name of resulting subbasin vector and raster map
#% answer: subbasins
#% gisprompt: new,cell,raster
#%end
#%Option
#% guisection: Output
#% key: catchments
#% type: string
#% required: no
#% multiple: no
#% key_desc: string
#% description: Name of resulting vector and raster of all stations catchments
#% answer: catchments
#% gisprompt: new,cell,raster
#%end
#%Option
#% guisection: Output
#% key: catchmentprefix
#% type: string
#% required: no
#% multiple: no
#% key_desc: string
#% description: Prefix of individual catchment vector for each station
#% answer: catchment_
#% gisprompt: new,cell,raster
#%end

#%Option
#% guisection: Output
#% key: slopesteepness
#% type: string
#% required: no
#% multiple: no
#% key_desc: string
#% description: USLE slopesteepness from r.watershed, useful for m.swim.substats
#% answer: slopesteepness
#% gisprompt: new,cell,raster
#%end

#%Option
#% guisection: Output
#% key: slopelength
#% type: string
#% required: no
#% multiple: no
#% key_desc: string
#% description: USLE slope length from r.watershed, useful for m.swim.substats
#% answer: slopelength
#% gisprompt: new,cell,raster
#%end

#%Option
#% guisection: Optional
#% key: streamcarve
#% type: string
#% required: no
#% multiple: no
#% key_desc: vector
#% description: Existing river network to be carved into the elevation raster
#% gisprompt: old,vector,vector
#%end

#%Option
#% guisection: Optional
#% key: predefined
#% type: string
#% required: no
#% multiple: no
#% key_desc: string
#% description: Raster of predefined units to include in the subbasin map
#% gisprompt: old,cell,raster
#%end

#%Option
#% guisection: Optional
#% key: rwatershedflags
#% type: string
#% required: no
#% multiple: no
#% key_desc: string
#% description: Flags parsed to r.watershed, check r.watershed --help
#% answer: s
#%end

#%Flag
#% guisection: Optional
#% key: s
#% label: Just print statistics of subbasins, --o must be set
#%end

#%Flag
#% guisection: Optional
#% key: k
#% label: Keep intermediat files (include __ in names)
#%end

import os, sys
import grass.script as grass
g_run=grass.run_command
gm   = grass.message
import numpy as np
import datetime as dt


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
        self.region['kmtocell'] = lambda km: int(round(np.mean(km)*10**6/(self.region['ewres']*self.region['nsres'])))
        self.region['celltokm'] = lambda c: c*(self.region['ewres']*self.region['nsres'])*1e-6
        
        # check if DEM to processed or if all inputs set
        if not ('accumulation' in self.options and 'drainage' in self.options and 'streams' in self.options):
            grass.fatal('Either of these not set: accumulation, drainage, streams.')
        
        # what to do with upthresh
        if 'upthreshcolumn' in self.options:
            gm('Will look for upper thresholds in the %s column.' %self.upthreshcolumn)
            # get thresholds from column in station vect
            try:
                threshs = grass.vector_db_select(self.stations,columns=self.upthreshcolumn)['values']
                self.upthresh = np.array(threshs.values(),dtype=float)[:,0]
            except:
                grass.fatal('Cant read the upper threshold from the column %s' %self.upthreshcolumn)
        
        # lothresh default
        if 'lothresh' not in self.options:
            self.lothresh = np.mean(self.upthresh)*0.05
            
        # streamthresh
        if 'streamthresh' in self.options:
            # convert to cells
            self.streamthresh = self.region['kmtocell'](self.streamthresh)
            # check if reasonable
            fract = float(self.streamthresh)/self.region['cells']
            if fract > 0.5 or fract < 0.01: grass.warning('streamthresh is %s percent of the region size!' %(fract*100))
        else:
            self.streamthresh = int(self.region['cells']*0.02)
         
        # if no r.watershed flags given
        if 'rwatershedflags' not in self.options: self.rwatershedflags='s'
        
        # check input for stats print
        if self.s:
            for o in ['streams','stations','catchmentprefix']:
                if o not in self.options: grass.fatal('%s needs to be set!')
            # get all catchments
            self.allcatchments=grass.mlist_strings('rast',self.catchmentprefix+'*')
            gm( 'Found these catchments %s' %self.allcatchments)
            # calculate station topology
            self.station_coor = snappoints(self.stations, self.streams)
            self.stationtopology = self.getTopology()
        
        # initialise subbasinsdone
        self.subbasinsdone={}

        return
        
    def procDEM(self):
        """Make drainage, accumulation and stream raster"""

        ######### decide on input arguments #######
        
        # from km2 to cells
        if type(self.upthresh) in [int,float]:
            uthresh=self.upthresh
        else: # take most common one in upthresh column
            uthresh = max(set(self.upthresh), key=list(self.upthresh).count)
            
        thresh = self.region['kmtocell'](uthresh)

        kwargs = {'elevation'   : self.elevation,
                  'threshold'   : thresh,
                  # Output
                  'accumulation': 'accum__float',
                  'drainage'    : self.drainage,
                  #'stream'      : self.streams, # created later with accum raster
                  'basin'       : 'standard__subbasins',
                  'slope_steepness': self.slopesteepness,
                  'length_slope': self.slopelength,
                  'flags'       : self.rwatershedflags}
        # check if depressions
        if 'depression' in self.options: kwargs['depression']=self.depression
        
        # carve streams
        if 'streamcarve' in self.options: kwargs['elevation'] = self.carvedelevation
            
        # run r.watershed command and show progress
        #environ=os.environ.copy()
        #environ['GRASS_MESSAGE_FORMAT'] = 'gui'
        grass.message(kwargs)
        g_run('r.watershed',overwrite=True,**kwargs) # the other keyword arguments
        
        # save subbasins in dictionary
        self.subbasinsdone[thresh] = 'standard__subbasins'
        
        # postprocess accumulation map
        grass.mapcalc("%s=int(if(accum__float <= 0,null(),accum__float))" %self.accumulation,
                      overwrite=True)          
        # make river network to vector
        grass.message('Making vector river network...')
        grass.mapcalc("{0}__thick = if({1} >= {2}, {1}, null())".format(self.streams,self.accumulation,self.streamthresh))
        g_run('r.thin',input=self.streams+'__thick',output=self.streams,quiet=True)
        g_run('r.to.vect', flags='s', input=self.streams, output=self.streams,type='line')
        
        return
    
    def carveStreams(self):
        '''Carve vector streams into the DEM, i.e. setting those cells =0'''
        # stream vector to raster cells
        streamrast = self.streamcarve.split('@')[0]+'__'
        g_run('v.to.rast',input=self.streamcarve,output=streamrast, type='line',
              use='val',val=1,quiet=True,overwrite=True)
        # carve
        self.carvedelevation = '%s__carved' %self.elevation.split('@')[0]
        grass.mapcalc("%s=if(isnull(%s),%s,0)" %(self.carvedelevation,streamrast,self.elevation),
                      overwrite=True)
        gm('Carved %s into the elevation %s' %(self.streamcarve,self.elevation))
        return
        
        
    def makeCatchments(self):
        '''Make catchment raster and if catchmentprefix is set also vectors
        for all stations'''
        self.allcatchments = []
        ####### snap stations to network #######
        gm('Snap stations to streams...')
        self.station_coor = snappoints(self.stations, self.streams)
            
        ######### WATERSHEDS #################
        gm('Creating station catchments...')
        # create watersheds for stations
        for i,s in enumerate(self.station_coor):
            if 'catchmentprefix' in self.options:
                name=self.catchmentprefix+'%s' %(i+1)
            else:
                name='watersheds__st%s' %(i+1)
            grass.message(('Calculating watershed for station %s' %(i+1)))
            g_run('r.water.outlet', input=self.drainage, overwrite=True,
                  output=name+'__all1', coordinates='%s,%s' %tuple(s))
            # give watershed number and put 0 to null()
            grass.mapcalc(name+' = if('+name+'__all1'+' == 1,%s,null())' %(i+1),
                          overwrite=True)
            # make vector of catchment as well
            if 'catchmentprefix' in self.options:
                g_run('r.to.vect', overwrite=True, quiet=True, flags='vs',
                      input=name, output=name, type='area')
            # update maps dictionary
            self.allcatchments +=[name]

        return self.allcatchments
    
    def getTopology(self):
        '''Return a dictionary with staions as keys and lists of watershed ids
        as values'''
        # get station topology
        stopo=rwhat(self.allcatchments, self.station_coor)
        # list of numpy arrays with indeces of nonzero cat values
        stopo=stopo.transpose()
        topo=[]
        for i,d in enumerate(stopo):
            d[i]=0
            s=np.nonzero(d)[0]
            topo+=[s]
        topology = dict(zip(range(1,len(topo)+1),topo))
        return topology
        
    def makeSubbasins(self):
        """Create subbasins with corrected stations shape file
        with maximum threshold in square km from the maps processed in procDEM
        """
    
        out={} # for output maps
        # calculate station topology from individual watershed rasters
        self.stationtopology = self.getTopology()
        ########### SUBBASINS ##########################
        # use upthresh for all if self.upthresh not already converted to list in __init__
        if type(self.upthresh) in [int,float]:
            self.upthresh = [self.upthresh]*len(self.stationtopology)
        # report station topology and make subbasins
        watersheds=self.allcatchments
        subareas=[]
        subbasins=[]
        gm('Station topology:')
        for i,included in enumerate(self.stationtopology.values()):
            gm( 'Calculating subbasins for %s'%watersheds[i])
            ####### SUBAREAS #######
            if len(included) > 0:
                gm('includes watersheds: ')
                for ii in included: gm(watersheds[ii])
                # subarea name
                subarea_name='subareas__'+watersheds[i]
                # make list of watersheds to exclude i.e. the ones that are included
                masked_waters=['isnull('+watersheds[ii]+')' for ii in included]
                # make subarea excluding the upstream watersheds
                # if primary watershed is not (~) null and (&) the masked_waters
                # are null, i.e. the area downstream of the included watersheds
                exp=subarea_name+'=if(~isnull('+watersheds[i]+') & '+' & '.join(masked_waters)+',%s,null())' %(i+1)
                grass.mapcalc(exp, overwrite=True, quiet=True)
                subareas+=[subarea_name]
                ### check if outlet, ie. if stations-1 are included no others are downstream
                if len(included)==len(self.stationtopology)-1:
                    self.outletcoor = self.station_coor[i]
            else:
                subareas+=[watersheds[i]]
            #### MASK SUBAREA #####
            #g_run('r.mask', rast=subareas[i], overwrite=True, quiet=True)
                              
            # prepare inputs for the subbasins
            subbasins_name='subbasins__'+subareas[i]
            # calculate threshold from sq km to cells
            thresh = int(round(self.upthresh[i]*1000**2/(self.region['ewres']*self.region['nsres'])))
            grass.message('Subbasin threshold: %s km2, %s cells' %(self.upthresh[i],thresh))
            
            # check if already calculated with that threshold
            if thresh in self.subbasinsdone:
                subbasins_uncut = self.subbasinsdone[thresh]
                gm('Using %s, already calculated.' %subbasins_uncut)
            else:
                subbasins_uncut = subbasins_name+'__uncut'
                kwargs ={'elevation': self.elevation,
                         'basin'    : subbasins_uncut,
                         'threshold': thresh,
                         'flags'    : self.rwatershedflags}
                # carved elevation
                if 'streamcarve' in self.options:
                    kwargs['elevation']=self.carvedelevation
                    
                ##### r.watershed            
                g_run('r.watershed', overwrite=True, quiet=True, **kwargs)
                
                # add to done subbasins list
                self.subbasinsdone[thresh] = subbasins_uncut
                
            # cut out subbasins for subarea
            exp = subbasins_name + '=if(isnull(%s), null(), %s)' %(subareas[i],subbasins_uncut)
            grass.mapcalc(exp)
           
            # g_run('r.mask',flags='r',quiet=True) #remove mask
            subbasins+=[subbasins_name]
            
        #update maps dictionary
        out['subbasins']= subbasins
        out['subareas'] = subareas
    
    
        ### Make sure no subbasins have the same cat
        lastmax = 0 # in case only 1 station is used
        for i,s in enumerate(out['subbasins']):
            # if more than one subarea add the last max to the current
            if i>0:
                grass.mapcalc('{0}__lastmax={0} + {1}'.format(s,lastmax),
                                      overwrite=True, quiet=True)
                out['subbasins'][i] = s+'__lastmax'
            # get classes and check if subbasins were produced
            classes = getclasses(out['subbasins'][i])
            if len(classes)==0:
                gm( '%s has no subbasins and will be omitted (station too close to others?)' %s)
                continue
            lastmax =classes.max()
        
        ### PREDEFINED
        if 'predefined' in self.options:
            gm('Including predefined subbasins %s' %self.predefined)
            gm('Catchment boundaries disregarded, doublecheck %s' %self.catchments)
            # avoid same numbers occur in subbasins
            predef = self.predefined.split('@')[0]+'__aboverange'
            grass.mapcalc('%s=%s+%s' %(predef,self.predefined,lastmax),overwrite=True)
            # add to beginning of out subbasins
            out['subbasins'] = [predef]+out['subbasins']
            
        ### PATCHING
        grass.message('Patching catchments and subbasins')
        # subbasins maps
        patchbasins(out['subbasins'], outname=self.subbasins)
        # watersheds/subareas
        patchbasins(out['subareas'], outname=self.catchments)
        # make vector
        g_run('r.to.vect', overwrite=True, quiet=True, flags='vs',
                      input=self.catchments, output=self.catchments,
                      type='area')
        ### clean subbasin raster and vector keeping the same name
        self.cleanSubbasins()
        ### make continuous subbasinIDs
        self.continousCats()
        
        gm('Uploading catchmentID, elevation statistics, centroid coordinates, size to the subbasin attribute table...')
        # assign catchment id
        g_run('v.db.addcolumn', map=self.subbasins, column='catchmentID int',quiet=True)
        g_run('v.what.rast', map=self.subbasins, raster=self.catchments,
              column='catchmentID',quiet=True,type='centroid')
        
        # mean,min,max and centroid elevation and subbasin size
        cols = ['%s_elevation double' %s for s in ['average','max','min','centroid']]
        cols += ['size double','centroid_x double','centroid_y double']
        g_run('v.db.addcolumn', map=self.subbasins,quiet=True, column=','.join(cols))
        g_run('v.what.rast', map=self.subbasins, raster=self.elevation,
                   column='centroid_elevation',type='centroid',quiet=True)
        for s in ['min','max','average']:
            g_run('r.stats.zonal', base=self.subbasins,cover=self.elevation, method=s,
                   output='%s__elevation' %s, overwrite=True)
            #grass.mapcalc("%s__elevation=int(%s__elevation)" %(s,s), overwrite=True)
            g_run('v.what.rast', map=self.subbasins, raster='%s__elevation' %s,
                   column='%s_elevation' %s, type='centroid', quiet=True)
        # size
        g_run('v.to.db',map=self.subbasins,option='area',units='kilometers',
              columns='size',quiet=True)
        
        # centroid x,y
        g_run('v.to.db',map=self.subbasins,option='coor',units='meters',
              columns='centroid_x,centroid_y',quiet=True)
        
        # random colormap
        g_run('r.colors', quiet=True, map=self.subbasins, color='random')
        g_run('r.colors', quiet=True, map=self.catchments, color='random')
        
        grass.message('''Created %s and %s as raster and vector maps.
        ''' %(self.subbasins,self.catchments))
        
        return
        
    def cleanSubbasins(self):
        '''Make vector and remove areas smaller than lothresh'''
        
        #### clean subbasins
        grass.message('Clean subbasin map...')
        # add little areas of watershed to subbasin map that arent covered
        exp=''''subbasins__0'=if(~isnull('{1}') & isnull({0}),9999,'{0}')'''
        grass.mapcalc(exp.format(self.subbasins,self.catchments), overwrite=True)
        g_run('g.remove', rast=self.subbasins,quiet=True)
        g_run('g.rename', rast='subbasins__0,%s' %self.subbasins,quiet=True)

        # convert subbasins to vector
        g_run('r.to.vect', overwrite=True, quiet=False, flags='',
                          input=self.subbasins,
                          output=self.subbasins+'__unclean',
                          type='area')
        # remove small subbasins smaller than a thenth of threshold (m2)
        prunedist = float(np.mean(self.upthresh)*3)
        g_run('v.clean',  overwrite=True, quiet=True,
                          input=self.subbasins+'__unclean',
                          output=self.subbasins, flags='c',
                          type='area', tool='rmarea,prune',
                          thresh='%s,%s' %(self.lothresh*1000**2,prunedist))
        gm("Don't worry about these warnings!")
        g_run('v.build',map=self.subbasins,overwrite=True, quiet=True)
        
        return
        
    def continousCats(self):
        '''Make continous vector subbasin numbering with the outlet as 1, and update raster,
        in the process, also assigns drainage areas to subbasin table'''
        
        # TODO, not optimal: make raster and then vector again to have continuous cats
        g_run('v.to.rast', input=self.subbasins,output=self.subbasins,type='area',
              use='cat', overwrite=True, quiet=True)
        g_run('r.to.vect', overwrite=True, quiet=True, flags='s',type='area',
                          input=self.subbasins, output=self.subbasins)
        # delete default label column
        g_run('v.db.dropcolumn',map=self.subbasins, column='value,label', quiet=True)
        # add separate subbasinID column
        g_run('v.db.addcolumn', map=self.subbasins, columns='subbasinID int', quiet=True)
        g_run('v.db.update', map=self.subbasins, column='subbasinID', qcol='cat', quiet=True)
        
        
        # get drainage area via accumulation map in sq km
        g_run('r.stats.zonal', base=self.subbasins,cover=self.accumulation, method='max',
             output='max__accum__cells', overwrite=True)
        cellareakm = self.region['nsres']*self.region['ewres']*10**-6
        grass.mapcalc("max__accum=max__accum__cells*%s" %cellareakm, overwrite=True)

        # upload to subbasin table
        g_run('v.db.addcolumn', map=self.subbasins, column='darea double',quiet=True)
        g_run('v.what.rast', map=self.subbasins,raster='max__accum',column='darea',
              type='centroid',quiet=True)
        
        # change subbasin with the greatest drainage area (=outlet subbasin) to 1
        tbl = getTable(self.subbasins, dtype=(int,float),columns='subbasinID,darea')
        # get max cat and old 1 cat
        catmax = np.argmax(tbl['darea'])+1
        
        # swap both values
        g_run('v.db.update',map=self.subbasins, column='subbasinID',
                      where='cat=%s' %catmax, value=1)
        g_run('v.db.update',map=self.subbasins, column='subbasinID',
                      where='cat=1', value=catmax)        
        # reclass to subbasinID via copy
        g_run('g.copy',vect=self.subbasins+',unreclassed__subbasins',quiet=True)
        g_run('v.reclass', input='unreclassed__subbasins', output=self.subbasins,
              column='subbasinID',overwrite=True,quiet=True)
        g_run('v.db.addtable', map=self.subbasins,key='subbasinID',quiet=True)
        g_run('v.db.join', map=self.subbasins,column='subbasinID',
              otable='unreclassed__subbasins', ocolumn='subbasinID',quiet=True)
        g_run('v.db.dropcolumn',map=self.subbasins,column='cat',quiet=True)
        # make raster again
        g_run('v.to.rast',input=self.subbasins,output=self.subbasins,
              use='cat', overwrite=True, quiet=True)
        return


    def statistics(self):
        '''Output some statistics of the subbasin and subcatchment map'''
        # subcatchments
        scs = getAreas(self.catchments)
        # subbasin sizes
        sbs = getTable(self.subbasins,dtype=(int,int,float),
                       columns='subbasinID,catchmentID,size')
        
        gm('-----------------------------------------------------------------')
        print( '''Catchment sizes :
ID  excl. upstream   incl. upstream  outlet subbasin  upstream stations''')        
        outletsb = rwhat([self.subbasins],self.station_coor)
        for i,a in enumerate(scs):
            upstsize = np.sum(scs[self.stationtopology[a[0]],1])+a[1]
            print( '%3i %14.2f %16.2f %16i  %s' %(a[0],a[1],upstsize,
                                outletsb[i],self.stationtopology[a[0]]+1))

        # compile nice rows with total in the first column (first initialise dict, then add a column for each station)
        sub = {'st':'%8s ' %'total','n':'%8i ' %len(sbs),'min':'%8.2f ' %sbs['size'].min(),
               'mean':'%8.2f ' %sbs['size'].mean(), 'max':'%8.2f ' %sbs['size'].max()}
        cols = np.unique(sbs['catchmentID'])
        for c in cols:
            subs = sbs['size'][sbs['catchmentID']==c]
            if len(subs)==0: continue # in case sb outside catchments
            sub['st']   += '%8i ' %c
            sub['n']    += '%8i ' %len(subs)
            sub['min']  += '%8.2f ' %np.min(subs)
            sub['mean'] += '%8.2f ' %np.mean(subs)
            sub['max']  += '%8.2f ' %np.max(subs)
        
        print( '''
Subbasin statistics (km2):
Station: {st}
  Count: {n}
    Min: {min}
   Mean: {mean}
    Max: {max}'''.format(**sub))
        
        print('-----------------------------------------------------------------')        
        return scs,sbs


def rreclass(in_raster, in_list, out_list, proper=True):
    """Reclass a GRASS raster map from via an in list and outlist \n
    Patches together the rules, writes it to file, performs r.reclass,
    deletes in_raster and rules file and renames the outraster"""
    
    # temporary rules file
    temp_rules=grass.tempfile()
    # put lists in easy writable numpy array
    rules=np.array((in_list,out_list)).transpose()
    #write rules to file
    np.savetxt(temp_rules, rules ,delimiter='=', fmt='%i')
    # reclass raster in grass
    g_run('r.reclass',input=in_raster,
                      overwrite=True, quiet=True,
                      output=in_raster+'__',
                      rules=temp_rules)
    # make reclassed raster a proper raster, remove in_rast and rename output
    if proper:
        grass.mapcalc('__temp='+in_raster+'__', overwrite=True, quiet=True)
        g_run('g.remove', rast=in_raster+'__,'+in_raster, quiet=True)
        g_run('g.rename', rast='__temp,'+in_raster, quiet=True)
    return
    
def rwhat(rasters,coordinates):
    '''Get point values of rasters [list] at the coordinates [list of tuple pairs]
    '''
    # string of coordinate pairs
    coor_pairs=['%s,%s' %tuple(cs) for cs in coordinates]
    what=grass.read_command('r.what',
                             map=','.join(rasters),
                             null=0,
                             coordinates=','.join(coor_pairs),
                             separator=',').split('\n')[:-1]
    # put category values into numpy array of integers
    what_array=np.array([map(int,l.split(',')[-len(rasters):]) for l in what])
    
    return what_array
    
def getPoints(vector):
    '''Get point coordinates from point vector'''
    
    if vector==False:
        print 'Stations not defined!'
        return -1
    # get station coordinates from shapes[stations]
    coords=grass.read_command('v.report',
                              quiet=True,
                              map=vector,
                              option='coor').split('\n')[1:-1]
    #put coordinates into list of double tuples
    station_coor=[]
    for cs in coords: 
        cs=cs.split('|')
        station_coor+=[(float(cs[-3]),float(cs[-2]))]
    return station_coor

def snappoints(points, lines):
    '''correct points to lines by snapping'''

    # get distances to rivernetwork and nearest x and y
    snapped_points=grass.read_command('v.distance',flags='p',quiet=True,
                         _from=points,to=lines,
                         from_type='point',to_type='line',
                         upload='dist,to_x,to_y',column='d,x,y').split()
    # format, report and reassign station_coor
    snapped_coor=[]
    for d in snapped_points[1:]:
        d=d.split('|')
        snapped_coor+=[(float(d[2]),float(d[3]))]
        gm('Station %s moved by %8.1fm to: %s' %(d[0],float(d[1]),snapped_coor[-1]))

    return snapped_coor
        
def patchbasins(rastlist, outname):
    #patch all subbs together if more than one station
    sb_len=len(rastlist)
    if sb_len == 1:
        g_run('g.rename', quiet=True, rast=rastlist[0]+','+outname)
    elif sb_len > 1:
        g_run('r.patch',  input=','.join(rastlist),
                          output=outname,
                          overwrite=True, quiet=True)
    else:
        grass.fatal('No maps in out[subbasins]')
    return

def getclasses(raster):
    # get classes as integers in list
    classes=grass.read_command('r.stats', quiet=True,
                               input=raster, flags='n').split()
    return np.array(classes,dtype=int)



def getAreas(raster):
    subb_sizes=grass.read_command('r.stats',flags='na',
                                  input=raster).split()
    subs=np.array(subb_sizes, dtype=float).reshape((-1,2))
    subs[:,1] = subs[:,1]*10**-6
    return subs    

def getTable(vector,dtype='S250',**kw):
    '''Get a vector table into a numpy field array, dtype can either be one 
    for all or a list for each column'''
    tbl = grass.vector_db_select(vector,**kw)
    cols = tbl['columns']
    values = [tuple(row) for row in tbl['values'].values()]
    dtypes = {}
    if type(dtype) not in [list,tuple]:
        dtypes.update(dict(zip(cols,[dtype]*len(tbl['columns']))))
    elif len(dtype)!=len(cols):
        raise IOError('count of dtype doesnt match the columns!')
    else:
        dtypes.update(dict(zip(cols,dtype)))
    
    # first check for empty entries
    tbl = np.array(values,dtype=zip(cols,['S250']*len(cols)))
    convertedvals = []
    for c in cols:
        i = tbl[c]==''
        if len(tbl[c][i])>0:
            print 'Column %s has %s empty cells, will be parsed as float.' %(c,len(tbl[c][i]))
            if dtypes[c] in [float,int]:
                dtypes[c]=float
                tbl[c][i]='nan'
        # actual type conversion
        convertedvals += [np.array(tbl[c],dtype=dtypes[c])]
    # now properly make it
    tbl = np.array(zip(*convertedvals),dtype=[(c,dtypes[c]) for c in cols])
    # now set nans
    #for c in ix: tbl[c][ix[c]]=np.nan
    return tbl
    
if __name__=='__main__':
    # start time
    st = dt.datetime.now()
    # get options and flags
    o, f = grass.parser()
    grass.message(('GIS Environment:',grass.gisenv()))       
    grass.message(('Parameters:',o,f))
    
    # send all to main
    keywords = o; keywords.update(f)
    main=main(**keywords)
    
    # execute
    if main.s:
        main.statistics()
        sys.exit()
        
    # carve streams
    if 'streamcarve' in main.options:
        gm('Carving streams...')
        main.carveStreams()
        
    # process DEM if need be
    if not main.d:
        gm('Will process DEM first to derive accumulation, drainage and streams.')
        main.procDEM()
    
    # make catchments
    main.makeCatchments()
    # make subbasins
    main.makeSubbasins()
    #### write some statistics
    substats = main.statistics()
    
    # clean
    if not main.k:
        grass.run_command('g.mremove',rast='*__*',vect='*__*',flags='fb',quiet=True)
    # report time it took
    delta = dt.datetime.now()-st
    grass.message('Execution took %s hh:mm:ss' %delta)
