#!/usr/bin/env python
#
############################################################################
#
# MODULE:      m.swim.subbasins v1.0
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
#% key: stations_snapped
#% type: string
#% required: no
#% multiple: no
#% key_desc: string
#% description: Stations as snapped to streams plus some additional info in the table
#% answer: stations_snapped
#% gisprompt: new,cell,vector
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

import sys
import numpy as np
import datetime as dt
from collections import OrderedDict
import grass.script as grass
grun = grass.run_command
gread = grass.read_command
gm = grass.message
gwarn = grass.warning
gdebug = grass.debug
gprogress = grass.core.percent


def interpret_options(optionsandflags):
    options = {}
    for o in optionsandflags:
        if optionsandflags[o] != '':
            try:
                options[o] = int(optionsandflags[o])             # int
            except ValueError:
                try:
                    options[o] = float(optionsandflags[o])       # float
                except ValueError:
                    options[o] = optionsandflags[o]  # str
    return options


class main:
    def __init__(self, **optionsandflags):
        '''Process all arguments and prepare processing'''
        # add all options and flags as attributes (only nonempty ones)
        self.options = interpret_options(optionsandflags)
        self.__dict__.update(self.options)

        # save region for convenience
        self.region = grass.region()
        self.region['kmtocell'] = lambda km: (int(round(np.mean(km)*10**6 /
                                  (self.region['ewres']*self.region['nsres']))))
        self.region['celltokm'] = lambda c: c*(self.region['ewres']*
                                               self.region['nsres'])*1e-6

        # check if DEM to processed or if all inputs set
        if not self.is_set('accumulation', 'drainage', 'streams'):
            grass.fatal('Either of these not set: accumulation, drainage, streams.')

        # what to do with upthresh
        if self.is_set('upthreshcolumn'):
            gm('Will look for upper thresholds in the %s column.' %
               self.upthreshcolumn)
            # get thresholds from column in station vect
            try:
                threshs = grass.vector_db_select(
                          self.stations, columns=self.upthreshcolumn)['values']
                self.upthresh = OrderedDict([(k, float(v[0]))
                                             for k, v in sorted(threshs.items())])
            except:
                grass.fatal('Cant read the upper threshold from the column %s'
                            % self.upthreshcolumn)

        # lothresh default
        if 'lothresh' not in self.options:
            self.lothresh = np.mean(self.upthresh.values())*0.05

        # streamthresh
        if 'streamthresh' in self.options:
            # convert to cells
            self.streamthresh = self.region['kmtocell'](self.streamthresh)
            # check if reasonable
            fract = float(self.streamthresh)/self.region['cells']
            if fract > 0.5 or fract < 0.01:
                grass.warning('streamthresh is %s percent of the region size!' %
                              (fract*100))
        else:
            self.streamthresh = int(self.region['cells']*0.02)

        # if no r.watershed flags given
        if 'rwatershedflags' not in self.options:
            self.rwatershedflags = 's'

        # check input for stats print
        if self.s:
            for o in ['streams', 'stations', 'catchmentprefix']:
                if not self.is_set(o):
                    grass.fatal('%s needs to be set!')
            # get all catchments
            rst = grass.list_strings('rast', self.catchmentprefix+'*')
            rolist = [(int(r.split('@')[0].replace(self.catchmentprefix, '')), r)
                      for r in sorted(rst) if '__' not in r]
            self.catchment_rasters = OrderedDict(rolist)
            gm('Found these catchments %s' % self.catchment_rasters)
            # calculate station topology
            self.snap_stations()
            self.get_stations_topology()

        # initialise subbasinsdone
        self.subbasinsdone = {}

        return

    def execute(self):
        # execute
        if self.s:
            self.print_statistics()
            sys.exit()

        # carve streams
        if 'streamcarve' in self.options:
            gm('Carving streams...')
            self.carve_streams()

        # process DEM if need be
        if not self.d:
            self.process_DEM()

        # always process
        self.snap_stations()
        self.make_catchments()
        self.get_stations_topology()
        self.make_subareas()
        self.make_subbasins()
        self.postprocess_catchments()
        self.postprocess_subbasins()
        self.write_stations_snapped()
        self.print_statistics()

        # clean
        if not self.k:
            grass.run_command('g.remove', type='raster,vector', pattern='*__*',
                              flags='fb', quiet=True)
        return

    def is_set(self, *options):
        return all([hasattr(self, i) for i in options])

    def process_DEM(self):
        """Make drainage, accumulation and stream raster"""
        gm('Processing DEM to derive accumulation, drainage and streams...')
        # decide on input arguments #######
        # from km2 to cells
        if type(self.upthresh) in [int, float]:
            uthresh = self.upthresh
        else:  # take most common one in upthresh column
            uthresh = max(set(self.upthresh.values()),
                          key=list(self.upthresh).count)

        thresh = self.region['kmtocell'](uthresh)

        kwargs = {'elevation': self.elevation,
                  'threshold': thresh,
                  # Output
                  'accumulation': 'accum__float',
                  'drainage': self.drainage,
                  'basin': 'standard__subbasins',
                  'slope_steepness': self.slopesteepness,
                  'length_slope': self.slopelength,
                  'flags': self.rwatershedflags}
        # check if depressions
        if self.is_set('depression'):
            kwargs['depression'] = self.depression
        # carve streams
        if self.is_set('streamcarve'):
            kwargs['elevation'] = self.carvedelevation

        grun('r.watershed', **kwargs)

        # save subbasins in dictionary
        self.subbasinsdone[thresh] = 'standard__subbasins'

        # postprocess accumulation map
        grass.mapcalc("%s=int(if(accum__float <= 0,null(),accum__float))" %
                      self.accumulation)
        # make river network to vector
        gm('Making vector river network...')
        # stream condition
        scon = '{0} >= {1}'.format(self.accumulation, self.streamthresh)
        # include those streams that were carved as well
        if 'streamcarve' in self.options:
            scon += ' | !isnull(%s)' % self.streamrastcarved
        # extract out of accumulation and make vector
        grass.mapcalc(self.streams+"__thick = if(%s, %s, null())" %
                      (scon, self.accumulation))
        grun('r.thin', input=self.streams+'__thick', output=self.streams,
             quiet=True)
        grun('r.to.vect', flags='s', input=self.streams, output=self.streams,
             type='line', quiet=True)

        return

    def carve_streams(self):
        '''Carve vector streams into the DEM, i.e. setting those cells =0'''
        gm('Carving %s into the elevation %s...' % (self.streamcarve, self.elevation))
        # stream vector to raster cells
        self.streamrastcarved = self.streamcarve.split('@')[0]+'__'
        grun('v.to.rast', input=self.streamcarve, output=self.streamrastcarved,
             type='line', use='val', val=1, quiet=True)
        # carve
        self.carvedelevation = '%s__carved' % self.elevation.split('@')[0]
        grass.mapcalc("%s=if(isnull(%s),%s,0)" % (self.carvedelevation,
                      self.streamrastcarved, self.elevation))
        return

    def snap_stations(self):
        '''Correct stations by snapping them to the streams vector.
        Snapped stations are written out to stations_snapped if given.
        '''
        warning_threshold = 1000  # m
        gm('Snapping stations to streams...')
        # types
        dtnames = ('stationID', 'distance', 'x', 'y')
        dtpy = (int, float, float, float)
        # get distances to rivernetwork and nearest x and y
        snapped_points = gread('v.distance', flags='p', quiet=True,
                               from_=self.stations, to=self.streams,
                               from_type='point', to_type='line',
                               upload='dist,to_x,to_y').split()
        # format, report and reassign stations_snapped_coor
        snapped_coor = np.array([tuple(d.split('|'))
                                 for d in snapped_points[1:]],
                                dtype=zip(dtnames, dtpy))
        # warn if above threshold
        snapped_over_thresh = snapped_coor[snapped_coor['distance'] >
                                           warning_threshold]
        if len(snapped_over_thresh) > 0:
            gwarn('These stations were moved further than %sm '
                  '(stationID: distance):' % warning_threshold)
            for i, d in snapped_over_thresh[['stationID', 'distance']]:
                gwarn('%i %1.0f' % (i, d))
        # save results
        lo = [(i, snapped_coor[i]) for i in dtnames]
        self.stations_snapped_columns = OrderedDict(lo)
        lo = [(i, (x, y)) for i, x, y in snapped_coor[['stationID', 'x', 'y']]]
        self.stations_snapped_coor = OrderedDict(lo)
        return

    def make_catchments(self):
        '''Make catchment raster and if catchmentprefix is set also vectors
        for all stations'''
        self.catchment_rasters = OrderedDict()

        nfmt = '%' + '0%ii' % len(str(max(self.stations_snapped_coor.keys())))

        gm('Creating catchments...')
        # create watersheds for stations
        for i, (si, (x, y)) in enumerate(self.stations_snapped_coor.items()):
            if self.is_set('catchmentprefix'):
                name = self.catchmentprefix + nfmt % si
            else:
                name = 'watersheds__st%s' % si
            gdebug(('station %s' % si))
            grun('r.water.outlet', input=self.drainage, output=name+'__all1',
                 coordinates='%s,%s' % (x, y))
            # give watershed number and put 0 to null()
            grass.mapcalc(name+' = if('+name+'__all1'+' == 1,%s,null())' % si,
                          quiet=True)
            # make vector of catchment as well
            if 'catchmentprefix' in self.options:
                grun('r.to.vect', quiet=True, flags='vs',
                     input=name, output=name, type='area')
            # update maps dictionary
            self.catchment_rasters[si] = name
            # report progress
            gprogress(i+1, len(self.stations_snapped_coor), 1)
        return

    def get_stations_topology(self):
        '''Return a dictionary with staions as keys and lists of watershed ids
        as values'''
        gm('Check station topology...')
        # get station topology
        stopo = rwhat(self.catchment_rasters.values(),
                      self.stations_snapped_coor.values())
        # list of numpy arrays with indeces of nonzero cat values
        stationid_array = np.array(self.stations_snapped_coor.keys())
        stopo = stopo.transpose()
        topo = []
        for i, d in enumerate(stopo):
            d[i] = 0
            s = np.nonzero(d)[0]
            topo += [stationid_array[s]]
        self.stations_topology = OrderedDict(zip(self.stations_snapped_coor.keys(), topo))

        # save downstream stationID
        ts = {k: len(v) for k, v in self.stations_topology.items() if len(v) > 0}
        ts = sorted(ts, key=ts.get)
        dsid = {}
        # find first occurence of id in length sorted downstream ids
        for i in sorted(self.stations_topology.keys()):
            for ii in ts:
                if i in self.stations_topology[ii] and i not in dsid:
                    dsid[i] = ii
                    break
            if i not in dsid:
                dsid[i] = -1
        dsid_ar = np.array([dsid[k] for k in sorted(dsid.keys())], dtype=int)
        self.stations_snapped_columns['ds_stationID'] = dsid_ar
        return

    def make_subareas(self):
        """Create catchment areas without headwater catchments."""
        gm('Creating catchment subareas...')
        self.subarea_rasters = OrderedDict()
        for i, (sid, included) in enumerate(self.stations_topology.items()):
            if len(included) > 0:
                # subarea name
                subarea_name = 'subareas__'+self.catchment_rasters[sid]
                # make list of watersheds to exclude i.e. the ones that are included
                masked_waters = ['isnull('+self.catchment_rasters[ii]+')'
                                 for ii in included]
                # make subarea excluding the upstream watersheds
                # if primary watershed is not (~) null and (&) the masked_waters
                # are null, i.e. the area downstream of the included watersheds
                exp = (subarea_name + '=if(~isnull(' +
                       self.catchment_rasters[sid] + ') & ' +
                       ' & '.join(masked_waters) + ', %s, null())' % sid)
                grass.mapcalc(exp, quiet=True)
                self.subarea_rasters[sid] = subarea_name
                # check if outlet, ie. if stations-1 are included no others are downstream
                if len(included) == len(self.stations_topology) - 1:
                    self.outletcoor = self.stations_snapped_coor[sid]
            else:
                self.subarea_rasters[sid] = self.catchment_rasters[sid]
            # report progress
            gprogress(i+1, len(self.stations_topology), 1)
        # path watersheds/subareas
        patch_basins(self.subarea_rasters.values(), outname=self.catchments)
        return

    def make_subbasins(self):
        """Create subbasins with corrected stations shape file
        with maximum threshold in square km from the maps processed in process_DEM
        """
        self.subbasins_rasters = OrderedDict()

        # use upthresh for all if self.upthresh not already converted to list in __init__
        if type(self.upthresh) in [int, float]:
            self.upthresh = OrderedDict([(i, self.upthresh)
                                         for i in self.stations_topology])

        for i, sid in enumerate(self.stations_topology.keys()):
            # prepare inputs for the subbasins
            subbasins_name = 'subbasins__' + self.subarea_rasters[sid]
            # calculate threshold from sq km to cells
            thresh = int(round(self.upthresh[sid] * 1000**2 /
                               (self.region['ewres']*self.region['nsres'])))
            gdebug('Subbasin threshold: %s km2, %s cells' %
                   (self.upthresh[sid], thresh))

            # check if already calculated with that threshold
            if thresh in self.subbasinsdone:
                subbasins_uncut = self.subbasinsdone[thresh]
                gdebug('Using %s, already calculated.' % subbasins_uncut)
            else:
                subbasins_uncut = subbasins_name+'__uncut'
                kwargs = {'elevation': self.elevation,
                          'basin'    : subbasins_uncut,
                          'threshold': thresh,
                          'flags'    : self.rwatershedflags}
                # carved elevation
                if 'streamcarve' in self.options:
                    kwargs['elevation'] = self.carvedelevation

                # r.watershed
                grun('r.watershed', overwrite=True, quiet=True, **kwargs)

                # add to done subbasins list
                self.subbasinsdone[thresh] = subbasins_uncut

            # cut out subbasins for subarea
            exp = (subbasins_name + '=if(isnull(%s), null(), %s)' %
                   (self.subarea_rasters[sid], subbasins_uncut))
            grass.mapcalc(exp)
            self.subbasins_rasters[sid] = subbasins_name
            # report progress
            gprogress(i+1, len(self.stations_topology), 1)

        # Make sure no subbasins have the same cat
        lastmax = 0  # in case only 1 station is used
        for i, (sid, srast) in enumerate(self.subbasins_rasters.items()):
            # if more than one subarea add the last max to the current
            if i > 0:
                grass.mapcalc('{0}__lastmax={0} + {1}'.format(srast, lastmax),
                              quiet=True)
                self.subbasins_rasters[sid] = srast + '__lastmax'
            # get classes and check if subbasins were produced
            classes = gread('r.stats', input=self.subbasins_rasters[sid],
                            quiet=True, flags='n').split()
            if len(classes) == 0:
                gwarn('%s has no subbasins and will be omitted'
                      ' (station too close to others?)' % srast)
                continue
            lastmax = max(classes)
            # report progress
            gprogress(i+1, len(self.subbasins_rasters), 1)

        if self.is_set('predefined'):
            gm('Including predefined subbasins %s...' % self.predefined)
            gwarn('Catchment boundaries are disregarded, doublecheck %s' %
                  self.catchments)
            # avoid same numbers occur in subbasins
            predef = self.predefined.split('@')[0]+'__aboverange'
            grass.mapcalc('%s=%s+%s' % (predef, self.predefined, lastmax))
            # add to beginning of subbasins_rasters
            self.subbasins_rasters = OrderedDict([('predefined', predef)] +
                                                 self.subbasins_rasters.items())

        # PATCHING subbasins maps
        patch_basins(self.subbasins_rasters.values(),
                     outname=self.subbasins)

        # clean subbasin raster and vector keeping the same name
        self.clean_subbasins()

        return

    def postprocess_catchments(self):
        gm('Creating catchments vector map...')
        # ## make vector from subbasins to match those
        # grun('v.dissolve', quiet=True, input=self.subbasins,
        #      output=self.catchments, column='catchmentID')
        # grun('v.what.rast', map=self.catchments, raster=self.catchments,
        #      column='catchmentID', type='centroid', quiet=True)
        grun('r.to.vect', input=self.catchments, output=self.catchments,
             type='area', flags='vst', quiet=True)
        grun('v.db.addtable', map=self.catchments, quiet=True,
             key='catchmentID', column='size double')
        grun('v.to.db', map=self.catchments, option='area', units='kilometers',
             columns='size', quiet=True)

        grun('r.colors', quiet=True, map=self.catchments, color='random')
        return

    def postprocess_subbasins(self):
        gm('Adding subbasin info to subbasin attribute table...')
        # add subbasinIDs to stations_snapped and write out stations_snapped
        ds_sbid = rwhat([self.subbasins],
                        self.stations_snapped_coor.values()).flatten()
        self.stations_snapped_columns['outlet_subbasinID'] = ds_sbid

        # assign catchment id
        grun('v.db.addcolumn', map=self.subbasins, column='catchmentID int',
             quiet=True)
        grun('v.what.rast', map=self.subbasins, raster=self.catchments,
             column='catchmentID', quiet=True, type='centroid')

        # mean,min,max and centroid elevation and subbasin size
        cols = ['%s_elevation double' % s
                for s in ['average', 'max', 'min', 'centroid']]
        cols += ['size double', 'centroid_x double', 'centroid_y double']
        grun('v.db.addcolumn', map=self.subbasins, quiet=True,
             column=','.join(cols))
        grun('v.what.rast', map=self.subbasins, raster=self.elevation,
             column='centroid_elevation', type='centroid', quiet=True)
        for s in ['min', 'max', 'average']:
            grun('r.stats.zonal', base=self.subbasins, cover=self.elevation,
                 method=s, output='%s__elevation' % s, quiet=True)
            grun('v.what.rast', map=self.subbasins, raster='%s__elevation' % s,
                 column='%s_elevation' % s, type='centroid', quiet=True)
        # size
        grun('v.to.db', map=self.subbasins, option='area', units='kilometers',
             columns='size', quiet=True)

        # centroid x,y
        grun('v.to.db', map=self.subbasins, option='coor', units='meters',
             columns='centroid_x,centroid_y', quiet=True)

        # random colormap
        grun('r.colors', quiet=True, map=self.subbasins, color='random')
        return

    def clean_subbasins(self):
        '''Make vector and remove areas smaller than lothresh and
        make continous vector subbasin numbering with the outlet as 1
        and update raster in the process.
        Also assigns drainage areas to subbasin table
        '''
        gm('Cleaning subbasin map...')
        tmp_subbasins = '%s__unclean' % self.subbasins
        # rename to ovoid overwrite
        grun('g.rename', raster=self.subbasins+','+tmp_subbasins, quiet=True)
        # add little areas of watershed to subbasin map that arent covered
        exp = "subbasins__0=if(~isnull('{1}') & isnull({0}),9999,'{0}')"
        grass.mapcalc(exp.format(tmp_subbasins, self.catchments))

        # convert subbasins to vector
        grun('r.to.vect', quiet=True, flags='', input='subbasins__0',
             output=self.subbasins + '__unclean', type='area')
        # remove small subbasins smaller than a thenth of threshold (m2)
        prunedist = float(np.mean(self.upthresh.values())*3)
        subbasins_cleaned = self.subbasins + '__cleaned'
        grun('v.clean', quiet=True, input=self.subbasins+'__unclean', flags='bc',
             output=subbasins_cleaned, type='area', tool='rmarea,prune',
             thresh='%s,%s' % (self.lothresh*1000**2, prunedist))
        grun('v.build', map=subbasins_cleaned, quiet=True)

        gm('Assigning continuous categories to subbasins map...')
        # TODO, not optimal: make raster and then vector again to have continuous cats
        subbasins_continuous = self.subbasins + '__continuous'
        grun('v.to.rast', input=subbasins_cleaned, output=subbasins_continuous,
             type='area', use='cat', quiet=True)
        grun('r.to.vect', quiet=True, flags='s', type='area',
             input=subbasins_continuous, output=subbasins_continuous)

        # delete default label column
        grun('v.db.dropcolumn', map=subbasins_continuous, column='value,label',
             quiet=True)
        # add separate subbasinID column
        grun('v.db.addcolumn', map=subbasins_continuous, columns='subbasinID int',
             quiet=True)
        grun('v.db.update', map=subbasins_continuous, column='subbasinID', qcol='cat',
             quiet=True)

        # get drainage area via accumulation map in sq km
        grun('r.stats.zonal', base=subbasins_continuous, cover=self.accumulation,
             method='max', output='max__accum__cells', quiet=True)
        cellareakm = self.region['nsres']*self.region['ewres']*10**-6
        grass.mapcalc("max__accum=max__accum__cells*%s" % cellareakm, quiet=True)

        # upload to subbasin table
        grun('v.db.addcolumn', map=subbasins_continuous, column='darea double',
             quiet=True)
        grun('v.what.rast', map=subbasins_continuous, raster='max__accum',
             column='darea', type='centroid', quiet=True)

        # change subbasin with the greatest drainage area (=outlet subbasin) to 1
        tbl = get_table(subbasins_continuous, dtype=(int, float),
                        columns='subbasinID,darea')
        # get max cat and old 1 cat
        catmax = np.argmax(tbl['darea'])+1

        # swap both values
        grun('v.db.update', map=subbasins_continuous, column='subbasinID',
             where='cat=%s' % catmax, value=1)
        grun('v.db.update', map=subbasins_continuous, column='subbasinID',
             where='cat=1', value=catmax)
        # reclass to subbasinID via copy
        grun('g.copy', vect=subbasins_continuous + ',unreclassed__subbasins',
             quiet=True)
        # TODO: here final self.subbasins vector is created
        grun('v.reclass', input='unreclassed__subbasins', output=self.subbasins,
             column='subbasinID', quiet=True)
        grun('v.db.addtable', map=self.subbasins, key='subbasinID', quiet=True)
        grun('v.db.join', map=self.subbasins, column='subbasinID',
             otable='unreclassed__subbasins', ocolumn='subbasinID', quiet=True)
        grun('v.db.dropcolumn', map=self.subbasins, column='cat', quiet=True)
        # make raster again
        # TODO: here final self.subbasins raster is created
        grun('v.to.rast', input=self.subbasins, output=self.subbasins,
             use='cat', quiet=True)
        return

    def write_stations_snapped(self):
        """Write out the stations_snapped to a vector."""
        types = {'i': 'int', 'f': 'double'}
        # columns
        cols = self.stations_snapped_columns
        cols_dt = [' '.join([i, types[cols[i].dtype.kind]]) for i in cols.keys()]
        cols_fmt = '|'.join(['%'+cols[i].dtype.kind for i in cols.keys()])
        data = np.column_stack(cols.values())
        # create vector if needed
        p = grass.feed_command('v.in.ascii', input='-', x=3, y=4, cat=1, quiet=True,
                               columns=cols_dt, output=self.stations_snapped)
        np.savetxt(p.stdin, data, delimiter='|', fmt=cols_fmt)
        p.stdin.close()
        p.wait()

    def print_statistics(self):
        '''Output some statistics of the subbasin and subcatchment map'''
        # subcatchments
        scs = get_table(self.catchments, columns='catchmentID,size',
                        dtype=(int, float))
        # subbasin sizes
        sbs = get_table(self.subbasins, dtype=(int, int, float),
                        columns='subbasinID,catchmentID,size')
        outletsb = rwhat([self.subbasins], self.stations_snapped_coor.values())
        gm('-----------------------------------------------------------------')
        print('''Catchment sizes :
ID  excl. upstream   incl. upstream  outlet subbasin  upstream stations''')
        for i, a in enumerate(scs):
            upix = [np.where(scs['catchmentID'] == c)[0][0]
                    for c in self.stations_topology[a[0]] if c in scs['catchmentID']]
            upstsize = np.sum(scs['catchmentID'][upix])+a[1]
            print('%3i %14.2f %16.2f %16i  %s' % (a[0], a[1], upstsize,
                                                  outletsb[i], self.stations_topology[a[0]]))

        # compile nice rows with total in the first column (first initialise dict, then add a column for each station)
        sub = {'st': '%8s ' % 'total',
               'n': '%8i ' % len(sbs),
               'min': '%8.2f ' % sbs['size'].min(),
               'mean': '%8.2f ' % sbs['size'].mean(),
               'max': '%8.2f ' % sbs['size'].max()}
        cols = np.unique(sbs['catchmentID'])
        for c in cols:
            subs = sbs['size'][sbs['catchmentID']==c]
            if len(subs) == 0:
                continue  # in case sb outside catchments
            sub['st']   += '%8i ' % c
            sub['n']    += '%8i ' % len(subs)
            sub['min']  += '%8.2f ' % np.min(subs)
            sub['mean'] += '%8.2f ' % np.mean(subs)
            sub['max']  += '%8.2f ' % np.max(subs)

        print('''
Subbasin statistics (km2):
Station: {st}
  Count: {n}
    Min: {min}
   Mean: {mean}
    Max: {max}'''.format(**sub))

        print('-----------------------------------------------------------------')
        return scs, sbs


def rreclass(in_raster, in_list, out_list, proper=True):
    """Reclass a GRASS raster map from via an in list and outlist \n
    Patches together the rules, writes it to file, performs r.reclass,
    deletes in_raster and rules file and renames the outraster"""

    # temporary rules file
    temp_rules = grass.tempfile()
    # put lists in easy writable numpy array
    rules = np.array((in_list, out_list)).transpose()
    # write rules to file
    np.savetxt(temp_rules, rules, delimiter='=', fmt='%i')
    # reclass raster in grass
    grun('r.reclass', input=in_raster,
         overwrite=True, quiet=True,
         output=in_raster + '__',
         rules=temp_rules)
    # make reclassed raster a proper raster, remove in_rast and rename output
    if proper:
        grass.mapcalc('__temp=' + in_raster + '__', quiet=True)
        grun('g.remove', type='rast', name=in_raster + '__,' + in_raster,
             flags='f', quiet=True)
        grun('g.rename', rast='__temp,' + in_raster, quiet=True)
    return


def rwhat(rasters, coordinates):
    '''Get point values of rasters [list] at the coordinates [list of tuple pairs]
    '''
    # string of coordinate pairs
    coor_pairs = ['%s,%s' % tuple(cs) for cs in coordinates]
    what = gread('r.what',
                 map=','.join(rasters),
                 null=0,
                 coordinates=','.join(coor_pairs),
                 separator=',').split('\n')[:-1]
    # put category values into numpy array of integers
    what_array = np.array(
        [map(int, l.split(',')[-len(rasters):]) for l in what])

    return what_array


def patch_basins(rastlist, outname):
    # patch all subbs together if more than one station
    sb_len = len(rastlist)
    if sb_len == 1:
        grun('g.rename', quiet=True, rast=rastlist[0] + ',' + outname)
    elif sb_len > 1:
        grun('r.patch',  input=','.join(rastlist),
             output=outname,
             overwrite=True, quiet=True)
    else:
        grass.fatal('No maps in out[subbasins]')
    return


def get_table(vector, dtype='S250', **kw):
    '''Get a vector table into a numpy field array, dtype can either be one
    for all or a list for each column'''
    tbl = grass.vector_db_select(vector, **kw)
    cols = tbl['columns']
    values = [tuple(row) for row in tbl['values'].values()]
    dtypes = {}
    if type(dtype) not in [list, tuple]:
        dtypes.update(dict(zip(cols, [dtype] * len(tbl['columns']))))
    elif len(dtype) != len(cols):
        raise IOError('count of dtype doesnt match the columns!')
    else:
        dtypes.update(dict(zip(cols, dtype)))

    # first check for empty entries
    tbl = np.array(values, dtype=zip(cols, ['S250'] * len(cols)))
    convertedvals = []
    for c in cols:
        i = tbl[c] == ''
        if len(tbl[c][i]) > 0:
            gm('Column %s has %s empty cells, will be parsed as float.' %
               (c, len(tbl[c][i])))
            if dtypes[c] in [float, int]:
                dtypes[c] = float
                tbl[c][i] = 'nan'
        # actual type conversion
        convertedvals += [np.array(tbl[c], dtype=dtypes[c])]
    # now properly make it
    tbl = np.array(zip(*convertedvals), dtype=[(c, dtypes[c]) for c in cols])
    tbl.sort()
    return tbl


if __name__ == '__main__':
    # start time
    st = dt.datetime.now()
    # get options and flags
    o, f = grass.parser()
    fmt = lambda d: '\n'.join(['%s: %s' % (k, v) for k, v in d.items()])+'\n'
    grass.message('GIS Environment:\n'+fmt(grass.gisenv()))
    grass.message('Parameters:\n'+fmt(o)+fmt(f))

    # warn if MASKED
    if 'MASK' in grass.list_grouped('rast')[grass.gisenv()['MAPSET']]:
        maskcells = gread('r.stats', input='MASK', flags='nc').split()[1]
        grass.message('!!! MASK active with %s cells, will only process those !!!'
                      % maskcells)

    # send all to main
    keywords = o
    keywords.update(f)
    main = main(**keywords)

    main.execute()

    # report time it took
    delta = dt.datetime.now() - st
    grass.message('Execution took %s hh:mm:ss' % delta)
