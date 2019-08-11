#!/usr/bin/env python
#
############################################################################
#
# MODULE:      m.swim.substats v1.2
# AUTHOR(S):   Michel Wortmann, wortmann@pik-potsdam.de
# PURPOSE:     Preprocessing suit for the Soil and Water Integrated Model (SWIM)
# COPYRIGHT:   (C) 2012-2019 by Wortmann/PIK
#
#              This program is free software under the GNU General Public
#              License (>=v2). Read the file COPYING that comes with GRASS
#              for details.
#
#############################################################################

#%Module
#% description: Soil and Water Integrated Model (SWIM) preprocessor: subbasin statistics
#%End

#%Option
#% guisection: Required
#% key: projectname
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% description: Name of project
#%end

#%Option
#% guisection: Required
#% key: projectpath
#% type: string
#% required: yes
#% multiple: no
#% key_desc: path
#% description: Path to project folder
#% gisprompt: old,dir,dir
#%end

#%Option
#% guisection: Required
#% key: subbasins
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% description: Subbasin vector map, statistics will be updated in table
#% gisprompt: old,vector,vector
#%end

#%Option
#% guisection: Subbasin
#% key: elev0
#% type: string
#% required: no
#% multiple: no
#% key_desc: raster or value
#% description: Reference elevation of the climate data to correct T and P.
#% answer: 0
#%end

#%Option
#% guisection: Subbasin
#% key: sdtsav
#% type: string
#% required: no
#% multiple: no
#% key_desc: raster or value
#% description: Initial water storage in the subbasin, m3
#% answer: 0
#%end

#%Option
#% guisection: Subbasin
#% key: sl
#% type: string
#% required: no
#% multiple: no
#% key_desc: raster or value
#% description: USLE slope length raster from r.watershed
#% answer: slopelength
#%end

#%Option
#% guisection: Subbasin
#% key: stp
#% type: string
#% required: no
#% multiple: no
#% key_desc: raster or value
#% description: USLE slope steepness raster from r.watershed
#% answer: slopesteepness
#%end

#%Option
#% guisection: Subbasin
#% key: salb
#% type: string
#% required: no
#% multiple: no
#% key_desc: raster or value
#% description: Soil albedo (default value or raster)
#% answer: 0.15
#%end

#%Option
#% guisection: Subbasin
#% key: sno
#% type: string
#% required: no
#% multiple: no
#% key_desc: raster or value
#% description: Initial snow water content, mm (default value or raster)
#% answer: 0
#%end

#%Option
#% guisection: Subbasin
#% key: ovn
#% type: string
#% required: no
#% multiple: no
#% key_desc: raster or value
#% description: Overland flow N value (default value or raster)
#% answer: 0.15
#%end

#%Option
#% guisection: Subbasin
#% key: rt
#% type: string
#% required: no
#% multiple: no
#% key_desc: raster or value
#% description: Return flow travel time, days (0=SWIM estimates, default value or raster)
#% answer: 0
#%end

#%Option
#% guisection: Subbasin
#% key: css
#% type: string
#% required: no
#% multiple: no
#% key_desc: raster or value
#% description: Channel slope, m/m (default value or raster)
#% answer: 0.5
#%end

#%Option
#% guisection: Subbasin
#% key: ecp
#% type: string
#% required: no
#% multiple: no
#% key_desc: raster or value
#% description: USLE erosion control practice factor P (default value or raster)
#% answer: 0.5
#%end

#%Option
#% guisection: Subbasin
#% key: suborder
#% type: string
#% required: no
#% multiple: no
#% key_desc: order
#% description: Order of values in .sub file (variables as given as arguments)
#% answer: salb,sno,chl,chs,chw,chk,chn,ovn,rt,css,ecp,sl,stp,lat,elev0,sdtsav
#%end

#%Option
#% guisection: Routing
#% key: mainstreams
#% type: string
#% required: no
#% multiple: no
#% key_desc: name
#% description: Main streams for each vector, made with m.swim.routing
#% gisprompt: old,raster,raster
#% answer: mainstreams
#%end

#%Option
#% guisection: Routing
#% key: elevation
#% type: string
#% required: no
#% multiple: no
#% key_desc: name
#% description: Elevation raster for slopes
#% gisprompt: old,cell,raster
#%end

#%Option
#% guisection: Routing
#% key: drainage
#% type: string
#% required: no
#% multiple: no
#% key_desc: name
#% description: drainage raster from m.swim.subbasins/r.watershed
#% gisprompt: old,cell,raster
#% answer: drainage
#%end

#%Option
#% guisection: Routing
#% key: accumulation
#% type: string
#% required: no
#% multiple: no
#% key_desc: name
#% description: accumulation raster from m.swim.subbasins/r.watershed
#% gisprompt: old,cell,raster
#% answer: accumulation
#%end

#%Option
#% guisection: Routing
#% key: chl
#% type: string
#% required: no
#% multiple: no
#% key_desc: raster or value
#% description: Main channel length, km (blank to calculate, otherwise value or raster)
#% answer:
#%end

#%Option
#% guisection: Routing
#% key: chs
#% type: string
#% required: no
#% multiple: no
#% key_desc: raster or value
#% description: Main channel slope, m/m (blank to calculate, value or raster)
#% answer:
#%end

#%Option
#% guisection: Routing
#% key: chw
#% type: string
#% required: no
#% multiple: no
#% key_desc: raster or value
#% description: Main channel width, m (blank to calculate, value or raster)
#% answer:
#%end

#%Option
#% guisection: Routing
#% key: chd
#% type: string
#% required: no
#% multiple: no
#% key_desc: raster or value
#% description: Main channel depth, m (blank to calculate, value or raster)
#% answer:
#%end

#%Option
#% guisection: Routing
#% key: chk
#% type: string
#% required: no
#% multiple: no
#% key_desc: raster or value
#% description: Effective hydraulic conductivity of main channel (default value or raster)
#% answer: 0.37
#%end

#%Option
#% guisection: Routing
#% key: chn
#% type: string
#% required: no
#% multiple: no
#% key_desc: raster or value
#% description: Channel N value (default value or raster)
#% answer: 0.075
#%end

#%Option
#% guisection: Routing
#% key: chxk
#% type: string
#% required: no
#% multiple: no
#% key_desc: raster or value
#% description: channel USLE K factor (default value or raster)
#% answer: 0.28
#%end

#%Option
#% guisection: Routing
#% key: chc
#% type: string
#% required: no
#% multiple: no
#% key_desc: raster or value
#% description: channel USLE C factor (default value or raster)
#% answer: 1
#%end

#%Option
#% guisection: Routing
#% key: rteorder
#% type: string
#% required: no
#% multiple: no
#% key_desc: order
#% description: Order of values in .rte file (see help for values)
#% answer: chw,chd,chs,chl,chn,chk,chxk,chc
#%end


#%Option
#% guisection: Groundwater
#% key: gwht
#% type: string
#% required: no
#% multiple: no
#% key_desc: raster or value
#% description: Initial groundwater height, m (default value or raster)
#% answer: 1
#%end

#%Option
#% guisection: Groundwater
#% key: gwq
#% type: string
#% required: no
#% multiple: no
#% key_desc: raster or value
#% description: Initial groundwater flow contribution to streamflow, mm (default value or raster)
#% answer: 0.5
#%end

#%Option
#% guisection: Groundwater
#% key: abf
#% type: string
#% required: no
#% multiple: no
#% key_desc: raster or value
#% description: Groundwater alpha factor (default value or raster)
#% answer: 0.048
#%end

#%Option
#% guisection: Groundwater
#% key: syld
#% type: string
#% required: no
#% multiple: no
#% key_desc: raster or value
#% description: Specific yield of aquifer (default value or raster)
#% answer: 0.003
#%end

#%Option
#% guisection: Groundwater
#% key: delay
#% type: string
#% required: no
#% multiple: no
#% key_desc: raster or value
#% description: Groundwater delay, days (default value or raster)
#% answer: 200
#%end

#%Option
#% guisection: Groundwater
#% key: revapc
#% type: string
#% required: no
#% multiple: no
#% key_desc: raster or value
#% description: Fraction of root zone percolation to revap (0-1 or raster)
#% answer: 0.2
#%end

#%Option
#% guisection: Groundwater
#% key: rchrgc
#% type: string
#% required: no
#% multiple: no
#% key_desc: raster or value
#% description: Fraction of root zone percolation into deep aquifer (0-1 or raster)
#% answer: 0.05
#%end

#%Option
#% guisection: Groundwater
#% key: revapmn
#% type: string
#% required: no
#% multiple: no
#% key_desc: raster or value
#% description: Initial revap storage, mm (default value or raster)
#% answer: 0
#%end

#%Option
#% guisection: Groundwater
#% key: gworder
#% type: string
#% required: no
#% multiple: no
#% key_desc: order
#% description: Order of values in .gw file (see help for values)
#% answer: gwht,gwq,abf,syld,delay,revapc,rchrgc,revapmn
#%end

#%Option
#% guisection: Optional
#% key: dummy
#% type: string
#% required: no
#% multiple: no
#% key_desc: raster or value
#% description: Dummy value as fill value and example
#% answer: 0
#%end

#%Flag
#% guisection: Optional
#% key: k
#% label: Keep intermediate files (those named *__*)
#%end

import grass.script as grass
import numpy as np
import os
import datetime as dt
grun = grass.run_command
gread= grass.read_command
gm   = grass.message

class main:
    def __init__(self,**optionsandflags):
        '''Check all input and prepare the processing'''
        # add all options and flags as attributes (only nonempty ones)
        self.options = {}
        for o in optionsandflags:
            if optionsandflags[o]!='': self.options[o] = optionsandflags[o]
        self.__dict__.update(self.options)

        ### parameters that have functions
        self.functions = {'chl':self.mainChannelLength, ### MAIN CHANNEL LENGTH
                          'chs':self.mainChannelSlope,  ### MAIN CHANNEL SLOPE
                          'chw':self.channelWidth,      ### MAIN CHANNEL WIDTH
                          'chd':self.channelDepth,      ### MAIN CHANNEL DEPTH
                          'lat':self.centroid_latitude}
        self.functionsrasts = ['elevation','drainage','accumulation','mainstreams']

        # get needed parameters for the three files
        self.orders = {'subbasin': self.suborder,
                       'groundwater': self.gworder,
                       'routing': self.rteorder}
        # check if set properly
        for f in self.orders:
            # if not set at all
            if len(self.orders[f])==0: grass.fatal('%sorder not set!' %f)
            # make lists
            self.orders[f] = self.orders[f].split(',')
            for p in self.orders[f]:
                # skip flu, subbasin area fractions
                if p=='flu': continue
                # not in options
                if p not in self.options:
                    if p in self.functions:
                        if any([r not in self.options for r in self.functionsrasts]):
                            grass.fatal('To calculate %s, %s must be given!' %(p,','.join(self.functionsrasts)))
                    else:
                        # if no function exists and not given as argument
                        grass.fatal('Not sure what to do with %s given in %sorder. No default value or raster given.' %(p,f))

        # check if project dir exists
        if not os.path.exists(self.projectpath):
                grass.fatal('projectpath %s doesnt exists.' %self.projectpath)

        ### .sub, .rte, .gw files location
        path=os.path.join(self.projectpath,'Sub')
        # check if Sub directory exists
        if not os.path.exists(path):
            grass.message('%s doesnt exist and will be created.' %path)
            os.mkdir(path)
        self.subpath = path

        # make rast from subbasin vector
        self.subbasinrast='subbasin__rast'
        grun('v.to.rast',input=self.subbasins,output=self.subbasinrast,
              use='cat', overwrite=True, quiet=True)

        # mask only subbasins
        grun('r.mask',rast=self.subbasinrast,quiet=True,overwrite=True)

        # decide if any parameter for the functions missing or if needed at all
        # channel width and depth
        if 'chw' not in self.options or 'chd' not in self.options:
            if 'accumulation' not in self.options:
                grass.fatal('If chw or chd is empty, accumulation needs to be set.')
        # channel length and slope
        if 'chl' not in self.options or 'chs' not in self.options:
            if 'mainstreams' not in self.options:
                grass.fatal('If chl or chs is empty, mainstreams needs to be set.')
            if 'chl' not in self.options and 'drainage' not in self.options:
                grass.fatal('If chl is empty, drainag needs to be set.')
            if 'chs' not in self.options and 'elevation' not in self.options:
                grass.fatal('If chs is empty, elevation needs to be set.')

            # only if chl and chs is empty, make proper mainstream rast
            self.mainstreamrast = self.makeMainStreamRast()

        # get current maps and environment
        self.maps = {t:grass.list_grouped(t) for t in ['rast','vect']}
        self.env  = grass.gisenv()
        return

    def subbasinStats(self):
        '''Calculate all subbasin stats for the Sub/ files, defaults and list of
        parameters given in 'swim_defaults.py that needs to be importable,
        all listed parameters can also be given as grass rasters, for each subbasin
        one value, e.g. a reclass of the subbasin map.
        Returns the dataframe compiled with named columns'''

        # get subbasin area fractions to initialise
        subbfractions   = self.areaFractions()
        self.nsubbasins = len(subbfractions)

        params = {'flu':subbfractions} # will be filled with all others
        for f in self.orders.keys():
            for p in self.orders[f]:
                # get each parameter
                if p not in params:
                    par = self.getParam(p)
                    # check if long enough
                    if len(par)!=self.nsubbasins:
                        grass.fatal("""
Can only find/calculate %s values for %s, but there are %s subbasins.""" %(len(par),p,self.nsubbasins))
                    # save for later writing
                    params[p] = par

        return params


    def getParam(self, param):
        '''Check if param is given as default, needs to be calculated or is a map'''

        # check if parameter given as input
        if param in self.options:
            # then its either a float or a map
            try: # default value given
                stats = float(self.options[param])*np.ones(self.nsubbasins)
                grass.message( 'Using default value for %s = %s' %(param,stats[0]))
            except: # map name given
                rast = self.options[param].split('@')
                rast,maps = {1: (rast[0],self.env['MAPSET']),2: rast}[len(rast)]
                if rast not in self.maps['rast'][maps]:
                    grass.fatal('%s not found.' %self.options[param])
                grass.message( 'Will use average subbasin values of %s for %s.' %(self.options[param],param))
                # upload mean values to subbasins table
                rastname = self.meanSubbasin(self.options[param])
                stats    = self.upload2Subbasins(rastname,
                                column=self.options[param].split('@')[0])

        else: # try to calculate it
            stats = self.functions[param]() # call the respective function
            if type(stats) == str:
                stats = self.upload2Subbasins(stats, column=stats)
            # correct channelLength
            if param=='chl': stats = self.correctChannelLength()
        # report statistics
        nn = ~np.isnan(stats)
        nnans = len(np.where(~nn)[0])
        gm("%s statistics:" % param)
        gm("min: %s mean: %s max: %s number of nans: %s" %
           (np.min(stats[nn]), np.mean(stats[nn]), np.max(stats[nn]), nnans))
        return stats

    def meanSubbasin(self,raster, method='average'):
        '''Upload mean value of raster to the table of self.subbasins and return
        values as array.
        '''
        base = self.subbasinrast
        tmpname='%s__%s__%s' %(base.split('@')[0],raster.split('@')[0],method)
        # average over the subbasins
        grun('r.stats.zonal',base=base,cover=raster,method=method,
             output=tmpname,overwrite=True,quiet=True)
        grass.message('%s of %s over each %s saved in %s' %(method,raster,base,tmpname))
        return tmpname

    def upload2Subbasins(self, raster, column):
        '''Upload the values given in raster at the centroids to the subbasins
        table and return them as array (sorted by subbasins cats/subbasinID'''
        # upload to vector
        grun('v.db.addcolumn',map=self.subbasins,columns=column+' double precision')
        grun('v.what.rast', map=self.subbasins, raster=raster, column=column,
             type='centroid',quiet=True)

        # get column out of table
        stats=getTable(self.subbasins,dtype=float,
                       columns='subbasinID,%s' %column)[column]
        return stats

    def makeMainStreamRast(self):
        # make stream rast with subbasin categories
        streamrast='%s__rast' %self.mainstreams.split('@')[0]
        grun('v.to.rast',input=self.mainstreams,type='line',output=streamrast+'__1',use='val',
             overwrite=True, quiet=True)
        grass.mapcalc("'{0}'=if(isnull('{1}'),null(),'{2}')".format(streamrast,streamrast+'__1','subbasin__rast'),
                      overwrite=True)

        # check if it has same number of categories than subbasins
        sbstats=rstats('subbasin__rast')
        strstats=rstats(streamrast)
        if len(sbstats)!=len(strstats):
            self.nostreams = [sb for sb in sbstats['id'] if sb not in strstats['id']]
            gm('''These subbasins dont seem to have a mainstream:''')
            gm(','.join(map(str,self.nostreams)))
            gm('''I will assume they have at least one stream cell.''')
        return streamrast

    def areaFractions(self):
        '''Calculate the fractions for each subbasin of the entire catchment'''
        subbcells=rstats(self.subbasinrast,flags='cn')
        fraction = np.float64(subbcells['value'])/np.sum(subbcells['value'])
        # report stats
        v = fraction
        grass.message('''Subbasin fraction statistics:
        min: %s mean: %s max: %s n: %s''' %(v.min(),v.mean(),v.max(),len(v)))

        return fraction

    def centroid_latitude(self):
        centroids = 'subbasin__centroids'
        grass.run_command('v.extract', input=self.subbasins, type='centroid',
                          output=centroids, flags='t', quiet=True, overwrite=True)
        tmpdir = grass.tempdir()
        tmploc = 'lonlat'
        grass.core.create_location(tmpdir, tmploc, epsg=4326)
        grass.run_command('g.mapset', mapset='PERMANENT', location=tmploc,
                          dbase=tmpdir, quiet=True)
        grass.run_command('v.proj', input=centroids, mapset=self.env['MAPSET'],
                          location=self.env['LOCATION_NAME'],
                          dbase=self.env['GISDBASE'], quiet=True)
        tbl = grass.read_command('v.report', map=centroids, option='coor')
        lat = np.array([l.split('|')[2] for l in tbl.split()[1:]], dtype=float)
        grass.run_command('g.mapset', mapset=self.env['MAPSET'],
                          location=self.env['LOCATION_NAME'],
                          dbase=self.env['GISDBASE'], quiet=True)
        return lat

    def mainChannelLength(self,rasterout='mainChannelLength'):
        '''Calculate the main channel length from the main channel vector
        for the subbasins'''
        grass.message('Calculating main channel length...')

        # make drainage direction weighted cell length
        streaminfo = rinfo(self.mainstreamrast)
        ewres, nsres = streaminfo['ewres'],streaminfo['nsres']
        exp = "'cell__len'=if(isnull('{streams}'),null(),"
        exp+= "if('{d}'==4 || '{d}'==8, {ew},0)+if('{d}' == 2 || '{d}' == 6,{ns},0)"
        exp+= "+if(%s,{dia},0))"%(" || ".join(["'{d}' == %s" %i for i in [1,3,5,7]]))
        exp = exp.format(d=self.drainage,streams=self.mainstreamrast,ew=ewres,ns=nsres,
                         dia=np.sqrt(float(ewres)**2+float(nsres)**2))
        grass.mapcalc(exp, overwrite=True)

        # report the sum of the cell length in the subbasins
        grun('r.stats.zonal',base=self.mainstreamrast,cover='cell__len',
             method='sum',output='cell__len__mainstreams__m',overwrite=True,quiet=True)

        # convert to km
        exp='cell__len__mainstreams=cell__len__mainstreams__m*0.001'
        grass.mapcalc(exp=exp, overwrite=True)
        # upload and get values
        lengthrast = self.meanSubbasin('cell__len__mainstreams')
        # change name to preserve raster
        grun('g.rename',rast=lengthrast+','+rasterout,quiet=True,overwrite=True)
        # set attributed for correction
        self.chl = rasterout

        return rasterout

    def correctChannelLength(self):
        '''Constrain channel length by subbasin perimeter and set nans to cell
        size, assuming that those without calculated chl have at least one cell as channel'''
        # get subbasin perimeter as an upper limit
        grun('v.db.addcolumn',map=self.subbasins,columns='perim__ double',quiet=True)
        grun('v.to.db',map=self.subbasins,columns='perim__',option='perimeter',
             units='kilometer',quiet=True)
        # get subbasin table
        tbl=getTable(self.subbasins,columns='subbasinID,%s,%s' %('perim__',self.chl),
                       dtype=[int,float,float])
        # where nan, set to grid size
        nolength = tbl[np.isnan(tbl[self.chl])]['subbasinID']
        if len(nolength)>0:
            where  = 'subbasinID IN (%s)' %(','.join(map(str,nolength)))
            res      = grass.region()['ewres']*1e-3 # m to km
            grun('v.db.update', map=self.subbasins, column=self.chl,
                 value=res, where=where,quiet=True)
            gm('%s subbasins have a minimum main channel length of %skm:' %(len(nolength),res))
            gm(','.join(map(str,nolength)))

        # too large
        toolong = tbl[tbl[self.chl]>tbl['perim__']]['subbasinID']
        if len(toolong)>0:
            gm('%s subbasins have a maximum main channel length equal to their perimeter:' %len(toolong))
            gm(','.join(map(str,toolong)))
            where  = 'subbasinID IN (%s)' %(','.join(map(str,toolong)))
            grun('v.db.update', map=self.subbasins, column=self.chl,
                 qcolumn='perim__', where=where,quiet=True)

            # remove perimeter column
            grun('v.db.dropcolumn',map=self.subbasins,column='perim__',quiet=True)
        # make raster again if needed
        if len(nolength)>0 or len(toolong)>0:
            grun('v.to.rast',input=self.subbasins,output=self.chl,use='attr',
                 attrcolumn=self.chl,type='area',overwrite=True,quiet=True)
        # return array of corrected lengths
        tbl=getTable(self.subbasins,dtype=[int,float],columns='subbasinID,%s' %self.chl)
        return tbl[self.chl]

    def mainChannelSlope(self,rasterout='mainChannelSlope'):
        '''Calculate main channel slope from the DEM and the streams'''
        grass.message('Calculating main channel slope...')

        # make slope
        sloperast='%s__slope' %self.elevation.split('@')[0]
        grun('r.slope.aspect',elevation=self.elevation,slope=sloperast,
             format='degrees',overwrite=True)

        # get mean values over main channels
        srast=sloperast+'__mainstreams'
        grun('r.stats.zonal',base=self.mainstreamrast,cover=sloperast,
             method='average',output=srast,overwrite=True,quiet=True)
        # convert degrees to m/m and minimum slope=0.0001 degrees
        grass.mapcalc(exp='%s=tan(max(%s,0.0001))' %('mainstream__slope__tan',srast),
                      overwrite=True)

        # mean over subbasins
        sloperastsubbasins = self.meanSubbasin('mainstream__slope__tan')
        # fill no data subbasins with minumum value
        exp='{0}=if(isnull({1}) & ~isnull({2}),0.0001,{1})'.format(rasterout,
             sloperastsubbasins,self.subbasinrast)
        grass.mapcalc(exp,overwrite=True)

        return rasterout

    def channelWidth(self,maxwidth = 3000,rasterout='channelWidth'):
        '''Calculate the mean channel width for each subbasin given with the drainage
        area / max accumulation according to an empirical approach:

        width = 1.29 * darea[km2] ^ 0.6
        maximum width = 3000m
        '''
        grass.message('Calculating main channel width...')
        # get the maximum accumulation from the accumulation map
        grun('r.stats.zonal',base=self.subbasinrast,cover=self.accumulation,
             method='max',output='max__accum',overwrite=True,quiet=True)
        res = rinfo(self.accumulation)['ewres']
        # calculate width
        exp=rasterout+'=min(1.29 * (max__accum * (%s^2)/1000000)^0.6, %s)' %(res,maxwidth)
        grass.mapcalc(exp, overwrite=True)

        return rasterout

    def channelDepth(self,maxdepth = 50, rasterout='channelDepth'):
        '''Calculate the mean channel depth for each subbasin given with the drainage
        area / max accumulation according to an empirical approach:

        depth = 0.13 * darea[km2] ^ 0.4 [m]
        maximum depth = 50m
        '''
        grass.message('Calculating main channel depth...')
        # get the maximum accumulation from the accumulation map
        grun('r.stats.zonal',base=self.subbasinrast,cover=self.accumulation,
             method='max',output='max__accum',overwrite=True,quiet=True)
        res = rinfo(self.accumulation)['ewres']
        # calculate width
        exp = rasterout+'=min(0.13 * (max__accum * (%s^2) / 1000000)^0.4, %s)' %(res,maxdepth)
        grass.mapcalc(exp,overwrite=True)

        return rasterout

    def writeSubFiles(self,data):
        '''Creates or overwrites files in the subpath with the
        .sub, .rte and .gw files from the data given and the structure given in parameters'''

        precision = '14.5'
        cfmt = '%-'+precision.split('.')[0]+'s'

        # build up formats
        grass.message('Writing subbasin.tab, routing.tab and '+
                      'groundwater.tab files to %s' %self.subpath)
        for p in self.orders:
            fname = os.path.join(self.subpath,'%s.tab' % p)
            tbl = np.column_stack([np.arange(1,self.nsubbasins+1)]+[data[c] for c in self.orders[p]])
            with open(fname,'w') as f:
                f.write(' '.join([cfmt%s for s in ['sub']+self.orders[p]])+'\n')
                fmts = '%14i'+ ' '.join(['%'+precision+'f']*len(self.orders[p]))
                np.savetxt(f, tbl, fmt=fmts)
        return

    def writeFileCio(self):
        '''Write the file.cio file for the SWIM input, subbasins should be
        array with subbasinIDs'''
        outname = os.path.join(self.projectpath, 'file.cio')
        # file.cio filename length in file list
        fnlen = 13
        # get subbasin cats
        subbasins = np.array(gread('r.stats',input='subbasin__rast',
                    flags='n', quiet=True).split(),dtype=int)
        # needed files / header
        named = ['.cod', '.fig', '.str', '.bsn', '.lut']
        headerlines = [  # number of empty lines at start
                       [self.projectname+ex for ex in named],
                       ('crop.dat','agman.dat','wgen.dat'),
                       ('soil.cio',),
                       ('runoff.dat',),
                       ('clim1.dat',),
                       ('clim2.dat',)] + [[]]*2  # number of empty lines at end
        # columns
        cols = [list(subbasins)] + [[0]*len(subbasins)]*4
        ilen = len(str(int(self.nsubbasins)))
        if ilen+len(self.projectname)+5 > fnlen:
            raise NameError('Name %s should be no longer than %s.' %(self.projectname,fnlen-ilen-5))
        for e in ['sub','rte','gw']:
            cols+= [[(self.projectname+'%0'+str(ilen)+'i.%s') %(n,e) for n in subbasins]]
        # write
        with open(outname, 'w') as f:
            # write header
            for l in headerlines:
                f.write(('%13s'*len(l)+'\n') % tuple(l))
                fmt = 5*'%4i'+3*(' %'+str(fnlen-1)+'s')+'\n'
            for l in zip(*cols):
                f.write(fmt % l)
        grass.message('Wrote %s' % outname)
        return



def rstats(rast,flags='n'):
    '''Return r.stats output as a sorted array'''
    values=grass.parse_command('r.stats',input=rast,flags=flags,separator='=')
    array = np.array(list(zip(values.keys(), values.values())),
                     dtype=[('id', int), ('value', float)])
    array.sort(order='id')
    return array

def rinfo(rast):
    return grass.parse_command('r.info',map=rast,flags='ge')
def runivar(rast):
    return grass.parse_command('r.univar',map=rast,flags='g')

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
    tbl = np.array(values, dtype=list(zip(cols, ['S250']*len(cols))))
    convertedvals = []
    for c in cols:
        i = tbl[c]==''
        if len(tbl[c][i]) > 0:
            grass.warning('Column %s has %s empty cells, will be parsed as '
                          'float.' % (c, len(tbl[c][i])))
            if dtypes[c] in [float,int]:
                dtypes[c]=float
                tbl[c][i]='nan'
        # actual type conversion
        convertedvals += [np.array(tbl[c],dtype=dtypes[c])]
    # now properly make it
    tbl = np.array(list(zip(*convertedvals)),
                   dtype=[(c, dtypes[c]) for c in cols])
    # now set nans
    #for c in ix: tbl[c][ix[c]]=np.nan
    return tbl


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

    # execute
    # calculate statistics and write .sub/.gw/.rte files
    data = main.subbasinStats()

    # Write to files
    main.writeSubFiles(data)

    ### file.cio
    main.writeFileCio()

    # clean
    grun('r.mask',flags='r')
    if not main.k:
        grass.run_command('g.remove',type='raster,vector', pattern='*__*',flags='fb',quiet=True)

    # report time it took
    delta = dt.datetime.now()-st
    grass.message('Execution took %s hh:mm:ss' %delta)
