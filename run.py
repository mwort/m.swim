#!/usr/bin/env python
############################################################################
#
# MODULE:       m.swim.run
# AUTHOR(S):    wortmann
# PURPOSE:      GRASS frontend for the SWIMpy python classes and functions
# COPYRIGHT:    (C) 2014 by M. Wortmann, and the GRASS Development Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
############################################################################

#%module
#% description: SWIM interface: running and postprocessing SWIM
#% keywords: SWIM, hydrological modelling
#%end


#%option
#% key: notes
#% type: string
#% required: no
#% multiple: no
#% key_desc: string
#% label: Notes to be saved with run
#% description: 
#% guisection: Run
#%end

#%option
#% key: purpose
#% type: string
#% required: no
#% multiple: no
#% key_desc: string
#% options: calibration,validation,autocalibration,testing
#% label: Purpose of run
#% description: Categories to later group runs better
#% answer: calibration
#% guisection: Run
#%end

#%option
#% key: interactive
#% type: string
#% required: no
#% multiple: no
#% key_desc: string
#% label: Station to observe interactively while running SWIM
#% description: Needs to be a valid station ID
#% guisection: Run
#%end

#%Flag
#% key: l
#% label: Send to llsubmit
#% guisection: Run
#%end

#%Flag
#% key: u
#% label: Preprocess SWIM files only (SWIM not started and no postprocessing)
#% description: e.g. start SWIM from commandline
#% guisection: Run
#%end

### TIME ##########################################################
#%option
#% key: iyr
#% type: integer
#% required: no
#% multiple: no
#% key_desc: year
#% options: 1900-3000
#% label: Start year
#% description: First year from which SWIM starts running
#% guisection: Time
#%end

#%option
#% key: nbyr
#% type: integer
#% required: no
#% multiple: no
#% key_desc: years
#% options: 0-500
#% label: Number of years
#% description: Run for that many years
#% guisection: Time
#%end

### BASIN PARAMETERS #########################################################
#%option
#% key: paramid
#% type: integer
#% required: no
#% multiple: no
#% key_desc: paramID
#% label: Parameter ID of previous run
#% description: 
#% guisection: Basin parameters
#%end
#%option
#% key: setsubcatch
#% type: string
#% required: no
#% multiple: yes
#% key_desc: stationID:paramID
#% label: Update the distributed subcatch parameters for each station (stationID:paramID,...)
#% description: Will update the subcatch.bsn file with the parameter values in paramID
#% guisection: Basin parameters
#%end
#%option
#% key: bff
#% type: double
#% required: no
#% multiple: no
#% key_desc: rate
#% label: Baseflow factor
#% description: 
#% guisection: Basin parameters
#%end

#%option
#% key: ecal
#% type: double
#% required: no
#% multiple: no
#% key_desc: rate
#% label: Evapotranspiration correction factor
#% description: 
#% guisection: Basin parameters
#%end

#%option
#% key: roc2
#% type: double
#% required: no
#% multiple: no
#% key_desc: rate
#% label: Routing factor for overland flow
#% description: 
#% guisection: Basin parameters
#%end

#%option
#% key: roc4
#% type: double
#% required: no
#% multiple: no
#% key_desc: rate
#% label: Routing factor for subsurface flow
#% description: 
#% guisection: Basin parameters
#%end

#%option
#% key: sccor
#% type: double
#% required: no
#% multiple: no
#% key_desc: rate
#% label: Saturated conductivity
#% description: 
#% guisection: Basin parameters
#%end

#%option
#% key: smrate
#% type: double
#% required: no
#% multiple: no
#% key_desc: rate
#% label: Snow melt rate, cm/deg/day
#% description: Intensity of snow melt in cm per degree per day
#% guisection: Basin parameters
#%end

#%option
#% key: gmrate
#% type: double
#% required: no
#% multiple: no
#% key_desc: rate
#% label: Glacier melt rate, mm/deg/day
#% description: Intensity of glacier melt in mm per degree per day
#% guisection: Basin parameters
#%end

### COMPONENTS #############################################################
#%Flag
#% key: s
#% label: Estimate saturated conductivity (isc)
#% description: instead of reading it from the soil files
#% guisection: Components
#%end
#%Flag
#% key: n
#% label: Use global curve numbers cnum1/cnum3 (icn)
#% description: Set cnum1 and cnum3 in .bsn file
#% guisection: Components
#%end
#%Flag
#% key: i
#% label: Consider interception storage (intercep)
#% guisection: Components
#%end

#%Flag
#% key: t
#% label: Use Turc-Ivanov evapotranspiration (iemeth)
#% description: Instead of Prestley-Taylor
#% guisection: Components
#%end
#%Flag
#% key: r
#% label: Use reservoir module (res_switch)
#% guisection: Components
#%end
#%Flag
#% key: c
#% label: Use distributed subcatchment parameters (subcatch)
#% description: subcatch.def/.bsn need to exist
#% guisection: Components
#%end
#%Flag
#% key: g
#% label: Use irrigation module (iirrig)
#% description: irrig.dat needs to exist
#% guisection: Components
#%end

components = {'s':'isc', 'n':'icn','i':'intercep','t':'iemeth','r':'res_switch',
              'c':'subcatch','g':'iirrig'}

### PEST #####################################################################

#%Flag
#% key: a
#% label: Run with PEST
#% guisection: Autocalibration
#%end

#%option
#% key: peststation
#% type: string
#% required: no
#% multiple: no
#% key_desc: stationID
#% label: Autocalibrate to observations of
#% description:
#% guisection: Autocalibration
#%end

#%option
#% key: pestinstructions
#% type: string
#% required: no
#% key_desc: path
#% label: Pest parameter file header 
#% description: Only upto the parameter data section, the rest is produced automatically. Comments are lines starting #
#% gisprompt: old,file,file
#% guisection: Autocalibration
#%end

### POSTPROCESSING ##########################################################
#%Flag
#% key: p
#% label: Postprocessing only
#% guisection: Postprocess
#%end

#%Flag
#% key: v
#% label: Save run statistics for stations (only valid if -p set)
#% guisection: Postprocess
#%end

#%option
#% key: savestations
#% type: string
#% required: no
#% multiple: yes
#% key_desc: stationIDs
#% label: Save stations runs (either 'all' or stationIDs)
#% description: Stations for which to save statistics and do the postprocessing
#% answer: all
#% guisection: Postprocess
#%end

#%option
#% key: showruns
#% type: string
#% required: no
#% multiple: yes
#% key_desc: stationIDs or runIDs
#% label: Show simulated vs observed Q for last run if stationID (or 'all') or for a particular runID
#% description: If multiple for same station, the runs will be plotted in same plot.
#% guisection: Postprocess
#%end

#%option
#% key: upsubbasins
#% type: string
#% required: no
#% multiple: yes
#% key_desc: files
#% label: Upload subbasin values to subbasin raster
#% description:
#% guisection: Postprocess
#%end

#%option
#% key: uphydrotopes
#% type: string
#% required: no
#% multiple: yes
#% key_desc: files
#% label: Upload hydrotope values to hydrotope raster
#% description:
#% guisection: Postprocess
#%end


### SETUP mySWIM ############################################################

#%Flag
#% key: 0
#% label: Set up a SWIM project or change any of the parameters below
#% guisection: Settings
#%end

#%option
#% key: prodir
#% type: string
#% required: no
#% key_desc: path
#% label: SWIM project folder (all other paths can be relative to this)
#% description: This folder should include all SWIM files/folders.
#% gisprompt: old,dir,dir
#% guisection: Settings
#%end

#%option
#% key: proname
#% type: string
#% required: no
#% key_desc: name
#% label: Name of project (keep it short)
#% description: Will try to find .bsn/.str/.cod files according to this name
#% guisection: Settings
#%end

#%option
#% key: simf
#% type: string
#% required: no
#% key_desc: path
#% label: Results file with all simulated subbasin discharge 
#% description: SubbasinIDs in the stations table should correspond to the columns
#% gisprompt: old,file,file
#% guisection: Settings
#%end

#%option
#% key: stationsvect
#% type: string
#% required: no
#% key_desc: raster
#% label: Stations point vector 
#% description: Should have a subbasinID column given in subbasinidcolumn or it will be uploaded if not given.
#% gisprompt: old,vector,vector
#% guisection: Settings
#%end

#%option
#% key: subbasinidcolumn
#% type: string
#% required: no
#% multiple: no
#% key_desc: name
#% label: SubbasinID column in stationsvect (corresponding to simf columns)
#% description: Leave blank to query subbasins and upload a subbasinID column to stationsvect
#% guisection: Settings
#%end

#%option
#% key: stationidcolumn
#% type: string
#% required: no
#% multiple: no
#% key_desc: name
#% label: Identifier column with stationIDs refereing to the columns in csvdata
#% description: Will be converted to 5 character string, if not given cats are used.
#% guisection: Settings
#%end

#%option
#% key: csvdata
#% type: string
#% required: no
#% multiple: no
#% key_desc: name
#% label: Observed discharge csv file with stationsvect idcolumn content as columns 
#% description: columns: (time [YYYY-MM-DD], stationID1, ..., stationIDn), (or pandas .pa DataFrame file)
#% gisprompt: old,file,file
#% guisection: Settings
#%end

#%option
#% key: subbasins
#% type: string
#% required: no
#% key_desc: raster
#% label: Subbasin raster
#% description: 
#% gisprompt: old,cell,raster
#% guisection: Settings
#%end

#%option
#% key: hydrotopes
#% type: string
#% required: no
#% key_desc: raster
#% label: Hydrotopes raster
#% description: 
#% gisprompt: old,cell,raster
#% guisection: Settings
#%end

#%option
#% key: llcmds
#% type: string
#% required: no
#% key_desc: specific
#% label: llsubmit commands
#% description: Will be written into the llsubmit jobfiles
#% guisection: Settings
#%end

import os, sys, fnmatch
import datetime as dt
import subprocess
import pandas as pa
import grass.script as grass
import pylab as pl
import numpy as np
grun = grass.run_command
gm = grass.message

try:
    import swimtxt as swim
except ImportError:
    print 'SWIMpy is not on your PYTHONPATH. Exiting...'
    sys.exit()

def loadDefs():
    if '--ui' in sys.argv: sys.argv.remove('--ui')
    # get names of module parameters
    options, flags = grass.parser()
    # get default values for the arguments
    defoptions = p.__dict__
    defoptions.update(p.readPar())
    
    defflags = []
    switches = dict(zip(components.values(),components.keys()))
    # update options
    for a in defoptions:
        # parameter in header and not already set
        if a in options and options[a]=='':
            options[a] = defoptions[a]
        # check flags
        if a in switches and (defoptions[a] != 0 or (a in flags and flags[a])):
            defflags += [switches[a]]
    
    # recreate commandline arguments
    sys.argv += ['-%s' %''.join(defflags)]
    sys.argv += ['%s=%s' %s for s in options.items()]
    sys.argv += ['--ui']
    return

def findFiles(filedict,findin='.',fatal=True):
    '''Find files that are given as values in filedict in the findin directory'''
    ffiles = filedict
    # dictionary of file identifiers with empty lists
    foundf = dict(zip(ffiles.keys(),[[] for i in range(len(ffiles))]))
    for root, dirnames, filenames in os.walk(findin):
        for f in ffiles:
            foundf[f]+=[os.path.join(root, fn) for fn in fnmatch.filter(filenames, ffiles[f])]
            
    # check if dublicated or none found
    for f in foundf.keys():
        if len(foundf[f]) != 1:
            if fatal: grass.fatal('Found no or more than one %s file: %r' %(ffiles[f],foundf[f]))
            else:
                grass.warning('Found no or more than one %s file: %r' %(ffiles[f],foundf[f]))
                foundf[f] = 'file not found'
        else:
            foundf[f] = foundf[f][0]
            
    return foundf

def writeTxtDB(path,fmt):
    # make nice header
    colnames, formats = [],[]
    for fm in fmt:
        length = max((len(fm[0]),len(fm[1]),abs(int(float(fm[1][:-1])))))
        colnames += [fm[0]+' '*(length-len(fm[0]))] # pad with space
        formats  += [' '*(length-len(fm[1]))+fm[1]] # pad with space
    # write files
    f=file(path,'w')
    f.write(' '.join(colnames)+'\n')
    f.write(' '.join(formats)+'\n')
    f.close()
    gm( 'Created DB in %s' %(path))
    return

def getStations(resourcedir='mySWIM',**datakwargs):
    ####### STATIONS ################################################    
    stationsinfo = grass.vector_info(options['stationsvect'])    
    # get stations table    
    try: stationstbl = grass.vector_db_select(options['stationsvect'])
    except: grass.fatal('Cant read the attribute table of %s, has it got one?' %options['stationsvect'])
    # make pandas df
    stationstbl = pa.DataFrame(stationstbl['values'].values(),columns=stationstbl['columns'])

    # check if it finds the stationidcolumn
    if len(options['stationidcolumn'])>0:
        if options['stationidcolumn'] not in stationstbl:
            grass.fatal('Cant find stationidcolumn %s' %options['stationidcolumn'])
        idcol = stationstbl[options['stationidcolumn']]
    elif 'stationID' in stationstbl:
        gm('Found stationID column in the subbasinsvect table.')
        idcol = stationstbl['stationID']
    else:
        ### stationID DESIGN #################################################
        gm('''Will use the order of stations in the stationvect as stationIDs and will assign these stationIDs:''')
        idcol = ['s%s' %s for s in range(1,len(stationstbl)+1)]
        gm(idcol)
        # make subbasinID column and upload them
        grun('v.db.addcolumn',map=options['stationsvect'],
                          columns='stationID varchar(5)',quiet=True)
        catcol = stationstbl.icol(1)
        for i,s in enumerate(catcol):
            grun('v.db.update',map=options['stationsvect'], column='stationID',
                 where='%s=%s' %(catcol.name,s), value=idcol[i],quiet=True)
        # make sure the csvdata file is read without header
        datakwargs.update({'names':['time']+idcol,'skiprows':1})
    # set as index
    stationstbl.index = idcol
    stationstbl.index.name='stationID'
    
    # check if it finds the subbasinidcolumn
    if len(options['subbasinidcolumn'])> 0:
        if options['subbasinidcolumn'] not in stationstbl:
            grass.fatal('Cant find subbasinidcolumn %s' %options['subbasinidcolumn'])
        # change name of column to subbasinID
        stationstbl['subbasinID'] = stationstbl[options['subbasinidcolumn']]
        stationstbl.drop(options['subbasinidcolumn'],axis=1,inplace=True)
        
    elif 'subbasinID' in stationstbl:
        gm('Found subbasinID column in the subbasinsvect table.')
        
    else:
        gm('Will upload subbasinIDs from the subbasins rast to the stationsvect table. Double check them and run the setup again if changed.')
        # first check if in same mapset then stationsvect
        if stationsinfo['mapset'] != grass.gisenv()['MAPSET']:
            grass.warning(' !!! Changing into mapset: %s !!!' %stationsinfo['mapset'])
            grun('g.mapset',mapset=stationsinfo['mapset'],quiet=True)
        grun('v.db.addcolumn',map=options['stationsvect'],
                          columns='subbasinID int',quiet=True)
        grun('v.what.rast',map=options['stationsvect'],
                          raster=options['subbasins'],column='subbasinID',quiet=True)
        subids = grass.vector_db_select(options['stationsvect'],columns='subbasinID')
        stationstbl['subbasinID'] = np.array(subids['values'].values()).flatten()
    # convert subbasinIDs to int
    for i,s in enumerate(stationstbl['subbasinID']):
        try: stationstbl['subbasinID'][i] = int(s)
        except: stationstbl['subbasinID'][i] = np.nan
        
    # report on found subbasinIDs
    gm('Found these subbasinIDs:\n%s' %stationstbl['subbasinID'])
    
    
    # read data
    gm('Attempting to read discharge data. I am expecting a header like this:')
    gm(','.join(['YYYY-MM-DD']+list(idcol)))
    data = pa.read_csv(options['csvdata'],parse_dates=0,index_col=0,**datakwargs)
    data.index.name = 'time'
    # check if read properlyq
    # check if all cats are in csvdata as columns
    for s in stationstbl.index:
        if s not in data:
            grass.warning('Cant find %s in the csvdata table header. This station wont have any data.' %s)
            gm(stationstbl.ix[s].to_string())
            data[s]=np.nan
    # get days per year covered
    peran = (~data.isnull()).astype(int).resample('a','sum')
    peran.index=peran.index.year
    # report
    gm('Overvations found (days per year):')
    for l in peran.to_string().split('\n'): gm(l)
#        # write individual files to resource folder
#        else:
#            path = os.path.join(resourcedir,s+'.pa')
#            sdat = pa.DataFrame({'Q':data[s]})
#            sdat.index.name = 'time'
#            sdat.to_pickle(path)
#            gm('Saved Q data for the station %s to %s in a python.pandas format.' %(s,path))
#            gm(sdat.head().to_string())
#            fpaths += [path]
#    # attach paths to stationstbl
#    stationstbl['data'] = fpaths
            
    # write out table to resourcedir
    path = os.path.join(resourcedir,'observations.csv')
    data.to_csv(path)
    # return as dictionary of stations
    stations = stationstbl.T.to_dict()
    return stations
        
def setupPro(resourcedir='mySWIM',parfile='mySWIM/myswim.py'):
    '''Set up all files needed for SWIM and mySWIM and write mySWIM folder and parameterfile'''
    # collect parameters in dict
    p = {}
    
    # check if needed params are set
    for e in ['proname','prodir','stationsvect','simf','subbasins','hydrotopes',
              'csvdata']:
        if len(options['proname'])<2: grass.fatal('%s must be set.' %e)
        p[e] = options[e]
    
    # change into prodir as all subsequent paths maybe relative to that
    os.chdir(p['prodir'])
        
    # find files in prodir
    ffiles = {'bsn':'.bsn','cod':'.cod','strf':'.str'}
    ffiles = dict([(e,options['proname']+ffiles[e]) for e in ffiles])
    foundf = findFiles(ffiles,findin=p['prodir'])
    # add to pro
    p.update(foundf)
    
    ### check for other swim files
    swimfiles = {'swim':'swim',
                 'soilcio':'soil.cio',
                 'filecio':'file.cio',
                 'figf':p['proname']+'.fig',
                 'runoff':'runoff.dat',
                 'simf':'rvaddQ_subbasin.prn',
                 'conf':'swim.conf',
                 'clim1':'clim1.dat',
                 'clim2':'clim2.dat',
                 'cntab':'cntab.dat',
                 'cropd':'crop.dat',
                 'subcatchf':'subcatch.def',
                 'subcatchbsn':'subcatch.bsn',
                 'wgen':'wgen.dat',
                 'wstor':'wstor.dat'}
    # simf
    if len(options['simf'])>0: swimfiles['simf']=os.path.basename(options['simf'])
    # check for files with warning if not exists
    swimfiles = findFiles(swimfiles,findin=p['prodir'],fatal=False)
    p.update(swimfiles) # only those that were found and not double    
    
    # check if mySWIM already exists and if not create dir
    resourcedirs = {'resourcedir':os.path.join(p['prodir'],resourcedir)}
    for n,d in [('Qdir','Q_files'),('clusterdir','lljobs'),('pestdir','pest')]:
        resourcedirs[n] = os.path.join(resourcedir,d)
    for n,d in resourcedirs.items():
        if os.path.exists(d):
            gm('Found the resource dir in: %s' %d)
        else:
            gm('Creating new resource dir in: %s' %d)
            os.makedirs(d)
    # attach to p
    p.update(resourcedirs)
            
    # get parameters from bsn file
    try: par,tmplt = swim.readBsnf(p['bsn'])
    except:
        grass.fatal('''Cant read the bsn file properly. Is it formatted like this:\n
        switch parameters
        value  param
        value  param
        ...
        _______________________________________________________________________
        basin, initialisation & calibration parameters
        param  param ... param description blabla
        value  value ... value
        ...
        _______________________________________________________________________
        CO2 EFFECT ON NET PHOTOSYNTHESIS (alpha) & TRANSPIRATION (beta)
        (ialpha,ibeta) = (1,0) OR (1,1) ONLY FOR SCENARIO PERIODS!
        ialpha    ibeta     C3C4crop  CO2-ref   CO2-scen
        0/1       0/1       3/4       346       406-436   OPTIONS & RANGES
        0         0         0         0         0
        _______________________________________________________________________''')
        
    gm('Found these parameters in the bsn file: %s' %','.join(par.keys()))
    gm('Those parameters will be saved in the future!')
    # make parameter save formats
    pfmts = [(e,{True:'14.3f',False:'7.4'}[par[e]>50]) for e in sorted(par.keys())]
    
    # check if resource files exists and if not create them
    rfiles = {'runsf': [('runID','04i'),('paramID','04i'),('climID','04i'),
                        ('runTIME','26s'),('station','5s'),('NSE','+7.3f'),
                        ('bias','+6.1f'),('start','10s'),('end','10s'),
                        ('purpose','-12s'),('notes','-64s')],
              'paramf': [('paramID','04i')]+pfmts,
              'climsf': [('climID','04i'), ('title','50s'),
                         ('start','10s'), ('end','10s'), ('notes','-64s')]}
    for f in rfiles:
        name = os.path.join(resourcedir,'%s.%s' %(p['proname'],f[:-1]))
        # exists?
        if os.path.exists(name):
            gm('%s file exists in: %s' %(f,name))
            p[f] = name
        else: # create
            writeTxtDB(name,rfiles[f])
            p[f] = name
    
    # get stations info and data
    p['stations'] = getStations(resourcedir=resourcedir)
    p['obsf']     = os.path.join(resourcedir,'observations.csv')
    
    # write all in p to a parameter file
    parfpath = os.path.join(resourcedir,'myswim.py')
    parf = file(parfpath,'w')
    parf.write('### mySWIM parameter file, saved on: %s\n' %dt.datetime.now())
    for e in sorted(p.keys()): parf.write('%s=%r\n' %(e,p[e]))
    
    gm('All set up! To proceed, uncheck the -0 flag and change into %s' %p['prodir'])
    return

def runinteractive(station):
    '''Interactive plotting while SWIM is running'''
    # open process and pipe it to python
    swimcmd = p.swim+' '+p.prodir+'/'
    ps = subprocess.Popen(swimcmd, shell=True, 
                      stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # initial simfile change time
    chtime = os.stat(p.simf).st_mtime
    # set up graph
    pl.ion()
    ax = pl.gca()
    ax.clear()
    Q = {'%s obs' %station: p.stations[station]['data']['Q'],
         '%s sim' %station: np.nan}
    Q = pa.DataFrame(Q, index=pa.DatetimeIndex(start=p.startdate,end=p.enddate,freq='D'))
    ax = Q.plot(style=['k--','r'],ax=ax)
    oline,sline = ax.lines
    pl.draw()
    
    while True:
        # print out model stdout
        nextline = ps.stdout.readline()
        if nextline == '' and ps.poll() != None:
            break
        gm(nextline)
        # check if simfile change time
        if os.stat(p.simf).st_mtime != chtime:
            # update graph
            Q['sim'] = swim.readSim(p.simf,simc=p.stations[station]['subbasinID'],
                               hotstart=False)['Q']
            sline.set_ydata(Q['sim'])
            pl.title(nextline)
            pl.draw()
    sys.stdout.flush()
    return
    
def checkInput():
    '''Check all user input when pre/pro/postprocessing, not when setting up'''
    # convert and empty options that are not needed
    ops={}    # local options
    for o in options:
        if options[o]!='':
            try: ops[o] = int(options[o])                      # is int
            except ValueError:
                try: ops[o] = float(options[o])                # is float
                except ValueError: ops[o] = options[o].strip() # is str
    # stations to save
    if 'savestations' in ops:
        if ops['savestations']=='all':
            ops['savestations'] = p.stations.keys()
        else:
            ops['savestations'] = ops['savestations'].split(',')
            for s in ops['savestations']:
                if s not in p.stations.keys(): grass.fatal('This is not a valid station ID: %s' %s)
            if len(ops['savestations'])==1: ops['savestations']=ops['savestations'][0]
    # setsubcatch
    if 'setsubcatch' in ops:
        dic = {}
        for s in ops['setsubcatch'].split(','):
            st_id=s.split(':')
            if len(st_id)!=2 or st_id[0] not in p.stations or st_id[1] not in p.readparams().index:
                grass.fatal('setsubcatch set incorrectly, need stationID:paramID pairs, got %s' %s)
            dic[st_id[0]]=st_id[1]
        ops=dic
    # peststation
    if flags['a']:    
        # check if peststation set properly
        if 'peststation' in ops and ops['peststation'] in p.stations:
            grass.fatal("Have you set peststation to a valid station ID?")
    # interactive
    if 'interactive' in ops:
        # check if peststation set properly
        if 'interactive' in ops and ops['interactive'] not in p.stations:
            grass.fatal("Have you set interactive to a valid station ID?")
        grass.warning('Interactive mode is experimental! The run statistics wont be saved!')
    
    # showruns
    if 'showruns' in ops:
        # split into list or if only one runID make list with one entry
        if type(ops['showruns'])==int: ops['showruns']=[ops['showruns']]
        else: ops['showruns'] = ops['showruns'].split(',')
        # check if valid stationIDs
        runids = p.readruns().index.tolist()
        for i,s in enumerate(ops['showruns']):
            try:
                s = int(s)
                ops['showruns'][i] = s
            except: continue
            # now check if either
            if s not in (p.stations.keys()+['all']):
                gm('%s is not a stationID' %s)
                if s not in runids:
                    grass.fatal('showruns seems to have an invalid runID or stationID: %s' %s)
    
    # decide what to parse to run/storeRuns
    # needed input        
    parse = {'stations':options['savestations']}
    # optional input
    oparse = ['notes','purpose']
    parse.update({n: options[n] for n in oparse if n in options})
    ops['parse'] = parse
    
    return ops
    
    
def preprocess():
    '''Check all input and write to swim files'''
                    
    ### COD FILE if needed #################################################
    cod = {}
    for i in ['iyr','nbyr']:
        if i in options: cod[i]=options[i]
    p.rwcod(**cod)
    
    ### parameters #########################################################
    if 'paramid' in options and options['paramid']!=0:
        params = options['paramid']
    else:
        params = options.copy()
        # translate from TRUE/FALSE flags to SWIM 1/0 params
        for f in components: params[components[f]] = int(flags[f])
    # write to .bsn file
    p.writePar(params)
    
    # update subcatch file if given
    if 'setsubcatch' in options:
        p.rwSubcatchPars(options['setsubcatch'])

    ### PEST ###############################################################
    if flags['a']:    
        # set up pest, might through exceptions
        options['pestname'] = p.setupPest('%s_pest' %options['peststation'],
                                          stations=options['peststation'])

    return
    
def process():
    '''Run swim or pest'''
    # set buffer environment variables to write output every line
    os.environ['GFORTRAN_UNBUFFERED_ALL'] = '1'
    
    if flags['a']:    
        # PEST
        # submit to cluster        
        ll = {}
        if flags['l']: ll='m.swim.run_pest'
        p.runPest(options['pestname'],**ll)
        
    elif 'interactive' in options:
        runinteractive(options['interactive'])
    else:
        if flags['l']: options['parse']['ll']='m.swim.run'
        # run        
        p.last = p.run(**options['parse'])
    
    return

def postprocess():
    '''Postprocessing stuff'''
    
    # store runs if postprocessing only
    if flags['p'] and flags['v']:
        stored = p.storeRuns(**options['parse'])
        
    ### PLOTTING Q ##########################################################
    if 'showruns' in options:
        # collect all runs in dictionary of stationIDs
        st = p.stations.keys()
        # empty stations dictionary
        runs = dict(zip(st,[[] for i in range(len(st))]))
        for s in options['showruns']:  # s is either all, any stationID or runID
            if s=='all':
                r = p[p.stations.keys()]
                for sta in r: runs[r[sta]['station']] += [r[sta]]
            else:
                r=p[s]
                runs[r['station']] += [r]
        # plot
        for s in runs:
            runslist = runs[s]
            # swim.run.plot() function
            if len(runslist)>0:
                junk = runslist[0].plot(runslist[1:])
            del runslist
        pl.show()
            
    return
    
    
if __name__ == "__main__":

#    # check if prodir to change into is given
#    for a in sys.argv:
#        if a.startswith('prodir='):
#            os.chdir(a.split('=')[1])

    # try to start up SWIMpy
    if '-0' not in sys.argv:
        try:
            p = swim.pro()
        except:
            grass.warning('''\
            No SWIM project found, make sure to be in the SWIM project folder or \
            set up a SWIM project in the Settings section and setting the -0 flag.''')
            # add -0 setup flag
            sys.argv += ['-0','--ui']
        
    # check if project defaults can/must be loaded
    if 'p' in vars() and ('--ui' in sys.argv or len(sys.argv)<2):
        loadDefs()
        
    # normal start
    options, flags = grass.parser()
    
    # setup project
    if flags['0']:
        setupPro()
    else:
        # Input checking
        options = checkInput()
    
        # report
        grass.message(('GIS Environment:',grass.gisenv()))
        grass.message(('Parameters:',options,flags))
        
        if not flags['p']:
            preprocess()
        
        if not flags['p'] and not flags['u']:
            process()
            
        if not flags['u']:
            postprocess()
        
    

# run.py -0 prodir=/data/sumario/wortmann/Tarim/testing proname=tari \
#             simf=/data/sumario/wortmann/Tarim/testing/Output/Res/rvaddQ_subbasin.prn \
#             stationsvect=Tarim_stations@PERMANENT stationidcolumn=id \
#             csvdata=/data/sumario/wortmann/Tarim/testing/mySWIM/observations.csv \
#             subbasins=subbasins@subbasins_new hydrotopes=hydrotopes@hydrotopes
