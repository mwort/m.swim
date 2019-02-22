#!/usr/bin/env python
#
############################################################################
#
# MODULE:      m.swim.routing v1.2
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
#% description: Soil and Water Integrated Model (SWIM) preprocessor: routing
#%End
#%Option
#% guisection: Required
#% key: subbasins
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% description: Subbasin vector, nextID and inletID will be uploaded to table
#% gisprompt: old,vector,vector
#%end
#%Option
#% guisection: Required
#% key: accumulation
#% type: string
#% required: yes
#% multiple: no
#% key_desc: name
#% description: CELL/integer accumulation raster, e.g. m.swim.subbasins
#% gisprompt: old,cell,raster
#% answer: accumulation
#%end

#%Flag
#% guisection: Output
#% key: r
#% label: Only remake routing network and .fig file (eg. after manual editing, outlets and inlets needed!)
#%end

#%Option
#% guisection: Output
#% key: routingnet
#% type: string
#% required: no
#% multiple: no
#% key_desc: name
#% description: Vector of routing network to be created
#% gisprompt: new,vector,vector
#% answer: routingnetwork
#%end

#%Option
#% guisection: Output
#% key: mainstreams
#% type: string
#% required: no
#% multiple: no
#% key_desc: name
#% description: Main streams vector and raster to be created (routingnet needed)
#% gisprompt: new,vector,vector
#% answer: mainstreams
#%end

#%Option
#% guisection: Output
#% key: outlets
#% type: string
#% required: no
#% multiple: no
#% key_desc: name
#% description: Outlets point vector
#% gisprompt: new,vector,vector
#% answer: subbasinoutlets
#%end

#%Option
#% guisection: Output
#% key: inlets
#% type: string
#% required: no
#% multiple: no
#% key_desc: name
#% description: Inlets point vector
#% gisprompt: new,vector,vector
#% answer: subbasininlets
#%end

#%Option
#% guisection: SWIM files
#% key: figpath
#% type: string
#% required: no
#% multiple: no
#% key_desc: path/myproj.fig
#% description: path to the .fig file
#% gisprompt: new,file,file
#%end

#%Option
#% guisection: Optional
#% key: fromto
#% type: string
#% required: no
#% multiple: yes
#% key_desc: from,to
#% description: Force routing manually (pairs of subbasinID,nextID, is slow)
#%end

#%Option
#% guisection: Optional
#% key: streams
#% type: string
#% required: no
#% multiple: no
#% key_desc: vect
#% description: Predefined streams to include in mainstreams output
#% gisprompt: old,vector,vector
#%end

#%Option
#% guisection: Optional
#% key: minmainstreams
#% type: double
#% required: no
#% multiple: no
#% label: Minimal drainage area headwater mainstreams, km2.
#% description: All headwater subbasins below this area have a mainstream of 10% of their drainage area.
#% answer: 50
#%end

#%Option
#% guisection: Optional
#% key: rivercourse
#% type: string
#% required: no
#% multiple: yes
#% key_desc: subbasinID
#% label: SubbasinID(s) to calculate river course for (to column: 'course_$id')
#% description: Uploads cumulative mainChannelLength for all subbasins in column: course_subbasinID
#%end

#%Flag
#% guisection: Optional
#% key: k
#% label: Keep intermediate files (those named *__*)
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
            if optionsandflags[o]!='': self.options[o] = optionsandflags[o]
        self.__dict__.update(self.options)
        #TODO: check if accumulation map is an int and has minus values within subbasin MASK

        # make rast from subbasin vector
        self.subbasinrast = 'subbasin__rast'
        self.nextsubb_col = 'nextID'
        self.subb_col     = 'subbasinID'

        # make sure subbasins are in current mapset
        VectMapset=self.subbasins.split('@')
        if len(VectMapset)==2 and VectMapset[1]!=grass.gisenv()['MAPSET']:
            grass.fatal('!!! %s needs to be in the current mapset because I will update its table !!!' %(self.subbasins))

        # check what columns to use
        cols = grass.vector_columns(self.subbasins)
        if 'subbasinID' not in cols and 'cat' in cols:
            grun('v.db.addcolumn',map=self.subbasins,columns=self.subb_col+' int',quiet=True)
            grun('v.db.update', map=self.subbasins, column=self.subb_col, qcol='cat')

        # make rast
        grun('v.to.rast',input=self.subbasins,output=self.subbasinrast,
                 use='attr', attrcolumn=self.subb_col,overwrite=True, quiet=True)
        # check fromto and prepare
        if 'fromto' in self.options:
            try:
                # make array
                a = np.array(self.fromto.split(','),dtype=int).reshape(-1,2)
                # add a 0 inlet column to it
                a = np.column_stack((a,np.zeros(len(a))))
                # bring into rec array
                self.fromto = np.array(zip(*a.T),dtype=[('subbasinID',int),
                                       ('nextID',int),('inletID',int)])
            except:
                grass.fatal('No integers or uneven number of values in fromto pairs. %s' %self.fromto)

        # check river course
        if 'rivercourse' in self.options:
            if type(self.rivercourse)==str:
                self.rivercourse = map(int,self.rivercourse.split(','))
            elif type(self.rivercourse)==int:
                self.rivercourse = [self.rivercourse]

        return

    def routing(self, searchradius=1.5, overwrite=True, quiet=True):
        '''Calculate the routing for the subbasins given as vector. A column will
        be add to the subbasin table called nextSubbasinID and overwrites the following:
        a vector called subbasinoutlets
        a vector called streams

        Search radius (in cells) determines how many cells around the outlet should
        be checked for inlet, should be greater than 1. The larger, the more likely
        are subbasins routed to the subsequent subbasin downstream.

        If vector=True output a vector with a nextID as a column in its table.
        If False output nextID as a raster.
        '''
        # some parameters
        kw = {'overwrite':overwrite, 'quiet':quiet}
        grun('r.mask',raster=self.subbasinrast,overwrite=True)

        # make outlet raster points, border subbasins with neg accumul. not correct
        grass.message('Searching outlets...')
        grun('r.stats.zonal', base=self.subbasinrast,cover=self.accumulation,
             method='max', output='maxaccum__', flags='r', **kw)
        exp = "outlets__=if('%s' == int(@maxaccum__),%s,null())" %(self.accumulation,self.subbasinrast)
        grass.mapcalc(exp, **kw)

        grass.message('Growing outlets...')
        # join outlets (set all to 1) and then grow outlets
        grass.mapcalc("outlets1__=if(isnull(outlets__),null(),1)",overwrite=True)
        # grow to all 8 neighbours (1.5 times cell size)
        grun('r.grow', input='outlets1__', output='outlets__grown', radius=searchradius,**kw)
        grun('r.clump', input='outlets__grown', output='outlets__grown__clumped',**kw)

        # make inlets
        grass.message('Searching inlets...')
        grun('r.stats.zonal', base='outlets__grown__clumped', method='max',
             cover=self.accumulation, output='maxoutlet__clumps', flags='r', **kw)
        grass.mapcalc("inlets__=if(%s == int(@maxoutlet__clumps),%s,null())" %(self.accumulation,self.subbasinrast), **kw)

        # transfer inlet subbasinID to clumps
        grun('r.stats.zonal', base='outlets__grown__clumped', cover='inlets__',
             method='max', output='clumps__subbasinID', **kw)

        # make outlets vector with nice columns
        grun('r.to.vect', input='outlets__', output=self.outlets,type='point',
             flags='v',**kw)
        grun('v.db.addcolumn', map=self.outlets, columns='subbasinID int', quiet=quiet)
        grun('v.db.update', map=self.outlets, column='subbasinID', qcol='cat', quiet=quiet)
        grun('v.db.dropcolumn',map=self.outlets, column='label',**kw)

        # make inlets point vector
        gm('Getting inletIDs...')
        grun('r.to.vect', input='inlets__', output=self.inlets, type='point',**kw)
        grun('v.db.renamecolumn', map=self.inlets, column='value,subbasinID',quiet=quiet)
        grun('v.db.dropcolumn', map=self.inlets, column='label',quiet=quiet)
        grun('v.db.addcolumn',map=self.inlets,columns='inletID int',**kw)
        grun('v.what.rast',map=self.inlets,raster='outlets__grown__clumped',
             column='inletID',quiet=quiet)

        # load inlet ID and subbasinID into outlets vector
        grun('v.db.addcolumn', map=self.outlets, columns='%s int,inletID int' %self.nextsubb_col, **kw)
        grun('v.what.rast', map=self.outlets, raster='clumps__subbasinID',
             column=self.nextsubb_col, quiet=quiet)
        grun('v.what.rast', map=self.outlets, raster='outlets__grown__clumped',
             column='inletID', quiet=quiet)

        # upload nextID into subbasin vector either presuming cats are subbasinIDs or the column exists
        grun('v.db.join',map=self.subbasins,column=self.subb_col,
             otable=self.outlets,ocolumn='subbasinID',scolumns=self.nextsubb_col+',inletID',
             quiet=True)
        # check outlets
        self.checkOutletAndFromto()

        gm('Routing successfully completed. nextID and inletID added to %s table.' %self.subbasins)

        return

    def checkOutletAndFromto(self):
        '''Check number of outlets and change to negative values in outlets and subbasins table'''
        # check basin outlets
        grass.message('''Find outlets (if other than ID 1, check if the accumulation
                     map has negative, ie. off-map flow)...''')
        # get fromto columns
        try: sboutin = readSubNxtID(self.subbasins)
        except ValueError: grass.fatal('Cant convert the subbasinID, nextID or inletID columns of %s to integers' %self.subbasins)
        outlets = sboutin[sboutin['subbasinID']==sboutin['nextID']]
        outlets['nextID'] = np.negative(outlets['nextID'])
        outlets['inletID']= 0
        change = outlets.copy()
        # manual change using fromto
        if 'fromto' in self.options:
            # append to change array
            change=np.append(change,self.fromto)
        # update subbasins and outlets table
        update = lambda sbb,map,col,val: grun('v.db.update',map=map,column=col,
                  where='subbasinID=%s' %sbb, value=val)
        # change for each row in change, in both vectors, both columns
        for s in change:
            for m in [self.subbasins,self.outlets]:
                for c in ['nextID','inletID']:
                    update(s['subbasinID'],m,c,s[c])

        # check outlets again for reporting
        sboutin = readSubNxtID(self.subbasins)
        outlets = sboutin[sboutin['subbasinID']==sboutin['nextID']]
        outlets = np.append(outlets,sboutin[sboutin['nextID']<=0])
        grass.message('Subbasin(s) %s is(are) outlet(s)' %outlets['subbasinID'])
        return

    def buildRoutingNet(self):
        '''Connect with subbasin centroids when subbasincentroids set to subbasins
        vector (x,y will be uploaded to table)'''
        # get inlets and outlets information
        oinfo = vreport(self.outlets,index='subbasinID')
        iinfo  = vreport(self.inlets,index='inletID')

        line = lambda fromxy,toxy: '%s,%s\n%s,%s\nNaN,NaN\n' %(fromxy['x'],fromxy['y'],
                                                               toxy['x'],toxy['y'])
        # loop over outlets
        lines=[]
        manual=[]
        for oID in oinfo:
            # inlet ID of outlet or downstream subbasin
            sbinletID=int(oinfo[oID]['inletID'])
            # for manually changed subbasins or outlets
            if sbinletID==0:
                if int(oinfo[oID]['nextID'])>0: manual+=[oinfo[oID]]
                continue
            # from outlets to inlets
            lines+=[line(oinfo[oID],iinfo[sbinletID])]
        gm('Connected %s outlets' %len(oinfo))

        # connect centroids with outlets and inlets with centroids if -c flag set
        # get subbasin centroids
        sbinfo = vreport(self.subbasins,index='subbasinID')
        # from centroids to outlets
        for o in oinfo:
            lines+=[line(sbinfo[o],oinfo[o])]
        # from inlets to centroids
        for i in iinfo:
            # get subbasinId of inlet
            insb = int(iinfo[i]['subbasinID'])
            lines+=[line(iinfo[i],sbinfo[insb])]

        # connect outlets with centroids of nextID for all manually changed subbasins
        for sb in manual:
            lines+=line(oinfo[int(sb['subbasinID'])],sbinfo[int(sb['nextID'])])
        # write to tmpfile and read as lines
        tf=grass.tempfile()
        f=file(tf,'w')
        f.writelines(lines)
        f.close()
        # make line vector
        grun('v.in.lines',input=tf,output=self.routingnet,separator=',',quiet=True)
        return


    def mkstreams(self):
        '''Create minimal stream network reaching all subbasins and with nice main streams'''
        # get max accumulation and cell count for each subbasin
        rsmaxaccum = gread('r.stats', input='maxaccum__', flags='lcn')
        rsmaxaccum = np.array(rsmaxaccum.split(), dtype=float).reshape((-1, 3)).astype(int)
        subbasinIDs = rsmaxaccum[:, 0]
        maxaccum = rsmaxaccum[:, 1]
        # cellcounts must not be larger than maxaccum (may happen in subbasins with more than 1 outlet)
        cellcounts = np.min([rsmaxaccum[:, 2], maxaccum], axis=0)
        # calculate optima accumulation for nice headwater mainstreams
        accr = grass.parse_command('r.info', map=self.accumulation, flags='g')
        minaccum = np.int32(round(float(self.minmainstreams)*1e6 /
                            (float(accr['nsres'])*float(accr['ewres']))))
        minaccum = np.ones_like(maxaccum) * minaccum
        optiaccum   = np.min([maxaccum-np.int32(cellcounts*0.1), minaccum], axis=0)
        # check incoming subbasins maxaccum and take the smallest accumulation to update optiaccum
        subnext     = readSubNxtID(self.outlets)
        optiaccum   = dict(zip(subbasinIDs,optiaccum))
        maxaccum    = dict(zip(subbasinIDs,maxaccum))
        for sn in np.unique(subnext['nextID']):
            if sn < 0: continue
            for sb in subnext[subnext['nextID']==sn]['subbasinID']:
                optiaccum[sn] = min(optiaccum[sn],maxaccum[sb]-1)
        # make raster
        tempf = grass.tempfile()
        np.savetxt(tempf,optiaccum.items(),fmt='%i=%i')
        grass.run_command('r.reclass',input='maxaccum__',output='optiaccum__',
                          rules=tempf,quiet=True)
        # get accumulation and make lines
        grass.mapcalc("{0}__unthin=if({1} > {2},{1},null())".format(self.mainstreams,
                      self.accumulation,'optiaccum__'),overwrite=True)

        # make sure all pixels between outlets and inlets are included
        #grun('v.to.rast',input=self.outletinletlines,output='routingnet__rast',
        #     use='val',type='line',quiet=True)
        #grass.mapcalc("{0}=if(isnull({0})&~isnull({1}),{2},{0})".format(self.mainstreams,
        #              'routingnet__rast',self.accumulation),overwrite=True)
        grun('r.thin', input=self.mainstreams+'__unthin', output=self.mainstreams, overwrite=True, quiet=True) # may exclude outlet/inlet points

        # use predefined streams for subbasins that have them
        if 'streams' in self.options:
            grun('v.to.rast',input=self.streams,output='stream__rast',
             use='val',type='line',val=1,quiet=True)
            grun('r.thin', input='stream__rast', output='stream__rast__thin',
                 overwrite=True, quiet=True)
            grun('r.stats.zonal',base=self.subbasinrast,cover='stream__rast__thin',
                 method='sum',output='n__streams__subbasins',quiet=True,overwrite=True)
            # get subbasins that have at least x cells of streams
            grass.mapcalc('{0}=if(n__streams__subbasins>10,{1},{0})'.format(self.mainstreams,'stream__rast__thin'),
                            overwrite=True)
        # make final vector
        grun('r.to.vect', input=self.mainstreams, output=self.mainstreams,
             type='line',quiet=True)
        return

    def fig_file(self):
        '''Write the .fig file needed for SWIM, subbasin vect needs to have a
        subbasinID column and a nextID column'''
        # get subbasinID and nextID or fromto array
        fromto = readSubNxtID(self.subbasins)
        # sort according to next and then subbasin ID
        fromto = np.sort(fromto, order=('subbasinID',))

        # lines list to be appended
        lines = []

        # write subbasin section
        for sbcat in fromto:
            sb = sbcat['subbasinID']
            lines += [['subbasin', 1, sb, sb, sb, '']]

        # ADD and ROUTE
        # calculate stream order for each subbasin
        order = subbasinorder(fromto)
        # upload to subbasin table
        otf = grass.tempfile()
        # make 1 indexed
        oar = np.array([order[s] for s in fromto['subbasinID']]) + 1
        np.savetxt(otf, np.column_stack([fromto['subbasinID'], oar]),
                   fmt='%i=%i')
        grass.run_command('r.reclass', input=self.subbasinrast,
                          output='sb__order', rules=otf, quiet=True)
        grun('v.db.addcolumn', map=self.subbasins,
             columns='strahler_order int')
        grun('v.what.rast', map=self.subbasins, raster='sb__order',
             type='centroid', column='strahler_order', quiet=True)
        # get order of nextID subbasin
        downstorder = fromto.copy()
        for i, sb in enumerate(fromto['nextID']):
            if sb > 0:
                downstorder['nextID'][i] = order[sb]
            else:
                # outlet order
                downstorder['nextID'][i] = max(order.values())

        # next storage location, i.e. non subbasin
        sID = fromto['subbasinID'].max()+1
        # loop over subbasin orders and get addroute lines
        maxorder = downstorder['nextID'].max()
        downstorderNsb = []
        for o in range(1, maxorder+1):
            subs = fromto[downstorder['nextID'] == o]  # nextID is next order
            # check if negative nextIDs/outlet is in there
            subs = subs[subs['nextID'] > 0]
            if len(subs) == 0:
                continue  # in case of the real outlet
            downstorderNsb += [(o, len(subs))]
            # get the addrout lines
            ls, routes = addroute(sID, subs)
            lines += ls
            # replace subbasinIDs by storageIDs from routes in fromto array
            for r in routes:
                ix = fromto['subbasinID'] == r['subbasinID']
                fromto['subbasinID'][ix] = r['storageID']
            # next storageID
            sID = routes['storageID'].max() + 1
        # add finish
        lines += [['finish', 0, '', '', '', '']]

        # write via numpy
        np.savetxt(self.figpath, np.array(lines), fmt='%-15s%1s'+'%6s'*4)
        gm('Wrote %s' % self.figpath)
        return {'order': order, 'downstorder': downstorder,
                'downstorderNsubbasins': downstorderNsb}

    def getCourse(self,headsb):
        '''Create river course from the headwater subbasin headsb to the outlet,
        that is reached when nextID<0.
        Subbasin vector needs to have subbasinID,nextID,mainChannelLength columns

        Uploads cumulative river lengths for the subbasins of the river course.
        '''
        gm('Calculating course for subbasinID %s...' %headsb)
        # get subbasinIDs,nextIDs,mainChannelLength
        subbasins = grass.vector_db_select(self.subbasins,
                    columns='subbasinID,nextID')['values'].values()
        nextids   = {int(s[0]):int(s[1]) for s in subbasins}

        # make course column
        ccol = 'course_%s' %headsb
        grun('v.db.addcolumn',map=self.subbasins,columns='%s double' %ccol)

        # find river course
        riverlength = 0
        riversb     = []
        sb = headsb
        while sb>0:
            riversb     += [sb]
            riverlength += 1
            grun('v.db.update',map=self.subbasins,column=ccol, value=riverlength,
                 where='subbasinID=%s' %sb)
            sb = nextids[sb]
        # report
        grass.message('''
        Uploaded cumulative river length from the %s to the outlet to the
        subbasin table in column %s (number of subbasins: %s). Check subbasin
        table and sort by %s
        To extract the subbasins use:
        v.extract map=%s where="%s!=''"
        ''' %(headsb,ccol,riverlength,ccol,self.subbasins,ccol))
        return 0

def vreport(vect,index):
    '''Return the v.report table with coordinates as a dictionary with subbasinID
    as index and column names in the entries for a point vector'''
    t=gread('v.report', map=vect, option='coor').split()
    t=[tuple(l.split('|')) for l in t]
    t=np.array(t[1:],dtype=[(i,object) for i in t[0]])
    lookup = {}
    for sb in t: lookup[int(sb[index])]=sb
    return lookup

def readSubNxtID(subbasinsvect,columns=('subbasinID','nextID','inletID')):
    '''Vector needs subbasinID, nextID and inletID column'''
    tbl=grass.vector_db_select(subbasinsvect,columns=','.join(columns))['values'].values()
    # check if empty cells
    tbl=np.array(tbl,dtype=np.unicode)
    for i,c in enumerate(columns):
        empty=tbl[tbl[:,i]==u'',i]
        if len(empty)>0: grass.fatal('The table %s has %s null values in column %s' %(subbasinsvect,len(empty),c))
    # convert to numpy rec array
    t = np.array(zip(*tbl.T),dtype=zip(columns,(int,)*len(columns)))
    return t


def subbasinorder(fromto):
    '''Calculate maximum number of upstream subbasin for each subbasin and
    return as dictionary.'''
    order = {}
    # find headwater sb without inflows to start from
    noinlets = np.array([i not in fromto['nextID'] for i in fromto['subbasinID']])
    hw = fromto[noinlets]
    orderX = hw
    i = 0
    # as long as there is one routed that isnt the largest outlet
    while len(orderX)>=1 and any(orderX['nextID']>0):
        # loop over each subbasin and correct order in order dict
        #sys.stdout.write('\rCalculating stream order %s' %i)
        for sb in orderX['subbasinID']:
            if sb in order:
                order[sb] = max([order[sb],i])
            else:
                order[sb] = i
        # get downstream subbasin info
        ins = np.unique(orderX['nextID']) # all downstream ids
        # get their info
        orderX = fromto[np.array([s in ins for s in fromto['subbasinID']])]
        # increase order
        i += 1
    #sys.stdout.write('\n')
    # assign remaining unique inlet, i.e. outlet the next higher order and for
    # its nextID (the same negative) another index higher
    order[int(orderX['subbasinID'])]=i
    order[-int(orderX['subbasinID'])]=i+1

    gm('Outlet subbasin order: %s' %i)
    return order

def addroute(sID,fromto):
    '''Write the add and route commands for the subbasinID and nextID columns
    given in fromto continuing from the storage location given in sID'''
    lines = []
    routesIDs = {}
    # write add and routes
    for inlet in np.unique(fromto['nextID']):
        # get subbasin that inflow into this sb
        upsb = fromto['subbasinID'][fromto['nextID']==inlet]
        # add all subbasins together if more than one
        if len(upsb)>1:
            for i,sb in enumerate(upsb[1:]):
                if i==0: lastsID=upsb[0]
                ### ADD
                lines += [['add', 5, sID, lastsID, sb, sb]] # 5 is the swim code for add
                #print 'add', sb, 'to', lastsID, '=', sID
                lastsID = sID
                sID += 1
        else:
            lastsID = upsb[0]
        ### ROUTE
        lines += [['route', 2, sID, inlet, lastsID, inlet]] # 2 is for route
        #print 'route',lastsID, 'through', inlet, '=',sID
        # add inlet station to routed storage location
        lastsID = sID
        sID += 1
        lines += [['add', 5, sID, lastsID, inlet, inlet]]
        #print 'add inlet',inlet, 'to routed', lastsID
        # save sID for further routing/adding
        routesIDs[inlet] = sID
        # increase storage location counter for next inlet
        sID += 1
    # format routes nicely
    routes = np.array(zip(routesIDs.keys(),routesIDs.values()),
                      dtype=[('subbasinID',int),('storageID',int)])
    return (lines,routes)


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

    # calculate routing
    if not main.r:
        grass.message('Will calculate routing for %s' %main.subbasins)
        main.routing()

    # if r set, check outlets and fromto again
    if main.r: main.checkOutletAndFromto()

    # calculate routing network and mainstreams if set
    if 'routingnet' in main.options:
        grass.message('Building routing vector...')
        main.buildRoutingNet()

    # routing network needs to be set to calculate mainstreams
    if not main.r and 'mainstreams' in main.options and 'routingnet' in main.options:
        grass.message('Creating mainstream network in %s...' %main.mainstreams)
        main.mkstreams()

    ### .fig file
    if 'figpath' in main.options:
        grass.message('Will write .fig file to %s!' % main.figpath)
        main.fig_file()

    # rivercouse
    if 'rivercourse' in main.options:
        for s in main.rivercourse: main.getCourse(s)

    # clean
    if not main.k:
        grass.run_command('g.remove',type='raster,vector', pattern='*__*',flags='fb',quiet=True)
        # remove mask
        grun('r.mask',flags='r',quiet=True)

    # report time it took
    delta = dt.datetime.now()-st
    grass.message('Execution took %s hh:mm:ss' %delta)
