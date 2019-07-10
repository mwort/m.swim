Bugs
====



Improvements
============
## m.swim.subbasins
- add lon/lat to subbasin table (centroids)
- tidy clean_subbasins method, maybe move some to postprocess_subbasins
- convert accumulation cells to km2 and remove all kmtocell/celltokm conversions,
  possible with routing?
- postprocess mainstreams
  ```
  # add vertices a bit short of the DEM resolution
  v.split input=mainstreams@subbasins output=mainstreams__split length=40
  # avoid unnessary nodes
  v.build.polylines input=mainstreams__split output=mainstreams__polylines type=line
  # add nodes at subbasin intersections and adds subbasinID to table
  v.overlay ainput=mainstreams__polylines binput=subbasins@subbasins output=mainstreams__subbasins op=and --o
  # clean table and columns
  g.copy vect=mainstreams__subbasins,mainstreams
  v.db.droptable mainstreams layer=1 -f
  v.db.addtable mainstreams
  # get subbasinID column from old vector
  v.db.join mainstreams column=cat other_table=mainstreams__subbasins other_column=cat subset='b_subbasinID'
  v.db.renamecolumn mainstreams column=b_subbasinID,subbasinID
  ```

## m.swim.substats
- add lon, lat of subbasin centroids to table for stat-outdat, possibly add
  reference elevation input existing column
```
v.to.points input=subbasins@subbasins type=centroid output=subbasin_centroids
g.mapset map=PERMANENT loc=Danube_latlon
v.proj input=subbasin_centroids map=subbasins loc=Danube_LAEA
v.db.addcolumn subbasin_centroids columns='lon double,lat double'
v.to.db subbasin_centroids option=coor columns=lon,lat
g.mapset map=subbasins loc=Danube_LAEA
v.proj input=subbasin_centroids map=PERMANENT loc=Danube_latlon
v.db.join map=subbasins column=subbasinID \
          other_table=subbasin_centroids_subbasin_centroids_1 \
          other_column=subbasinID subset_columns=lon,lat
v.db.select map=subbasins columns=subbasinID,lon,lat,average_elevation \
            separator=space file=input/stat-outdat.prn --o
```
- provide additional files, eg. subcatch, bsn etc.

## m.swim.hydrotopes
- soils need to be in sequential order
