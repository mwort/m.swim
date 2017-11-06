Bugs
====



Improvements
============
## m.swim.subbasins
- add lon/lat to subbasin table (centroids)
- tidy clean_subbasins method, maybe move some to postprocess_subbasins

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
