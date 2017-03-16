Bugs
====



Improvements
============
## m.swim.subbasins
- loop over input station cats, rather than just range(len(stations))
- output a snapped stations vector with terminating
  subbasinID, downstream stationID  and other station related infos, synchronised with catchments table
```
v.distance from=stations_germany_danube to=streams -p upload=to_x,to_y,dist | v.in.ascii output=stations_snapped skip=1 columns='stationID int,x double,y double,distance double' cat=1 x=2 y=3 input=-
v.what.rast map=stations_snapped raster=subbasins column=subbasinID
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
