# The GRASS7 Soil and Water Integrated Model (SWIM) preprocessor

Version v1.1

## Version history

*Version v1.1* - 2017-01-27:
Some adjustments made to align output with current SWIM version (but not yet fully synchronised) and compatibility with g.extension (GRASS 7.2) for convenient installation.

*Version v1.0* - 2016-09-18:
First version develped between 2012-2016 producing a fully functional SWIM project but not sychronised with the development of the SWIM code.


## Documentation

The documentation is still located here:
http://www.pik-potsdam.de/~wortmann/m.swim/doc


## Testing

There are tests for each grass module based on the Blankenstein catchment (upper Saale (Elbe). All required input maps are in the grassdb/utm32n/PERMANENT mapset, all tests should therefore be executed in:
```
cd test
grass70 grassdb/utm32n/PERMANENT
```
Run all tests:
```
./test_all.grass
```
Or test modules individually like:
```
./test_subbasins.grass
./test_routing.grass
```

## TODO
- adapt columns in structure file to include management, wetland, glaciers and elevations (not only when -d)
- adapt header of file.cio
- allow writeout of 3 subfiles only in m.swim.substats
- soils need to be in sequential order
- add lon/lat to subbasin table (centroids)
- add simple grass output for additional files to tests
