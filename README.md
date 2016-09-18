# The GRASS7 Soil and Water Integrated Model (SWIM) preprocessor

Version v1.0

## Version history

*Version v1.0*
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
