# The GRASS7 Soil and Water Integrated Model (SWIM) preprocessor

Version v1.1

## Version history

*Version v1.1* - 2017-01-27:
Some adjustments made to align output with current SWIM version (but not yet fully synchronised) including a small test suite and compatibility with g.extension (GRASS 7.2) for convenient installation.

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
make
```
Or test modules individually like:
```
make subbasins
make routing
```
The non-grass output of tests is verified by a sha1sum comparison.
Hash and paths are stored and committed to git in the file output.sha1 .
To verify that the output is the same run:
```
make checkoutput
```
Checking the diff of the output.sha1 file will indicate output files that have changed (if any).
