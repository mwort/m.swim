# The GRASS7 Soil and Water Integrated Model (SWIM) preprocessor

Version v1.2

## Documentation

The documentation is located here:
https://github.com/mwort/m.swim/wiki

## Version history

*upcomming*
- `streamthresh` is now a required parameter to `m.swim.subbasins`
- fixed hydrotopes counting starting from 1
- predefined subbasins will only be included if they are within the catchment

*Version v1.2* - 2019-01-09
- test project in sync with SWIM develop branch
- output `stations_snapped` vector in `m.swim.subbasins` with topology
  catchment info
- docs migrated to markdown and as Github Wiki
- `minmainstreams` argument to `m.swim.routing`
- add subbasin order in `m.swim.routing` to subbasin table
- `m.swim.substats` outputs just three files instead of three for each subbasin
- `rwatershedmemory` argument to `m.swim.subbasins`
- `m.swim.subbasins` output catchment names get station category postfix
  (instead of running counter)
- numerous bug fixes and refactoring

*Version v1.1* - 2017-01-27
- Some adjustments made to align output with current SWIM version (but not yet fully synchronised) including a small test suite and compatibility with g.extension (GRASS 7.2) for convenient installation.

*Version v1.0* - 2016-09-18
- First version develped between 2012-2016 producing a fully functional SWIM project but not sychronised with the development of the SWIM code.


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


## Releasing
- change version in `README.md` and in header of all module files.
- add change log in `README.md`
- `git commit -a -m 'Bump version to vX.X.'`
- `git tag vX.X`
- `git push & git push --tags`
