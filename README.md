# The GRASS7 Soil and Water Integrated Model (SWIM) preprocessor

Version v1.3

## Documentation

The documentation is located here:
https://github.com/mwort/m.swim/wiki

## Version history

*upcoming release*
- `m.swim.routing`: topologically correct mainstreams vector output
  - downstream directions, nodes at intersections and subbasin boundaries
  - subbasinID and accumulation (min, mean, max) columns in table
  - network nodes in layer 2 (for network analysis)
  - optional `minmainstreams` parameter, by default headwater subbasins have
    1-cell mainstream
  - `-d` flag allows `minmainstreams` to also apply to downstream subbasins
    (mainly intended to be used for pluvial flood investigations)
  - new required parameter `drainage` instead of `streams` which is set to the
    default output of `m.swim.subbasins`
- `m.swim.subbasin`: grided subbasin support
  - `-g` enables grid subbasins in current locations projects
  - `-l` enables grid subbasins in lon-lat projection
- subbasinoutlets vector no longer has category IDs equal to the subbasinID


*Version v1.3* - 2019-08-11
- `streamthresh` is now a required parameter to `m.swim.subbasins`
- fixed hydrotopes counting starting from 1
- predefined subbasins will only be included if they are within the catchment
- `m.swim.substats` is also writing out centroid latitude, reference elevation and initial water storage.
- Python3-proofed to run with grass 7.8


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
