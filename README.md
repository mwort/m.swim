# The GRASS7 Soil and Water Integrated Model (SWIM) preprocessor

Version v2.0

## Documentation

The documentation is located here:
http://www.pik-potsdam.de/~wortmann/m.swim/

## Version history

*Version v2.0* - 2022-03-22
- first version to support SWIM redesign csv input files
- a long-term support version for the old input file formats is available on the `v1_lts` branch
- several argument name changes to be more consistent with common grass modules and the corresponding SWIM variable names (see below).
- relies on the `mswim` library for the first time
- `m.swim.hydrotopes`:
  - `subbasins, landuse, soil, management` arguments renamed to `subbasin_id, landuse_id, soil_id, crop_management_id` in-line with new column names of the `hydrotope.csv` file
  - `strfilepath` renamed to `output`
  - optional columns with default values (e.g. wetland, crop_management_id etc.) are no longer written to file unless a map is given to populate them
- `m.swim.routing`:
  - `figpath` argument renamed to `output`
- `m.swim.substats`:
  - `projectpath` renamed to `output` for the `subbasin.csv` file
  - optional columns with default values are no longer written to file unless a map is given to populate them. If no-standard columns are needed, they will need to be added to the `gworder, suborder, routeorder` arguments.


*Version v1.6* - 2022-01-19
- last version supporting legacy SWIM input
- added the `m.swim.glaciers` module to create SWIM-G input files
- `m.swim.hydrotopes`
  - `countoursrast` is now only output file, custom contours raster is passed as string
    to `contours`
  - `-c` is not needed anymore and is deprecated
- several bug fixes

*Version v1.5* - 2020-07-28
- alpha implementation of the `mswim` python package to contain abstracted/shared
  functionality between all modules, mainly to test if it's working
- `-v` to show the version and install date of modules and grass itself
- `m.swim.climate`:
  - support for R-based climate interpolation input dropped
  - point-defined grid definition enabled via `grid=`, `-d` and optionally
    `lon_column, lat_column`, e.g. for rotated grids
  - `ncinfopath` argument renamed to `gridfilepath` with deprecation notice

*Version v1.4* - 2019-12-02
- `m.swim.routing`: topologically correct mainstreams vector output
  - downstream directions, nodes at intersections and subbasin boundaries
  - subbasinID and accumulation (min, mean, max) columns in table
  - network nodes in layer 2 (for network analysis)
  - optional `minmainstreams` parameter, by default headwater subbasins have
    1-cell mainstream
  - new required parameter `drainage` instead of `streams` which is set to the
    default output of `m.swim.subbasins`
- `m.swim.subbasin`: grided subbasin support
  - `-g` enables grid subbasins in current locations projects
  - `-l` enables grid subbasins in lon-lat projection
- documentation moved to http://www.pik-potsdam.de/~wortmann/m.swim
- subbasinoutlets vector no longer has category IDs equal to the subbasinID
- dropped `m.swim.run` as succeeded by [swimpy](http://www.pik-potsdam.de/~wortmann/swimpy/)


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

There are tests for each grass module based on the Blankenstein catchment
(upper Saale (Elbe)). All required input maps are in the
`grassdb/utm32n/PERMANENT` mapset, all tests should therefore be executed in:
```
cd test
grass grassdb/utm32n/PERMANENT
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
- publish docs: `cd doc; make [URL=...] publish`
