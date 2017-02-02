# all tests depend on the following maps in grassdb/utm32n/PERMANENT
# vector: stations, rivers
# raster: elevation, landuse, soil

GRASSLOC=grassdb/utm32n

MODULES = subbasins routing substats hydrotopes

default: all

all: $(MODULES)

clean:
	rm -r project/input/* $(GRASSLOC)/swim



$(MODULES):
	GRASS_BATCH_JOB=./test_$@.grass GRASS_OVERWRITE=1 grass $(GRASSLOC)/PERMANENT
