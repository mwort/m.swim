# all tests depend on the following maps in grassdb/utm32n/PERMANENT
# vector: stations, rivers
# raster: elevation, landuse, soil

GRASSLOC=grassdb/utm32n

MODULES = subbasins routing substats hydrotopes

OUTPUTDIRS = project
GRASSOUTPUT = $(GRASSLOC)/swim

default: all

all: $(MODULES) checkoutput

checkoutput:
	find $(OUTPUTDIRS) -type f -exec sha1sum {} + > output.sha1

clean:
	rm -r $(OUTPUTDIRS) $(GRASSOUTPUT)



$(MODULES):
	GRASS_BATCH_JOB=./test_$@.grass GRASS_OVERWRITE=1 grass $(GRASSLOC)/PERMANENT
	echo Run make checkoutput to see if output differs!
