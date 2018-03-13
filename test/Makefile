# all tests depend on the following maps in grassdb/utm32n/PERMANENT
# vector: stations, rivers
# raster: elevation, landuse, soil
GRASS ?= grass

GRASSTESTLOC=../test/grassdb/utm32n

MODULES = subbasins routing substats hydrotopes climate glaciers

TESTOUTPUTDIRS = project
GRASSTESTOUTPUT = $(GRASSTESTLOC)/swim $(GRASSTESTLOC)/glaciers

default: all

all: $(MODULES) checkoutput

checkoutput:
	find $(TESTOUTPUTDIRS) -type f -exec sha1sum {} + > output.sha1

clean:
	rm -r $(TESTOUTPUTDIRS) $(GRASSTESTOUTPUT)



$(MODULES):
	GRASS_BATCH_JOB=./test_$@.grass GRASS_OVERWRITE=1 $(GRASS) $(GRASSTESTLOC)/PERMANENT
	@echo Run make checkoutput to see if output differs!
