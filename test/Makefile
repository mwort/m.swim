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

# get shasum for each output file in pathname sorted list
checkoutput:
	find $(TESTOUTPUTDIRS) -type f | sort -f | xargs shasum > output.sha1

clean:
	rm -r $(TESTOUTPUTDIRS) $(GRASSTESTOUTPUT)



$(MODULES):
	GRASS_OVERWRITE=1 $(GRASS) $(GRASSTESTLOC)/PERMANENT --exec bash ./test_$@.grass
	@echo Run make checkoutput to see if output differs!
