MODULE_TOPDIR =../..

PGM = m.swim

SUBDIRS = m.swim.subbasins \
          m.swim.hydrotopes \
          m.swim.routing \
          m.swim.substats \
          m.swim.climate 
#          m.swim.run   # not running yet

include $(MODULE_TOPDIR)/include/Make/Dir.make

default: parsubdirs htmldir

install: installsubdirs
	$(INSTALL_DATA) $(PGM).html $(INST_DIR)/docs/html/