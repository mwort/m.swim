MODULE_TOPDIR =../..

PGM = m.swim

SUBDIRS = mswim \
          m.swim.subbasins \
          m.swim.hydrotopes \
          m.swim.routing \
          m.swim.substats \
          m.swim.climate 

include $(MODULE_TOPDIR)/include/Make/Dir.make

default: parsubdirs htmldir

install: installsubdirs
	$(INSTALL_DATA) $(PGM).html $(INST_DIR)/docs/html/
