MODULE_TOPDIR = ../../..

include $(MODULE_TOPDIR)/include/Make/Other.make
include $(MODULE_TOPDIR)/include/Make/Python.make

PKG = mswim
TOOLBOX = m.swim
MODULES = __init__ utils inout

ETCDIR = $(ETC)/$(TOOLBOX)/$(PKG)
ETCINSTALLDIR = $(INST_DIR)/etc/$(TOOLBOX)

PYFILES := $(patsubst %,$(ETCDIR)/%.py,$(MODULES))
PYCFILES := $(patsubst %,$(ETCDIR)/%.pyc,$(MODULES))

default: $(PYFILES) $(PYCFILES)

$(ETCDIR):
	$(MKDIR) $@

$(ETCDIR)/%: % | $(ETCDIR)
	$(INSTALL_DATA) $< $@

install:
	$(MKDIR) $(ETCINSTALLDIR)
	cp -r $(ETCDIR) $(ETCINSTALLDIR)
