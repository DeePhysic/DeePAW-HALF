VASPML_BUILD_DIR = build
VASPML_SRC_DIR = src
VASPML_RES_DIR = res

# Include compiler/preprocessor/linker information.
include makefile.include

# Common variables (e.g. for colored output)
include $(VASPML_RES_DIR)/makefile.common

TARGETS = libvaspml applications external
CLEAN_TARGETS = $(patsubst %, clean-%, $(TARGETS))

.PHONY: all \
        logo \
        pre \
        $(TARGETS) \
        $(CLEAN_TARGETS) \
        test \
        clean \
        clean-external-configure \
        clean-external-build \
        clean-test

all: $(filter-out external,$(TARGETS))

ifneq (,$(findstring no-logo,$(MAKE_OPTIONS)))
logo: ;
else
logo:
	@printf "$(C_VPURPLE) _    _____   _____ ____            __\n"
	@printf "$(C_VBLUE)| |  / /   | / ___// __ \____ ___  / /\n"
	@printf "$(C_VCYAN)| | / / /| | \__ \/ /_/ / __ \`__ \/ / \n"
	@printf "$(C_VGREEN)| |/ / ___ |___/ / ____/ / / / / / /  \n"
	@printf "$(C_VRED)|___/_/  |_/____/_/   /_/ /_/ /_/_/$(C_OFF)\n\n"
endif

pre: logo

libvaspml: pre
	@printf "Building VASPml C++ library in directory $(VASPML_BUILD_DIR)/$@...\n"
	if [ ! -d $(VASPML_BUILD_DIR)/$@ ] ; then mkdir -p $(VASPML_BUILD_DIR)/$@ ; fi
	cp $(VASPML_SRC_DIR)/$@/makefile $(VASPML_BUILD_DIR)/$@
	$(MAKE) -C $(VASPML_BUILD_DIR)/$@ depend
	$(MAKE) -C $(VASPML_BUILD_DIR)/$@ build
	@printf -- "$(LOG_SEPARATOR)\n"

applications: libvaspml
	@printf "Building VASPml applications in directory $(VASPML_BUILD_DIR)/$@...\n"
	if [ ! -d $(VASPML_BUILD_DIR)/$@ ] ; then mkdir -p $(VASPML_BUILD_DIR)/$@ ; fi
	cp $(VASPML_SRC_DIR)/$@/makefile $(VASPML_BUILD_DIR)/$@
	$(MAKE) -C $(VASPML_BUILD_DIR)/$@ build
	@printf -- "$(LOG_SEPARATOR)\n"

external: libvaspml
	@printf "Building third-party software with VASPml support in directory $(VASPML_BUILD_DIR)/$@...\n"
	if [ ! -d $(VASPML_BUILD_DIR)/$@ ] ; then mkdir -p $(VASPML_BUILD_DIR)/$@ ; fi
	cp $(VASPML_SRC_DIR)/$@/makefile $(VASPML_BUILD_DIR)/$@
	$(MAKE) -C $(VASPML_BUILD_DIR)/$@ build
	@printf -- "$(LOG_SEPARATOR)\n"

test:
	@printf "Building and running tests in directory ./test...\n"
	$(MAKE) -C test

clean: $(CLEAN_TARGETS)

clean-external-configure:
	@printf "Removing configure settings and build of third-party software...\n"
	$(MAKE) -C $(VASPML_BUILD_DIR)/external clean-lammps-configure

clean-external-build:
	@printf "Removing build of third-party software...\n"
	$(MAKE) -C $(VASPML_BUILD_DIR)/external clean-lammps-build

clean-test:
	@printf "Removing build directory for tests...\n"
	$(MAKE) -C test clean

$(CLEAN_TARGETS):
	@printf "Removing build directory for target:   \"$(subst clean-,,$@)\".\n"
	$(RM) -r $(VASPML_BUILD_DIR)/$(subst clean-,,$@)
	$(RM) bin/*
	$(RM) lib/*
	@printf -- "$(LOG_SEPARATOR)\n"

$(VERBOSE).SILENT:
