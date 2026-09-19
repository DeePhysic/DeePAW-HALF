# Location of source files.
BASE_DIR = ..
RES_DIR = $(BASE_DIR)/res
BUILD_DIR = build
SRC_DIR = src
COV_DIR = $(BASE_DIR)/

# Include compiler/preprocessor/linker information.
include $(BASE_DIR)/makefile.include

# Common variables (e.g. for colored output)
include $(RES_DIR)/makefile.common

TARGETS = libvaspml
CLEAN_TARGETS = $(patsubst %, clean-%, $(TARGETS))

.PHONY: all \
	$(TARGETS) \
	coverage \
	clean \
	$(CLEAN_TARGETS) \
	clean-coverage

all: $(TARGETS)

main-libvaspml:
	$(MAKE) -C $(BASE_DIR) libvaspml

$(TARGETS): common
ifdef SELECT_TEST
	@printf "Building and running test part: \"$(@)\" (selected test: $(SELECT_TEST)) in directory $(BUILD_DIR)/$@...\n"
else
	@printf "Building and running test part: \"$(@)\" in directory $(BUILD_DIR)/$@...\n"
endif
	if [ ! -d $(BUILD_DIR)/$@ ] ; then mkdir -p $(BUILD_DIR)/$@ ; fi
	cp $(SRC_DIR)/$@/makefile $(BUILD_DIR)/$@
	$(MAKE) -C $(BUILD_DIR)/$@ depend
	$(MAKE) -C $(BUILD_DIR)/$@ run
	@printf -- "$(LOG_SEPARATOR)\n"

common: main-libvaspml
	@printf "Building common test infrastructure in directory $(BUILD_DIR)/$@...\n"
	if [ ! -d $(BUILD_DIR)/$@ ] ; then mkdir -p $(BUILD_DIR)/$@ ; fi
	cp $(SRC_DIR)/$@/makefile $(BUILD_DIR)/$@
	$(MAKE) -C $(BUILD_DIR)/$@ depend
	$(MAKE) -C $(BUILD_DIR)/$@ build
	@printf -- "$(LOG_SEPARATOR)\n"

coverage:
	@printf "$(C_OFF)Generating coverage report in directory cov-html...\n$(C_VGREY)"
	gcovr --xml-pretty --exclude-unreachable-branches -e "build/common/*" -e "build/libvaspml/*" -e "../build/benchmarks/*" -e "../build/applications/test_*" --print-summary -o cov-xml/coverage.xml --html-details cov-html/index.html -r $(COV_DIR)
	@printf -- "$(C_OFF)$(LOG_SEPARATOR)\n"

clean: $(CLEAN_TARGETS) clean-common clean-coverage

clean-main-libvaspml:
	$(MAKE) -C $(BASE_DIR) clean

clean-all: clean-main-libvaspml clean

$(CLEAN_TARGETS):
	@printf "Removing build directory for test part:   \"$(subst clean-,,$@)\".\n"
	$(RM) -r $(BUILD_DIR)/$(subst clean-,,$@)
	@printf -- "$(LOG_SEPARATOR)\n"

clean-common:
	@printf "Removing build directory for common test infrastructure.\n"
	$(RM) -r $(BUILD_DIR)/$(subst clean-,,$@)
	$(RM) include/*
	$(RM) lib/*
	@printf -- "$(LOG_SEPARATOR)\n"

clean-coverage:
	@printf "Removing coverage data.\n"
	$(RM) cov-html/*.html cov-xml/*.xml
	$(RM) $(BASE_DIR)/build/libvaspml/*.gcda
	$(RM) $(BUILD_DIR)/common/*.gcda
	$(RM) $(BUILD_DIR)/libvaspml/*.gcda
	@printf -- "$(LOG_SEPARATOR)\n"

# Silence compiler output (re-enable with "VERBOSE=1" make argument).
$(VERBOSE).SILENT:
