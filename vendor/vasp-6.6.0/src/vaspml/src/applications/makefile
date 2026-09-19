.PHONY: \
	all \
	build \
	clean

# Names of all applications.
APP = mlab2json \
      mlab2bin \
      mlab2xyz \
      mlff2json \
      mlffdiff \
      incar2json \
      vaspml-refit \
      vaspml-predict \
      vaspml-predict-grace \
      vaspml-concat \
      vaspml-inference \
      vaspml-select

# Location of source files.
BASE_DIR = ../..
SRC_DIR = $(BASE_DIR)/src/applications
RES_DIR = $(BASE_DIR)/res
LIB_DIR = $(BASE_DIR)/lib
INC_DIR = $(BASE_DIR)/build/libvaspml
BIN_DIR = $(BASE_DIR)/bin

LIB_PREREQ = $(LIB_DIR)/libvaspml.a
LIB_LINK = -lvaspml

# Include compiler/preprocessor/linker information.
include $(BASE_DIR)/makefile.include

# Common variables (e.g. for colored output)
include $(RES_DIR)/makefile.common

# Derive source files from application names.
SRC = $(addsuffix .cpp,$(APP))

# Static pattern rule to fetch all changed source files.
$(SRC): %: $(SRC_DIR)/%
	@printf $(LOG_INOUT) "$(C_VGREY)Collecting source $< " $(LOG_ARROW) "$(C_VPURPLE)$(@)$(C_OFF)\n"
	rsync -u $< .

# Default compile rule.
$(APP): %: %.cpp $(LIB_PREREQ)
	@printf $(LOG_INOUT) "$(C_VGREY)Compiling $< " $(LOG_ARROW) "$(C_VCYAN)$(@)$(C_OFF)\n"
	$(CXX) $(CXXFLAGS) $(INCLUDE) -I$(INC_DIR) $(CPP_OPTIONS) -c $<
	$(LD) $(CXXFLAGS) $(INCLUDE) -I$(INC_DIR) $(CPP_OPTIONS) -o $@ $@.o -L$(LIB_DIR) $(LIB_LINK) $(LDFLAGS)
	@printf $(LOG_INOUT) "$(C_VGREY)Copying executable $(@)" $(LOG_ARROW) "$(C_VRED)$(BIN_DIR)/$(@)$(C_OFF)\n";
	cp $@ $(BIN_DIR)

# Build all applications.
build: $(APP)

# Remove all sources, objects, temporary and resulting files.
ifneq ($(realpath $(SRC_DIR)),$(CURDIR))
clean:
	$(RM) *.cpp *.o $(APP)
else
clean:
	echo "ERROR: Invoking clean target in source folder, this is dangerous, aborting..."
	exit 1
endif

# Silence compiler output (re-enable with "VERBOSE=1" make argument).
$(VERBOSE).SILENT:
