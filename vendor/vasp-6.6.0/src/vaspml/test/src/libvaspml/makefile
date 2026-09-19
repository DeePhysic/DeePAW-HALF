.PHONY: \
       	run \
       	build \
        depend \
        clean

# Name of test part.
TEST_PART = libvaspml

# Location of source files.
BASE_DIR = ../../..
SRC_DIR = $(BASE_DIR)/test/src/$(TEST_PART)
RES_DIR = $(BASE_DIR)/res
LIB_DIR = $(BASE_DIR)/lib
INC_DIR = $(BASE_DIR)/build/libvaspml
BIN_DIR = $(BASE_DIR)/bin
TEST_INC_DIR = $(BASE_DIR)/test/build/common
TEST_LIB_DIR = $(BASE_DIR)/test/lib

# Default number of ranks for MPI tests.
TEST_MPI_RANKS ?= 4

# Include compiler/preprocessor/linker information.
include $(BASE_DIR)/makefile.include

# Common variables (e.g. for colored output)
include $(RES_DIR)/makefile.common

# Source files with .cpp extension and prefix "test".
SRC = $(notdir $(wildcard $(SRC_DIR)/test*.cpp))
ifdef SELECT_TEST
SRC := $(filter test_$(SELECT_TEST).cpp,$(SRC))
endif
# Header files with .hpp extension.
HEADERS = $(notdir $(wildcard $(SRC_DIR)/*.hpp))
# Test executables created from source files.
TEST_ALL = $(SRC:.cpp=.x)
TEST_MPI = $(filter test_mpi%, $(TEST_ALL))
TEST = $(filter-out test_mpi%, $(TEST_ALL))
# For each source file a dependency file needs to be generated.
DEP = $(SRC:.cpp=.d)

LIB_PREREQ = $(LIB_DIR)/libvaspml.a
LIB_LINK = -lvaspml

# Run script needs to sum up all error codes to properly forward a combined exit status to make.
run: build
	err=0; \
	for x in $(TEST); do \
	printf $(LOG_INOUT) "$(C_VGREY)Running test executable " $(LOG_ARROW) "$(C_VRED)$${x}$(C_OFF)\n"; \
	export BOOST_TEST_LOGGER="$${BOOST_TEST_LOGGER:-HRF,test_suite}"; \
	./$${x}; \
	err=$$(($${err} + $$?)); \
	done; \
	for x in $(TEST_MPI); do \
	printf $(LOG_INOUT) "$(C_VGREY)Running test executable " $(LOG_ARROW) "$(C_VRED)$${x}$(C_OFF)\n"; \
	export BOOST_TEST_LOGGER="$${BOOST_TEST_LOGGER:-HRF,test_suite}"; \
	mpirun -np $(TEST_MPI_RANKS) ./$${x}; \
	err=$$(($${err} + $$?)); \
	done; \
	exit $${err};

build: $(TEST_ALL)

# Static pattern rule to fetch all changed source files.
$(SRC): %: $(SRC_DIR)/%
	@printf $(LOG_INOUT) "$(C_VGREY)Collecting source $< " $(LOG_ARROW) "$(C_VPURPLE)$(@)$(C_OFF)\n"
	rsync -u $< .

# Static pattern rule to fetch all changed header files.
$(HEADERS): %: $(SRC_DIR)/%
	@printf $(LOG_INOUT) "$(C_VGREY)Collecting header $< " $(LOG_ARROW) "$(C_VPURPLE)$(@)$(C_OFF)\n"
	rsync -u $< .

# Create dependency files for a single source file.
%.d: %.cpp $(HEADERS)
	@printf $(LOG_INOUT) "$(C_VGREY)Creating dependency for $<" $(LOG_ARROW) "$(C_VBLUE)$(@)$(C_OFF)\n"
	$(CPP_DEP) $(CXXFLAGS) $(INCLUDE) $(BOOST_INCLUDE) -I$(INC_DIR) -I$(TEST_INC_DIR) $(CPP_OPTIONS) -MM $< -MF $@

# Create dependency files for all source files.
depend: $(DEP)

# Default compile rule.
%.o: %.cpp
	@printf $(LOG_INOUT) "$(C_VGREY)Compiling $< " $(LOG_ARROW) "$(C_VCYAN)$(@)$(C_OFF)\n"
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(BOOST_INCLUDE) -I$(INC_DIR) -I$(TEST_INC_DIR) $(CPP_OPTIONS) -o $@ -c $<

# Linking test executable
$(TEST_ALL): %.x: %.o $(LIB_PREREQ) $(TEST_LIB_DIR)/libvaspml_test.a
	@printf $(LOG_INOUT) "$(C_VGREY)Linking $< " $(LOG_ARROW) "$(C_VGREEN)$(@)$(C_OFF)\n"
	$(CXX) $(CXXFLAGS) $(CPP_OPTIONS) -o $@ $< -L$(LIB_DIR) -L$(TEST_LIB_DIR) -lvaspml_test $(LIB_LINK) $(BOOST_LIB) $(LDFLAGS)

# Include dependency files created during "depend" stage.
ifeq (,$(filter clean,$(MAKECMDGOALS)))
-include $(DEP)
endif

# Remove all sources, objects, temporary and resulting files.
ifneq ($(realpath $(SRC_DIR)),$(CURDIR))
clean: clean-depend
	$(RM) *.o *.cpp *.hpp *.gcno *.gcda $(TEST_ALL)
else
clean:
	echo "ERROR: Invoking clean target in source folder, this is dangerous, aborting..."
	exit 1
endif

# Remove only dependencies helper file.
clean-depend:
	$(RM) *.d

# Silence compiler output (re-enable with "VERBOSE=1" make argument).
$(VERBOSE).SILENT:
