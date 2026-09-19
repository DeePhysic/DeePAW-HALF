.PHONY: \
	all \
	build \
	depend \
	link \
	clean

# Name of output library.
LIB_NAME = libvaspml

# Location of source files.
BASE_DIR = ../..
SRC_DIR = $(BASE_DIR)/src/$(LIB_NAME)
RES_DIR = $(BASE_DIR)/res
LIB_DIR = $(BASE_DIR)/lib

# Include compiler/preprocessor/linker information.
include $(BASE_DIR)/makefile.include

# Common variables (e.g. for colored output)
include $(RES_DIR)/makefile.common

# Source files with .cpp extension.
SRC = $(notdir $(wildcard $(SRC_DIR)/*.cpp))
# Header files with .hpp extension.
HEADERS = $(notdir $(wildcard $(SRC_DIR)/*.hpp))
# List of object files generated from source files list.
OBJ = $(SRC:.cpp=.o)
# For each source file a dependency file needs to be generated.
DEP = $(SRC:.cpp=.d)

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
	$(CPP_DEP) $(CXXFLAGS) $(INCLUDE) $(CPP_OPTIONS) -MM $< -MF $@

# Create dependency files for all source files.
depend: $(DEP)

# Default compile rule.
%.o: %.cpp
	@printf $(LOG_INOUT) "$(C_VGREY)Compiling $< " $(LOG_ARROW) "$(C_VCYAN)$(@)$(C_OFF)\n"
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(CPP_OPTIONS) -o $@ -c $<

# Link library.
link: $(OBJ)
	@printf $(LOG_INOUT) "$(C_VGREY)Linking library" $(LOG_ARROW) "$(C_VGREEN)$(LIB_NAME).a$(C_OFF)\n";
	$(AR) $(ARFLAGS) $(LIB_NAME).a $+

# Copy library to lib path.
build: link
	@printf $(LOG_INOUT) "$(C_VGREY)Copying library $(LIB_NAME).a" $(LOG_ARROW) "$(C_VRED)$(LIB_DIR)/$(LIB_NAME).a$(C_OFF)\n";
	cp $(LIB_NAME).a $(LIB_DIR)

# Include dependency files created during "depend" stage.
ifeq (,$(filter clean,$(MAKECMDGOALS)))
-include $(DEP)
endif

# Remove all sources, objects, temporary and resulting files.
ifneq ($(realpath $(SRC_DIR)),$(CURDIR))
clean: clean-depend
	$(RM) *.cpp *.hpp *.o $(LIB_NAME).a $(LIB_DIR)/$(LIB_NAME).a
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
