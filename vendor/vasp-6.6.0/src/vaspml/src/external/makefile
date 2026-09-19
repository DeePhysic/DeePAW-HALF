.PHONY: \
	build \
	lammps \
	lammps-fetch \
	lammps-switch-branch \
	lammps-configure \
	lammps-build \
	lammps-copy \
	clean-lammps-fetch \
	clean-lammps-configure \
	clean-lammps-build

# Location of source files.
BASE_DIR = ../..
SRC_DIR = $(BASE_DIR)/src/applications
RES_DIR = $(BASE_DIR)/res
LIB_DIR = $(BASE_DIR)/lib
INC_DIR = $(BASE_DIR)/build/libvaspml
BIN_DIR = $(BASE_DIR)/bin

LIB_PREREQ = $(LIB_DIR)/libvaspml.a
LIB_LINK = -lvaspml

LAMMPS_BRANCH ?= vasp-mlff-6.6.0
CMAKE_BUILD_EXTRA ?= --parallel

# Include compiler/preprocessor/linker information.
include $(BASE_DIR)/makefile.include

# Common variables (e.g. for colored output)
include $(RES_DIR)/makefile.common

# Build all applications.
build: lammps

lammps: $(LIB_PREREQ) lammps-copy

lammps-fetch:
	@printf $(LOG_INOUT) "$(C_VGREY)Fetching LAMMPS from vasp-dev/lammps repository " $(LOG_ARROW) "$(C_VPURPLE)lammps$(C_OFF)\n"
	if [ ! -e lammps ]; then \
	  git clone https://github.com/vasp-dev/lammps.git; \
	fi

lammps-switch-branch: lammps-fetch
	@printf $(LOG_INOUT) "$(C_VGREY)Switching LAMMPS repository branch " $(LOG_ARROW) "$(C_VBLUE)$(LAMMPS_BRANCH)$(C_OFF)\n"
	cd lammps && \
	git checkout $(LAMMPS_BRANCH);

lammps-configure: lammps-switch-branch
	@printf $(LOG_INOUT) "$(C_VGREY)Configuring LAMMPS CMAKE build " $(LOG_ARROW) "$(C_VBLUE)lammps/build$(C_OFF)\n"
	cd lammps && \
	cmake -S cmake -B build -DCMAKE_CXX_COMPILER="$(CXX)" -DCMAKE_CXX_FLAGS="$(CXXFLAGS)" -DCMAKE_EXE_LINKER_FLAGS="$(LDFLAGS)" -DPKG_ML-VASP=yes -DVASPML_DIR=../../../../

lammps-build: lammps-configure
	@printf $(LOG_INOUT) "$(C_VGREY)Building LAMMPS via CMAKE " $(LOG_ARROW) "$(C_VCYAN)lammps/build$(C_OFF)\n"
	cd lammps && \
	cmake --build build $(CMAKE_BUILD_EXTRA)

lammps-copy: lammps-build
	@printf $(LOG_INOUT) "$(C_VGREY)Copying LAMMPS executable " $(LOG_ARROW) "$(C_VRED)bin/lmp$(C_OFF)\n"
	cp lammps/build/lmp $(BIN_DIR)

clean-lammps-fetch:
	$(RM) -r lammps

clean-lammps-configure:
	$(RM) -r lammps/build

clean-lammps-build:
	cd lammps && \
	cmake --build build --target clean

# Silence compiler output (re-enable with "VERBOSE=1" make argument).
$(VERBOSE).SILENT:
