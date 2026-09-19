.SUFFIXES:
.SUFFIXES: .mod .o .f90
.PRECIOUS: %.f90
# .SECONDEXPANSION:

BUILD_INFO=1

FUFFIX=.F
SUFFIX=.f90

OFLAG_1=-O1
OFLAG_2=-O2
OFLAG_3=-O3
OFLAG_4=

OFLAG=$(OFLAG_2)
OFLAG_IN=$(OFLAG)

LIB=lib parser
LLIB=-Llib -ldmy -Lparser -lparser

SRCDIR=../../src
BINDIR=../../bin

include .objects
include makefile.include

# Add vaspml to libs when -Dlibvaspml
VASPML=$(findstring Dlibvaspml,$(CPP))
ifeq ($(VASPML),Dlibvaspml)
LIB+=vaspml
LLIB+=-Lvaspml/lib -lvaspml
endif

# FPP=$(CPP) $*$(FUFFIX) $*$(SUFFIX) $(CPP_OPTIONS)

FPP=$(filter-out -DwNGZhalf -DNGZhalf -DwNGXhalf -DNGXhalf,$(CPP))

DMPI=$(findstring DMPI,$(subst DMPI_,,$(CPP)))

LIB+=$(filter-out lib,$(LIBS))
LLIB+=$(filter-out -Llib -ldmy,$(LLIBS))

PLUGINS=$(findstring DPLUGINS,$(CPP))
ifeq ($(PLUGINS),DPLUGINS)
    LIB+=plugins
    LLIB+=-Lplugins -lplugins
    INCS+=-Iplugins
endif

OBJS=$(SOURCE) $(OBJECTS)

ifndef VERSION
   ifeq ($(DMPI),DMPI)
      FPP+=-DNGZhalf
   else
      FPP+=-DNGXhalf
   endif
   EXE=vasp_std
else
ifeq ($(VERSION),std)
   ifeq ($(DMPI),DMPI)
      FPP+=-DNGZhalf
   else
      FPP+=-DNGXhalf
   endif
   EXE=vasp_std
else ifeq ($(VERSION),gam)
   ifeq ($(DMPI),DMPI)
      FPP+=-DNGZhalf -DwNGZhalf
   else
      FPP+=-DNGXhalf -DwNGXhalf
   endif
   EXE=vasp_gam
else ifeq ($(VERSION),ncl)
   EXE=vasp_ncl
else
   $(error $(VERSION) is not a valid target)
endif
endif

DUSELIBXC=$(findstring DUSELIBXC,$(CPP))
ifeq ($(DUSELIBXC),DUSELIBXC)
   FPP+=$(INCS)
endif

OBJS_O1=$(filter-out $(OBJECTS_O1) $(OBJECTS_O2) $(OBJECTS_O3) $(OBJECTS_IN),$(SOURCE_O1)) $(OBJECTS_O1)
OBJS_O2=$(filter-out $(OBJECTS_O1) $(OBJECTS_O2) $(OBJECTS_O3) $(OBJECTS_IN),$(SOURCE_O2)) $(OBJECTS_O2)
OBJS_O3=$(filter-out $(OBJECTS_O1) $(OBJECTS_O2) $(OBJECTS_O3) $(OBJECTS_IN),$(SOURCE_O3)) $(OBJECTS_O3)
OBJS_IN=$(filter-out $(OBJECTS_O1) $(OBJECTS_O2) $(OBJECTS_O3) $(OBJECTS_IN),$(SOURCE_IN)) $(OBJECTS_IN)

F90SRC_O1=$(OBJS_O1:.o=$(SUFFIX))
F90SRC_O2=$(OBJS_O2:.o=$(SUFFIX))
F90SRC_O3=$(OBJS_O3:.o=$(SUFFIX))
F90SRC_IN=$(OBJS_IN:.o=$(SUFFIX))

ifeq ($(FFLAGS_1),)
   FFLAGS_1=$(FFLAGS)
endif
ifeq ($(FFLAGS_2),)
   FFLAGS_2=$(FFLAGS)
endif
ifeq ($(FFLAGS_3),)
   FFLAGS_3=$(FFLAGS)
endif
ifeq ($(FFLAGS_IN),)
   FFLAGS_IN=$(FFLAGS)
endif

ifeq ($(OFLAG_1),)
   OFLAG_1=$(OFLAG)
endif
ifeq ($(OFLAG_2),)
   OFLAG_2=$(OFLAG)
endif
ifeq ($(OFLAG_3),)
   OFLAG_3=$(OFLAG)
endif

ifeq ($(INCS_1),)
   INCS_1=$(INCS)
endif
ifeq ($(INCS_2),)
   INCS_2=$(INCS)
endif
ifeq ($(INCS_3),)
   INCS_3=$(INCS)
endif
ifeq ($(INCS_IN),)
   INCS_IN=$(INCS)
endif


OBJCTS=$(filter %.o, $(OBJS) $(FFT3D))
#OBJCTS_f90=$(filter-out getshmem.o, $(OBJCTS))
F90SRC=$(OBJCTS:.o=$(SUFFIX))
F90SRC+=main.f90

ifdef MODS
INCMODS=-Imodfiles
endif

-include .depend

# export

.PHONY: all cleanall clean sources objects modfiles dependencies depend libs $(LIB)

all: libs sources
	rm -f vasp ; $(MAKE) vasp ; cp vasp $(BINDIR)/$(EXE)

vasp: $(OBJS) $(FFT3D) $(INC) main.o
	$(FCL) -o vasp $(OBJS) main.o $(FFT3D) $(LLIB) $(LINK)

objects: sources $(LIB) $(OBJS) $(FFT3D) $(INC) main.o

sources:
	rsync -u $(SRCDIR)/*.F $(SRCDIR)/*.inc .
ifeq ($(BUILD_INFO),1)
	$(MAKE) build_info
endif
#	-rsync -u $(SRCDIR)/*.f .

libs: $(LIB)

$(filter-out vaspml,$(LIB)):
	rsync -ru $(SRCDIR)/$@ .
	cp makefile.include $@
	$(MAKE) -C $@ -j1
#	$(MAKE) -C $@

vaspml:
	rsync -ru $(SRCDIR)/$@ .
	cp $@/arch/makefile.include.vasp $@/makefile.include
	$(MAKE) -C $@ libvaspml

dependencies: sources
	$(MAKE) depend

depend: $(F90SRC)
	awk -f makedeps.awk $(F90SRC) > .depend
#	makedepf90 $(F90SRC) > .depend

modfiles: sources
	if [ ! -d modfiles ] ; then mkdir modfiles ; fi
	$(MAKE) -j1 mods ; cp *.mod modfiles

mods: $(F90SRC)
	$(FC) $(FREE) $(FFLAGS) $(INCS) $(FFSYNTAX) -c $(F90SRC)

%$(SUFFIX): %$(FUFFIX)
	$(FPP)

%.o: %$(SUFFIX)
	$(FC) $(FREE) $(FFLAGS) $(OFLAG) $(INCS) $(INCMODS) -c $*$(SUFFIX)

# $(OBJS_O1): $(F90SRC_O1)
$(OBJS_O1): %.o: %$(SUFFIX)
	$(FC) $(FREE) $(FFLAGS_1) $(OFLAG_1) $(INCS_1) $(INCMODS) -c $*$(SUFFIX)

# $(OBJS_O2): $(F90SRC_O2)
$(OBJS_O2): %.o: %$(SUFFIX)
	$(FC) $(FREE) $(FFLAGS_2) $(OFLAG_2) $(INCS_2) $(INCMODS) -c $*$(SUFFIX)

# $(OBJS_O3): $(F90SRC_O3)
$(OBJS_O3): %.o: %$(SUFFIX)
	$(FC) $(FREE) $(FFLAGS_3) $(OFLAG_3) $(INCS_3) $(INCMODS) -c $*$(SUFFIX)

# $(OBJS_IN): $(F90SRC_IN)
$(OBJS_IN): %.o: %$(SUFFIX)
	$(FC) $(FREE) $(FFLAGS_IN) $(OFLAG_IN) $(INCS_IN) $(INCMODS) -c $*$(SUFFIX)

cleandependencies:
	rm -f .depend

cleanall: clean
	rm -rf $(LIB) *.F test.x

clean:
	rm -f *.o *.mod *.f90

#-----------------------------------------------------------------------
# build info
#-----------------------------------------------------------------------
ifeq ($(BUILD_INFO),1)
# small helper function to quote strings
PQUOTE =$(subst ",\",$(1))

CPP_OPTIONS_STRING =$(patsubst %,&\n&%,$(CPP_OPTIONS))
LINK_LINE_STRING   =$(patsubst %,&\n&%,$(call PQUOTE,$(FCL) $(FFT3D) $(LLIB) $(LINK)))

FC_STRING     =$(patsubst %,&\n&%,$(call PQUOTE,$(FC)))
FCL_STRING    =$(patsubst %,&\n&%,$(call PQUOTE,$(FCL)))
FFLAGS_STRING =$(patsubst %,&\n&%,$(call PQUOTE,$(FFLAGS)))
LLIBS_STRING  =$(patsubst %,&\n&%,$(call PQUOTE,$(LLIBS)))
INCS_STRING   =$(patsubst %,&\n&%,$(call PQUOTE,$(INCS)))

build_info:
	# CPP_OPTIONS_STRING may have spaces in variables, so we format it differently
	printf "    character(len=*), parameter :: cpp_options = '&\n" >  build_info.inc
	for part in $(CPP_OPTIONS); do printf "&%s &\n" $$part; done >> build_info.inc
	printf "    &'\n" >> build_info.inc
	printf "    character(len=*), parameter :: link_line   = '$(LINK_LINE_STRING)'\n"   >> build_info.inc
	printf "    character(len=*), parameter :: fc     = '$(FC_STRING)'\n"     >>  build_info.inc
	printf "    character(len=*), parameter :: fcl    = '$(FCL_STRING)'\n"    >> build_info.inc
	printf "    character(len=*), parameter :: fflags = '$(FFLAGS_STRING)'\n" >> build_info.inc
	printf "    character(len=*), parameter :: llibs  = '$(LLIBS_STRING)'\n"  >>  build_info.inc
	printf "    character(len=*), parameter :: incs   = '$(INCS_STRING)'\n"   >> build_info.inc

##CPP_OPTIONS_STRING=$(patsubst %,\&\n\&%,$(CPP_OPTIONS))
##LINK_LINE_STRING  =$(patsubst %,\&\n\&%,$(FCL) $(FFT3D) $(LLIB) $(LINK))
##FC_STRING     =$(patsubst %,\&\n\&%,$(FC))
##FCL_STRING    =$(patsubst %,\&\n\&%,$(FCL))
##FFLAGS_STRING =$(patsubst %,\&\n\&%,$(FFLAGS))
##LLIBS_STRING  =$(patsubst %,\&\n\&%,$(LLIBS))
##INCS_STRING   =$(patsubst %,\&\n\&%,$(INCS))
##
##build_info:
##	sed -i.bak -E "s?(cpp_options .*= .*)'Unknown'?\1'$(CPP_OPTIONS_STRING)'?" build_info.inc
##	sed -i.bak -E "s?(link_line .*= .*)'Unknown'?\1'$(LINK_LINE_STRING)'?" build_info.inc
##	sed -i.bak -E "s?(fc .*= .*)'Unknown'?\1'$(FC_STRING)'?" build_info.inc
##	sed -i.bak -E "s?(fcl .*= .*)'Unknown'?\1'$(FCL_STRING)'?" build_info.inc
##	sed -i.bak -E "s?(fflags .*= .*)'Unknown'?\1'$(FFLAGS_STRING)'?" build_info.inc
##	sed -i.bak -E "s?(llibs .*= .*)'Unknown'?\1'$(LLIBS_STRING)'?" build_info.inc
##	sed -i.bak -E "s?(incs .*= .*)'Unknown'?\1'$(INCS_STRING)'?" build_info.inc
endif
#-----------------------------------------------------------------------
# sanity checks
#-----------------------------------------------------------------------

DOFFLOAD=$(findstring DACC_OFFLOAD,$(CPP))$(findstring DOMP_OFFLOAD,$(CPP))

check:
# check if DOFFLOAD is empty, then we can compile with -DUSELIBXC:
ifeq ($(DOFFLOAD),)
# if DOFFLOAD is not empty we cannot compile with -DUSELIBXC:
else
ifeq ($(DUSELIBXC),DUSELIBXC)
	@echo "CONFLICT: the GPU port (ACC_OFFLOAD or OMP_OFFLOAD) of VASP does not support compilation with -DUSELIBXC"
	exit 1
endif
endif
	exit 0

#-----------------------------------------------------------------------
# special rules
#-----------------------------------------------------------------------

main.o : main$(SUFFIX)
	$(FC) $(FREE) $(FFLAGS) $(DEBUG) $(INCS) $(INCMODS) -c $*$(SUFFIX)
