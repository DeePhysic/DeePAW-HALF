include makefile.include

.PHONY: all clean

FUFFIX=.F
SUFFIX=.f90

all: libplugins.a

libplugins.a: interface_fortran.o interface_cpp.o
	rm -f $@
	ar vq $@ $^

%$(SUFFIX): %$(FUFFIX)
	$(CPP)

interface_fortran.o : interface_fortran.f90
	$(FC) $(FREE) $(FFLAGS) $(OFLAG) $(INCS) $(INCMODS) -c $<

interface_cpp.o : interface_cpp.cpp interface_cpp.h
	$(CXX_PARS) $(CXX_FLAGS) -c $<

clean:
	rm *$(SUFFIX) *.o *.a
