CXX_ZE      ?= icpx
##CXXFLAGS_ZE ?= -fiopenmp -fopenmp-targets=spir64 -lze_loader
CXXFLAGS_ZE ?= -fiopenmp -fopenmp-targets=spir64

include makefile.include

TARGET = libZeInterface.a
HOBJS  = ze_device_get_uuid.o

.PHONY: all clean

all: ${TARGET}

%.o: %.cpp
	${CXX_ZE} ${CXXFLAGS_ZE} -o $@ -c $<

${TARGET}: $(HOBJS)
	rm -f ${TARGET}
	ar rucv ${TARGET} ${HOBJS}

clean:
	rm -rf *.o *.a
