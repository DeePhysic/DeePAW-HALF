ROCM_PATH ?= /opt/rocm
HIPCC     ?= ${ROCM_PATH}/bin/hipcc

include makefile.include

TARGET = libHipInterface.a
HOBJS  = runtime_interface.o \
         rocblas_interface.o \
         rocfft_interface.o \
         rccl_interface.o

.PHONY: all clean

all: ${TARGET}

%.o: %.cpp
	${HIPCC} ${HIPFLAGS} ${OFLAG} ${ROCM_INCS} -o $@ -c $<

${TARGET}: $(HOBJS)
	rm -f ${TARGET}
	ar rucv ${TARGET} ${HOBJS}

clean:
	rm -rf *.o *.a
