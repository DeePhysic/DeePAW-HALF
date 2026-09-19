#include "rccl/rccl.h"
#include <hip_runtime.h>
#include <iostream>
#include <cstdlib>

inline void rccl_check(ncclResult_t err) {
    if (err != ncclSuccess) {
        switch (err) {
            case ncclUnhandledCudaError:
                std::cerr << "ncclUnhandledCudaError" << std::endl; break;
            case ncclSystemError:
                std::cerr << "ncclSystemError" << std::endl; break;
            case ncclInternalError:
                std::cerr << "ncclInternalError" << std::endl; break;
            case ncclInvalidArgument:
                std::cerr << "ncclInvalidArgument" << std::endl; break;
            case ncclInvalidUsage:
                std::cerr << "ncclInvalidUsage" << std::endl; break;
            case ncclNumResults:
                std::cerr << "ncclNumResults" << std::endl; break;
            default:
                std::cerr << "Unknown rccl_status" << std::endl; break;
        }
        exit(-1);
    }
}

extern "C" {
    void nccl_groupstart() {
        rccl_check(ncclGroupStart());
    }

    void nccl_groupend() {
        rccl_check(ncclGroupEnd());
    }

    // ncclRes must be passed by reference (pointer), so the caller sees the return code
    void nccl_getuniqueid(char* ncclUID,
                          ncclResult_t* ncclRes) {

        // local variables
        int i;
        ncclUniqueId ncclUID2;

        *ncclRes=ncclGetUniqueId(&ncclUID2);

        // copy result as chunk directly
		memcpy(ncclUID, &ncclUID2, sizeof(ncclUniqueId));

        rccl_check(*ncclRes);
    }

    void nccl_comminitrank(void* nccl_comm,
                           int ncpu,
                           char ncclUID[],
                           int node_me,
                           ncclResult_t* ncclRes) {

        // local variables
        int i;
        ncclUniqueId ncclUID2;
        ncclComm_t* nccl_comm2=(ncclComm_t*) nccl_comm;

        // Copy UniqueId from caller into local struct
        memcpy(&ncclUID2, ncclUID, sizeof(ncclUniqueId));

        // Initialize communicator
        *ncclRes=ncclCommInitRank(nccl_comm2, ncpu, ncclUID2, node_me);

        rccl_check(*ncclRes); // abort if error
    }

    void nccl_bcast(void* vec,
                    int counts,
                    ncclDataType_t datatype,
                    int root,
                    void* comm,
                    void* stream,
                    ncclResult_t *ncclRes) {

        // local variables
        double* ptr=(double*) vec;
        ncclComm_t* nccl_comm= (ncclComm_t*) comm;
        hipStream_t* hipstream= (hipStream_t*) stream;

        *ncclRes=ncclBcast(ptr,counts,datatype,root,*nccl_comm,*hipstream);

        rccl_check(*ncclRes); // abort if error
    }

    void nccl_reduce(void* vec1,
                     void* vec2,
                     int counts,
                     ncclDataType_t datatype,
                     ncclRedOp_t op,
                     int root,
                     void* comm,
                     void* stream,
                     ncclResult_t *ncclRes) {

        // local variables
        double* ptr1=(double*) vec1;
        double* ptr2=(double*) vec2;
        ncclComm_t* nccl_comm= (ncclComm_t*) comm;
        hipStream_t* hipstream= (hipStream_t*) stream;

        *ncclRes=ncclReduce(ptr1,ptr2,counts,datatype,op,root,*nccl_comm,*hipstream);

        rccl_check(*ncclRes); // abort if error
    }

    void nccl_send(void* vec,
                   int counts,
                   ncclDataType_t datatype,
                   int peer,
                   void* comm,
                   void* stream,
                   ncclResult_t *ncclRes) {

        // local variables
        double* ptr=(double*) vec;
        ncclComm_t* nccl_comm= (ncclComm_t*) comm;
        hipStream_t* hipstream= (hipStream_t*) stream;

        *ncclRes=ncclSend(ptr,counts,datatype,peer,*nccl_comm,*hipstream);

        rccl_check(*ncclRes); // abort if error
    }

    void nccl_recv(void* vec,
                   int counts,
                   ncclDataType_t datatype,
                   int peer,
                   void* comm,
                   void* stream,
                   ncclResult_t *ncclRes) {

        // local variables
        double* ptr=(double*) vec;
        ncclComm_t* nccl_comm= (ncclComm_t*) comm;
        hipStream_t* hipstream= (hipStream_t*) stream;

        *ncclRes=ncclRecv(ptr,counts,datatype,peer,*nccl_comm,*hipstream);

        rccl_check(*ncclRes); // abort if error
    }

    void nccl_allgather(void* vec1,
                        void* vec2,
                        int counts,
                        ncclDataType_t datatype,
                        void* comm,
                        void* stream,
                        ncclResult_t *ncclRes) {

        // local variables
        double* ptr1=(double*) vec1;
        double* ptr2=(double*) vec2;
        ncclComm_t* nccl_comm= (ncclComm_t*) comm;
        hipStream_t* hipstream= (hipStream_t*) stream;

        *ncclRes=ncclAllGather(ptr1,ptr2,counts,datatype,*nccl_comm,*hipstream);

        rccl_check(*ncclRes); // abort if error
    }

    void nccl_allreduce(void* vec1,
                        void* vec2,
                        int counts,
                        ncclDataType_t datatype,
                        ncclRedOp_t op,
                        void* comm,
                        void* stream,
                        ncclResult_t *ncclRes) {

        // local variables
        double* ptr1=(double*) vec1;
        double* ptr2=(double*) vec2;
        ncclComm_t* nccl_comm= (ncclComm_t*) comm;
        hipStream_t* hipstream= (hipStream_t*) stream;

        *ncclRes=ncclAllReduce(ptr1,ptr2,counts,datatype,op,*nccl_comm,*hipstream);

        rccl_check(*ncclRes); // abort if error
    }
}
