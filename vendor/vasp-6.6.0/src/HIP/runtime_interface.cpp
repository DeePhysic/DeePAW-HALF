#include <hip_runtime.h>
#include <iostream>

#define HIP_ERR(stat)                                             \
    if(stat != hipSuccess)                                        \
    {                                                             \
        std::cerr << "HIP error: " << hipGetErrorString(stat) <<  \
        " in file " << __FILE__ << ":" << __LINE__ << std::endl;  \
        exit(-1);                                                 \
    }                                                             \

extern "C" {
    void hip_init(int device_num){
        HIP_ERR(hipInit(0));
        HIP_ERR(hipSetDevice(device_num));
    }
    void hip_set_device(int devid) {
        HIP_ERR(hipSetDevice(devid));
    }
    void hip_device_synchronize() {
        HIP_ERR(hipDeviceSynchronize());
    }
    void hip_device_get_uuid(char *id, int device_num){
        hipUUID_t uuid;
        hipDevice_t device = device_num;
        HIP_ERR(hipDeviceGetUuid(&uuid, device));
        // copy data into id memory (only a pointer)
        // otherwise the values are not visible to fortran
        memcpy(id, uuid.bytes, 16);
    }
    void hip_mem_get_info(size_t *free, size_t *total, int *err){
        hipError_t ierr = hipMemGetInfo(free, total);
        *err = ierr;
    }
    void hip_stream_create(void *stream) {
        hipStream_t* hipstream= (hipStream_t*) stream;
        HIP_ERR(hipStreamCreate(hipstream));
    }
    void hip_stream_destroy(void *stream) {
        hipStream_t* hipstream= (hipStream_t*) stream;
        HIP_ERR(hipStreamDestroy(*hipstream));
    }
    void hip_stream_synchronize(void *stream) {
        hipStream_t* hipstream= (hipStream_t*) stream;
        HIP_ERR(hipStreamSynchronize(*hipstream));
    }
    void hip_wait(){
        // wait for all operations to finish on the current device
        HIP_ERR(hipDeviceSynchronize());
    }
}
