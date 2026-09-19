#include <cstdlib>
//#include <iostream>
#include <level_zero/ze_api.h>

extern "C" {
    void ze_device_get_uuid(char *id, int device_num, int *ierr){
        int i, d, n;
//        // Initialize the driver
//        zeInit(0);

        // initialize ierr = -1: if no GPU "device_num" is found ierr = -1 on return.
        *ierr = -1;

        // Discover all the driver instances
        uint32_t driverCount = 0;
        zeDriverGet(&driverCount, nullptr);

        // if no driver is found set ierr = 1 and return
        if(driverCount == 0) {
            *ierr = 1;
            return;
        }

        ze_driver_handle_t* allDrivers = (ze_driver_handle_t*)std::malloc(driverCount * sizeof(ze_driver_handle_t));
        zeDriverGet(&driverCount, allDrivers);

//        std::cout << "driverCount: " << driverCount << std::endl;

        for(i = 0; i < driverCount; ++i) {
            uint32_t deviceCount = 0;
            zeDeviceGet(allDrivers[i], &deviceCount, nullptr);

            ze_device_handle_t* allDevices =(ze_device_handle_t*)std::malloc(deviceCount * sizeof(ze_device_handle_t));
            zeDeviceGet(allDrivers[i], &deviceCount, allDevices);

//            std::cout << "deviceCount: " << deviceCount << std::endl;

            for(d = 0; d < deviceCount; ++d) {
                ze_device_properties_t device_properties {};
                device_properties.stype = ZE_STRUCTURE_TYPE_DEVICE_PROPERTIES;
                zeDeviceGetProperties(allDevices[d], &device_properties);

                if(ZE_DEVICE_TYPE_GPU == device_properties.type) {
//                    std::cout << "i: " << i << " d: " << d << " name: " << device_properties.name;
//                    std::cout << " uuid:";
//                    for(const auto& e : device_properties.uuid.id) {
//                        std::cout << " " << unsigned(e);
//                    }
//                    std::cout << std::endl;

                    if(d == device_num) {
//                        std::cout << " uuid:";
//                        std::cout << "i: " << i << " d: " << d << " name: " << device_properties.name;
//                        for(const auto& e : device_properties.uuid.id) {
//                            std::cout << " " << unsigned(e);
//                        }
//                        std::cout << std::endl; 

                        // if there is more than one driver that provides a GPU "device_num"
                        // set ierr = 3 and return
                        if(*ierr == 0) {
                            *ierr = 2;
                            return;
                        }

                        // copy device_properties.uuid.id to id and set ierr = 0
                        for(n = 0; n < ZE_MAX_DEVICE_UUID_SIZE; n++) {
                            id[n] = device_properties.uuid.id[n];
                        }
                        *ierr = 0;
                    }
                }
            }
            std::free(allDevices);
        }
        free(allDrivers);
    }
}
