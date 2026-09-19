#include "ParallelEnvironment.hpp"

#include "Record.hpp"
#include "Tutor.hpp"

#include <algorithm> // for find
#include <stdexcept> // for runtime_error

using namespace vaspml;

MapString vaspml::dataMapParallelEnvironment()
{
    return MapString{
        {"palgo",  "String"    },
        {"flavor", "String"    },
        {"palgos", "Vec1String"}
    };
}

ParallelEnvironment::ParallelEnvironment(ShRec record) :
    data(assignOrMakeRecord(record, dataMapParallelEnvironment())),
    initialized(false),
    palgo(data->get<String>("palgo")),
    flavor(data->get<String>("flavor")),
    palgos(data->get<Vec1String>("palgos"))
#ifdef VASPML_PALGO_GPU_NVIDIA
    ,
    handle(nullptr)
#endif
{
    if (!palgos.empty()) return;

    palgo = "off";
    palgos.push_back("off");
    palgos.push_back("serial");
#ifdef VASPML_PALGO_THREADS
    palgos.push_back("threads");
#endif
#ifdef VASPML_PALGO_GPU
    palgos.push_back("gpu");
    flavor = "" VASPML_PALGO_GPU;
#endif
}

void ParallelEnvironment::init(const String palgoRequested)
{
    if (initialized)
    {
        global::tutor.warning("Reinitializing already initialized ParallelEnvironment instance.");
    }

    Vec1String palgoOff{"off", "0", "f", ".false.", "no"};
    Vec1String palgoOn{"on", "1", "t", ".true.", "yes"};
    // If requested palgo = "off", use "classic" code path (no parallel algorithms).
    if (std::find(palgoOff.begin(), palgoOff.end(), palgoRequested) != palgoOff.end())
    {
        palgo = "off";
    }
    // If requested palgo = "on", pick "best" parallelization option (last one in list).
    else if (std::find(palgoOn.begin(), palgoOn.end(), palgoRequested) != palgoOn.end())
    {
        palgo = palgos.back();
    }
    else if (std::find(palgos.begin(), palgos.end(), palgoRequested) != palgos.end())
    {
        palgo = palgoRequested;
    }
    else
    {
        throw std::runtime_error("Unable to select parallelization algorithm \"" + palgoRequested
                                 + "\", this build supports: \"" + supported() + "\". ");
    }

#ifdef VASPML_PALGO_GPU
    if (palgo == "gpu")
    {
#ifdef VASPML_PALGO_GPU_NVIDIA
        if (handle != nullptr)
        {
            global::tutor.warning("Attempting to create cuBLAS handle which was already created "
                                  "before. Old handle will be destroyed and a new one created.");
            destroyHandle();
        }
        cublasStatus_t status = cublasCreate(&handle);
        if (status != CUBLAS_STATUS_SUCCESS)
        {
            throw std::runtime_error("Could not create cuBLAS library context handle, reason: \""
                                     + convert(status) + "\".");
        }
#endif
    }
#endif // VASPML_PALGO_GPU

    initialized = true;

    return;
}

ParallelEnvironment::~ParallelEnvironment()
{
#ifdef VASPML_PALGO_GPU_NVIDIA
    if (palgo == "gpu") { destroyHandle(); }
#endif // VASPML_PALGO_GPU_NVIDIA
}

bool ParallelEnvironment::isInitialized() const
{
    return initialized;
}

bool ParallelEnvironment::off() const
{
    return palgo == "off";
}

bool ParallelEnvironment::gpu() const
{
    return palgo == "gpu";
}

String ParallelEnvironment::selected() const
{
    return palgo;
}

String ParallelEnvironment::supported() const
{
    String result;
    for (Vec1String::const_iterator p = palgos.begin(); p != palgos.end(); ++p)
    {
        result += *p;
        if (p + 1 != palgos.end()) result += ",";
    }

    return result;
}

#ifdef VASPML_PALGO_GPU_NVIDIA

String ParallelEnvironment::convert(cublasStatus_t status) const
{
    String result = "UNKNOWN";

    if (status == CUBLAS_STATUS_SUCCESS) result = "CUBLAS_STATUS_SUCCESS";
    else if (status == CUBLAS_STATUS_NOT_INITIALIZED) result = "CUBLAS_STATUS_NOT_INITIALIZED";
    else if (status == CUBLAS_STATUS_ALLOC_FAILED) result = "CUBLAS_STATUS_ALLOC_FAILED";
    else if (status == CUBLAS_STATUS_INVALID_VALUE) result = "CUBLAS_STATUS_INVALID_VALUE";
    else if (status == CUBLAS_STATUS_ARCH_MISMATCH) result = "CUBLAS_STATUS_ARCH_MISMATCH";
    else if (status == CUBLAS_STATUS_MAPPING_ERROR) result = "CUBLAS_STATUS_MAPPING_ERROR";
    else if (status == CUBLAS_STATUS_EXECUTION_FAILED) result = "CUBLAS_STATUS_EXECUTION_FAILED";
    else if (status == CUBLAS_STATUS_INTERNAL_ERROR) result = "CUBLAS_STATUS_INTERNAL_ERROR";
    else if (status == CUBLAS_STATUS_NOT_SUPPORTED) result = "CUBLAS_STATUS_NOT_SUPPORTED";
    else if (status == CUBLAS_STATUS_LICENSE_ERROR) result = "CUBLAS_STATUS_LICENSE_ERROR";

    return result;
}

const cublasHandle_t& ParallelEnvironment::get_handle() const
{
    return handle;
}

void ParallelEnvironment::destroyHandle()
{
    if (handle == nullptr)
    {
        global::tutor.warning("Attempting to destroy cuBLAS handle which was not created before.");
        return;
    }

    cublasStatus_t status = cublasDestroy(handle);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        global::tutor.warning("Got status \"" + convert(status)
                              + "\" upon attempt to destroy cuBLAS handle.");
    }
    handle = nullptr;

    return;
}

#endif // VASPML_PALGO_GPU_NVIDIA
