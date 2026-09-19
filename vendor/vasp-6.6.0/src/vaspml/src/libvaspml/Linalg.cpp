#include "Linalg.hpp"

using namespace vaspml;
using namespace vaspml::linalg;

linalg::TransposeHost linalg::toHost(Transpose transpose)
{
    if (transpose == Transpose::NOTRANS) return noTransHost;
    else return transHost;
}

#ifdef VASPML_PALGO_GPU_NVIDIA
linalg::TransposeDevice linalg::toDevice(Transpose transpose)
{
    if (transpose == Transpose::NOTRANS) return noTransDevice;
    else return transDevice;
}
#endif // Variants of VASPML_PALGO_GPU
