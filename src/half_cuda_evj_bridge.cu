#include <cuda_runtime_api.h>
#include <cusolverDn.h>

// CUDA Fortran currently does not ship a binding for cusolverDnZhegvj.
// This C ABI bridge intentionally accepts already-device-resident buffers.
extern "C" int half_cuda_zhegvj(int n, void* a, void* b, void* w,
                                 double tolerance, int max_sweeps,
                                 double* residual, int* sweeps) {
  cusolverDnHandle_t handle = nullptr;
  syevjInfo_t params = nullptr;
  cuDoubleComplex* workspace = nullptr;
  int* device_info = nullptr;
  int lwork = 0, host_info = 0;
  cusolverStatus_t status = cusolverDnCreate(&handle);
  if (status != CUSOLVER_STATUS_SUCCESS) return static_cast<int>(status);
  status = cusolverDnCreateSyevjInfo(&params);
  if (status == CUSOLVER_STATUS_SUCCESS)
    status = cusolverDnXsyevjSetTolerance(params, tolerance);
  if (status == CUSOLVER_STATUS_SUCCESS)
    status = cusolverDnXsyevjSetMaxSweeps(params, max_sweeps);
  if (status == CUSOLVER_STATUS_SUCCESS)
    status = cusolverDnZhegvj_bufferSize(handle, CUSOLVER_EIG_TYPE_1,
      CUSOLVER_EIG_MODE_NOVECTOR, CUBLAS_FILL_MODE_UPPER, n,
      static_cast<cuDoubleComplex*>(a), n, static_cast<cuDoubleComplex*>(b), n,
      static_cast<double*>(w), &lwork, params);
  if (status == CUSOLVER_STATUS_SUCCESS &&
      cudaMalloc(&workspace, sizeof(cuDoubleComplex) * static_cast<size_t>(lwork)) != cudaSuccess)
    status = CUSOLVER_STATUS_ALLOC_FAILED;
  if (status == CUSOLVER_STATUS_SUCCESS && cudaMalloc(&device_info, sizeof(int)) != cudaSuccess)
    status = CUSOLVER_STATUS_ALLOC_FAILED;
  if (status == CUSOLVER_STATUS_SUCCESS)
    status = cusolverDnZhegvj(handle, CUSOLVER_EIG_TYPE_1,
      CUSOLVER_EIG_MODE_NOVECTOR, CUBLAS_FILL_MODE_UPPER, n,
      static_cast<cuDoubleComplex*>(a), n, static_cast<cuDoubleComplex*>(b), n,
      static_cast<double*>(w), workspace, lwork, device_info, params);
  if (status == CUSOLVER_STATUS_SUCCESS &&
      cudaMemcpy(&host_info, device_info, sizeof(int), cudaMemcpyDeviceToHost) != cudaSuccess)
    status = CUSOLVER_STATUS_EXECUTION_FAILED;
  if (status == CUSOLVER_STATUS_SUCCESS && residual)
    status = cusolverDnXsyevjGetResidual(handle, params, residual);
  if (status == CUSOLVER_STATUS_SUCCESS && sweeps)
    status = cusolverDnXsyevjGetSweeps(handle, params, sweeps);
  if (workspace) cudaFree(workspace);
  if (device_info) cudaFree(device_info);
  if (params) cusolverDnDestroySyevjInfo(params);
  if (handle) cusolverDnDestroy(handle);
  if (status != CUSOLVER_STATUS_SUCCESS) return static_cast<int>(status);
  return host_info == 0 ? 0 : 10000 + host_info;
}
