#!/usr/bin/env bash
set -euo pipefail

remote_root="${HALF_REMOTE_ROOT:-/share/home/limusen/deepaw-half}"
mamba="${HALF_MAMBA:-/share/home/limusen/app/anaconda/anaconda3/bin/mamba}"
cuda_env="${HALF_CUDA13_ENV:-/share/home/limusen/app/anaconda/anaconda3/envs/half-cuda13}"
nvhpc_root="${HALF_NVHPC25_ROOT:-/share/home/limusen/app/nvhpc-25.9-root}"
nvhpc_rpm_url="${HALF_NVHPC25_URL:-https://developer.download.nvidia.com/hpc-sdk/rhel/x86_64/nvhpc-25-9-25.9-1.x86_64.rpm}"
nvhpc_rpm="${HALF_NVHPC25_RPM:-/share/home/limusen/app/installers/nvhpc-25.9.rpm}"
nvhpc_rpm_sha256="${HALF_NVHPC25_SHA256:-71e5b2ee972000a445898027066b0bce7e7155447f6d26c39ef38d7de499c4a7}"
rpm_db="${remote_root}/toolchains/rpmdb"
key_dir="${remote_root}/toolchains/keys"

if [[ ! -x "${mamba}" ]]; then
  echo "mamba was not found at ${mamba}" >&2
  exit 2
fi

if [[ ! -x "${cuda_env}/bin/nvcc" ]]; then
  "${mamba}" create -y -p "${cuda_env}" -c nvidia -c conda-forge \
    "cuda-version=13.0" "cuda-compiler=13.0.3" \
    "cuda-cudart-dev=13.0.*" "cuda-driver-dev=13.0.*" \
    "gfortran_linux-64=11.*" cmake ninja
fi

if [[ ! -f "${nvhpc_root}/.half-install-complete" ]]; then
  mkdir -p "$(dirname "${nvhpc_root}")"
  mkdir -p "$(dirname "${nvhpc_rpm}")" "${key_dir}"
  curl --fail --location --continue-at - --output "${nvhpc_rpm}" "${nvhpc_rpm_url}"
  printf '%s  %s\n' "${nvhpc_rpm_sha256}" "${nvhpc_rpm}" | sha256sum --check -
  mkdir -p "${rpm_db}"
  curl --fail --location --output "${key_dir}/RPM-GPG-KEY-NVIDIA-HPC-SDK" \
    https://developer.download.nvidia.com/hpc-sdk/rhel/RPM-GPG-KEY-NVIDIA-HPC-SDK
  curl --fail --location --output "${key_dir}/RPM-GPG-KEY-NVIDIA-HPC-SDK-2022" \
    https://developer.download.nvidia.com/hpc-sdk/rhel/RPM-GPG-KEY-NVIDIA-HPC-SDK-2022
  rpm --dbpath "${rpm_db}" --initdb
  rpm --dbpath "${rpm_db}" --import \
    "${key_dir}/RPM-GPG-KEY-NVIDIA-HPC-SDK" "${key_dir}/RPM-GPG-KEY-NVIDIA-HPC-SDK-2022"
  rpm --dbpath "${rpm_db}" --checksig "${nvhpc_rpm}"
  mkdir -p "${nvhpc_root}"
  rpm2cpio "${nvhpc_rpm}" | (cd "${nvhpc_root}" && cpio --extract --make-directories --preserve-modification-time --quiet)
  touch "${nvhpc_root}/.half-install-complete"
fi

nvhpc_sdk="${nvhpc_root}/opt/nvidia/hpc_sdk/Linux_x86_64/25.9"
nvfortran="${nvhpc_sdk}/compilers/bin/nvfortran"
if [[ ! -x "${nvfortran}" ]]; then
  echo "nvfortran 25.9 was not found after extraction: ${nvfortran}" >&2
  exit 3
fi

export PATH="${cuda_env}/bin:${nvhpc_sdk}/compilers/bin:${PATH}"
export NVHPC_CUDA_HOME="${cuda_env}"

"${cuda_env}/bin/nvcc" --version
"${nvfortran}" --version
printf 'HALF_CUDA13_ENV=%s\n' "${cuda_env}"
printf 'HALF_NVHPC25_SDK=%s\n' "${nvhpc_sdk}"
