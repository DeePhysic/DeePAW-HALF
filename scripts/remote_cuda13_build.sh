#!/usr/bin/env bash
set -euo pipefail

remote_host="${1:-hz.icqms.group}"
remote_port="${2:-8122}"
iterations="${3:-5000}"
gpu_cc="${4:-auto}"
remote_root="${HALF_REMOTE_ROOT:-/share/home/limusen/deepaw-half}"
project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
happy_root="$(cd "${project_root}/../happy" && pwd)"

ssh -p "${remote_port}" "${remote_host}" "mkdir -p '${remote_root}/half' '${remote_root}/fixtures'"
rsync -az --delete --exclude build --exclude .pytest_cache \
  -e "ssh -p ${remote_port}" "${project_root}/" "${remote_host}:${remote_root}/half/"
rsync -az -e "ssh -p ${remote_port}" \
  "${happy_root}/examples/si/CHGCAR.smooth" \
  "${remote_host}:${remote_root}/fixtures/Si.CHGCAR.smooth"

if [[ "${gpu_cc}" == "auto" ]]; then
  gpu_cc="$(ssh -p "${remote_port}" "${remote_host}" \
    "nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d '.'")"
fi
if [[ ! "${gpu_cc}" =~ ^[0-9]+$ ]]; then
  echo "Invalid GPU compute capability: ${gpu_cc}" >&2
  exit 2
fi

ssh -p "${remote_port}" "${remote_host}" \
  "set -euo pipefail
   export HALF_REMOTE_ROOT='${remote_root}'
   '${remote_root}/half/scripts/install_cuda13_toolchain.sh'
   cuda_env='/share/home/limusen/app/anaconda/anaconda3/envs/half-cuda13'
   nvhpc_sdk='/share/home/limusen/app/nvhpc-25.9-root/opt/nvidia/hpc_sdk/Linux_x86_64/25.9'
   export PATH=\"\${cuda_env}/bin:\${nvhpc_sdk}/compilers/bin:\${PATH}\"
   export NVHPC_CUDA_HOME=\"\${cuda_env}\"
   mkdir -p '${remote_root}/toolchain-25.9'
   \"\${nvhpc_sdk}/compilers/bin/makelocalrc\" \"\${nvhpc_sdk}/compilers/bin\" \
     -x -d '${remote_root}/toolchain-25.9' -gcc /usr/bin/gcc -gpp /usr/bin/g++ -cuda 13.0
   export GCCLOCALRC='${remote_root}/toolchain-25.9/localrc'
   cd '${remote_root}/half'
   build_dir='build/cuda13-cc${gpu_cc}-release'
   echo 'HALF_GPU_CC=${gpu_cc}'
   cmake --fresh --preset cuda13-release -B \"\${build_dir}\" -DHALF_GPU_CC='${gpu_cc}'
   cmake --build \"\${build_dir}\" --parallel
   ctest --test-dir \"\${build_dir}\" --output-on-failure
   \"\${build_dir}/half-inspect\" ../fixtures/Si.CHGCAR.smooth 400
   run_dir='/tmp/limusen-half-cuda13-cc${gpu_cc}-benchmark'
   mkdir -p \"\${run_dir}\"
   cp \"\${build_dir}/half-benchmark\" \"\${run_dir}/\"
   cp ../fixtures/Si.CHGCAR.smooth \"\${run_dir}/\"
   cd \"\${run_dir}\"
   ./half-benchmark Si.CHGCAR.smooth '${iterations}'"
