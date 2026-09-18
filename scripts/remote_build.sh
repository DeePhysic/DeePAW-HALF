#!/usr/bin/env bash
set -euo pipefail

remote_host="${1:-hz.icqms.group}"
remote_port="${2:-8122}"
build_preset="${3:-cuda12-release}"
remote_root="${HALF_REMOTE_ROOT:-/share/home/limusen/deepaw-half}"
remote_sdk="${HALF_REMOTE_SDK:-/share/app/nvidia/hpc_sdk/Linux_x86_64/24.5}"
remote_cmake="${HALF_REMOTE_CMAKE:-/share/app/general/cmake-3.31.0-linux-x86_64/bin}"
project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
happy_root="$(cd "${project_root}/../happy" && pwd)"

case "${build_preset}" in
  cpu-release|cuda12-release) localrc_cuda="12.4" ;;
  cuda13-release) localrc_cuda="13.0" ;;
  *) echo "Unsupported preset: ${build_preset}" >&2; exit 2 ;;
esac

ssh -p "${remote_port}" "${remote_host}" "mkdir -p '${remote_root}/half' '${remote_root}/fixtures'"
rsync -az --delete \
  --exclude build --exclude .pytest_cache \
  -e "ssh -p ${remote_port}" \
  "${project_root}/" "${remote_host}:${remote_root}/half/"
rsync -az -e "ssh -p ${remote_port}" \
  "${happy_root}/examples/si/CHGCAR.smooth" \
  "${remote_host}:${remote_root}/fixtures/Si.CHGCAR.smooth"

ssh -t -p "${remote_port}" "${remote_host}" \
  "cd '${remote_root}/half' \
    && export PATH='${remote_cmake}:${remote_sdk}/compilers/bin':\"\$PATH\" \
    && mkdir -p '${remote_root}/toolchain' \
    && '${remote_sdk}/compilers/bin/makelocalrc' '${remote_sdk}/compilers/bin' \
       -x -d '${remote_root}/toolchain' -gcc /usr/bin/gcc -gpp /usr/bin/g++ -cuda '${localrc_cuda}' \
    && export GCCLOCALRC='${remote_root}/toolchain/localrc' \
    && cmake --fresh --preset '${build_preset}' \
    && cmake --build --preset '${build_preset}' --parallel \
    && ctest --preset '${build_preset}' \
    && ./build/'${build_preset}'/half inspect ../fixtures/Si.CHGCAR.smooth --encut 400"
