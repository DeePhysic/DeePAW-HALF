#!/usr/bin/env bash
set -euo pipefail

remote_host="${1:-hz.icqms.group}"
remote_port="${2:-8122}"
iterations="${3:-5000}"
project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

"${project_root}/scripts/remote_build.sh" "${remote_host}" "${remote_port}" cpu-release
"${project_root}/scripts/remote_build.sh" "${remote_host}" "${remote_port}" cuda12-release

ssh -p "${remote_port}" "${remote_host}" \
  "/share/home/limusen/deepaw-half/half/build/cuda12-release/half-benchmark \
   /share/home/limusen/deepaw-half/fixtures/Si.CHGCAR.smooth '${iterations}'"
