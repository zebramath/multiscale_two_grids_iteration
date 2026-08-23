#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="${project_dir}/build"

cmake -S "${project_dir}" -B "${build_dir}" \
    -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build "${build_dir}" --target experiments unit_core -j
ctest --test-dir "${build_dir}" --output-on-failure

for experiment in 1 2 3 4 5; do
    "${build_dir}/experiment${experiment}" "$@"
done
