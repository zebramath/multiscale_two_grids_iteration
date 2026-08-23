#!/usr/bin/env bash
set -euo pipefail

mode="${1:-quick}"
if [[ "$mode" != "quick" && "$mode" != "full" ]]; then
    echo "usage: $0 [quick|full]" >&2
    exit 2
fi

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="$project_dir/build"

cmake -S "$project_dir" -B "$build_dir" -DCMAKE_BUILD_TYPE=Release
cmake --build "$build_dir" -j
ctest --test-dir "$build_dir" --output-on-failure

cd "$project_dir"
if [[ "$mode" == "quick" ]]; then
    "$build_dir/experiment5" --threads=4 --max-cycles=4000
    "$build_dir/experiment6" --threads=4 --spectral-iters=30 --max-cycles=4000
    "$build_dir/experiment7" --threads=4 --seeds=1 --contrast=1e4 \
        --max-cycles=3000
else
    "$build_dir/generate_test_fields"
    for number in 1 2 3 4 5; do
        "$build_dir/experiment$number"
    done
    "$build_dir/experiment6" --threads=4 --spectral-iters=80 \
        --max-cycles=20000
    "$build_dir/experiment7" --threads=4 --seeds=3 \
        --max-cycles=10000
fi
