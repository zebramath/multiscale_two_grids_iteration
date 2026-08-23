#!/usr/bin/env sh
set -eu

mode="${1:-quick}"
threads="${TGI_THREADS:-4}"
build_dir="${TGI_BUILD_DIR:-build}"

cmake -S . -B "$build_dir" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build "$build_dir" --parallel "$threads"
ctest --test-dir "$build_dir" --output-on-failure

if [ "$mode" = "quick" ]; then
    "$build_dir/experiment7" --fine=32 --coarse=8 --threads="$threads" \
        --contrast=1e4 --max-cycles=12000
    "$build_dir/experiment8" --fine=32 --coarse=8 --threads="$threads" \
        --contrast=1e4 --max-cycles=12000
    "$build_dir/experiment9" --quick --threads="$threads"
elif [ "$mode" = "full" ]; then
    "$build_dir/experiment7" --fine=128 --coarse=16 --threads="$threads" \
        --contrast=1e4 --max-cycles=40000
    "$build_dir/experiment8" --fine=64 --coarse=8 --threads="$threads" \
        --contrast=1e4 --max-cycles=20000
    "$build_dir/experiment9" --threads="$threads"
else
    echo "usage: $0 [quick|full]" >&2
    exit 2
fi
