#!/usr/bin/env sh
set -eu

mode="${1:-quick}"
threads="${TGI_THREADS:-4}"
build_dir="${TGI_BUILD_DIR:-build}"
step_timeout="${TGI_STEP_TIMEOUT_SECONDS:-900}"

run_step() {
    label="$1"
    shift
    echo "[start] $label"
    if command -v timeout >/dev/null 2>&1; then
        timeout --foreground "$step_timeout" "$@"
    else
        "$@"
    fi
    echo "[done]  $label"
}

run_step configure cmake -S . -B "$build_dir" \
    -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
run_step build cmake --build "$build_dir" --parallel "$threads"
run_step tests ctest --test-dir "$build_dir" --output-on-failure

if [ "$mode" = "quick" ]; then
    run_step experiment7 "$build_dir/experiment7" \
        --fine=32 --coarse=8 --threads="$threads" \
        --contrast=1e4 --max-cycles=12000
    run_step experiment9 "$build_dir/experiment9" --quick --threads="$threads"
    run_step experiment11 "$build_dir/experiment11" --threads="$threads"
elif [ "$mode" = "full" ]; then
    run_step experiment7 "$build_dir/experiment7" \
        --fine=128 --coarse=16 --threads="$threads" \
        --contrast=1e4 --max-cycles=40000
    run_step experiment9 "$build_dir/experiment9" --threads="$threads"
    run_step experiment11 "$build_dir/experiment11" --threads="$threads"
    run_step experiment13 "$build_dir/experiment13" --threads="$threads"
    run_step experiment14 "$build_dir/experiment14" --threads="$threads"
else
    echo "usage: $0 [quick|full]" >&2
    exit 2
fi

