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
    run_step experiment1_two_grid_comparison \
        "$build_dir/experiment1_two_grid_comparison" \
        --quick --threads="$threads"
    run_step experiment4_multilevel_comparison \
        "$build_dir/experiment4_multilevel_comparison" \
        --quick --threads="$threads"
elif [ "$mode" = "full" ]; then
    run_step experiment1_two_grid_comparison \
        "$build_dir/experiment1_two_grid_comparison" --threads="$threads"
    run_step experiment2_finite_path \
        "$build_dir/experiment2_finite_path" --threads="$threads"
    run_step experiment3_oracle_validation \
        "$build_dir/experiment3_oracle_validation" --threads="$threads"
    run_step experiment4_multilevel_comparison \
        "$build_dir/experiment4_multilevel_comparison" --threads="$threads"
else
    echo "usage: $0 [quick|full]" >&2
    exit 2
fi
