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
    run_step experiment3_adaptive_oracle \
        "$build_dir/experiment3_adaptive_oracle" \
        --fine=32 --coarse=8 --threads="$threads" \
        --contrast=1e4 --max-cycles=12000
    run_step experiment4_robustness \
        "$build_dir/experiment4_robustness" --quick --threads="$threads"
elif [ "$mode" = "full" ]; then
    run_step experiment1_localization \
        "$build_dir/experiment1_localization" --threads="$threads"
    run_step experiment2_pcg_path \
        "$build_dir/experiment2_pcg_path" --threads="$threads"
    run_step experiment3_adaptive_oracle \
        "$build_dir/experiment3_adaptive_oracle" --threads="$threads"
    run_step experiment4_robustness \
        "$build_dir/experiment4_robustness" --threads="$threads"
    run_step experiment5_diagnostics \
        "$build_dir/experiment5_diagnostics" --threads="$threads"
    run_step experiment6_workload \
        "$build_dir/experiment6_workload" --threads="$threads"
else
    echo "usage: $0 [quick|full]" >&2
    exit 2
fi
