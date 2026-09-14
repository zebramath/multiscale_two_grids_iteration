#!/usr/bin/env sh
set -eu

mode="${1:-quick}"
threads="${TGI_THREADS:-4}"
build_dir="${TGI_BUILD_DIR:-build}"
step_timeout="${TGI_STEP_TIMEOUT_SECONDS:-10800}"
case "$mode" in
    quick|full) ;;
    *)
        echo "usage: $0 [quick|full]" >&2
        exit 2
        ;;
esac
case "$threads" in
    ''|*[!0-9]*|0)
        echo "TGI_THREADS must be a positive integer" >&2
        exit 2
        ;;
esac
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$script_dir/.."

if [ "$mode" = "quick" ]; then
    results_dir="${TGI_QUICK_RESULTS_DIR:-$build_dir/quick-results}"
else
    results_dir="${TGI_RESULTS_DIR:-results}"
fi

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

build_direct() {
    cxx="${CXX:-c++}"
    mkdir -p "$build_dir"
    common="-std=c++17 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror -pthread -I src"
    if "$cxx" -fopenmp -x c++ -E /dev/null >/dev/null 2>&1; then
        common="$common -fopenmp"
    fi
    for source in experiments/*.cpp; do
        program=$(basename "$source" .cpp)
        "$cxx" $common -DTGI_RESULTS_DIR=\"results\" \
            "$source" -o "$build_dir/$program"
    done
}

if command -v cmake >/dev/null 2>&1; then
    run_step configure cmake -S . -B "$build_dir" -DCMAKE_BUILD_TYPE=Release
    run_step build cmake --build "$build_dir" --parallel "$threads"
else
    echo "[info] cmake unavailable; using direct C++17 build"
    build_direct
fi

if [ "$mode" = "quick" ]; then
    run_step experiment1 env TGI_RESULTS_DIR="$results_dir" \
        "$build_dir/experiment1_finite_path_comparison" \
        --quick --threads="$threads"
    run_step experiment2 env TGI_RESULTS_DIR="$results_dir" \
        "$build_dir/experiment2_spectral_path" \
        --quick --spectral-iterations=80 --threads="$threads"
    run_step experiment3 env TGI_RESULTS_DIR="$results_dir" \
        "$build_dir/experiment3_scaling_propagation" \
        --quick --spectral-iterations=80 --threads="$threads"
    run_step experiment4-export env TGI_RESULTS_DIR="$results_dir" \
        "$build_dir/experiment4_local_diagnostic_export" \
        --maximum-steps=40 --threads="$threads"
else
    run_step experiment1 env TGI_RESULTS_DIR="$results_dir" \
        "$build_dir/experiment1_finite_path_comparison" --threads="$threads"
    run_step experiment2 env TGI_RESULTS_DIR="$results_dir" \
        "$build_dir/experiment2_spectral_path" --threads="$threads"
    run_step experiment3 env TGI_RESULTS_DIR="$results_dir" \
        "$build_dir/experiment3_scaling_propagation" --threads="$threads"
    run_step experiment4-export env TGI_RESULTS_DIR="$results_dir" \
        "$build_dir/experiment4_local_diagnostic_export" --threads="$threads"
fi

if [ "$mode" = "full" ]; then
    run_step experiment2-endpoints env TGI_RESULTS_DIR="$results_dir" \
        "$build_dir/experiment2_spectral_path" --endpoint-only \
        --spectral-iterations=200 --threads="$threads"
    run_step experiment2-summary python3 scripts/summarize_spectral_paths.py \
        --results "$results_dir"
fi

run_step experiment4-analysis env MPLCONFIGDIR="$build_dir/matplotlib" \
    python3 scripts/analyze_local_diagnostic.py --results "$results_dir"
run_step cross-spectral-plot env MPLCONFIGDIR="$build_dir/matplotlib" \
    python3 scripts/plot_spectral_path.py \
    "$results_dir/experiment2_cross_channel_spectral_path.csv" \
    "$results_dir/experiment2_cross_channel_spectral_path.png" \
    "Cross-channel spectral path"
run_step ring-spectral-plot env MPLCONFIGDIR="$build_dir/matplotlib" \
    python3 scripts/plot_spectral_path.py \
    "$results_dir/experiment2_winding_ring_spectral_path.csv" \
    "$results_dir/experiment2_winding_ring_spectral_path.png" \
    "Winding-ring spectral path"

if [ "$mode" = "full" ]; then
    run_step experiment5 env TGI_RESULTS_DIR="$results_dir" \
        "$build_dir/experiment5_robustness" --threads="$threads"
    run_step experiment6 env TGI_RESULTS_DIR="$results_dir" \
        "$build_dir/experiment6_endpoint_comparison" --threads="$threads"
    run_step experiment7 env TGI_RESULTS_DIR="$results_dir" \
        "$build_dir/experiment7_multilevel_pilot" --threads="$threads"
fi
echo "[info] results directory: $results_dir"
