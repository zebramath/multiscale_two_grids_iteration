#!/usr/bin/env sh
set -eu

mode="${1:-quick}"
threads="${TGI_THREADS:-4}"
build_dir="${TGI_BUILD_DIR:-build}"
step_timeout="${TGI_STEP_TIMEOUT_SECONDS:-900}"
case "$mode" in
    quick|supplemental|multilevel|full) ;;
    *)
        echo "usage: $0 [quick|supplemental|multilevel|full]" >&2
        exit 2
        ;;
esac
case "$threads" in
    ''|*[!0-9]*|0)
        echo "TGI_THREADS must be a positive integer" >&2
        exit 2
        ;;
esac
case "$step_timeout" in
    ''|*[!0-9]*|0)
        echo "TGI_STEP_TIMEOUT_SECONDS must be a positive integer" >&2
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
    for warning in -Wshadow -Wconversion -Wsign-conversion \
                   -Wduplicated-cond -Wduplicated-branches -Wlogical-op \
                   -Wnull-dereference -Wformat=2; do
        if "$cxx" "$warning" -x c++ -fsyntax-only /dev/null \
                >/dev/null 2>&1; then
            common="$common $warning"
        fi
    done
    if "$cxx" -fopenmp -x c++ -E /dev/null >/dev/null 2>&1; then
        common="$common -fopenmp"
    fi
    "$cxx" $common -DTGI_RESULTS_DIR=\"results\" \
        experiments/experiment1_two_grid_comparison.cpp \
        -o "$build_dir/experiment1_two_grid_comparison"
    "$cxx" $common -DTGI_RESULTS_DIR=\"results\" \
        experiments/experiment2_step_scan.cpp \
        -o "$build_dir/experiment2_step_scan"
    "$cxx" $common -DTGI_RESULTS_DIR=\"results\" \
        experiments/experiment3_oracle_validation.cpp \
        -o "$build_dir/experiment3_oracle_validation"
    "$cxx" $common -DTGI_RESULTS_DIR=\"results\" \
        experiments/experiment4_submission_robustness.cpp \
        -o "$build_dir/experiment4_submission_robustness"
    "$cxx" $common -DTGI_RESULTS_DIR=\"results\" \
        experiments/experiment5_stopping_ablation.cpp \
        -o "$build_dir/experiment5_stopping_ablation"
    "$cxx" $common -DTGI_RESULTS_DIR=\"results\" \
        experiments/experiment6_fixed_physical_refinement.cpp \
        -o "$build_dir/experiment6_fixed_physical_refinement"
    "$cxx" $common -DTGI_RESULTS_DIR=\"results\" \
        experiments/experiment7_multilevel_pilot.cpp \
        -o "$build_dir/experiment7_multilevel_pilot"
}

if command -v cmake >/dev/null 2>&1; then
    run_step configure cmake -S . -B "$build_dir" \
        -DCMAKE_BUILD_TYPE=Release
    run_step build cmake --build "$build_dir" --parallel "$threads"
else
    echo "[info] cmake unavailable; using the direct C++17 build"
    echo "[start] build-direct"
    build_direct
    echo "[done]  build-direct"
fi

if [ "$mode" = "quick" ]; then
    run_step experiment1-two-grid \
        env TGI_RESULTS_DIR="$results_dir" \
        "$build_dir/experiment1_two_grid_comparison" \
        --quick --threads="$threads"
elif [ "$mode" = "supplemental" ]; then
    run_step experiment5-stopping-ablation \
        env TGI_RESULTS_DIR="$results_dir" \
        "$build_dir/experiment5_stopping_ablation" --threads="$threads"
    run_step experiment6-fixed-physical \
        env TGI_RESULTS_DIR="$results_dir" \
        "$build_dir/experiment6_fixed_physical_refinement" --threads="$threads"
elif [ "$mode" = "multilevel" ]; then
    run_step experiment7-multilevel \
        env TGI_RESULTS_DIR="$results_dir" \
        "$build_dir/experiment7_multilevel_pilot" --threads="$threads"
else
    run_step experiment1-two-grid \
        env TGI_RESULTS_DIR="$results_dir" \
        "$build_dir/experiment1_two_grid_comparison" --threads="$threads"
    run_step experiment2-step-scan \
        env TGI_RESULTS_DIR="$results_dir" \
        "$build_dir/experiment2_step_scan" --threads="$threads"
    run_step experiment2-cross-plot \
        env MPLCONFIGDIR="$build_dir/matplotlib" \
        python3 scripts/plot_path_scan.py \
        "$results_dir/experiment2_cross_channel_path.csv" \
        "$results_dir/experiment2_cross_channel_path.png" \
        "Cross-channel: 128/16, contrast 1e4" 43
    run_step experiment2-ring-plot \
        env MPLCONFIGDIR="$build_dir/matplotlib" \
        python3 scripts/plot_path_scan.py \
        "$results_dir/experiment2_winding_ring_path.csv" \
        "$results_dir/experiment2_winding_ring_path.png" \
        "Winding-ring: 128/16, contrast 1e4" 43
    run_step experiment3-oracle \
        env TGI_RESULTS_DIR="$results_dir" \
        "$build_dir/experiment3_oracle_validation" --threads="$threads"
    run_step experiment4-submission \
        env TGI_RESULTS_DIR="$results_dir" \
        "$build_dir/experiment4_submission_robustness" --threads="$threads"
    run_step experiment5-stopping-ablation \
        env TGI_RESULTS_DIR="$results_dir" \
        "$build_dir/experiment5_stopping_ablation" --threads="$threads"
    run_step experiment6-fixed-physical \
        env TGI_RESULTS_DIR="$results_dir" \
        "$build_dir/experiment6_fixed_physical_refinement" --threads="$threads"
    run_step experiment7-multilevel \
        env TGI_RESULTS_DIR="$results_dir" \
        "$build_dir/experiment7_multilevel_pilot" --threads="$threads"
fi
echo "[info] results directory: $results_dir"
