#!/usr/bin/env sh
set -eu

mode="${1:-quick}"
threads="${TGI_THREADS:-4}"
build_dir="${TGI_BUILD_DIR:-build}"
step_timeout="${TGI_STEP_TIMEOUT_SECONDS:-900}"
case "$mode" in
    quick|full) ;;
    *)
        echo "usage: $0 [quick|full]" >&2
        exit 2
        ;;
esac
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$script_dir/.."

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
    # shellcheck disable=SC2086
    "$cxx" $common tests/unit_core.cpp -o "$build_dir/unit_core"
    # shellcheck disable=SC2086
    "$cxx" $common tests/regression_v56.cpp -o "$build_dir/regression_v56"
    # shellcheck disable=SC2086
    "$cxx" $common -DTGI_RESULTS_DIR=\"results\" \
        experiments/experiment1_two_grid_comparison.cpp \
        -o "$build_dir/experiment1_two_grid_comparison"
    # shellcheck disable=SC2086
    "$cxx" $common -DTGI_RESULTS_DIR=\"results\" \
        experiments/experiment2_finite_path.cpp \
        -o "$build_dir/experiment2_finite_path"
    # shellcheck disable=SC2086
    "$cxx" $common -DTGI_RESULTS_DIR=\"results\" \
        experiments/experiment3_oracle_validation.cpp \
        -o "$build_dir/experiment3_oracle_validation"
    # shellcheck disable=SC2086
    "$cxx" $common -DTGI_RESULTS_DIR=\"results\" \
        experiments/experiment4_submission_robustness.cpp \
        -o "$build_dir/experiment4_submission_robustness"
}

if command -v cmake >/dev/null 2>&1; then
    run_step configure cmake -S . -B "$build_dir" \
        -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
    run_step build cmake --build "$build_dir" --parallel "$threads"
    run_step tests ctest --test-dir "$build_dir" --output-on-failure
else
    echo "[info] cmake unavailable; using the direct C++17 build"
    echo "[start] build-direct"
    build_direct
    echo "[done]  build-direct"
    run_step unit-core "$build_dir/unit_core"
    run_step regression-v56 "$build_dir/regression_v56"
fi

if [ "$mode" = "quick" ]; then
    run_step experiment1-two-grid \
        "$build_dir/experiment1_two_grid_comparison" \
        --quick --threads="$threads"
else
    run_step experiment1-two-grid \
        "$build_dir/experiment1_two_grid_comparison" --threads="$threads"
    run_step experiment2-finite-path \
        "$build_dir/experiment2_finite_path" --threads="$threads"
    run_step experiment3-oracle \
        "$build_dir/experiment3_oracle_validation" --threads="$threads"
    run_step experiment4-submission \
        "$build_dir/experiment4_submission_robustness" --threads="$threads"
fi
