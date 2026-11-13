"""Interleave both actual C++ kernels within one process, measuring CPU time."""
from pathlib import Path
from ctypes import CDLL, POINTER, c_double, c_size_t
import json, math, statistics, time
ROOT = Path(__file__).resolve().parent
kernels = {}
for edition in ("baseline", "candidate"):
    lib = CDLL(str(ROOT / (edition + ".so")))
    lib.run_logsf.argtypes = [POINTER(c_double), c_size_t, POINTER(c_double)]
    lib.run_logsf.restype = None
    kernels[edition] = (lib, lib.run_logsf)
first = [(1 + i % 9) / 128 for i in range(32)]
workloads = [
    ("issue-32", first),
    ("issue-256", [1 / 128] * 256),
    ("issue-tiny-16", [1e-20] * 16),
    ("issue-32-with-deterministic-trials", first + [0.0] * 496 + [1.0] * 496),
    ("all-zero", [0.0] * 1024),
    ("all-one", [1.0] * 1024),
    ("empty", []),
]
compiler_record = ROOT / "compiler-version.txt"
report = {"clock": "process_time_ns",
          "compiler": compiler_record.read_text().splitlines()[0] if compiler_record.exists() else None,
          "build_flags": "-std=c++17 -O3 -shared -fPIC -fvisibility=hidden (run_native.sh)",
          "original_xsf_pr": "https://github.com/scipy/xsf/pull/495",
          "original_head": "0cfb4a6b1598a4d06e83a21cdcdfeadd0b8e9f7d", "cases": []}
for name, probabilities in workloads:
    data = (c_double * len(probabilities))(*probabilities)
    outputs = {edition: (c_double * (len(probabilities) + 1))() for edition in kernels}
    samples = {edition: [] for edition in kernels}
    # Actual work is identical; ABBA ordering balances short-term host changes.
    for edition in ("baseline", "candidate"):
        kernels[edition][1](data, len(probabilities), outputs[edition])
    max_difference = 0.0
    for original, candidate in zip(outputs["baseline"], outputs["candidate"]):
        if math.isfinite(original) and math.isfinite(candidate):
            max_difference = max(max_difference, abs(original - candidate))
        elif original != candidate:
            raise RuntimeError(f"{name}: nonfinite output disagreement")
    if max_difference > 2e-10:
        raise RuntimeError(f"{name}: output changed by {max_difference}")
    repeats = 20 if len(probabilities) >= 128 else 200
    for edition in ("baseline", "candidate", "candidate", "baseline") * 2:
        start = time.process_time_ns()
        for _ in range(repeats):
            kernels[edition][1](data, len(probabilities), outputs[edition])
        elapsed = time.process_time_ns() - start
        samples[edition].append(elapsed / repeats / 1000)
    medians = {edition: statistics.median(values) for edition, values in samples.items()}
    report["cases"].append({
        "name": name, "n": len(probabilities), "repeats": repeats,
        "cpu_microseconds": medians, "measurements": samples,
        "ratio_baseline_over_candidate": medians["baseline"] / medians["candidate"],
        "max_logsf_difference": max_difference,
    })
(ROOT / "native-comparison.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
