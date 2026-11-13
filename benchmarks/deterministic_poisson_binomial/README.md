# Deterministic Poisson-binomial support

This supplement is based on FazeelUsmani's XSF PR #495 at `0cfb4a6b1598a4d06e83a21cdcdfeadd0b8e9f7d`. The accurate-tail implementation remains that contribution, addressing ckrfeng's SciPy #26409, with edward-lcl's XSF #336 as prior work.

The source-only commit is `09f0909a06116d1a6454061307c5d208db4602d5`. It skips certain failures and successes in the common PMF/log-PMF recurrences, translates the result in the caller's buffer, and bounds SF/log-SF tail work by feasible support. With m variable trials, recurrence work becomes O(m²+n). Fully stochastic inputs retain the original recurrence after an O(n) detection scan.

## Recorded native execution

GCC 14.2.0, C++17/-O3, dependency include checkout `0495378922713867d41659b748b05e7955a6fa9e`. Both actual headers were compiled. Tables were compared with long-double convolution for double/float, non-contiguous input, empty inputs and deterministic boundaries. Maximum observed double log-tail error was 3.64e-12. The original difficult cases retained SF 1.196353171051496e-29, log-SF -1186.4216412494752 and tiny-probability SF 1.5999999999999994e-19.

`native-comparison.json` preserves alternating ABBA process-CPU blocks including identical ctypes overhead. For the declared 1,024-trial input (32 original variable probabilities plus 496 zeros and 496 ones), median call time changed from 11,272.95 to 12.81 microseconds. This deterministic-heavy workload does not establish a general speedup for ordinary stochastic inputs. No GPU execution is claimed.

Run `bash run_native.sh /absolute/path/to/xsf/include` from this directory to compile the retained headers and rerun the native workloads. The script records the actual compiler. Existing results are from the execution described above.
