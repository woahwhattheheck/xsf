#!/usr/bin/env bash
set -euo pipefail
B220_XSF_INCLUDE="${1:?Pass the path to an XSF include directory}"
B220_CXX="${B220_CXX:-g++}"
B220_SOURCE_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$B220_SOURCE_DIR"
test -f "$B220_XSF_INCLUDE/xsf/config.h"
"$B220_CXX" --version > compiler-version.txt
for edition in baseline candidate; do
  test -f "$edition/include/xsf/stats.h"
  "$B220_CXX" -std=c++17 -O3 -I "$edition/include" -I "$B220_XSF_INCLUDE" -I "$B220_XSF_INCLUDE/xsf" evaluate.cpp -o "$edition-run"
  "./$edition-run" > "$edition.csv" 2> "$edition.log"
  "$B220_CXX" -std=c++17 -O3 -shared -fPIC -fvisibility=hidden -I "$edition/include" -I "$B220_XSF_INCLUDE" -I "$B220_XSF_INCLUDE/xsf" bridge.cpp -o "$edition.so"
done
python compare.py

