#!/usr/bin/env bash
set -euo pipefail

while IFS= read -r -d '' file; do
    printf '%s\n' "$file"
    clang-format -i "$file"
done < <(git ls-files -z -- '*.cpp' '*.hpp' '*.c' '*.h')
