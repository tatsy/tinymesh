#!/bin/bash

EXTS=".cpp .hpp .c .h"
FILES=""
for e in $EXTS; do
    NEW=`git ls-files | grep $e\$`
    FILES="${FILES} ${NEW}"
done

for f in $FILES; do
    echo $f
    clang-format -i $f
done

