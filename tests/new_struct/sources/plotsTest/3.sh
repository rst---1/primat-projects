#!/bin/bash

# полный путь до скрипта
ABSOLUTE_FILENAME=`readlink -e "$0"`
# каталог в котором лежит скрипт
DIRECTORY=`dirname "$ABSOLUTE_FILENAME"`


"$DIRECTORY/Sections/4.py"
echo "4.py - done========================================="
"$DIRECTORY/Sections/5.sh"
echo "5.sh - done========================================="

