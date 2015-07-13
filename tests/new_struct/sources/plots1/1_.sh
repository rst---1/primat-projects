#!/bin/bash

# полный путь до скрипта
ABSOLUTE_FILENAME=`readlink -e "$0"`
# каталог в котором лежит скрипт
DIRECTORY=`dirname "$ABSOLUTE_FILENAME"`


"$DIRECTORY/Sections/2.py"
echo "2.py - done========================================="
"$DIRECTORY/Sections/3.sh"
echo "3.sh - done========================================="

