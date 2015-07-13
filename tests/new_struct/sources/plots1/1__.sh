#!/bin/bash

# полный путь до скрипта
ABSOLUTE_FILENAME=`readlink -e "$0"`
# каталог в котором лежит скрипт
DIRECTORY=`dirname "$ABSOLUTE_FILENAME"`


"$DIRECTORY/Sections/7.py"
echo "7.py - done========================================="
"$DIRECTORY/Sections/8.sh"
echo "8.sh - done========================================="
