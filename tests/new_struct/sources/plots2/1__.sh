#!/bin/bash

# полный путь до скрипта
ABSOLUTE_FILENAME=`readlink -e "$0"`
# каталог в котором лежит скрипт
DIRECTORY=`dirname "$ABSOLUTE_FILENAME"`


"$DIRECTORY/Sections/7.py"
echo "7.py - done========================================="
"$DIRECTORY/Sections/8.sh"
echo "8.sh - done========================================="
"$DIRECTORY/Sections/9_U_x_gradX_x.py"
echo "9_U_x_gradX_x.py - done============================="
"$DIRECTORY/Sections/9_U_x_gradX_y.py"
echo "9_U_x_gradX_y.py - done============================="
"$DIRECTORY/Sections/9_U_x_gradY_x.py"
echo "9_U_x_gradX_x.py - done============================="
"$DIRECTORY/Sections/9_U_x_gradY_y.py"
echo "9_U_x_gradX_y.py - done============================="
"$DIRECTORY/Sections/9_U_x_x.py"
echo "9_U_x_x.py - done==================================="
"$DIRECTORY/Sections/9_U_x_y.py"
echo "9_U_x_x.py - done==================================="
"$DIRECTORY/Sections/9_U_y_gradX_x.py"
echo "9_U_y_gradX_x.py - done============================="
"$DIRECTORY/Sections/9_U_y_gradX_y.py"
echo "9_U_y_gradX_y.py - done============================="
"$DIRECTORY/Sections/9_U_y_gradY_x.py"
echo "9_U_y_gradX_x.py - done============================="
"$DIRECTORY/Sections/9_U_y_gradY_y.py"
echo "9_U_y_gradX_y.py - done============================="
"$DIRECTORY/Sections/9_U_y_x.py"
echo "9_U_y_x.py - done==================================="
"$DIRECTORY/Sections/9_U_y_y.py"
echo "9_U_y_x.py - done==================================="
"$DIRECTORY/Sections/10.sh"
echo "10.sh - done========================================"

