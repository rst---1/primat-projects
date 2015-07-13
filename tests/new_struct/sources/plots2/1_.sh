#!/bin/bash

# полный путь до скрипта
ABSOLUTE_FILENAME=`readlink -e "$0"`
# каталог в котором лежит скрипт
DIRECTORY=`dirname "$ABSOLUTE_FILENAME"`


"$DIRECTORY/Sections/2.py"
echo "2.py - done========================================="
"$DIRECTORY/Sections/3.sh"
echo "3.sh - done========================================="
"$DIRECTORY/Sections/4_tau_xx_x.py"
echo "4_tau_xx_x.py - done================================"
"$DIRECTORY/Sections/4_tau_xx_y.py"
echo "4_tau_xx_y.py - done================================"
"$DIRECTORY/Sections/4_tau_xy_x.py"
echo "4_tau_xy_x.py - done================================"
"$DIRECTORY/Sections/4_tau_xy_y.py"
echo "4_tau_xy_y.py - done================================"
"$DIRECTORY/Sections/4_tau_yy_x.py"
echo "4_tau_yy_x.py - done================================"
"$DIRECTORY/Sections/4_tau_yy_y.py"
echo "4_tau_yy_y.py - done================================"
"$DIRECTORY/Sections/4_tau_zz_x.py"
echo "4_tau_zz_x.py - done================================"
"$DIRECTORY/Sections/4_tau_zz_y.py"
echo "4_tau_zz_y.py - done================================"
"$DIRECTORY/Sections/5.sh"
echo "5.sh - done========================================="

