#!/bin/bash
#Работает с файлами, где есть координаты, функция
# полный путь до скрипта
ABSOLUTE_FILENAME=`readlink -e "$0"`
# каталог в котором лежит скрипт
DIRECTORY=`dirname "$ABSOLUTE_FILENAME"`

gnuplot -e "set xlabel 'X'; set ylabel 'Y'; splot '$1' u 2:3:4,  '$2' u 2:3:4" -p
