#!/bin/bash
#Работает с файлами, где есть координаты, функция, производная и производная
# полный путь до скрипта
ABSOLUTE_FILENAME=`readlink -e "$0"`
# каталог в котором лежит скрипт
DIRECTORY=`dirname "$ABSOLUTE_FILENAME"`

gnuplot -e "set xlabel 'X'; set ylabel 'Y'; splot '$1' u 2:3:6" -p
