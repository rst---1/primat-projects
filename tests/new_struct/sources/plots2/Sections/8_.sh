#!/bin/bash

# полный путь до скрипта
ABSOLUTE_FILENAME=`readlink -e "$0"`
# каталог в котором лежит скрипт
DIRECTORY=`dirname "$ABSOLUTE_FILENAME"`

stringEnd=".gpl"
stringS1="$DIRECTORY/$1$stringEnd"
stringS2="$DIRECTORY/$2$stringEnd"

if [ "$3" = "1" ]
then
	gnuplot << EOF
	set term png enhanced size 1600,1200
	set style line 1 lt 2 lw 2 lc rgb 'red'
	set style line 2 lt 2 lw 2 lc rgb 'blue'
#	set style line 3 lt 2 lw 2 lc rgb 'blue'
#	set style line 4 lt 2 lw 4 lc rgb 'red'
	set output "$DIRECTORY/png/$1.png"
	set size 1.0, 1.0				#соотношение сторон общего графика
	set origin 0.0, 0.0				#изменит начало координат графика. x и y даны в "экранных" координатах.

		#  plot №1
		set size 1.0,0.8
#		set yrange [ -0.125 : 0.125 ]
		set yrange [ $5 : $6 ]
		set origin 0.0,0.2
		set xlabel " Y" font "Helvetica,18"
		set bmargin 4
#		set ylabel "U_z" font "Helvetica,18" offset 3
		set key top right
		set grid
		plot "$stringS1" u 5:6 title "$1", "$stringS2" u 5:6 title "$2"

#==============================================================
EOF
elif [ "$3" = "2" ]
then
	gnuplot << EOF
	set term png enhanced size 1600,1200
	set style line 1 lt 2 lw 2 lc rgb 'red'
	set style line 2 lt 2 lw 2 lc rgb 'blue'
#	set style line 3 lt 2 lw 2 lc rgb 'blue'
#	set style line 4 lt 2 lw 4 lc rgb 'red'
	set output "$DIRECTORY/png/$1.png"
	set size 1.0, 1.0				#соотношение сторон общего графика
	set origin 0.0, 0.0				#изменит начало координат графика. x и y даны в "экранных" координатах.

		#  plot №1
		set size 1.0,0.8
#		set yrange [ -0.125 : 0.125 ]
		set yrange [ $5 : $6 ]
		set origin 0.0,0.2
		set xlabel " X" font "Helvetica,18"
		set bmargin 4
#		set ylabel "U_z" font "Helvetica,18" offset 3
		set key top right
		set grid
		plot "$stringS1" u 4:6 title "$1", "$stringS2" u 4:6 title "$2"

#==============================================================
EOF
fi
