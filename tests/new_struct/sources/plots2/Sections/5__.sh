#!/bin/bash

# полный путь до скрипта
ABSOLUTE_FILENAME=`readlink -e "$0"`
# каталог в котором лежит скрипт
DIRECTORY=`dirname "$ABSOLUTE_FILENAME"`

stringEnd=".gpl__OUT"
stringS1="$DIRECTORY/$1$stringEnd"

if [ "$2" = "1" ]
then
gnuplot << EOF
	set term png enhanced size 1600,1200
	set style line 1 lt 2 lw 3 lc rgb 'red'
	set style line 2 lt 2 lw 3 lc rgb 'blue'
	set style line 3 lt 2 lw 2 lc rgb 'blue'
	set style line 4 lt 2 lw 4 lc rgb 'red'
	set output "$DIRECTORY/errors/$1.png"
	set size 1.0, 1.0				#соотношение сторон общего графика
	set origin 0.0, 0.0				#изменит начало координат графика. x и y даны в "экранных" координатах.

	set multiplot					#вход в режим multiplot
		#  plot №1
		set size 1.0,0.8
#		set yrange [ 0.0 : 1.0 ]
		set origin 0.0,0.2
		set title "abs"
		set xlabel " Y" font "Helvetica,18"
		set bmargin 4
		set ylabel "Abs" font "Helvetica,18" offset 1,0
		set key top right
		set grid
		plot "$stringS1" u 2:5 ls 1 title "abs"
		
		#  plot №2
		set size 0.3,0.3
#		set yrange [ 0.0 : 400.0 ]
		set origin 0.35,0.55
		set title "%"
		set xlabel " Y"
		set ylabel "Otn" offset 1,0
		set key top right
		set grid
		plot "$stringS1" u 2:6 ls 2 title "otn"
	
	unset multiplot
#==============================================================
EOF
elif [ "$2" = "2" ]
then
gnuplot << EOF
	set term png enhanced size 1600,1200
	set style line 1 lt 2 lw 3 lc rgb 'red'
	set style line 2 lt 2 lw 3 lc rgb 'blue'
	set style line 3 lt 2 lw 2 lc rgb 'blue'
	set style line 4 lt 2 lw 4 lc rgb 'red'
	set output "$DIRECTORY/errors/$1.png"
	set size 1.0, 1.0				#соотношение сторон общего графика
	set origin 0.0, 0.0				#изменит начало координат графика. x и y даны в "экранных" координатах.

	set multiplot					#вход в режим multiplot
		#  plot №1
		set size 1.0,0.8
#		set yrange [ 0.0 : 1.0 ]
		set origin 0.0,0.2
		set title "abs"
		set xlabel " X" font "Helvetica,18"
		set bmargin 4
		set ylabel "Abs" font "Helvetica,18" offset 1,0
		set key top right
		set grid
		plot "$stringS1" u 1:5 ls 1 title "abs"
		
		#  plot №2
		set size 0.3,0.3
#		set yrange [ 0.0 : 350.0 ]
		set origin 0.35,0.55
		set title "%"
		set xlabel " X"
		set ylabel "Otn" offset 1,0
		set key top right
		set grid
		plot "$stringS1" u 1:6 ls 2 title "otn"
	
	unset multiplot
#==============================================================
EOF
elif [ "$2" = "3" ]
then
gnuplot << EOF
	set term png enhanced size 1600,1200
	set style line 1 lt 2 lw 3 lc rgb 'red'
	set style line 2 lt 2 lw 3 lc rgb 'blue'
	set style line 3 lt 2 lw 2 lc rgb 'blue'
	set style line 4 lt 2 lw 4 lc rgb 'red'
	set output "$DIRECTORY/errors/$1.png"
	set size 1.0, 1.0				#соотношение сторон общего графика
	set origin 0.0, 0.0				#изменит начало координат графика. x и y даны в "экранных" координатах.

	set multiplot					#вход в режим multiplot
		#  plot №1
		set size 1.0,0.8
#		set yrange [ 0.0 : 1.0 ]
		set origin 0.0,0.2
		set title "abs"
		set xlabel " Y" font "Helvetica,18"
		set bmargin 4
		set ylabel "Abs" font "Helvetica,18" offset 1,0
		set key top right
		set grid
		plot "$stringS1" u 2:5 ls 1 title "abs"
		
		#  plot №2
		set size 0.3,0.3
#		set yrange [ 0.0 : 800.0 ]
		set origin 0.35,0.55
		set title "%"
		set xlabel " Y"
		set ylabel "Otn" offset 1,0
		set key top right
		set grid
		plot "$stringS1" u 2:6 ls 2 title "otn"
	
	unset multiplot
#==============================================================
EOF
elif [ "$2" = "4" ]
then
gnuplot << EOF
	set term png enhanced size 1600,1200
	set style line 1 lt 2 lw 3 lc rgb 'red'
	set style line 2 lt 2 lw 3 lc rgb 'blue'
	set style line 3 lt 2 lw 2 lc rgb 'blue'
	set style line 4 lt 2 lw 4 lc rgb 'red'
	set output "$DIRECTORY/errors/$1.png"
	set size 1.0, 1.0				#соотношение сторон общего графика
	set origin 0.0, 0.0				#изменит начало координат графика. x и y даны в "экранных" координатах.

	set multiplot					#вход в режим multiplot
		#  plot №1
		set size 1.0,0.8
#		set yrange [ 0.0 : 1.0 ]
		set origin 0.0,0.2
		set title "abs"
		set xlabel " X" font "Helvetica,18"
		set bmargin 4
		set ylabel "Abs" font "Helvetica,18" offset 1,0
		set key top right
		set grid
		plot "$stringS1" u 1:5 ls 1 title "abs"
		
		#  plot №2
		set size 0.3,0.3
#		set yrange [ 0.0 : 400.0 ]
		set origin 0.35,0.55
		set title "%"
		set xlabel " X"
		set ylabel "Otn" offset 1,0
		set key top right
		set grid
		plot "$stringS1" u 1:6 ls 2 title "otn"
	
	unset multiplot
#==============================================================
EOF
fi
