#!/bin/bash

# полный путь до скрипта
ABSOLUTE_FILENAME=`readlink -e "$0"`
# каталог в котором лежит скрипт
DIRECTORY=`dirname "$ABSOLUTE_FILENAME"`

while read myline
do
#-----------------------------------------------------Find Directory
	if [ "${myline:0:30}" = "1_pvpython.py_DIRECTORY_OF_OUT" ]
	then
		NEWDIRECTORY="${myline:31:100}"
	fi
#-----------------------------------------------------
done < "$DIRECTORY/0SETTINGS.txt"

stringEnd=".gpl"
stringS2="$NEWDIRECTORY/$2$stringEnd"
stringS3="$NEWDIRECTORY/$3$stringEnd"
stringS4="$NEWDIRECTORY/$4$stringEnd"
stringS5="$NEWDIRECTORY/$5$stringEnd"
stringS6="$NEWDIRECTORY/$6$stringEnd"

if [ "$1" = "1" ]
then
	gnuplot << EOF
	set term png enhanced size 1600,1200
	set style line 1 lt 2 lw 2 lc rgb 'red'
	set style line 2 lt 2 lw 2 lc rgb 'blue'
#	set style line 3 lt 2 lw 2 lc rgb 'blue'
#	set style line 4 lt 2 lw 4 lc rgb 'red'
	set output "$NEWDIRECTORY/$2.png"
	set size 1.0, 1.0				#соотношение сторон общего графика
	set origin 0.0, 0.0				#изменит начало координат графика. x и y даны в "экранных" координатах.

		#  plot №1
		set size 1.0,0.8
#		set yrange [ $5 : $6 ]
		set origin 0.0,0.2
		set xlabel " Y" font "Helvetica,18"
		set bmargin 4
#		set ylabel "U_z" font "Helvetica,18" offset 3
		set key top right
		set grid
		plot "$stringS2" u 5:6 lw 3 title "$2", "$stringS3" u 5:6 lw 3 title "$3", "$stringS4" u 5:6 lw 3 title "$4" , "$stringS5" u 5:6 lw 3 title "$5", "$stringS6" u 5:6 lw 3 title "$6"

#==============================================================
EOF
elif [ "$1" = "2" ]
then
	gnuplot << EOF
	set term png enhanced size 1600,1200
	set style line 1 lt 2 lw 2 lc rgb 'red'
	set style line 2 lt 2 lw 2 lc rgb 'blue'
#	set style line 3 lt 2 lw 2 lc rgb 'blue'
#	set style line 4 lt 2 lw 4 lc rgb 'red'
	set output "$NEWDIRECTORY/$2.png"
	set size 1.0, 1.0				#соотношение сторон общего графика
	set origin 0.0, 0.0				#изменит начало координат графика. x и y даны в "экранных" координатах.

		#  plot №1
		set size 1.0,0.8
#		set yrange [ $5 : $6 ]
		set origin 0.0,0.2
		set xlabel " X" font "Helvetica,18"
		set bmargin 4
#		set ylabel "U_z" font "Helvetica,18" offset 3
		set key top right
		set grid
		plot "$stringS2" u 4:6 lw 3 title "$2", "$stringS3" u 4:6 lw 3 title "$3", "$stringS4" u 4:6 lw 3 title "$4" , "$stringS5" u 4:6 lw 3 title "$5", "$stringS6" u 4:6 lw 3 title "$6"

#==============================================================
EOF
fi
