#!/bin/bash

# полный путь до скрипта
ABSOLUTE_FILENAME=`readlink -e "$0"`
# каталог в котором лежит скрипт
DIRECTORY=`dirname "$ABSOLUTE_FILENAME"`

while read myline
do
#-----------------------------------------------------Find Directory
	if [ "${myline:0:16}" = "DIRECTORY_OF_OUT" ]
	then
		NEWDIRECTORY="${myline:17:100}"
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
	set style line 1 lt 2 lw 3 lc rgb 'red'
	set style line 2 lt 2 lw 3 lc rgb 'blue'
#	set style line 3 lt 2 lw 3 lc rgb 'blue'
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
#		set ylabel "Tau_z_x" font "Helvetica,18" offset 3
		set key top right
		set grid
		plot "$stringS2" u 3:4 lw 3 title "$2", "$stringS3" u 3:4 lw 3 title "$3", "$stringS4" u 3:4 lw 3 title "$4" , "$stringS5" u 3:4 lw 3 title "$5", "$stringS6" u 3:4 lw 3 title "$6"

#==============================================================
EOF
elif [ "$1" = "2" ]
then
	gnuplot << EOF
	set term png enhanced size 1600,1200
	set style line 1 lt 2 lw 3 lc rgb 'red'
	set style line 2 lt 2 lw 3 lc rgb 'blue'
#	set style line 3 lt 2 lw 3 lc rgb 'blue'
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
#		set ylabel "Tau_z_x" font "Helvetica,18" offset 3
		set key top right
		set grid
		plot "$stringS2" u 2:4 lw 3 title "$2", "$stringS3" u 2:4 lw 3 title "$3", "$stringS4" u 2:4 lw 3 title "$4", "$stringS5" u 2:4 lw 3 title "$5", "$stringS6" u 2:4 lw 3 title "$6"

#==============================================================
EOF
elif [ "$1" = "3" ]
then
	gnuplot << EOF
	set term png enhanced size 1600,1200
	set style line 1 lt 2 lw 3 lc rgb 'red'
	set style line 2 lt 2 lw 3 lc rgb 'blue'
#	set style line 3 lt 2 lw 3 lc rgb 'blue'
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
#		set ylabel "Tau_z_y" font "Helvetica,18" offset 3
		set key top right
		set grid
		plot "$stringS2" u 3:4 lw 3 title "$2", "$stringS3" u 3:4 lw 3 title "$3", "$stringS4" u 3:4 lw 3 title "$4", "$stringS5" u 3:4 lw 3 title "$5", "$stringS6" u 3:4 lw 3 title "$6"

#==============================================================
EOF
elif [ "$1" = "4" ]
then
	gnuplot << EOF
	set term png enhanced size 1600,1200
	set style line 1 lt 2 lw 3 lc rgb 'red'
	set style line 2 lt 2 lw 3 lc rgb 'blue'
#	set style line 3 lt 2 lw 3 lc rgb 'blue'
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
#		set ylabel "Tau_z_y" font "Helvetica,18" offset 3
		set key top right
		set grid
		plot "$stringS2" u 2:4 lw 3 title "$2", "$stringS3" u 2:4 lw 3 title "$3", "$stringS4" u 2:4 lw 3 title "$4", "$stringS5" u 2:4 lw 3 title "$5", "$stringS6" u 2:4 lw 3 title "$6"
	
#==============================================================
EOF
fi
