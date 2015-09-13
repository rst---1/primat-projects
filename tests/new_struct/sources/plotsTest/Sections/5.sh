#!/bin/bash

# полный путь до скрипта
ABSOLUTE_FILENAME=`readlink -e "$0"`
# каталог в котором лежит скрипт
DIRECTORY=`dirname "$ABSOLUTE_FILENAME"`

number_x=0
number_y=0
while read myline
do
#-----------------------------------------------------Find Slice_x
	if [ "${myline:0:7}" = "Slice_x" ]
	then
		Slice_x_[$number_x]="${myline:8:8}"
		let "number_x += 1"
	fi
#-----------------------------------------------------Find Slice_y
	if [ "${myline:0:7}" = "Slice_y" ]
	then
		Slice_y_[$number_y]="${myline:8:8}"
		let "number_y += 1"
	fi
#-----------------------------------------------------Find Directory
	if [ "${myline:0:16}" = "DIRECTORY_OF_OUT" ]
	then
		NEWDIRECTORY="${myline:17:150}"
	fi
#-----------------------------------------------------Find Name of files
	if [ "${myline:0:14}" = "CSV_FILENAME_1" ]
	then
		CSV_FILENAME_1="${myline:15:50}"
	fi
#-----------------------------------------------------intervals
	if [ "${myline:0:29}" = "3.interval_UP___tau_xx_when_x" ]
	then
		UP___tau_xx_when_x="${myline:30:50}"
	fi
	if [ "${myline:0:29}" = "3.interval_DOWN_tau_xx_when_x" ]
	then
		DOWN_tau_xx_when_x="${myline:30:50}"
	fi

	if [ "${myline:0:29}" = "3.interval_UP___tau_xx_when_y" ]
	then
		UP___tau_xx_when_y="${myline:30:50}"
	fi
	if [ "${myline:0:29}" = "3.interval_DOWN_tau_xx_when_y" ]
	then
		DOWN_tau_xx_when_y="${myline:30:50}"
	fi


	if [ "${myline:0:29}" = "3.interval_UP___tau_yy_when_x" ]
	then
		UP___tau_yy_when_x="${myline:30:50}"
	fi
	if [ "${myline:0:29}" = "3.interval_DOWN_tau_yy_when_x" ]
	then
		DOWN_tau_yy_when_x="${myline:30:50}"
	fi

	if [ "${myline:0:29}" = "3.interval_UP___tau_yy_when_y" ]
	then
		UP___tau_yy_when_y="${myline:30:50}"
	fi
	if [ "${myline:0:29}" = "3.interval_DOWN_tau_yy_when_y" ]
	then
		DOWN_tau_yy_when_y="${myline:30:50}"
	fi


	if [ "${myline:0:29}" = "3.interval_UP___tau_xy_when_x" ]
	then
		UP___tau_xy_when_x="${myline:30:50}"
	fi
	if [ "${myline:0:29}" = "3.interval_DOWN_tau_xy_when_x" ]
	then
		DOWN_tau_xy_when_x="${myline:30:50}"
	fi

	if [ "${myline:0:29}" = "3.interval_UP___tau_xy_when_y" ]
	then
		UP___tau_xy_when_y="${myline:30:50}"
	fi
	if [ "${myline:0:29}" = "3.interval_DOWN_tau_xy_when_y" ]
	then
		DOWN_tau_xy_when_y="${myline:30:50}"
	fi


	if [ "${myline:0:29}" = "3.interval_UP___tau_zz_when_x" ]
	then
		UP___tau_zz_when_x="${myline:30:50}"
	fi
	if [ "${myline:0:29}" = "3.interval_DOWN_tau_zz_when_x" ]
	then
		DOWN_tau_zz_when_x="${myline:30:50}"
	fi

	if [ "${myline:0:29}" = "3.interval_UP___tau_zz_when_y" ]
	then
		UP___tau_zz_when_y="${myline:30:50}"
	fi
	if [ "${myline:0:29}" = "3.interval_DOWN_tau_zz_when_y" ]
	then
		DOWN_tau_zz_when_y="${myline:30:50}"
	fi
#-----------------------------------------------------
done < "$DIRECTORY/0SETTINGS.txt"
echo "number_x = "$number_x
echo "number_y = "$number_y
echo "CSV_FILENAME_1 = "$CSV_FILENAME_1
echo "NEWDIRECTORY = "$NEWDIRECTORY

#=========================================================================================


#координаты сечений функций Tau_xx по x
number=0
while (let "number<number_x")
do
	A_Tau_xx_x_[$number]="$CSV_FILENAME_1 tau_xx when x(${Slice_x_[$number]})"
	echo A_Tau_xx_x_[$number] = ${A_Tau_xx_x_[$number]}
	let "number += 1"
done

#координаты сечений функций Tau_xx по y
number=0
while (let "number<number_y")
do
	A_Tau_xx_y_[$number]="$CSV_FILENAME_1 tau_xx when y(${Slice_y_[$number]})"
	let "number += 1"
done



#координаты сечений функций Tau_yy по x
number=0
while (let "number<number_x")
do
	A_Tau_yy_x_[$number]="$CSV_FILENAME_1 tau_yy when x(${Slice_x_[$number]})"
	let "number += 1"
done

#координаты сечений функций Tau_yy по y
number=0
while (let "number<number_y")
do
	A_Tau_yy_y_[$number]="$CSV_FILENAME_1 tau_yy when y(${Slice_y_[$number]})"
	let "number += 1"
done



#координаты сечений функций Tau_xy по x
number=0
while (let "number<number_x")
do
	A_Tau_xy_x_[$number]="$CSV_FILENAME_1 tau_xy when x(${Slice_x_[$number]})"
	let "number += 1"
done

#координаты сечений функций Tau_xy по y
number=0
while (let "number<number_y")
do
	A_Tau_xy_y_[$number]="$CSV_FILENAME_1 tau_xy when y(${Slice_y_[$number]})"
	let "number += 1"
done



#координаты сечений функций Tau_zz по x
number=0
while (let "number<number_x")
do
	A_Tau_zz_x_[$number]="$CSV_FILENAME_1 tau_zz when x(${Slice_x_[$number]})"
	let "number += 1"
done

#координаты сечений функций Tau_zz по y
number=0
while (let "number<number_y")
do
	A_Tau_zz_y_[$number]="$CSV_FILENAME_1 tau_zz when y(${Slice_y_[$number]})"
	let "number += 1"
done


#=========================================================================================
#Вызовы скриптов для создания графиков в gnuplot'е
#-------------------------------------------------------------------------
#создание графика, содержащего сечения функции Tau_zx на разных координатах x
#"адрес скрипта" параметр1 параметр2 параметр3
#	параметр1 - передача файла 1 вида "CurveGrid tau_zx x=0.00000000.gpl" - сечение из deal.II
#	параметр2 - передача файла 4 вида "CurveGrid_Analytic tau_zx x(0.0).gpl" - сечение из deal.II
#	параметр3 - параметр, описывающий класс сечения
#			параметр3 = "1" - для tau_zx x
#			параметр3 = "2" - для tau_zx y
#			параметр3 = "3" - для tau_zy x
#			параметр3 = "4" - для tau_zy y
#			параметр4 = передаёт координату по которой строится сечение


"$DIRECTORY/5_.sh" "1" "${A_Tau_xx_x_[0]}" "${A_Tau_xx_x_[1]}" "${A_Tau_xx_x_[2]}" "${A_Tau_xx_x_[3]}" "${A_Tau_xx_x_[4]}"

"$DIRECTORY/5_.sh" "2" "${A_Tau_xx_y_[0]}" "${A_Tau_xx_y_[1]}" "${A_Tau_xx_y_[2]}" "${A_Tau_xx_y_[3]}" "${A_Tau_xx_y_[4]}"



"$DIRECTORY/5_.sh" "3" "${A_Tau_yy_x_[0]}" "${A_Tau_yy_x_[1]}" "${A_Tau_yy_x_[2]}" "${A_Tau_yy_x_[3]}" "${A_Tau_yy_x_[4]}"

"$DIRECTORY/5_.sh" "4" "${A_Tau_yy_y_[0]}" "${A_Tau_yy_y_[1]}" "${A_Tau_yy_y_[2]}" "${A_Tau_yy_y_[3]}" "${A_Tau_yy_y_[4]}"







"$DIRECTORY/5_.sh" "1" "${A_Tau_xy_x_[0]}" "${A_Tau_xy_x_[1]}" "${A_Tau_xy_x_[2]}" "${A_Tau_xy_x_[3]}" "${A_Tau_xy_x_[4]}"

"$DIRECTORY/5_.sh" "2" "${A_Tau_xy_y_[0]}" "${A_Tau_xy_y_[1]}" "${A_Tau_xy_y_[2]}" "${A_Tau_xy_y_[3]}" "${A_Tau_xy_y_[4]}"



"$DIRECTORY/5_.sh" "3" "${A_Tau_zz_x_[0]}" "${A_Tau_zz_x_[1]}" "${A_Tau_zz_x_[2]}" "${A_Tau_zz_x_[3]}" "${A_Tau_zz_x_[4]}"

"$DIRECTORY/5_.sh" "4" "${A_Tau_zz_y_[0]}" "${A_Tau_zz_y_[1]}" "${A_Tau_zz_y_[2]}" "${A_Tau_zz_y_[3]}" "${A_Tau_zz_y_[4]}"











