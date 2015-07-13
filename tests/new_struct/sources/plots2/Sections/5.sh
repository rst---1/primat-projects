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
#-----------------------------------------------------Find Name of files
	if [ "${myline:0:31}" = "4_tau_zx_x.py_ComparingFile_OUT" ]
	then
		ComparingFile_OUT="${myline:32:50}"
	fi
#-----------------------------------------------------
done < "$DIRECTORY/0SETTINGS.txt"
echo "number_x = "$number_x
echo "number_y = "$number_y



#=========================================================================================


#координаты сечений функций Tau_xx по x
number=0
while (let "number<number_x")
do
	A_Tau_xx_x_[$number]="$ComparingFile_OUT tau_xx when x(${Slice_x_[$number]})"
	let "number += 1"
done

#координаты сечений функций Tau_xx по y
number=0
while (let "number<number_y")
do
	A_Tau_xx_y_[$number]="$ComparingFile_OUT tau_xx when y(${Slice_y_[$number]})"
	let "number += 1"
done

#координаты сечений функций Tau_yy по x
number=0
while (let "number<number_x")
do
	A_Tau_yy_x_[$number]="$ComparingFile_OUT tau_yy when x(${Slice_x_[$number]})"
	let "number += 1"
done

#координаты сечений функций Tau_yy по y
number=0
while (let "number<number_y")
do
	A_Tau_yy_y_[$number]="$ComparingFile_OUT tau_yy when y(${Slice_y_[$number]})"
	let "number += 1"
done

#координаты сечений функций Tau_xy по x
number=0
while (let "number<number_x")
do
	A_Tau_xy_x_[$number]="$ComparingFile_OUT tau_xy when x(${Slice_x_[$number]})"
	let "number += 1"
done

#координаты сечений функций Tau_xy по y
number=0
while (let "number<number_y")
do
	A_Tau_xy_y_[$number]="$ComparingFile_OUT tau_xy when y(${Slice_y_[$number]})"
	let "number += 1"
done

#координаты сечений функций Tau_zz по x
number=0
while (let "number<number_x")
do
	A_Tau_zz_x_[$number]="$ComparingFile_OUT tau_zz when x(${Slice_x_[$number]})"
	let "number += 1"
done

#координаты сечений функций Tau_zz по y
number=0
while (let "number<number_y")
do
	A_Tau_zz_y_[$number]="$ComparingFile_OUT tau_zz when y(${Slice_y_[$number]})"
	let "number += 1"
done



#=========================================================================================
#Вызовы скриптов для создания графиков в gnuplot'е
#-------------------------------------------------------------------------
#создание графика, содержащего сечения функции Tau_zx на разных координатах x
#"адрес скрипта" параметр1 параметр2 параметр3
#	параметр1 - передача файла 1 вида "CurveGrid_Precision tau_zy x(0.00000000).gpl" - сечение из deal.II
#	параметр2 - параметр, описывающий класс сечения
#			параметр2 = "1" - для tau_zx x
#			параметр2 = "2" - для tau_zx y
#			параметр2 = "3" - для tau_zy x
#			параметр2 = "4" - для tau_zy y

number=0
while (let "number<number_x")
do
	"$DIRECTORY/5_.py" "${A_Tau_xx_x_[$number]}"
	"$DIRECTORY/5__.sh" "${A_Tau_xx_x_[$number]}" "1" ${Slice_x_[$number]}
	let "number += 1"
done

number=0
while (let "number<number_y")
do
	"$DIRECTORY/5_.py" "${A_Tau_xx_y_[$number]}"
	"$DIRECTORY/5__.sh" "${A_Tau_xx_y_[$number]}" "2" ${Slice_y_[$number]}
	let "number += 1"
done

number=0
while (let "number<number_x")
do
	"$DIRECTORY/5_.py" "${A_Tau_yy_x_[$number]}"
	"$DIRECTORY/5__.sh" "${A_Tau_yy_x_[$number]}" "3" ${Slice_x_[$number]}
	let "number += 1"
done

number=0
while (let "number<number_y")
do
	"$DIRECTORY/5_.py" "${A_Tau_yy_y_[$number]}"
	"$DIRECTORY/5__.sh" "${A_Tau_yy_y_[$number]}" "4" ${Slice_y_[$number]}
	let "number += 1"
done
number=0
while (let "number<number_x")
do
	"$DIRECTORY/5_.py" "${A_Tau_xy_x_[$number]}"
	"$DIRECTORY/5__.sh" "${A_Tau_xy_x_[$number]}" "1" ${Slice_x_[$number]}
	let "number += 1"
done

number=0
while (let "number<number_y")
do
	"$DIRECTORY/5_.py" "${A_Tau_xy_y_[$number]}"
	"$DIRECTORY/5__.sh" "${A_Tau_xy_y_[$number]}" "2" ${Slice_y_[$number]}
	let "number += 1"
done

number=0
while (let "number<number_x")
do
	"$DIRECTORY/5_.py" "${A_Tau_zz_x_[$number]}"
	"$DIRECTORY/5__.sh" "${A_Tau_zz_x_[$number]}" "3" ${Slice_x_[$number]}
	let "number += 1"
done

number=0
while (let "number<number_y")
do
	"$DIRECTORY/5_.py" "${A_Tau_zz_y_[$number]}"
	"$DIRECTORY/5__.sh" "${A_Tau_zz_y_[$number]}" "4" ${Slice_y_[$number]}
	let "number += 1"
done


cd "$DIRECTORY"
rm *.gpl__OUT





