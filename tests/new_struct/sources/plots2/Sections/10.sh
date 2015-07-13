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
	if [ "${myline:0:24}" = "9_U.py_ComparingFile_OUT" ]
	then
		ComparingFile_OUT="${myline:25:50}"
	fi
#-----------------------------------------------------
done < "$DIRECTORY/0SETTINGS.txt"
echo "number_x = "$number_x
echo "number_y = "$number_y


#=========================================================================================


#координаты сечений функций U_x по x
number=0
while (let "number<number_x")
do
	A_U_x_x[$number]="$ComparingFile_OUT U_x when x(${Slice_x_[$number]})"
	let "number += 1"
done

#координаты сечений функций U_x по y
number=0
while (let "number<number_y")
do
	A_U_x_y[$number]="$ComparingFile_OUT U_x when y(${Slice_y_[$number]})"
	let "number += 1"
done

#координаты сечений функций U_x_gradX по x
number=0
while (let "number<number_x")
do
	A_U_x_gradX_x[$number]="$ComparingFile_OUT U_x_gradX when x(${Slice_x_[$number]})"
	let "number += 1"
done

#координаты сечений функций U_x_gradX по y
number=0
while (let "number<number_y")
do
	A_U_x_gradX_y[$number]="$ComparingFile_OUT U_x_gradX when y(${Slice_y_[$number]})"
	let "number += 1"
done

#координаты сечений функций U_x_gradY по x
number=0
while (let "number<number_x")
do
	A_U_x_gradY_x[$number]="$ComparingFile_OUT U_x_gradY when x(${Slice_x_[$number]})"
	let "number += 1"
done

#координаты сечений функций U_x_gradY по y
number=0
while (let "number<number_y")
do
	A_U_x_gradY_y[$number]="$ComparingFile_OUT U_x_gradY when y(${Slice_y_[$number]})"
	let "number += 1"
done

#координаты сечений функций U_y по x
number=0
while (let "number<number_x")
do
	A_U_y_x[$number]="$ComparingFile_OUT U_y when x(${Slice_x_[$number]})"
	let "number += 1"
done

#координаты сечений функций U_y по y
number=0
while (let "number<number_y")
do
	A_U_y_y[$number]="$ComparingFile_OUT U_y when y(${Slice_y_[$number]})"
	let "number += 1"
done

#координаты сечений функций U_y_gradY по x
number=0
while (let "number<number_x")
do
	A_U_y_gradX_x[$number]="$ComparingFile_OUT U_y_gradY when x(${Slice_x_[$number]})"
	let "number += 1"
done

#координаты сечений функций U_y_gradY по y
number=0
while (let "number<number_y")
do
	A_U_y_gradX_y[$number]="$ComparingFile_OUT U_y_gradY when y(${Slice_y_[$number]})"
	let "number += 1"
done

#координаты сечений функций U_y_gradX по x
number=0
while (let "number<number_x")
do
	A_U_y_gradY_x[$number]="$ComparingFile_OUT U_y_gradY when x(${Slice_x_[$number]})"
	let "number += 1"
done

#координаты сечений функций U_y_gradX по y
number=0
while (let "number<number_y")
do
	A_U_y_gradY_y[$number]="$ComparingFile_OUT U_y_gradY when y(${Slice_y_[$number]})"
	let "number += 1"
done


#-------------------------------------------------------------------------
#создание графика, содержащего сечения функции U_z на разных координатах x
#"адрес скрипта" параметр1 параметр2 параметр3
#	параметр1 - передача файла 1 вида "CurveGrid_Precision U_z x(0.00000000).gpl" - сечение из deal.II
#	параметр2 - параметр, описывающий класс сечения
#			параметр2 = "1" - для U_* x
#			параметр2 = "2" - для U_* y


number=0
while (let "number<number_x")
do
	"$DIRECTORY/10_.py" "${A_U_x_x[$number]}"
	"$DIRECTORY/10__.sh" "${A_U_x_x[$number]}" "1" ${Slice_x_[$number]}
	let "number += 1"
done

number=0
while (let "number<number_y")
do
	"$DIRECTORY/10_.py" "${A_U_x_y[$number]}"
	"$DIRECTORY/10__.sh" "${A_U_x_y[$number]}" "2" ${Slice_y_[$number]}
	let "number += 1"
done

number=0
while (let "number<number_x")
do
	"$DIRECTORY/10_.py" "${A_U_x_gradX_x[$number]}"
	"$DIRECTORY/10__.sh" "${A_U_x_gradX_x[$number]}" "1" ${Slice_x_[$number]}
	let "number += 1"
done

number=0
while (let "number<number_y")
do
	"$DIRECTORY/10_.py" "${A_U_x_gradX_y[$number]}"
	"$DIRECTORY/10__.sh" "${A_U_x_gradX_y[$number]}" "2" ${Slice_y_[$number]}
	let "number += 1"
done

number=0
while (let "number<number_x")
do
	"$DIRECTORY/10_.py" "${A_U_x_gradY_x[$number]}"
	"$DIRECTORY/10__.sh" "${A_U_x_gradY_x[$number]}" "1" ${Slice_x_[$number]}
	let "number += 1"
done

number=0
while (let "number<number_y")
do
	"$DIRECTORY/10_.py" "${A_U_x_gradY_y[$number]}"
	"$DIRECTORY/10__.sh" "${A_U_x_gradY_y[$number]}" "2" ${Slice_y_[$number]}
	let "number += 1"
done

number=0
while (let "number<number_x")
do
	"$DIRECTORY/10_.py" "${A_U_y_x[$number]}"
	"$DIRECTORY/10__.sh" "${A_U_y_x[$number]}" "1" ${Slice_x_[$number]}
	let "number += 1"
done

number=0
while (let "number<number_y")
do
	"$DIRECTORY/10_.py" "${A_U_y_y[$number]}"
	"$DIRECTORY/10__.sh" "${A_U_y_y[$number]}" "2" ${Slice_y_[$number]}
	let "number += 1"
done

number=0
while (let "number<number_x")
do
	"$DIRECTORY/10_.py" "${A_U_y_gradX_x[$number]}"
	"$DIRECTORY/10__.sh" "${A_U_y_gradX_x[$number]}" "1" ${Slice_x_[$number]}
	let "number += 1"
done

number=0
while (let "number<number_y")
do
	"$DIRECTORY/10_.py" "${A_U_y_gradX_y[$number]}"
	"$DIRECTORY/10__.sh" "${A_U_y_gradX_y[$number]}" "2" ${Slice_y_[$number]}
	let "number += 1"
done

number=0
while (let "number<number_x")
do
	"$DIRECTORY/10_.py" "${A_U_y_gradY_x[$number]}"
	"$DIRECTORY/10__.sh" "${A_U_y_gradY_x[$number]}" "1" ${Slice_x_[$number]}
	let "number += 1"
done

number=0
while (let "number<number_y")
do
	"$DIRECTORY/10_.py" "${A_U_y_gradY_y[$number]}"
	"$DIRECTORY/10__.sh" "${A_U_y_gradY_y[$number]}" "2" ${Slice_y_[$number]}
	let "number += 1"
done


cd "$DIRECTORY"
rm *.gpl__OUT





