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
	if [ "${myline:0:16}" = "CSV_U_FILENAME_1" ]
	then
		CSV_U_FILENAME_1="${myline:17:50}"
	fi
#-----------------------------------------------------Find Name of files
	if [ "${myline:0:16}" = "CSV_U_FILENAME_2" ]
	then
		CSV_U_FILENAME_2="${myline:17:50}"
	fi
#-----------------------------------------------------intervals
	if [ "${myline:0:26}" = "8.interval_UP___U_x_when_x" ]
	then
		UP___U_x_when_x="${myline:27:50}"
	fi
	if [ "${myline:0:26}" = "8.interval_DOWN_U_x_when_x" ]
	then
		DOWN_U_x_when_x="${myline:27:50}"
	fi

	if [ "${myline:0:26}" = "8.interval_UP___U_x_when_y" ]
	then
		UP___U_x_when_y="${myline:27:50}"
	fi
	if [ "${myline:0:26}" = "8.interval_DOWN_U_x_when_y" ]
	then
		DOWN_U_x_when_y="${myline:27:50}"
	fi

	if [ "${myline:0:32}" = "8.interval_UP___U_x_gradX_when_x" ]
	then
		UP___U_x_gradX_when_x="${myline:33:50}"
	fi
	if [ "${myline:0:32}" = "8.interval_DOWN_U_x_gradX_when_x" ]
	then
		DOWN_U_x_gradX_when_x="${myline:33:50}"
	fi

	if [ "${myline:0:32}" = "8.interval_UP___U_x_gradX_when_y" ]
	then
		UP___U_x_gradX_when_y="${myline:33:50}"
	fi
	if [ "${myline:0:32}" = "8.interval_DOWN_U_x_gradX_when_y" ]
	then
		DOWN_U_x_gradX_when_y="${myline:33:50}"
	fi

	if [ "${myline:0:32}" = "8.interval_UP___U_x_gradY_when_x" ]
	then
		UP___U_x_gradY_when_x="${myline:33:50}"
	fi
	if [ "${myline:0:32}" = "8.interval_DOWN_U_x_gradY_when_x" ]
	then
		DOWN_U_x_gradY_when_x="${myline:33:50}"
	fi

	if [ "${myline:0:32}" = "8.interval_UP___U_x_gradY_when_y" ]
	then
		UP___U_x_gradY_when_y="${myline:33:50}"
	fi
	if [ "${myline:0:32}" = "8.interval_DOWN_U_x_gradY_when_y" ]
	then
		DOWN_U_x_gradY_when_y="${myline:33:50}"
	fi


	if [ "${myline:0:26}" = "8.interval_UP___U_y_when_x" ]
	then
		UP___U_y_when_x="${myline:27:50}"
	fi
	if [ "${myline:0:26}" = "8.interval_DOWN_U_y_when_x" ]
	then
		DOWN_U_y_when_x="${myline:27:50}"
	fi

	if [ "${myline:0:26}" = "8.interval_UP___U_y_when_y" ]
	then
		UP___U_y_when_y="${myline:27:50}"
	fi
	if [ "${myline:0:26}" = "8.interval_DOWN_U_y_when_y" ]
	then
		DOWN_U_y_when_y="${myline:27:50}"
	fi

	if [ "${myline:0:32}" = "8.interval_UP___U_y_gradX_when_x" ]
	then
		UP___U_y_gradX_when_x="${myline:33:50}"
	fi
	if [ "${myline:0:32}" = "8.interval_DOWN_U_y_gradX_when_x" ]
	then
		DOWN_U_y_gradX_when_x="${myline:33:50}"
	fi

	if [ "${myline:0:32}" = "8.interval_UP___U_y_gradX_when_y" ]
	then
		UP___U_y_gradX_when_y="${myline:33:50}"
	fi
	if [ "${myline:0:32}" = "8.interval_DOWN_U_y_gradX_when_y" ]
	then
		DOWN_U_y_gradX_when_y="${myline:33:50}"
	fi

	if [ "${myline:0:32}" = "8.interval_UP___U_y_gradY_when_x" ]
	then
		UP___U_y_gradY_when_x="${myline:33:50}"
	fi
	if [ "${myline:0:32}" = "8.interval_DOWN_U_y_gradY_when_x" ]
	then
		DOWN_U_y_gradY_when_x="${myline:33:50}"
	fi

	if [ "${myline:0:32}" = "8.interval_UP___U_y_gradY_when_y" ]
	then
		UP___U_y_gradY_when_y="${myline:33:50}"
	fi
	if [ "${myline:0:32}" = "8.interval_DOWN_U_y_gradY_when_y" ]
	then
		DOWN_U_y_gradY_when_y="${myline:33:50}"
	fi
#-----------------------------------------------------
done < "$DIRECTORY/0SETTINGS.txt"
echo "number_x = "$number_x
echo "number_y = "$number_y
echo "CSV_U_FILENAME_1 = "$CSV_U_FILENAME_1
echo "CSV_U_FILENAME_2 = "$CSV_U_FILENAME_2

#=========================================================================================

#координаты сечений функций U_x по x
number=0
while (let "number<number_x")
do
	A_U_x_x[$number]="$CSV_U_FILENAME_1 U_x when x(${Slice_x_[$number]})"
	D_U_x_x[$number]="$CSV_U_FILENAME_2 U_x when x(${Slice_x_[$number]})"
	let "number += 1"
done

#координаты сечений функций U_x по y
number=0
while (let "number<number_y")
do
	A_U_x_y[$number]="$CSV_U_FILENAME_1 U_x when y(${Slice_y_[$number]})"
	D_U_x_y[$number]="$CSV_U_FILENAME_2 U_x when y(${Slice_y_[$number]})"
	let "number += 1"
done

#координаты сечений функций U_x_gradX по x
number=0
while (let "number<number_x")
do
	A_U_x_gradX_x[$number]="$CSV_U_FILENAME_1 U_x_gradX when x(${Slice_x_[$number]})"
	D_U_x_gradX_x[$number]="$CSV_U_FILENAME_2 U_x_gradX when x(${Slice_x_[$number]})"
	let "number += 1"
done

#координаты сечений функций U_x_gradX по y
number=0
while (let "number<number_y")
do
	A_U_x_gradX_y[$number]="$CSV_U_FILENAME_1 U_x_gradX when y(${Slice_y_[$number]})"
	D_U_x_gradX_y[$number]="$CSV_U_FILENAME_2 U_x_gradX when y(${Slice_y_[$number]})"
	let "number += 1"
done

#координаты сечений функций U_x_gradY по x
number=0
while (let "number<number_x")
do
	A_U_x_gradY_x[$number]="$CSV_U_FILENAME_1 U_x_gradY when x(${Slice_x_[$number]})"
	D_U_x_gradY_x[$number]="$CSV_U_FILENAME_2 U_x_gradY when x(${Slice_x_[$number]})"
	let "number += 1"
done

#координаты сечений функций U_x_gradY по y
number=0
while (let "number<number_y")
do
	A_U_x_gradY_y[$number]="$CSV_U_FILENAME_1 U_x_gradY when y(${Slice_y_[$number]})"
	D_U_x_gradY_y[$number]="$CSV_U_FILENAME_2 U_x_gradY when y(${Slice_y_[$number]})"
	let "number += 1"
done

#координаты сечений функций U_y по x
number=0
while (let "number<number_x")
do
	A_U_y_x[$number]="$CSV_U_FILENAME_1 U_y when x(${Slice_x_[$number]})"
	D_U_y_x[$number]="$CSV_U_FILENAME_2 U_y when x(${Slice_x_[$number]})"
	let "number += 1"
done

#координаты сечений функций U_y по y
number=0
while (let "number<number_y")
do
	A_U_y_y[$number]="$CSV_U_FILENAME_1 U_y when y(${Slice_y_[$number]})"
	D_U_y_y[$number]="$CSV_U_FILENAME_2 U_y when y(${Slice_y_[$number]})"
	let "number += 1"
done

#координаты сечений функций U_y_gradX по x
number=0
while (let "number<number_x")
do
	A_U_y_gradX_x[$number]="$CSV_U_FILENAME_1 U_y_gradX when x(${Slice_x_[$number]})"
	D_U_y_gradX_x[$number]="$CSV_U_FILENAME_2 U_y_gradX when x(${Slice_x_[$number]})"
	let "number += 1"
done

#координаты сечений функций U_y_gradX по y
number=0
while (let "number<number_y")
do
	A_U_y_gradX_y[$number]="$CSV_U_FILENAME_1 U_y_gradX when y(${Slice_y_[$number]})"
	D_U_y_gradX_y[$number]="$CSV_U_FILENAME_2 U_y_gradX when y(${Slice_y_[$number]})"
	let "number += 1"
done

#координаты сечений функций U_y_gradY по x
number=0
while (let "number<number_x")
do
	A_U_y_gradY_x[$number]="$CSV_U_FILENAME_1 U_y_gradY when x(${Slice_x_[$number]})"
	D_U_y_gradY_x[$number]="$CSV_U_FILENAME_2 U_y_gradY when x(${Slice_x_[$number]})"
	let "number += 1"
done

#координаты сечений функций U_y_gradY по y
number=0
while (let "number<number_y")
do
	A_U_y_gradY_y[$number]="$CSV_U_FILENAME_1 U_y_gradY when y(${Slice_y_[$number]})"
	D_U_y_gradY_y[$number]="$CSV_U_FILENAME_2 U_y_gradY when y(${Slice_y_[$number]})"
	let "number += 1"
done


#=========================================================================================
#Вызовы скриптов для создания графиков в gnuplot'е
#-------------------------------------------------------------------------
#создание графика, содержащего сечения функции U_ на разных координатах x
#"адрес скрипта" параметр1 параметр2 параметр3
#	параметр1 - передача файла 1 вида "CurveGrid_U_* x=0.00000000.gpl" - сечение из deal.II
#	параметр2 - передача файла 4 вида "CurveGrid_U_*(Analytic) x(0.0).gpl" - сечение из deal.II
#	параметр3 - параметр, описывающий класс сечения
#			параметр3 = "1" - для U_* x
#			параметр3 = "2" - для U_* y

number=0
while (let "number<number_x")
do
	"$DIRECTORY/8_.sh" "${A_U_x_x[$number]}" "${D_U_x_x[$number]}" "1" ${Slice_x_[$number]} "$DOWN_U_x_when_x" "$UP___U_x_when_x"
	let "number += 1"
done

number=0
while (let "number<number_y")
do
	"$DIRECTORY/8_.sh" "${A_U_x_y[$number]}" "${D_U_x_y[$number]}" "2" ${Slice_y_[$number]} "$DOWN_U_x_when_y" "$UP___U_x_when_y"
	let "number += 1"
done

number=0
while (let "number<number_x")
do
	"$DIRECTORY/8_.sh" "${A_U_x_gradX_x[$number]}" "${D_U_x_gradX_x[$number]}" "1" ${Slice_x_[$number]} "$DOWN_U_x_gradX_when_x" "$UP___U_x_gradX_when_x"
	let "number += 1"
done

number=0
while (let "number<number_y")
do
	"$DIRECTORY/8_.sh" "${A_U_x_gradX_y[$number]}" "${D_U_x_gradX_y[$number]}" "2" ${Slice_y_[$number]} "$DOWN_U_x_gradX_when_y" "$UP___U_x_gradX_when_y"
	let "number += 1"
done

number=0
while (let "number<number_x")
do
	"$DIRECTORY/8_.sh" "${A_U_x_gradY_x[$number]}" "${D_U_x_gradY_x[$number]}" "1" ${Slice_x_[$number]} "$DOWN_U_x_gradY_when_x" "$UP___U_x_gradY_when_x"
	let "number += 1"
done

number=0
while (let "number<number_y")
do
	"$DIRECTORY/8_.sh" "${A_U_x_gradY_y[$number]}" "${D_U_x_gradY_y[$number]}" "2" ${Slice_y_[$number]} "$DOWN_U_x_gradY_when_y" "$UP___U_x_gradY_when_y"
	let "number += 1"
done

number=0
while (let "number<number_x")
do
	"$DIRECTORY/8_.sh" "${A_U_y_x[$number]}" "${D_U_y_x[$number]}" "1" ${Slice_x_[$number]} "$DOWN_U_y_when_x" "$UP___U_y_when_x"
	let "number += 1"
done

number=0
while (let "number<number_y")
do
	"$DIRECTORY/8_.sh" "${A_U_y_y[$number]}" "${D_U_y_y[$number]}" "2" ${Slice_y_[$number]} "$DOWN_U_y_when_y" "$UP___U_y_when_y"
	let "number += 1"
done

number=0
while (let "number<number_x")
do
	"$DIRECTORY/8_.sh" "${A_U_y_gradX_x[$number]}" "${D_U_y_gradX_x[$number]}" "1" ${Slice_x_[$number]} "$DOWN_U_y_gradX_when_x" "$UP___U_y_gradX_when_x"
	let "number += 1"
done

number=0
while (let "number<number_y")
do
	"$DIRECTORY/8_.sh" "${A_U_y_gradX_y[$number]}" "${D_U_y_gradX_y[$number]}" "2" ${Slice_y_[$number]} "$DOWN_U_y_gradX_when_y" "$UP___U_y_gradX_when_y"
	let "number += 1"
done

number=0
while (let "number<number_x")
do
	"$DIRECTORY/8_.sh" "${A_U_y_gradY_x[$number]}" "${D_U_y_gradY_x[$number]}" "1" ${Slice_x_[$number]} "$DOWN_U_y_gradY_when_x" "$UP___U_y_gradY_when_x"
	let "number += 1"
done

number=0
while (let "number<number_y")
do
	"$DIRECTORY/8_.sh" "${A_U_y_gradY_y[$number]}" "${D_U_y_gradY_y[$number]}" "2" ${Slice_y_[$number]} "$DOWN_U_y_gradY_when_y" "$UP___U_y_gradY_when_y"
	let "number += 1"
done








