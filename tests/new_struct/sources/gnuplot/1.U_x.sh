#!/bin/bash
# полный путь до скрипта
ABSOLUTE_FILENAME=`readlink -e "$0"`
# каталог в котором лежит скрипт
DIRECTORY=`dirname "$ABSOLUTE_FILENAME"`


DIRECTORY1="$DIRECTORY/../out/FU_x_grad.gpl"
DIRECTORY2="$DIRECTORY/../out_analytic/FAU_x_grad.gpl"
"$DIRECTORY/U_.sh" "$DIRECTORY1" "$DIRECTORY2"

DIRECTORY1="$DIRECTORY/../out_errors/ErrFU_x_gradX.gpl"
"$DIRECTORY/err_.sh" "$DIRECTORY1"

DIRECTORY1="$DIRECTORY/../out_errors/ErrFU_x_gradY.gpl"
"$DIRECTORY/err_.sh" "$DIRECTORY1"
