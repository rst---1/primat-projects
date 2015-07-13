#!/bin/bash
# полный путь до скрипта
ABSOLUTE_FILENAME=`readlink -e "$0"`
# каталог в котором лежит скрипт
DIRECTORY=`dirname "$ABSOLUTE_FILENAME"`


DIRECTORY1="$DIRECTORY/../out_errors/ErrFU_x.gpl"
"$DIRECTORY/err_.sh" "$DIRECTORY1"

DIRECTORY1="$DIRECTORY/../out_errors/ErrFU_y.gpl"
"$DIRECTORY/err_.sh" "$DIRECTORY1"

DIRECTORY1="$DIRECTORY/../out_errors/ErrFU_z.gpl"
"$DIRECTORY/err_.sh" "$DIRECTORY1"
