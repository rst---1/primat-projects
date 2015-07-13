#!/bin/bash
# полный путь до скрипта
ABSOLUTE_FILENAME=`readlink -e "$0"`
# каталог в котором лежит скрипт
DIRECTORY=`dirname "$ABSOLUTE_FILENAME"`


DIRECTORY1="$DIRECTORY/../out_errors/ErrFtau_xx.gpl"
"$DIRECTORY/err_.sh" "$DIRECTORY1"

DIRECTORY1="$DIRECTORY/../out_errors/ErrFtau_yy.gpl"
"$DIRECTORY/err_.sh" "$DIRECTORY1"

DIRECTORY1="$DIRECTORY/../out_errors/ErrFtau_xy.gpl"
"$DIRECTORY/err_.sh" "$DIRECTORY1"

DIRECTORY1="$DIRECTORY/../out_errors/ErrFtau_zz.gpl"
"$DIRECTORY/err_.sh" "$DIRECTORY1"

DIRECTORY1="$DIRECTORY/../out_errors/ErrFtau_zy.gpl"
"$DIRECTORY/err_.sh" "$DIRECTORY1"

DIRECTORY1="$DIRECTORY/../out_errors/ErrFtau_zx.gpl"
"$DIRECTORY/err_.sh" "$DIRECTORY1"
