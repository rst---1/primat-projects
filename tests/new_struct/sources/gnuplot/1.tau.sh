#!/bin/bash
# полный путь до скрипта
ABSOLUTE_FILENAME=`readlink -e "$0"`
# каталог в котором лежит скрипт
DIRECTORY=`dirname "$ABSOLUTE_FILENAME"`


DIRECTORY1="$DIRECTORY/../out/Ftau_xx.gpl"
DIRECTORY2="$DIRECTORY/../out_analytic/FAtau_xx.gpl"
"$DIRECTORY/tau_.sh" "$DIRECTORY1" "$DIRECTORY2"

DIRECTORY1="$DIRECTORY/../out/Ftau_yy.gpl"
DIRECTORY2="$DIRECTORY/../out_analytic/FAtau_yy.gpl"
"$DIRECTORY/tau_.sh" "$DIRECTORY1" "$DIRECTORY2"

DIRECTORY1="$DIRECTORY/../out/Ftau_xy.gpl"
DIRECTORY2="$DIRECTORY/../out_analytic/FAtau_xy.gpl"
"$DIRECTORY/tau_.sh" "$DIRECTORY1" "$DIRECTORY2"

DIRECTORY1="$DIRECTORY/../out/Ftau_zz.gpl"
DIRECTORY2="$DIRECTORY/../out_analytic/FAtau_zz.gpl"
"$DIRECTORY/tau_.sh" "$DIRECTORY1" "$DIRECTORY2"

DIRECTORY1="$DIRECTORY/../out/Ftau_zx.gpl"
DIRECTORY2="$DIRECTORY/../out_analytic/FAtau_zx.gpl"
"$DIRECTORY/tau_.sh" "$DIRECTORY1" "$DIRECTORY2"

DIRECTORY1="$DIRECTORY/../out/Ftau_zy.gpl"
DIRECTORY2="$DIRECTORY/../out_analytic/FAtau_zy.gpl"
"$DIRECTORY/tau_.sh" "$DIRECTORY1" "$DIRECTORY2"

