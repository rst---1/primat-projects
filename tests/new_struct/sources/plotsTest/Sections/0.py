Slice_x=-0.491
Slice_x=-0.333
Slice_x=-0.25
Slice_x=-0.166
Slice_x=0.0
Slice_x=0.166
Slice_x=0.25
Slice_x=0.333
Slice_x=0.491

Slice_y=-3.37
Slice_y=-2.94318
Slice_y=-2.52272
Slice_y=-2.10227
Slice_y=-1.68182
Slice_y=-1.26136
Slice_y=-0.84091
Slice_y=-0.42045
Slice_y=0.0
Slice_y=0.42045
Slice_y=0.84091
Slice_y=1.26136
Slice_y=1.68182
Slice_y=2.10227
Slice_y=2.52272
Slice_y=2.94318
Slice_y=3.37

#показывает путь, куда поместить файлы сечений
1_pvpython.py_DIRECTORY_OF_OUT=/home/rst/deal.II/rjkz_3/test_of_grid/Sections

#показывает точность построениея триангуляции
1_pvpython.py_TOLERANCE_OF_DELAUNAY2D=0.0

#Путь к исходному CSV-файлу поверхности, из которой нужно сделать графики сечений
1_pvpython.py_DIRECTORY_OF_IN_1=/home/rst/deal.II/rjkz_3/test_of_grid/Ftau_3696el.csv



#начало для имени файла,который содержит данные сечения вида 1 (поверхность tau)
CSV_FILENAME_1=CurveGrid_Dig

#начало для имени файла,который содержит данные сечения вида 2 (поверхность tau)
CSV_FILENAME_2=CurveGrid_Analytic



#начало для имени 1го файла для сравнения
4_tau_zx_x.py_ComparingFile_1=Curve_Analytic

#начало для имени 2го файла для сравнения
4_tau_zx_x.py_ComparingFile_2=Curve_Dig

#начало для имени файла для сохранения результата сравнения
4_tau_zx_x.py_ComparingFile_OUT=Curve_Precision




#показывает путь, куда поместить файлы сечений
6_pvpython.py_DIRECTORY_OF_OUT=/home/rst/deal.II/rjkz_3/test_of_grid/Sections

#показывает точность построениея триангуляции
6_pvpython.py_TOLERANCE_OF_DELAUNAY2D=0.0

#Путь к исходному CSV-файлу поверхности, из которой нужно сделать графики сечений
6_pvpython.py_DIRECTORY_OF_IN=/home/rst/deal.II/rjkz_3/test_of_grid/FU_z.csv




#начало для имени файла,который содержит данные сечения вида 1 (поверхность U_z)
CSV_U_z_FILENAME_1=Curve_Dig

#начало для имени файла,который содержит данные сечения вида 2 (поверхность U_z)
CSV_U_z_FILENAME_2=Curve_Analytic





#начало для имени 1го файла для сравнения
9_U_z_x.py_ComparingFile_1=Curve_Analytic

#начало для имени 2го файла для сравнения
9_U_z_x.py_ComparingFile_2=Curve_Dig

#начало для имени файла для сохранения результата сравнения
9_U_z_x.py_ComparingFile_OUT=Curve_Precision





#5__.sh
stringEnd=".gpl__OUT"
*Precision*







