#!/usr/bin/pvpython
#open file
from paraview.simple import *
reader_1 = CSVReader(FileName = "/home/rst/primat-projects/tests/new_struct/sources/out/T2.2/1x1/Fcsv.csv")
TableToPoints_1 = TableToPoints()
TableToPoints_1.XColumn = "XPoints:0"
TableToPoints_1.YColumn = "YPoints:1"
TableToPoints_1.ZColumn = "ZPoints:2"
Show()
Render()


from paraview.simple import *
reader_2 = CSVReader(FileName = "/home/rst/primat-projects/tests/new_struct/sources/out/T2.2/1x1_slayer3_curve_A0.1_T1/Fcsv.csv")
TableToPoints_2 = TableToPoints()
TableToPoints_2.XColumn = "XPoints:0"
TableToPoints_2.YColumn = "YPoints:1"
TableToPoints_2.ZColumn = "ZPoints:2"
Show()
Render()
