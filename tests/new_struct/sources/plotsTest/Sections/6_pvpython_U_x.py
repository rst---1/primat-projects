#!/usr/bin/pvpython
SETTINGS_FILE = "/home/rst/primat-projects/tests/new_struct/sources/plotsTest/Sections/0SETTINGS.txt"
#---------------------------------------------------------
Slice_x_ = []										#SLICES
Slice_y_ = []
with open (SETTINGS_FILE) as fd:					#open file
	t = fd.read()									#read file
	for line in t.splitlines():						#read each line
		if (line[:7] == "Slice_x"):
			Slice_x_.append( line[8:] )
		if (line[:7] == "Slice_y"):
			Slice_y_.append( line[8:] )
		if (line[:30] == "6_pvpython.py_DIRECTORY_OF_OUT"):
			DIRECTORY_OF_OUT = line[31:]
		if (line[:37] == "6_pvpython.py_TOLERANCE_OF_DELAUNAY2D"):
			TOLERANCE_OF_DELAUNAY2D = float(line[38:])
		if (line[:36] == "6_pvpython.py_DIRECTORY_OF_FU_x_grad"):
			DIRECTORY_OF_IN_U_x = line[37:]
		if (line[:37] == "6_pvpython.py_DIRECTORY_OF_FAU_x_grad"):
			DIRECTORY_OF_IN_U_x_A = line[38:]
		if (line[:36] == "6_pvpython.py_DIRECTORY_OF_FU_y_grad"):
			DIRECTORY_OF_IN_U_y = line[37:]
		if (line[:37] == "6_pvpython.py_DIRECTORY_OF_FAU_y_grad"):
			DIRECTORY_OF_IN_U_y_A = line[38:]
		if (line[:16] == "CSV_U_FILENAME_1"):
			CSV_U_FILENAME_1 = line[17:]
		if (line[:16] == "CSV_U_FILENAME_2"):
			CSV_U_FILENAME_2 = line[17:]
#---------------------------------------------------------
#DIRECTORY_OF_OUT = '/home/rst/deal.II/rjkz_3/test_of_grid/Sections'
#TOLERANCE_OF_DELAUNAY2D = 0.0
#DIRECTORY_OF_IN_1 = "/home/rst/deal.II/rjkz_3/test_of_grid/FU_z.csv"


#open file
from paraview.simple import *
reader_1 = CSVReader(FileName = [DIRECTORY_OF_IN_U_x])

#filter "Table To Points"
TableToPoints_1 = TableToPoints()
TableToPoints_1.XColumn = "Points:1"
TableToPoints_1.YColumn = "Points:2"
TableToPoints_1.ZColumn = "Points:3"
#Render()

#filter "Delaunay2D"
SetActiveSource(TableToPoints_1)
Delaunay2D_1 = Delaunay2D()
Delaunay2D_1.Tolerance = TOLERANCE_OF_DELAUNAY2D
Delaunay2D_1.ProjectionPlaneMode = "Best-Fitting Plane"
#Show()
#Render()

#filter "Clip"
#ClipFilter_1 = Clip()
#ClipFilter_1.ClipType = "Box"
#ClipFilter_1.ClipType.Scale = [ 0.414, 0.414, 2.0 ]
#Render()

#filter "Slice"
sliceFilter_1 = Slice()
Show(sliceFilter_1)
#Render()
#---------------------------------------------------------
for i_Files in range(len(Slice_x_)):
	sliceFilter_1.SliceType.Origin = [float(Slice_x_[i_Files]),0.0,0.0]
	sliceFilter_1.SliceType.Normal = [1.0,0.0,0.0]
	#Render()
	writer = CreateWriter(DIRECTORY_OF_OUT + "/" + CSV_U_FILENAME_1 + " U_x when x("+str(Slice_x_[i_Files])+").csv")
	writer.WriteAllTimeSteps = 1
	writer.FieldAssociation = "Points"
	writer.UpdatePipeline()
#---------------------------------------------------------
for i_Files in range(len(Slice_y_)):
	sliceFilter_1.SliceType.Origin = [0.0,float(Slice_y_[i_Files]),0.0]
	sliceFilter_1.SliceType.Normal = [0.0,1.0,0.0]
	#Render()
	writer = CreateWriter(DIRECTORY_OF_OUT + "/" + CSV_U_FILENAME_1 + " U_x when y("+str(Slice_y_[i_Files])+").csv")
	writer.WriteAllTimeSteps = 1
	writer.FieldAssociation = "Points"
	writer.UpdatePipeline()
#---------------------------------------------------------
#Delete filters1
Delete(sliceFilter_1)
#Delete(ClipFilter_1)
Delete(Delaunay2D_1)
Delete(TableToPoints_1)
Delete(reader_1)
























#open file
from paraview.simple import *
reader_3 = CSVReader(FileName = [DIRECTORY_OF_IN_U_x])

#filter "Table To Points"
TableToPoints_3 = TableToPoints()
TableToPoints_3.XColumn = "Points:1"
TableToPoints_3.YColumn = "Points:2"
TableToPoints_3.ZColumn = "Points:4"
#Render()

#filter "Delaunay2D"
SetActiveSource(TableToPoints_3)
Delaunay2D_3 = Delaunay2D()
Delaunay2D_3.Tolerance = TOLERANCE_OF_DELAUNAY2D
Delaunay2D_3.ProjectionPlaneMode = "Best-Fitting Plane"
#Show()
#Render()


#filter "Slice"
sliceFilter_3 = Slice()
Show(sliceFilter_3)
Render
#---------------------------------------------------------
for i_Files in range(len(Slice_x_)):
	sliceFilter_3.SliceType.Origin = [float(Slice_x_[i_Files]),0.0,0.0]
	sliceFilter_3.SliceType.Normal = [1.0,0.0,0.0]
	#Render()
	writer = CreateWriter(DIRECTORY_OF_OUT + "/" + CSV_U_FILENAME_1 + " U_x_gradX when x("+str(Slice_x_[i_Files])+").csv")
	writer.WriteAllTimeSteps = 1
	writer.FieldAssociation = "Points"
	writer.UpdatePipeline()
#---------------------------------------------------------
for i_Files in range(len(Slice_y_)):
	sliceFilter_3.SliceType.Origin = [0.0,float(Slice_y_[i_Files]),0.0]
	sliceFilter_3.SliceType.Normal = [0.0,1.0,0.0]
	#Render()
	writer = CreateWriter(DIRECTORY_OF_OUT + "/" + CSV_U_FILENAME_1 + " U_x_gradX when y("+str(Slice_y_[i_Files])+").csv")
	writer.WriteAllTimeSteps = 1
	writer.FieldAssociation = "Points"
	writer.UpdatePipeline()
#---------------------------------------------------------
#Delete filters2
Delete(sliceFilter_3)
#Delete(ClipFilter_3)
Delete(Delaunay2D_3)
Delete(TableToPoints_3)
Delete(reader_3)
























#open file
from paraview.simple import *
reader_5 = CSVReader(FileName = [DIRECTORY_OF_IN_U_x])

#filter "Table To Points"
TableToPoints_5 = TableToPoints()
TableToPoints_5.XColumn = "Points:1"
TableToPoints_5.YColumn = "Points:2"
TableToPoints_5.ZColumn = "Points:5"
#Render()

#filter "Delaunay2D"
SetActiveSource(TableToPoints_5)
Delaunay2D_5 = Delaunay2D()
Delaunay2D_5.Tolerance = TOLERANCE_OF_DELAUNAY2D
Delaunay2D_5.ProjectionPlaneMode = "Best-Fitting Plane"
#Show()
#Render()


#filter "Slice"
sliceFilter_5 = Slice()
Show(sliceFilter_5)
Render
#---------------------------------------------------------
for i_Files in range(len(Slice_x_)):
	sliceFilter_5.SliceType.Origin = [float(Slice_x_[i_Files]),0.0,0.0]
	sliceFilter_5.SliceType.Normal = [1.0,0.0,0.0]
	#Render()
	writer = CreateWriter(DIRECTORY_OF_OUT + "/" + CSV_U_FILENAME_1 + " U_x_gradY when x("+str(Slice_x_[i_Files])+").csv")
	writer.WriteAllTimeSteps = 1
	writer.FieldAssociation = "Points"
	writer.UpdatePipeline()
#---------------------------------------------------------
for i_Files in range(len(Slice_y_)):
	sliceFilter_5.SliceType.Origin = [0.0,float(Slice_y_[i_Files]),0.0]
	sliceFilter_5.SliceType.Normal = [0.0,1.0,0.0]
	#Render()
	writer = CreateWriter(DIRECTORY_OF_OUT + "/" + CSV_U_FILENAME_1 + " U_x_gradY when y("+str(Slice_y_[i_Files])+").csv")
	writer.WriteAllTimeSteps = 1
	writer.FieldAssociation = "Points"
	writer.UpdatePipeline()
#---------------------------------------------------------
#Delete filters2
Delete(sliceFilter_5)
#Delete(ClipFilter_5)
Delete(Delaunay2D_5)
Delete(TableToPoints_5)
Delete(reader_5)




















