#!/usr/bin/pvpython
SETTINGS_FILE = "/home/rst/primat-projects/tests/new_struct/sources/plots/Sections/0SETTINGS.txt"
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
Render()

#filter "Delaunay2D"
SetActiveSource(TableToPoints_1)
Delaunay2D_1 = Delaunay2D()
Delaunay2D_1.Tolerance = TOLERANCE_OF_DELAUNAY2D
Show()
Render()

#filter "Clip"
#ClipFilter_1 = Clip()
#ClipFilter_1.ClipType = "Box"
#ClipFilter_1.ClipType.Scale = [ 0.414, 0.414, 2.0 ]
#Render()

#filter "Slice"
sliceFilter_1 = Slice()
Show(sliceFilter_1)
Render()
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
reader_2 = CSVReader(FileName = [DIRECTORY_OF_IN_U_x_A])

#filter "Table To Points"
TableToPoints_2 = TableToPoints()
TableToPoints_2.XColumn = "Points:1"
TableToPoints_2.YColumn = "Points:2"
TableToPoints_2.ZColumn = "Points:3"
Render()

#filter "Delaunay2D"
SetActiveSource(TableToPoints_2)
Delaunay2D_2 = Delaunay2D()
Delaunay2D_2.Tolerance = TOLERANCE_OF_DELAUNAY2D
Show()
Render()


#filter "Slice"
sliceFilter_2 = Slice()
Show(sliceFilter_2)
Render
#---------------------------------------------------------
for i_Files in range(len(Slice_x_)):
	sliceFilter_2.SliceType.Origin = [float(Slice_x_[i_Files]),0.0,0.0]
	sliceFilter_2.SliceType.Normal = [1.0,0.0,0.0]
	#Render()
	writer = CreateWriter(DIRECTORY_OF_OUT + "/" + CSV_U_FILENAME_2 + " U_x when x("+str(Slice_x_[i_Files])+").csv")
	writer.WriteAllTimeSteps = 1
	writer.FieldAssociation = "Points"
	writer.UpdatePipeline()
#---------------------------------------------------------
for i_Files in range(len(Slice_y_)):
	sliceFilter_2.SliceType.Origin = [0.0,float(Slice_y_[i_Files]),0.0]
	sliceFilter_2.SliceType.Normal = [0.0,1.0,0.0]
	#Render()
	writer = CreateWriter(DIRECTORY_OF_OUT + "/" + CSV_U_FILENAME_2 + " U_x when y("+str(Slice_y_[i_Files])+").csv")
	writer.WriteAllTimeSteps = 1
	writer.FieldAssociation = "Points"
	writer.UpdatePipeline()
#---------------------------------------------------------
#Delete filters2
Delete(sliceFilter_2)
#Delete(ClipFilter_2)
Delete(Delaunay2D_2)
Delete(TableToPoints_2)
Delete(reader_2)







#open file
from paraview.simple import *
reader_3 = CSVReader(FileName = [DIRECTORY_OF_IN_U_x])

#filter "Table To Points"
TableToPoints_3 = TableToPoints()
TableToPoints_3.XColumn = "Points:1"
TableToPoints_3.YColumn = "Points:2"
TableToPoints_3.ZColumn = "Points:4"
Render()

#filter "Delaunay2D"
SetActiveSource(TableToPoints_3)
Delaunay2D_3 = Delaunay2D()
Delaunay2D_3.Tolerance = TOLERANCE_OF_DELAUNAY2D
Show()
Render()


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
reader_4 = CSVReader(FileName = [DIRECTORY_OF_IN_U_x_A])

#filter "Table To Points"
TableToPoints_4 = TableToPoints()
TableToPoints_4.XColumn = "Points:1"
TableToPoints_4.YColumn = "Points:2"
TableToPoints_4.ZColumn = "Points:4"
Render()

#filter "Delaunay2D"
SetActiveSource(TableToPoints_4)
Delaunay2D_4 = Delaunay2D()
Delaunay2D_4.Tolerance = TOLERANCE_OF_DELAUNAY2D
Show()
Render()


#filter "Slice"
sliceFilter_4 = Slice()
Show(sliceFilter_4)
Render
#---------------------------------------------------------
for i_Files in range(len(Slice_x_)):
	sliceFilter_4.SliceType.Origin = [float(Slice_x_[i_Files]),0.0,0.0]
	sliceFilter_4.SliceType.Normal = [1.0,0.0,0.0]
	#Render()
	writer = CreateWriter(DIRECTORY_OF_OUT + "/" + CSV_U_FILENAME_2 + " U_x_gradX when x("+str(Slice_x_[i_Files])+").csv")
	writer.WriteAllTimeSteps = 1
	writer.FieldAssociation = "Points"
	writer.UpdatePipeline()
#---------------------------------------------------------
for i_Files in range(len(Slice_y_)):
	sliceFilter_4.SliceType.Origin = [0.0,float(Slice_y_[i_Files]),0.0]
	sliceFilter_4.SliceType.Normal = [0.0,1.0,0.0]
	#Render()
	writer = CreateWriter(DIRECTORY_OF_OUT + "/" + CSV_U_FILENAME_2 + " U_x_gradX when y("+str(Slice_y_[i_Files])+").csv")
	writer.WriteAllTimeSteps = 1
	writer.FieldAssociation = "Points"
	writer.UpdatePipeline()
#---------------------------------------------------------
#Delete filters2
Delete(sliceFilter_4)
#Delete(ClipFilter_4)
Delete(Delaunay2D_4)
Delete(TableToPoints_4)
Delete(reader_4)







#open file
from paraview.simple import *
reader_5 = CSVReader(FileName = [DIRECTORY_OF_IN_U_x])

#filter "Table To Points"
TableToPoints_5 = TableToPoints()
TableToPoints_5.XColumn = "Points:1"
TableToPoints_5.YColumn = "Points:2"
TableToPoints_5.ZColumn = "Points:5"
Render()

#filter "Delaunay2D"
SetActiveSource(TableToPoints_5)
Delaunay2D_5 = Delaunay2D()
Delaunay2D_5.Tolerance = TOLERANCE_OF_DELAUNAY2D
Show()
Render()


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







#open file
from paraview.simple import *
reader_6 = CSVReader(FileName = [DIRECTORY_OF_IN_U_x_A])

#filter "Table To Points"
TableToPoints_6 = TableToPoints()
TableToPoints_6.XColumn = "Points:1"
TableToPoints_6.YColumn = "Points:2"
TableToPoints_6.ZColumn = "Points:5"
Render()

#filter "Delaunay2D"
SetActiveSource(TableToPoints_6)
Delaunay2D_6 = Delaunay2D()
Delaunay2D_6.Tolerance = TOLERANCE_OF_DELAUNAY2D
Show()
Render()


#filter "Slice"
sliceFilter_6 = Slice()
Show(sliceFilter_6)
Render
#---------------------------------------------------------
for i_Files in range(len(Slice_x_)):
	sliceFilter_6.SliceType.Origin = [float(Slice_x_[i_Files]),0.0,0.0]
	sliceFilter_6.SliceType.Normal = [1.0,0.0,0.0]
	#Render()
	writer = CreateWriter(DIRECTORY_OF_OUT + "/" + CSV_U_FILENAME_2 + " U_x_gradY when x("+str(Slice_x_[i_Files])+").csv")
	writer.WriteAllTimeSteps = 1
	writer.FieldAssociation = "Points"
	writer.UpdatePipeline()
#---------------------------------------------------------
for i_Files in range(len(Slice_y_)):
	sliceFilter_6.SliceType.Origin = [0.0,float(Slice_y_[i_Files]),0.0]
	sliceFilter_6.SliceType.Normal = [0.0,1.0,0.0]
	#Render()
	writer = CreateWriter(DIRECTORY_OF_OUT + "/" + CSV_U_FILENAME_2 + " U_x_gradY when y("+str(Slice_y_[i_Files])+").csv")
	writer.WriteAllTimeSteps = 1
	writer.FieldAssociation = "Points"
	writer.UpdatePipeline()
#---------------------------------------------------------
#Delete filters2
Delete(sliceFilter_6)
#Delete(ClipFilter_6)
Delete(Delaunay2D_6)
Delete(TableToPoints_6)
Delete(reader_6)






