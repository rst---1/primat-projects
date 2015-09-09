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
		if (line[:16] == "DIRECTORY_OF_OUT"):
			DIRECTORY_OF_OUT = line[17:]
		if (line[:20] == "1_pvpython.py_DIRECTORY_OF_Ftau_xx"):
			DIRECTORY_OF_IN_XX = line[21:]
		if (line[:20] == "1_pvpython.py_DIRECTORY_OF_Ftau_yy"):
			DIRECTORY_OF_IN_YY = line[21:]
		if (line[:20] == "1_pvpython.py_DIRECTORY_OF_Ftau_xy"):
			DIRECTORY_OF_IN_XY = line[21:]
		if (line[:20] == "1_pvpython.py_DIRECTORY_OF_Ftau_zz"):
			DIRECTORY_OF_IN_ZZ = line[21:]
		if (line[:14] == "CSV_FILENAME_1"):
			CSV_FILENAME_1 = line[15:]
#---------------------------------------------------------

#open file
from paraview.simple import *
reader_1 = CSVReader(FileName = [DIRECTORY_OF_IN_XX])

#filter "Table To Points"
TableToPoints_1 = TableToPoints()
TableToPoints_1.XColumn = "Points:1"
TableToPoints_1.YColumn = "Points:2"
TableToPoints_1.ZColumn = "Points:3"
#Render()
Render()

#filter "Delaunay2D"
SetActiveSource(TableToPoints_1)
Delaunay2D_1 = Delaunay2D()
Delaunay2D_1.Tolerance = 0.0
Delaunay2D_1.ProjectionPlaneMode = "Best-Fitting Plane"
#Show()
#Render()

#filter "Clip"
#ClipFilter_1_1 = Clip()
#ClipFilter_1_1.ClipType = "Sphere"
#ClipFilter_1_1.ClipType.Center = [ x_01, y_01, 0.2 ]
#ClipFilter_1_1.ClipType.Radius = RADIUS_OF_BIG_SPHERE


#==================================================================================


#filter "Slice"
sliceFilter_1 = Slice()
Show(sliceFilter_1)
#Render()

#---------------------------------------------------------
for i_Files in range(len(Slice_x_)):
	sliceFilter_1.SliceType.Origin = [float(Slice_x_[i_Files]),0.0,0.0]
	sliceFilter_1.SliceType.Normal = [1.0,0.0,0.0]
#	Render()
	writer = CreateWriter(DIRECTORY_OF_OUT + "/" + CSV_FILENAME_1 + " tau_xx when x("+str(Slice_x_[i_Files])+").csv")
	writer.WriteAllTimeSteps = 1
	writer.FieldAssociation = "Points"
	writer.UpdatePipeline()
#---------------------------------------------------------
for i_Files in range(len(Slice_y_)):
	sliceFilter_1.SliceType.Origin = [0.0,float(Slice_y_[i_Files]),0.0]
	sliceFilter_1.SliceType.Normal = [0.0,1.0,0.0]
	Render()
	writer = CreateWriter(DIRECTORY_OF_OUT + "/" + CSV_FILENAME_1 + " tau_xx when y("+str(Slice_y_[i_Files])+").csv")
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



