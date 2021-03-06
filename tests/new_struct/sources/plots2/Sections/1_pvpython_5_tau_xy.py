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
		if (line[:30] == "1_pvpython.py_DIRECTORY_OF_OUT"):
			DIRECTORY_OF_OUT = line[31:]
		if (line[:37] == "1_pvpython.py_TOLERANCE_OF_DELAUNAY2D"):
			TOLERANCE_OF_DELAUNAY2D = float(line[38:])
		if (line[:34] == "1_pvpython.py_DIRECTORY_OF_Ftau_xx"):
			DIRECTORY_OF_IN_XX = line[35:]
		if (line[:35] == "1_pvpython.py_DIRECTORY_OF_FAtau_xx"):
			DIRECTORY_OF_IN_XX_A = line[36:]
		if (line[:34] == "1_pvpython.py_DIRECTORY_OF_Ftau_yy"):
			DIRECTORY_OF_IN_YY = line[35:]
		if (line[:35] == "1_pvpython.py_DIRECTORY_OF_FAtau_yy"):
			DIRECTORY_OF_IN_YY_A = line[36:]
		if (line[:34] == "1_pvpython.py_DIRECTORY_OF_Ftau_xy"):
			DIRECTORY_OF_IN_XY = line[35:]
		if (line[:35] == "1_pvpython.py_DIRECTORY_OF_FAtau_xy"):
			DIRECTORY_OF_IN_XY_A = line[36:]
		if (line[:34] == "1_pvpython.py_DIRECTORY_OF_Ftau_zz"):
			DIRECTORY_OF_IN_ZZ = line[35:]
		if (line[:35] == "1_pvpython.py_DIRECTORY_OF_FAtau_zz"):
			DIRECTORY_OF_IN_ZZ_A = line[36:]
		if (line[:14] == "CSV_FILENAME_1"):
			CSV_FILENAME_1 = line[15:]
		if (line[:14] == "CSV_FILENAME_2"):
			CSV_FILENAME_2 = line[15:]
#---------------------------------------------------------

#open file
from paraview.simple import *
reader_5 = CSVReader(FileName = [DIRECTORY_OF_IN_XY])

#filter "Table To Points"
TableToPoints_5 = TableToPoints()
TableToPoints_5.XColumn = "Points:1"
TableToPoints_5.YColumn = "Points:2"
TableToPoints_5.ZColumn = "Points:3"
Render()

#filter "Delaunay2D"
SetActiveSource(TableToPoints_5)
Delaunay2D_5 = Delaunay2D()
Delaunay2D_5.Tolerance = TOLERANCE_OF_DELAUNAY2D
Show()
Render()

#=======================================================================


#filter "Slice"
sliceFilter_5 = Slice()
Show(sliceFilter_5)
Render()

#---------------------------------------------------------
for i_Files in range(len(Slice_x_)):
	sliceFilter_5.SliceType.Origin = [float(Slice_x_[i_Files]),0.0,0.0]
	sliceFilter_5.SliceType.Normal = [1.0,0.0,0.0]
	Render()
	writer = CreateWriter(DIRECTORY_OF_OUT + "/" + CSV_FILENAME_1 + " tau_xy when x("+str(Slice_x_[i_Files])+").csv")
	writer.WriteAllTimeSteps = 1
	writer.FieldAssociation = "Points"
	writer.UpdatePipeline()
#---------------------------------------------------------
for i_Files in range(len(Slice_y_)):
	sliceFilter_5.SliceType.Origin = [0.0,float(Slice_y_[i_Files]),0.0]
	sliceFilter_5.SliceType.Normal = [0.0,1.0,0.0]
	Render()
	writer = CreateWriter(DIRECTORY_OF_OUT + "/" + CSV_FILENAME_1 + " tau_xy when y("+str(Slice_y_[i_Files])+").csv")
	writer.WriteAllTimeSteps = 1
	writer.FieldAssociation = "Points"
	writer.UpdatePipeline()
#---------------------------------------------------------
#Delete filters5
Delete(sliceFilter_5)
Delete(Delaunay2D_5)
Delete(TableToPoints_5)
Delete(reader_5)










