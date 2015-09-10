#!/usr/bin/python
SETTINGS_FILE = "/home/rst/primat-projects/tests/new_struct/sources/plotsTest/Sections/0SETTINGS.txt"
#---------------------------------------------------------
Slice_x_ = []										#SLICES
Slice_y_ = []
DIRECTORY_OF_FUNCTION = []
NAME_OF_FUNCTION = []
NAME_OF_MESH = []
with open (SETTINGS_FILE) as fd:					#open file
	t = fd.read()									#read file
	for line in t.splitlines():						#read each line
		if (line[:7] == "Slice_x"):
			Slice_x_.append( line[8:] )
		if (line[:7] == "Slice_y"):
			Slice_y_.append( line[8:] )
		if (line[:16] == "DIRECTORY_OF_OUT"):
			DIRECTORY_OF_OUT = line[17:]
		if (line[:20] == "DIRECTORY_OF_Ftau_xx"):
			DIRECTORY_OF_IN_XX = line[21:]
		if (line[:20] == "DIRECTORY_OF_Ftau_yy"):
			DIRECTORY_OF_IN_YY = line[21:]
		if (line[:20] == "DIRECTORY_OF_Ftau_xy"):
			DIRECTORY_OF_IN_XY = line[21:]
		if (line[:20] == "DIRECTORY_OF_Ftau_zz"):
			DIRECTORY_OF_IN_ZZ = line[21:]
		if (line[:14] == "CSV_FILENAME_1"):
			CSV_FILENAME_1 = line[15:]

		if (line[:21] == "DIRECTORY_OF_MOUTTTTT"):
			DIRECTORY_OF_MOUTTTTT = line[22:]
		if (line[:21] == "DIRECTORY_OF_FUNCTION"):
			DIRECTORY_OF_FUNCTION.append( line[22:] )
			NAME_OF_FUNCTION.append( line[-10:-4] )
			NAME_OF_MESH.append( line[-41:-12] )
#			NAME_OF_PROBLEM.append( line[-57:-43] )



print
print

print 
print NAME_OF_FUNCTION

print
print NAME_OF_MESH





