#!/usr/bin/python
SETTINGS_FILE = "/home/rst/primat-projects/tests/new_struct/sources/plotsTest/Sections/0SETTINGS.txt"
#---------------------------------------------------------
Slice_x_ = []										#SLICES
Slice_y_ = []
DIRECTORY_OF_FUNCTION = []
NAME_OF_FUNCTION = []
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

		if (line[:21] == "DIRECTORY_OF_FUNCTION"):
			DIRECTORY_OF_FUNCTION.append( line[22:] )
			NAME_OF_FUNCTION.append( line[-10:-4] )

print CSV_FILENAME_1
print DIRECTORY_OF_IN_XX
print DIRECTORY_OF_IN_YY
print DIRECTORY_OF_IN_XY
print DIRECTORY_OF_IN_ZZ
print
print DIRECTORY_OF_FUNCTION
print
print NAME_OF_FUNCTION
print


for i_Files in range(len(Slice_x_)):
	print DIRECTORY_OF_OUT + "/" + CSV_FILENAME_1 + " tau_xx when x("+str(Slice_x_[i_Files])+").csv"

print
for i_Files in range(len(Slice_x_)):
	print DIRECTORY_OF_OUT + "/" + CSV_FILENAME_1 + " " + NAME_OF_FUNCTION[0] + " when x("+str(Slice_x_[i_Files])+").csv"

print
for i_Files in range(len(Slice_x_)):
	print DIRECTORY_OF_OUT + "/" + CSV_FILENAME_1 + " " + NAME_OF_FUNCTION[1] + " when x("+str(Slice_x_[i_Files])+").csv"

print
for i_Files in range(len(Slice_x_)):
	print DIRECTORY_OF_OUT + "/" + CSV_FILENAME_1 + " " + NAME_OF_FUNCTION[2] + " when x("+str(Slice_x_[i_Files])+").csv"
