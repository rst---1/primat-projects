#!/usr/bin/python

import sys
DIRECTORY = sys.path[0]

#---------------------------------------------------------
with open (DIRECTORY + "/Sections/0SETTINGS.txt") as fd:		#open file
	t = fd.read()									#read file
	for line in t.splitlines():						#read each line
		if (line[:30] == "1_pvpython.py_DIRECTORY_OF_OUT"):
			NEWDIRECTORY = line[31:]
			print NEWDIRECTORY
#---------------------------------------------------------


FileName = sys.argv[1]
String1 = "0.330000"
String2 = "0.660000"

#============delete strings of file========================================
print '.	.	.	.	.	.	delete strings of file . .'

rst = []
with open (NEWDIRECTORY + "/" + FileName) as fd:					#open file
	t = fd.read()							#read the file
	for line in t.splitlines():				#read each line
		subline = line.split()
#		print subline[1][:8]
		if (subline[1][:8] != String1) and (subline[1][:8] != String2):
			rst.append(line)
		else:
			line1 = ""
			line1 = "    " + subline[0] + "\t\t" + subline[1] + "\t\t" + subline[2] + "\t\t" + "0.0"
			rst.append(line1)
with open(NEWDIRECTORY + "/" + FileName, 'w') as fd:
	fd.write('\n'.join(rst))
rst = []

print '.	.	.	.	.	.	SUCCESS!'
