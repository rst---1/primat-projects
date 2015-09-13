#!/usr/bin/python
import sys
import os
import subprocess

DIRECTORY = sys.path[0]

#---------------------------------------------------------

with open (DIRECTORY + "/Sections/0SETTINGS.txt") as fd:		#open file
	t = fd.read()									#read file
	for line in t.splitlines():						#read each line
		if (line[:21] == "DIRECTORY_OF_MOUTTTTT"):
			NEWDIRECTORY = line[22:]
			print NEWDIRECTORY
		if (line[:21] == "DIRECTORY_OF_MOUTTTTT"):
			DIRECTORY_OF_MOUTTTTT = line[22:]

#---------------------------------------------------------
print "\n\n\n ======================================== \n\n\n"
#---------------------------------------------------------

for dirname, dirnames, filenames in os.walk(DIRECTORY_OF_MOUTTTTT):
    for filename in filenames:
		if (filename[-4:] == '.gpl'):
			PathOfFile = os.path.join(dirname, filename)
			print PathOfFile
			os.system(DIRECTORY + "/deleter_of_strings.py " + PathOfFile)
#---------------------------------------------------------
print "\n\n\n ======================================== \n\n\n"
#---------------------------------------------------------

for dirname, dirnames, filenames in os.walk(DIRECTORY_OF_MOUTTTTT):
    for filename in filenames:
		if (filename[-4:] == '.gpl'):
			PathOfFile = os.path.join(dirname, filename)
			print PathOfFile
			os.system(DIRECTORY + "/converter_gpl_to_csv.py " + PathOfFile)



