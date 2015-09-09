#!/usr/bin/python
import sys
import os
import subprocess

DIRECTORY = sys.path[0]

#---------------------------------------------------------

with open (DIRECTORY + "/Sections/0SETTINGS.txt") as fd:		#open file
	t = fd.read()									#read file
	for line in t.splitlines():						#read each line
		if (line[:16] == "DIRECTORY_OF_OUT"):
			NEWDIRECTORY = line[17:]
			print NEWDIRECTORY

#---------------------------------------------------------
print "\n\n\n ======================================== \n\n\n"
#---------------------------------------------------------

os.chdir(NEWDIRECTORY)
print os.getcwd()
for filename in os.listdir('.'):					#view each file in directory
	if (filename[-4:] == '.gpl'):
		os.system(DIRECTORY + "/deleter_of_strings.py " + filename)
		print filename

#---------------------------------------------------------
print "\n\n\n ======================================== \n\n\n"
#---------------------------------------------------------

os.chdir(NEWDIRECTORY)
print os.getcwd()
for filename in os.listdir('.'):					#view each file in directory
	if (filename[-4:] == '.gpl'):
		os.system(DIRECTORY + "/converter_gpl_to_csv.py " + filename)
		print filename

