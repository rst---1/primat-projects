#!/usr/bin/python

import sys
DIRECTORY = sys.path[0]

#---------------------------------------------------------

with open (DIRECTORY + "/0SETTINGS.txt") as fd:		#open file
	t = fd.read()									#read file
	for line in t.splitlines():						#read each line
		if (line[:30] == "1_pvpython.py_DIRECTORY_OF_OUT"):
			NEWDIRECTORY = line[31:]
			print NEWDIRECTORY
		if (line[:14] == "CSV_FILENAME_1"):
			CSV_FILENAME_1 = line[15:]
			LENGTH_OF_CSV_FILENAME_1 = len(CSV_FILENAME_1) + 5
			print CSV_FILENAME_1
			print "LENGTH_OF_CSV_FILENAME_1 = " + str(LENGTH_OF_CSV_FILENAME_1)
		if (line[:14] == "CSV_FILENAME_2"):
			CSV_FILENAME_2 = line[15:]
			LENGTH_OF_CSV_FILENAME_2 = len(CSV_FILENAME_2) + 5
			print CSV_FILENAME_2
			print "LENGTH_OF_CSV_FILENAME_2 = " + str(LENGTH_OF_CSV_FILENAME_2)

#============delete fist string of file========================================
print '.	.	.	.	.	.	delete fist string of file . .'

rst = []
#FirstString1 = '"Points:3","Points:4","Points:5","Points:6","Points:7","Points:8","Points:9","Points:0","Points:1","Points:2"'
#FirstString2 = '"Points:2","Points:4","Points:5","Points:6","Points:7","Points:8","Points:9","Points:0","Points:1","Points:2"'
#FirstString3 = '"Points:2","Points:3","Points:4","Points:6","Points:7","Points:8","Points:9","Points:0","Points:1","Points:2"'
#FirstString4 = '"Points:2","Points:3","Points:5","Points:6","Points:7","Points:8","Points:9","Points:0","Points:1","Points:2"'

#------------------------------------------------------------------------------

import os
print '.	.	.	.	.	.	this directory:'
print os.getcwd()
os.chdir(NEWDIRECTORY)

#------------------------------------------------------------------------------

for filename in os.listdir("."):					#view each file in directory
	if (((filename[-4:] == '.csv') and (filename[:LENGTH_OF_CSV_FILENAME_1] == CSV_FILENAME_1 + " tau_")) or ((filename[-4:] == '.csv') and (filename[:LENGTH_OF_CSV_FILENAME_2] == CSV_FILENAME_2 + " tau_"))):
		print filename
		with open (filename) as fd:					#open file
			t = fd.read()							#read the file
			IS_FirstString = 1						#Flag for first string of a file
			for line in t.splitlines():				#read each line
				if (IS_FirstString == 0):
					rst.append(line)
				else:
					IS_FirstString = 0
		with open(filename, 'w') as fd:
			fd.write('\n'.join(rst))
	rst = []

print '.	.	.	.	.	.	SUCCESS!'

#============replacement commas by tabs========================================
print '.	.	.	.	.	.	replacement commas by tabs . .'

rst = []
for filename in os.listdir("."):					#view each file in directory
	if (((filename[-4:] == '.csv') and (filename[:LENGTH_OF_CSV_FILENAME_1] == CSV_FILENAME_1 + " tau_")) or ((filename[-4:] == '.csv') and (filename[:LENGTH_OF_CSV_FILENAME_2] == CSV_FILENAME_2 + " tau_"))):
		print filename
		with open (filename) as fd:					#open file
			t = fd.read()							#read the file
			for line in t.splitlines():				#read each line
				for letter in line:
					if letter == ',':
						letter ='\t'
					rst.append(letter)
				rst.append('\n')
		with open(filename, 'w') as fd:
			fd.write(''.join(rst))
	rst = []

print '.	.	.	.	.	.	SUCCESS!'

#============rename a file extension===========================================
print '.	.	.	.	.	.	rename a file extension . .'

for filename in os.listdir("."):
	if (((filename[-4:] == '.csv') and (filename[:LENGTH_OF_CSV_FILENAME_1] == CSV_FILENAME_1 + " tau_")) or ((filename[-4:] == '.csv') and (filename[:LENGTH_OF_CSV_FILENAME_2] == CSV_FILENAME_2 + " tau_"))):
		os.rename(filename, filename[:-6] + ".gpl")

print '.	.	.	.	.	.	SUCCESS!'

