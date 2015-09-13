#!/usr/bin/python

import sys
DIRECTORY = sys.path[0]

#---------------------------------------------------------
with open (DIRECTORY + "/0SETTINGS.txt") as fd:		#open file
	t = fd.read()									#read file
	for line in t.splitlines():						#read each line
		if (line[:16] == "DIRECTORY_OF_OUT"):
			NEWDIRECTORY = line[17:]
			print NEWDIRECTORY
		if (line[:14] == "CSV_FILENAME_1"):
			CSV_FILENAME_1 = line[15:]
			LENGTH_OF_CSV_FILENAME_1 = len(CSV_FILENAME_1) + 5
			print CSV_FILENAME_1
			print "LENGTH_OF_CSV_FILENAME_1 = " + str(LENGTH_OF_CSV_FILENAME_1)

		if (line[:21] == "DIRECTORY_OF_MOUTTTTT"):
			DIRECTORY_OF_MOUTTTTT = line[22:]
			print NEWDIRECTORY
#============delete fist string of file========================================
print '.	.	.	.	.	.	delete fist string of file . .'

rst = []

#------------------------------------------------------------------------------

import os
print '.	.	.	.	.	.	this directory:'
print os.getcwd()
os.chdir(NEWDIRECTORY)

#------------------------------------------------------------------------------

for dirname, dirnames, filenames in os.walk(DIRECTORY_OF_MOUTTTTT):
    for filename in filenames:
		if (filename[-6:] == ".0.csv"):
			PathOfFile = os.path.join(dirname, filename)
			print PathOfFile
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
for dirname, dirnames, filenames in os.walk(DIRECTORY_OF_MOUTTTTT):
    for filename in filenames:
		if (filename[-6:] == ".0.csv"):
			PathOfFile = os.path.join(dirname, filename)
			print PathOfFile
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

for dirname, dirnames, filenames in os.walk(DIRECTORY_OF_MOUTTTTT):
    for filename in filenames:
		if (filename[-6:] == ".0.csv"):
			PathOfFile = os.path.join(dirname, filename)
			print PathOfFile
			os.rename(filename, filename[:-6] + ".gpl")

print '.	.	.	.	.	.	SUCCESS!'

