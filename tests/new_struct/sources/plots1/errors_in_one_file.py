#!/usr/bin/python

import sys
DIRECTORY = sys.path[0]
FileNameOUT = DIRECTORY + "/" + sys.argv[1]
FirstString1 = 'Max_Box_abs'
FirstString1 = 'Max_Box_otn'
FirstString1 = 'Ave_Box_abs'
FirstString1 = 'Ave_Box_otn'
FirstString1 = 'Ave_Box'

#-----------------------------------------------------------------------

import os
print os.getcwd()
os.chdir(DIRECTORY + "/Sections")

rst = []
for filename in os.listdir('.'):					#view each file in directory
	if (filename[:len(sys.argv[1])] == sys.argv[1]) and (filename[-4:] == '.gpl'):
		print filename
		rst.append('\n' + filename + ':\n\n')
		lines = []
		with open (filename) as fd:
			lines = fd.readlines()
			for a in range(7):
				print lines[len(lines) - 7 + a]
				rst.append(lines[len(lines) - 7 + a])
				#print "number of lines = " + str(len(lines))
			rst.append('\n//==============================')
fd.close()

#-----------------------------------------------------------------------

with open(DIRECTORY + "/" + sys.argv[1], 'w') as fd:
	fd.write(''.join(rst))

rst = []


print '.	.	.	.	.	.	SUCCESS!'
