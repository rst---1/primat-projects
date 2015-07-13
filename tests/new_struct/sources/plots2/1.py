#!/usr/bin/python
import sys
import os
import subprocess

DIRECTORY = sys.path[0]

os.chdir(DIRECTORY)
print os.getcwd()
for filename in os.listdir('.'):					#view each file in directory
	if (filename[-4:] == '.gpl'):
		os.system(DIRECTORY + "/converter_gpl_to_csv.py " + filename)
		print filename

#a=subprocess.call(["cd "+DIRECTORY+"/Sections/png | ls"], shell=True)

os.chdir(DIRECTORY+"/Sections")
#print os.getcwd()
for filename in os.listdir('.'):					#view each file in directory
	if ((filename[-4:] == '.gpl')or(filename[-4:] == '.csv')or(filename[-1:] == '~')):
		#print filename
		os.remove(filename)
print "*.gpl   *.csv   *.*~ files are deleted in "+os.getcwd()

os.chdir(DIRECTORY+"/Sections/png")
#print os.getcwd()
for filename in os.listdir('.'):					#view each file in directory
	if ((filename[-4:] == '.png')or(filename[-1:] == '~')):
		#print filename
		os.remove(filename)
print "*.png   *.*~ files are deleted in "+os.getcwd()

os.chdir(DIRECTORY+"/Sections/errors")
#print os.getcwd()
for filename in os.listdir('.'):					#view each file in directory
	if ((filename[-4:] == '.png')or(filename[-1:] == '~')):
		#print filename
		os.remove(filename)
print "*.png   *.*~ files are deleted in "+os.getcwd()


