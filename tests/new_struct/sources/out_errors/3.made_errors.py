#!/usr/bin/python
#for files: U_x_grad.gpl, U_y_grad.gpl, U_z_grad.gpl
#made files with their errros, where:
#	i,	x,	y,	Abs,	Otn
import sys
import numpy
DIRECTORY = sys.path[0]
import os
#print '.	.	.	.	.	.	this directory:'
#print os.getcwd()
os.chdir(DIRECTORY)


#========================read files=============================================
filename1 = sys.argv[1]
filename2 = sys.argv[2]
filenameErrX = sys.argv[3]
filenameErrY = sys.argv[4]


rst1 = []
rst2 = []
filename1EXIST = 0
filename2EXIST = 0
os.chdir(DIRECTORY + "/../out")
for filename in os.listdir('.'):					#view each file in directory
	if (filename == filename1):
#		print filename
		filename1EXIST = 1
		with open (filename) as opfile1:					#open file
			readfile1 = opfile1.read()						#read the file
			Len1 = readfile1.splitlines()					#read each line
			for j in range(len(Len1)):
				line1 = readfile1.splitlines()[j]
				str10 = line1.split()[0]
				str11 = line1.split()[1]
				str12 = line1.split()[2]
				str13 = line1.split()[3]
				str14 = line1.split()[4]
				str15 = line1.split()[5]
#				rst1.append(str10)
				rst1.append(str11)
				rst1.append(str12)
				rst1.append(str13)
				rst1.append(str14)
				rst1.append(str15)
#		rst1.append('\n')
#		print rst1
#		print '\nwwwww\n'
#------------------
os.chdir(DIRECTORY + "/../out_analytic")
for filename in os.listdir('.'):					#view each file in directory
	if (filename == filename2):
#		print filename
		filename2EXIST = 1
		with open (filename) as opfile2:					#open file
			readfile2 = opfile2.read()						#read the file
			Len2 = readfile2.splitlines()					#read each line
			for j in range(len(Len2)):
				line2 = readfile2.splitlines()[j]
				str20 = line2.split()[0]
				str21 = line2.split()[1]
				str22 = line2.split()[2]
				str23 = line2.split()[3]
				str24 = line2.split()[4]
				str25 = line2.split()[5]
#				rst2.append(str20)
				rst2.append(str21)
				rst2.append(str22)
				rst2.append(str23)
				rst2.append(str24)
				rst2.append(str25)
#		rst2.append('\n')
#		print rst2
#		print '\nzzzz\n'

#====================checking files in directory================================
if (filename1EXIST == 0):
#	print "filename1EXIST"
	sys.exit(1)
if (filename2EXIST == 0):
#	print "filename2EXIST"
	sys.exit(1)
#====================create new file============================================

os.chdir(DIRECTORY + "/../out_errors")
rst1 = numpy.array(rst1)
rst1 = rst1.astype('float')
rst2 = numpy.array(rst2)
rst2 = rst2.astype('float')

fdErrX = open (filenameErrX, 'w')
fdErrY = open (filenameErrY, 'w')

absError = -100.0
otnError = -100.0
num = 0
i = 0
for i in range(len(Len1)):					#number of strings in FilleBox
	x = rst1[i*5+0]							#x
	y = rst1[i*5+1]							#y
	Fun = rst1[i*5+2]						#Function
	Fun_X = rst1[i*5+3]						#Function_gradX
	Fun_Y = rst1[i*5+4]						#Function_gradY
	FunA_X = rst2[i*5+3]					#Analytic Function_gradX
	FunA_Y = rst2[i*5+4]					#Analytic Function_gradX
#------------------------
	absError1 = abs(Fun_X-FunA_X)
	if (FunA_X != 0.0):
		otnError1 = abs(absError1 / FunA_X * 100.0)
	else:
		otnError1 = 0.0
	fdErrX.write( '{0:8d}\t\t{1:.35f}\t\t{2:.35f}\t\t{3:.35f}\t\t{4:.35f}\t\t{5:.35f}\t\t{6:.35f}\n'.format(num,x,y,Fun_X,FunA_X,absError1,otnError1) )
#------------------------
	absError2 = abs(Fun_Y-FunA_Y)
	if (FunA_Y != 0.0):
		otnError2 = abs(absError2 / FunA_Y * 100.0)
	else:
		otnError2 = 0.0
	fdErrY.write( '{0:8d}\t\t{1:.35f}\t\t{2:.35f}\t\t{3:.35f}\t\t{4:.35f}\t\t{5:.35f}\t\t{6:.35f}\n'.format(num,x,y,Fun_Y,FunA_Y,absError2,otnError2) )
#------------------------
	num += 1

fdErrX.close()
fdErrY.close()
rst1 = []
rst2 = []

#===============================================================================









