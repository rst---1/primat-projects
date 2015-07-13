#!/usr/bin/python
#Delete files in directories: out, out_analytic, out_errors
import sys
import os
import subprocess

DIRECTORY = sys.path[0]

Files = ["FU_x.gpl", "FU_y.gpl",
		"FU_z.gpl", "FU_x_grad.gpl", 
		"FU_y_grad.gpl", "FU_z_grad.gpl", 
		"Ftau_xx.gpl", "Ftau_yy.gpl",
		"Ftau_xy.gpl", "Ftau_zz.gpl",
		"Ftau_zx.gpl", "Ftau_zy.gpl"]

FilesA = ["FAU_x.gpl", "FAU_y.gpl",
		"FAU_z.gpl", "FAU_x_grad.gpl", 
		"FAU_y_grad.gpl", "FAU_z_grad.gpl", 
		"FAtau_xx.gpl", "FAtau_yy.gpl",
		"FAtau_xy.gpl", "FAtau_zz.gpl",
		"FAtau_zx.gpl", "FAtau_zy.gpl"]

FilesErr = ["ErrFU_x.gpl", "ErrFU_y.gpl",
		"ErrFU_z.gpl", "ErrFU_x_grad.gpl", 
		"ErrFU_y_grad.gpl", "ErrFU_z_grad.gpl", 
		"ErrFtau_xx.gpl", "ErrFtau_yy.gpl",
		"ErrFtau_xy.gpl", "ErrFtau_zz.gpl",
		"ErrFtau_zx.gpl", "ErrFtau_zy.gpl",
		"ErrFU_x_gradX.gpl", "ErrFU_x_gradY.gpl",
		"ErrFU_y_gradX.gpl", "ErrFU_y_gradY.gpl",
		"ErrFU_z_gradX.gpl", "ErrFU_z_gradY.gpl"]

#print Files

os.chdir(DIRECTORY+"/out_analytic")
#print os.getcwd()
for ii in range(len(FilesA)):					#view each file in directory
	for filename in os.listdir('.'):
		if (filename == FilesA[ii]):
			print filename
			os.remove(filename)

print "\n\n"

os.chdir(DIRECTORY+"/out_errors")
#print os.getcwd()
for ii in range(len(FilesErr)):
	for filename in os.listdir('.'):					#view each file in directory
		if (filename == FilesErr[ii]):
			print filename
			os.remove(filename)

print "\n\n"

for dirname, dirnames, filenames in os.walk(DIRECTORY+"/out"):
    for filename in filenames:
		if(  (filename[-4:] == ".gpl")or(filename[-1:] == "~")or(filename[-4:] == ".csv")or(filename[-4:] == ".png")or(filename[-4:] == ".vtk")or(filename[-4:] == ".txt")  ):
			print(os.path.join(dirname, filename))
			os.remove(os.path.join(dirname, filename))
