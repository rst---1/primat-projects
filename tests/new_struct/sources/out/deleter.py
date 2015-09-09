#!/usr/bin/python
import sys
import os

DIRECTORY = sys.path[0]


for dirname, dirnames, filenames in os.walk(DIRECTORY):
    for filename in filenames:
		if(  (filename[-4:] == ".gpl")or(filename[-1:] == "~")or(filename[-4:] == ".csv")or(filename[-4:] == ".png")or(filename[-4:] == ".vtk")or(filename[-4:] == ".txt")  ):
			print(os.path.join(dirname, filename))
			os.remove(os.path.join(dirname, filename))



