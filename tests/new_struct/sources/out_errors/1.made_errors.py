#!/usr/bin/python

import sys
import numpy
DIRECTORY = sys.path[0]
import os
os.chdir(DIRECTORY)

result=0

result = os.system(DIRECTORY+"/2.made_errors.py FU_x.gpl FAU_x.gpl ErrFU_x.gpl")
if(result==0):
	print DIRECTORY
else:
	print "\nERROR==============ERROR==============ERROR==============ERROR"
	print "AAAAAAAAAAAAAAAAAAAAAAAAAAAAWTFAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
	print "ERROR==============ERROR==============ERROR==============ERROR\n"


result = os.system(DIRECTORY+"/2.made_errors.py FU_y.gpl FAU_y.gpl ErrFU_y.gpl")
if(result==0):
	print DIRECTORY
else:
	print "\nERROR==============ERROR==============ERROR==============ERROR"
	print "AAAAAAAAAAAAAAAAAAAAAAAAAAAAWTFAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
	print "ERROR==============ERROR==============ERROR==============ERROR\n"


result = os.system(DIRECTORY+"/2.made_errors.py FU_z.gpl FAU_z.gpl ErrFU_z.gpl")
if(result==0):
	print DIRECTORY
else:
	print "\nERROR==============ERROR==============ERROR==============ERROR"
	print "AAAAAAAAAAAAAAAAAAAAAAAAAAAAWTFAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
	print "ERROR==============ERROR==============ERROR==============ERROR\n"





result = os.system(DIRECTORY+"/3.made_errors.py FU_x_grad.gpl FAU_x_grad.gpl ErrFU_x_gradX.gpl ErrFU_x_gradY.gpl")
if(result==0):
	print DIRECTORY
else:
	print "\nERROR==============ERROR==============ERROR==============ERROR"
	print "AAAAAAAAAAAAAAAAAAAAAAAAAAAAWTFAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
	print "ERROR==============ERROR==============ERROR==============ERROR\n"


result = os.system(DIRECTORY+"/3.made_errors.py FU_y_grad.gpl FAU_y_grad.gpl ErrFU_y_gradX.gpl ErrFU_y_gradY.gpl")
if(result==0):
	print DIRECTORY
else:
	print "\nERROR==============ERROR==============ERROR==============ERROR"
	print "AAAAAAAAAAAAAAAAAAAAAAAAAAAAWTFAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
	print "ERROR==============ERROR==============ERROR==============ERROR\n"


result = os.system(DIRECTORY+"/3.made_errors.py FU_z_grad.gpl FAU_z_grad.gpl ErrFU_z_gradX.gpl ErrFU_z_gradY.gpl")
if(result==0):
	print DIRECTORY
else:
	print "\nERROR==============ERROR==============ERROR==============ERROR"
	print "AAAAAAAAAAAAAAAAAAAAAAAAAAAAWTFAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
	print "ERROR==============ERROR==============ERROR==============ERROR\n"





result = os.system(DIRECTORY+"/2.made_errors.py Ftau_xx.gpl FAtau_xx.gpl ErrFtau_xx.gpl")
if(result==0):
	print DIRECTORY
else:
	print "\nERROR==============ERROR==============ERROR==============ERROR"
	print "AAAAAAAAAAAAAAAAAAAAAAAAAAAAWTFAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
	print "ERROR==============ERROR==============ERROR==============ERROR\n"


result = os.system(DIRECTORY+"/2.made_errors.py Ftau_yy.gpl FAtau_yy.gpl ErrFtau_yy.gpl")
if(result==0):
	print DIRECTORY
else:
	print "\nERROR==============ERROR==============ERROR==============ERROR"
	print "AAAAAAAAAAAAAAAAAAAAAAAAAAAAWTFAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
	print "ERROR==============ERROR==============ERROR==============ERROR\n"


result = os.system(DIRECTORY+"/2.made_errors.py Ftau_xy.gpl FAtau_xy.gpl ErrFtau_xy.gpl")
if(result==0):
	print DIRECTORY
else:
	print "\nERROR==============ERROR==============ERROR==============ERROR"
	print "AAAAAAAAAAAAAAAAAAAAAAAAAAAAWTFAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
	print "ERROR==============ERROR==============ERROR==============ERROR\n"


result = os.system(DIRECTORY+"/2.made_errors.py Ftau_zz.gpl FAtau_zz.gpl ErrFtau_zz.gpl")
if(result==0):
	print DIRECTORY
else:
	print "\nERROR==============ERROR==============ERROR==============ERROR"
	print "AAAAAAAAAAAAAAAAAAAAAAAAAAAAWTFAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
	print "ERROR==============ERROR==============ERROR==============ERROR\n"





result = os.system(DIRECTORY+"/2.made_errors.py Ftau_zx.gpl FAtau_zx.gpl ErrFtau_zx.gpl")
if(result==0):
	print DIRECTORY
else:
	print "\nERROR==============ERROR==============ERROR==============ERROR"
	print "AAAAAAAAAAAAAAAAAAAAAAAAAAAAWTFAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
	print "ERROR==============ERROR==============ERROR==============ERROR\n"


result = os.system(DIRECTORY+"/2.made_errors.py Ftau_zy.gpl FAtau_zy.gpl ErrFtau_zy.gpl")
if(result==0):
	print DIRECTORY
else:
	print "\nERROR==============ERROR==============ERROR==============ERROR"
	print "AAAAAAAAAAAAAAAAAAAAAAAAAAAAWTFAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
	print "ERROR==============ERROR==============ERROR==============ERROR\n"
