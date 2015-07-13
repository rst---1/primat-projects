#!/usr/bin/python
#X
import sys
import numpy

DIRECTORY = sys.path[0]

Slice_y_ = []										#SLICES
with open (DIRECTORY + "/0SETTINGS.txt") as fd:		#open file
	t = fd.read()									#read file
	for line in t.splitlines():						#read each line
		if (line[:7] == "Slice_y"):
			Slice_y_.append( line[8:] )
		if (line[:29] == "4_tau_zx_x.py_ComparingFile_1"):
			ComparingFile_1 = line[30:]
		if (line[:29] == "4_tau_zx_x.py_ComparingFile_2"):
			ComparingFile_2 = line[30:]
		if (line[:31] == "4_tau_zx_x.py_ComparingFile_OUT"):
			ComparingFile_OUT = line[32:]

#---------------------------------------------------------------
#===============================================================
#---------------------------------------------------------------

#Files = [Slice_y_1, Slice_y_2, Slice_y_3, Slice_y_4, Slice_y_5, Slice_y_6, Slice_y_7, Slice_y_8, Slice_y_9]

for i_Files in range(len(Slice_y_)):					#Big Loop
	FileName0 = DIRECTORY + "/" + ComparingFile_1 + " tau_xy when y(" + str(Slice_y_[i_Files]) + ").gpl"
	FileName1 = DIRECTORY + "/" + ComparingFile_2 + " tau_xy when y(" + str(Slice_y_[i_Files]) + ").gpl"
	FileNameOut = DIRECTORY + "/" + ComparingFile_OUT + " tau_xy when y(" + str(Slice_y_[i_Files]) + ").gpl"
	
	fd1 = open (FileName1, 'r')						#open BoxFile
	t1 = fd1.read()									#read BoxFile
	fd0 = open (FileName0, 'r')						#open AnalFile
	t0 = fd0.read()		
	
#-------------------creating arrayes----------------------------
	
	
	rst0 = []
	Len0 = t0.splitlines()
	for j in range(len(Len0)):
		line0 = t0.splitlines()[j]
		str00 = line0.split()[1]
		str01 = line0.split()[2]
		str02 = line0.split()[3]
		rst0.append(str00)
		rst0.append(str01)
		rst0.append(str02)
	rst0 = numpy.array(rst0)
	rst0 = rst0.astype('float')
	
	
	rst1 = []
	Len1 = t1.splitlines()
	for j in range(len(Len1)):
		line1 = t1.splitlines()[j]
		str00 = line1.split()[1]
		str01 = line1.split()[3]
		rst1.append(str00)
		rst1.append(str01)
	rst1 = numpy.array(rst1)
	rst1 = rst1.astype('float')
	
		
	fdOut = open (FileNameOut, 'w')
#---------------------------------------------------------------
	
	Max_Box_abs = -100.0							#varnames of Ave and Max errors
	Max_Box_otn = -100.0
	
	Ave_Box_abs = 0.0
	Ave_Box_otn = 0.0

	Ave_Box = 0.0

	i = 0
	for i in range(len(Len0)):					#number of strings in FilleBox
		x = rst0[i*3+0]							#x
		y = rst0[i*3+1]							#y
		Box = rst0[i*3+2]						#Box
		j = 0
		for j in range (len(Len1)):
			xx = rst1[j*2+0]						#xx
			Curve = rst1[j*2+1]						#Curve
			if ( round(x,2) == round(xx,2)):
				Box_abs = abs(Curve - Box)				#calculating of errors
				if (Box != 0.0):
					Box_otn = abs(Box_abs / Box * 100)
				else:
					Box_otn = 0.0
					
				if(Max_Box_abs < Box_abs ):
					Max_Box_abs = Box_abs
				if(Max_Box_otn < Box_otn ):
					Max_Box_otn = Box_otn
					
				Ave_Box_abs = (i * Ave_Box_abs + Box_abs) / (i+1)
				Ave_Box_otn = (i * Ave_Box_otn + Box_otn) / (i+1)
		
				Ave_Box = (i * Ave_Box + Curve) / (i+1)
			
				#precision FILE: "x","y","Box","Curve","Box_abs","Box_otn"
				fdOut.write( str(x) + "\t" + str(y) + "\t" + str(Box) + "\t" + str(Curve) + "\t" + str(Box_abs) +  "\t" + str(Box_otn) + "\n" )
				break
	
#---------------------------------------------------------------
	
	fdOut.write("\n")
	fdOut.write("Max_Box_abs = " + str(Max_Box_abs) + "\n")
	fdOut.write("Max_Box_otn = " + str(Max_Box_otn) + "\n")
	fdOut.write("\n")
	fdOut.write("Ave_Box_abs = " + str(Ave_Box_abs) + "\n")
	fdOut.write("Ave_Box_otn = " + str(Ave_Box_otn) + "\n")
	fdOut.write("Ave_Box = " + str(Ave_Box) + "\n")
	fdOut.write("\n")
	
	fd1.close()
	fd0.close()
	fdOut.close()
	print FileNameOut
	
#---------------------------------------------------------------
#===============================================================
#---------------------------------------------------------------































