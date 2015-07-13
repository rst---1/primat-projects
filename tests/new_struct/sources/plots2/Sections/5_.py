#!/usr/bin/python

import sys
DIRECTORY = sys.path[0]
FileNameIN = DIRECTORY + "/" + sys.argv[1] + ".gpl"
FileNameOUT = DIRECTORY + "/" + sys.argv[1] + ".gpl__OUT"

rst = []	
fd = open (FileNameIN, 'r')				#open file
t = fd.read()							#read the file

i = 0
#for i in range(255):					#read only 255 lines
A = len(t.splitlines())
for i in range(A-8):
	rst.append(t.splitlines()[i])
fd.close()

with open(FileNameOUT, 'w') as fd:
	fd.write('\n'.join(rst))
rst = []

print FileNameOUT
print FileNameIN
print '.	.	.	.	.	.	SUCCESS!'


