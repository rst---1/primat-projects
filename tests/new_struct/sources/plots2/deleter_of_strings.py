#!/usr/bin/python

import sys
DIRECTORY = sys.path[0]
FileName = sys.argv[1]
String1 = "-0.17187"

#============delete strings of file========================================
print '.	.	.	.	.	.	delete strings of file . .'

rst = []
with open (DIRECTORY + "/" + FileName) as fd:					#open file
	t = fd.read()							#read the file
	for line in t.splitlines():				#read each line
		if (line[:8] != String1):
			rst.append(line)
with open(DIRECTORY + "/" + FileName, 'w') as fd:
	fd.write('\n'.join(rst))
rst = []

print '.	.	.	.	.	.	SUCCESS!'
