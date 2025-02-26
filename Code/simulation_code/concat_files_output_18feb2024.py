#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f','--list_files',dest="list_files",required=True) #input file with the list of files to combine
parser.add_argument('-o','--output',dest="output",required=True) #output file
args = parser.parse_args()

input_file=open(args.list_files, 'r')

output_file=args.output

list_files=[]

for line in input_file:
	list_files.append(line.rstrip('\n'))
input_file.close()


first_file=True

output_string=''

for temp_file in list_files:
	file=open(temp_file, 'r')

	first_line=True

	for line in file:
		if(first_file):
			header=line.rstrip('\n')
			output_string+=line
			first_file=False
			first_line=False
		else:
			if(first_line):
				if(line.rstrip('\n')!=header):
					print('Error headers do not match!\nProgram execution stopped.')
					file.close()
					exit()
				first_line=False
			else:
				output_string+=line


output=open(output_file, 'w')
output.write(output_string)
output.close()

