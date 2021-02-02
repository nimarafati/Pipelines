#!/usr/bin/python

########################################################
# NBIS 2020 - Sweden                                   #  
# nimarafati@gmail.com	                               #
# Please cite the script by referencing to github      #
# repository 	    				       #
########################################################
import argparse
import re 
import sys

parser = argparse.ArgumentParser(description = 'Parse bbduk results')
parser.add_argument('-file', help = 'bbduk summary results', action = 'store')
args = parser.parse_args()

fh = open(args.file)
sample = args.file
sample = sample.replace('_bbduk.txt', '')
print('Sample\tReads_percentage\tBases_percentage')
for line in fh:
	cntr = 0
	if not line.strip().startswith('#'):
		if line.strip().startswith('Total Remove'):
#			sys.stdin.readline()
#			res = re.findall(r'\(\d+\)', line)
			res = re.findall(r'\((.*?)%\)', line)
			print(sample + '\t' + str(res[0]) + '\t' + str(res[1]))
#		Total Removed:          	1250616 reads (3.15%) 	159889945 bases (3.12%)
	
