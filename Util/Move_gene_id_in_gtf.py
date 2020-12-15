#!/usr/bin/python

########################################################
# NBIS 2020 - Sweden                                   #  
# nimarafati@gmail.com	                               #
# Please cite the script by referencing to github      #
# repository 	                           				       #
########################################################
import argparse
import re 
import sys

parser = argparse.ArgumentParser(description = 'This script parse gtf file to move gene_id to begining of the record')
parser.add_argument('-file', help = 'Provide a gtf file', action = 'store')
args = parser.parse_args()

fh = open(args.file)
sample = args.file
for line in fh:
	cntr = 0
	if not line.strip().startswith('#'):
		line_list = line.strip().split('\t')
		attribute = line_list[8].split(';')
		for index, val in enumerate(attribute):
		    if "gene_id" in val:
		        attribute[index], attribute[0] = attribute[0], val
		attribute_txt = '; '.join(map(str,attribute))
		gtf = '\t'.join(map(str,line_list[0:7]))
		print(gtf + '\t' + attribute_txt)
