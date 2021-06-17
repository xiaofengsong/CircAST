#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# This program is designed for extracting circular transctipt sequence outputted in CircAST_result.txt
import sys
if sys.version_info[0] != 2 or sys.version_info[1] != 7:
	print >>sys.stderr, "\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " CPAT needs python2.7!\n"
	sys.exit()
import os
import optparse
parse = optparse.OptionParser()
parse.add_option('-F','--fafolder',dest = 'fafolder',action = 'store',metavar = 'fa file folder',help = 'Please enter your fa file folder, which include all the decompressed file of chromFa.tar.gz')
parse.add_option('-I','--input',dest = 'inputfile',action = 'store',metavar = 'input file',help = 'Please enter CircAST_result.txt or some lines of this file')
parse.add_option('-O','--output',dest = 'outputfile',action = 'store',metavar = 'output file',help = 'Please enter output file name,such as *.txt or *.fa')


(options,args) = parse.parse_args()
#check input files
for file in ([options.fafolder,options.inputfile]):
	if not (file):
		print >>sys.stderr,"\nError: Lack of input file!\n"
		parse.print_help()
		sys.exit(0)

inputCircASTFile = options.inputfile
inputFaFolder = options.fafolder
outputTFile = options.outputfile
if inputFaFolder[-1] == "/":
	inputFaFolder = inputFaFolder[:-1]


chr_set = set()
in_list = open(inputCircASTFile)
for line in in_list:
	line = line.strip('\n')
	array = line.split("\t")
	chr_set.add(array[0])
in_list.close()

for i in chr_set:
	if os.path.exists("temp_CircAST_" + i + ".txt"):
		os.remove("temp_CircAST_" + i + ".txt")


if os.path.exists(outputTFile):
	os.remove(outputTFile)

in_list = open(inputCircASTFile)
for line in in_list:
	line = line.strip('\n')
	array = line.split("\t")
	out_list = open("temp_CircAST_" + array[0] + ".txt","a")
	out_list.write(line + "\n")
	out_list.close()
in_list.close()

line = ""
out_list = open(outputTFile,"a")
for i in chr_set:
	fa_name = inputFaFolder + "/" + i + ".fa"
	if not os.path.exists(fa_name):
		out_list.write("Error: " + i + ".fa" + " is not in the folder: " + inputFaFolder + "\n")
out_list.close()

out_list = open(outputTFile,"a")
for i in chr_set:
	fa_name = inputFaFolder + "/" + i + ".fa"
	if os.path.exists(fa_name):
		in_fa = open(fa_name,"r")
		in_fa.readline()
		chr_fa_string = in_fa.read()
		chr_fa_string = chr_fa_string.replace("\n","")
		in_fa.close()
		in_list = open("temp_CircAST_" + i + ".txt","r")
		for line in in_list:
			line = line.strip('\n')
			array = line.split("\t")
			start = array[1]
			end = array[2]
			exon_start = array[8].split(",")
			exon_len = array[7].split(",")
			exon_end = [0] * len(exon_len)
			tr_seq = ""
			strain = array[5]
			for t in range(0,len(exon_len)):
				temp_start = int(exon_start[t]) + int(start)
				temp_end = temp_start + int(exon_len[t]) - 1
				tr_seq = tr_seq + chr_fa_string[temp_start - 1:temp_end]
			tr_seq = tr_seq.upper()
			if strain == "-":
				tr_seq = tr_seq[::-1]
				tr_seq = tr_seq.replace("A","X")
				tr_seq = tr_seq.replace("T","A")
				tr_seq = tr_seq.replace("X","T")
				tr_seq = tr_seq.replace("C","Y")
				tr_seq = tr_seq.replace("G","C")
				tr_seq = tr_seq.replace("Y","G")
			out_list.write(">" + array[3] + "\n")
			out_list.write(tr_seq + "\n")
		in_list.close()
out_list.close()
		
for i in chr_set:
	os.remove("temp_CircAST_" + i + ".txt")

	
