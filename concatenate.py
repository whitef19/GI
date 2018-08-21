# -*- coding: utf-8 -*-
# python 3.5
import os
import sys
import csv
import datetime
import argparse
import pandas as pd
from pathlib import Path
import numpy as np

def timestamp(name):
	print("{0} : Timestamp: {1:%Y-%m-%d %H:%M:%S}".format(name, datetime.datetime.now()))

def argsparse():
	parser=argparse.ArgumentParser(description='Concatenate tables from different database and keep information for island in common between source.')
	parser.add_argument('-i', action='store', dest='input', metavar='file', help='file of concatenated tables',required=True)
	parser.add_argument('-o', action='store', dest='output', help='output file database (default=database.tsv)',default='database.tsv')
	args=parser.parse_args()
	return args

def uniq(file):
	Islands=[]
	IDs=[]
	sources=['ICEberg','islander','PAIDB','Dimob','islandviewer','Islander','Sigi','Islandpick']
	islands=[{} for i in sources]
	with open (file,'r') as f:
		for line in f:
			if not line.startswith("#"):
				line=line.replace('\n','').split('\t')
				if ('' not in line[:4]) and (len(line)==8):
					source=line[5]
					ID=line[0]+'-'+line[2]+'-'+line[3]
					islands[sources.index(source)][ID]=(line)
					IDs.append(ID)
	for ID in list(set(IDs)):
		table=[source[ID] for source in islands if ID in source]
		if len(table)>1:
			colums=[(list(set([row[i] for row in table]))) for i in range(8) ]
			line=[','.join(col) if idx!=7 else col[0] for idx,col in enumerate(colums)]
			Islands.append(line)
		else: 
			Islands.append(table[0])
	return Islands

def writing(islands,output):
	out=open(output,'w')
	out.write('#ACCESSION'+'\t'+'ORGANISM'+'\t'+'START'+'\t'+'END'+'\t'+'INSERTION'+'\t'+'DETECTION'+'\t'+'REFERENCE'+'\t'+'SEQUENCE')
	for line in islands:
		out.write('\n'+'\t'.join(line))
	out.close()

def main():
	args=argsparse()	
	islands=uniq(args.input)
	writing(islands,args.output)

if __name__ == "__main__":
	timestamp("STARTING")
	main()
	timestamp("DONE")