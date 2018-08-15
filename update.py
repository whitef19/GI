# -*- coding: utf-8 -*-
# python 3.5
import os
import sys
import csv
import datetime
import argparse
import pandas as pd
import _pickle as cpickle
from pathlib import Path
import numpy as np

def timestamp(name):
	print("{0} : Timestamp: {1:%Y-%m-%d %H:%M:%S}".format(name, datetime.datetime.now()))

def pickle(name,objects) :
	with open(name,'wb') as file:
		pickler = cpickle.Pickler(file)
		pickler.dump(objects)

def unpickle(name):
	with open(name,'rb') as f:
		pickler=cpickle.Unpickler(f)
		objects=pickler.load()
	return objects

def argsparse():
	parser=argparse.ArgumentParser(description='Parse genomic islands databases')
	parser.add_argument('-i', action='store', dest='input', metavar='file', help='file of concatenated tables',required=True)
	parser.add_argument('-o', action='store', dest='output', help='output file database (default=database.txt)',default='database.txt')
	parser.add_argument('-seq', action='store_true', dest='seq', help='add sequence from downloaded fasta')
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
				if ('' not in line[:4]) and len(line)==8:
					source=line[5]
					ID=line[0]+'-'+line[2]+'-'+line[3]
					islands[sources.index(source)][ID]=(line)
					IDs.append(ID)
	for ID in list(set(IDs)):
		table=[source[ID] for source in islands if ID in source]
		if len(table)>1:
			colums=[(list(set([row[i] for row in table]))) for i in range(8) ]
			line=[col[0] if idx!=5 else ','.join(col) for idx,col in enumerate(colums)]
			Islands.append(line)
		else: 
			Islands.append(table[0])
	return Islands

def sequences(file):
	Islands=[]
	errors=[]
	with open (file,'r') as f:
		for island in f:
			if not island.startswith("#"):
				island=island.replace('\n','').split('\t')
				ID=island[0]+'-'+island[2]+'-'+island[3]

				header=[]
				error_fasta=False
				if os.path.isfile('sequences/island.'+ID+'.fa'):
					with open('sequences/island.'+ID+'.fa','r') as f:
						for line in f:
							if line.startswith('>') :
								header.append(line.replace('\n',''))
								if len(header)!=1:
									errors.append(ID)
									error_fasta=True

							elif (not error_fasta) and (island[2]!='1'): 
								if island[7]=='':
									island[7]=line.replace('\n','')
									Islands.append(island)
							else: 
								Islands.append(island)
				else:
					Islands.append(island)
					if island[7]!='':
						o=open('sequences/island.'+ID+'.fa','w')
						o.write('>'+island[0]+' '+island[1]+' db:'+island[5]+' '+ID)
						o.write('\n'+island[7])
						o.close()

	return Islands,errors

def writing(islands,output):
	out=open(output,'w')
	out.write('#ACCESSION'+'\t'+'ORGANISM'+'\t'+'START'+'\t'+'END'+'\t'+'INSERTION'+'\t'+'DETECTION'+'\t'+'REFERENCE'+'\t'+'SEQUENCE')
	for line in islands:
		out.write('\n'+'\t'.join(line))
	out.close()

def main():
	args=argsparse()
	file=args.input
	output=args.output
	
	if args.seq :
		islands,errors=sequences(file)
		writing(islands,output)
		o=open('errors.out','w')
		for e in errors:
			o.write('\n'+e)
		o.close()
	else:
		islands=uniq(file)
		writing(islands,output)

	
if __name__ == "__main__":
	timestamp("STARTING")
	main()
	timestamp("DONE")