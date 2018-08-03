# -*- coding: utf-8 -*-
# python 3.5
import os
import sys
import csv
import argparse
import datetime
import pandas as pd
import _pickle as cpickle

class Ilot:
	def __init__(self,line,col,desired):
		self.col=col
		self.ID=line[col.index('ACCESSION')]+'_'+str(line[col.index('START')])+'_'+str(line[col.index('END')])
		self.line=line
		self.positif=False
		if 'REFERENCE' in col:
			self.positif=True 
		
	def get_ajusted(self,desired):
		self.ajusted=[str(self.line[self.col.index(c)]) if c in self.col else '' for c in desired]
		return self.ajusted
	def __eq__ (self, other):
		return self.ID==other.ID
	def __eq__ (self, ID):
		return self.ID==ID

def timestamp(name):
	print("{0} : Timestamp: {1:%Y-%m-%d %H:%M:%S}".format(name, datetime.datetime.now()))

def argsparse():
	parser=argparse.ArgumentParser(description='Parsing of genomic island database')
	parser.add_argument('-bp','--pad', action='store', dest='pad',type=int, help='number of base pair to add', default=50)
	parser.add_argument('-p','--pickle', action='store', dest='pickle', help='pickle file of genomic island', required=True)
	parser.add_argument('-sh', action='store', dest='sh', help='output sh script (default=sequences.sh)', default='sequences.sh')
#	parser.add_argument('-o','--output', action='store', dest='output', help='output prefix (default=output.out)', default='output.out')
#	parser.add_argument('-key', action='store', nargs='*', dest='key', help='keyword for research')
	args=parser.parse_args()
	return args

def pickle(name,objects) :
	with open(name,'wb') as file:
		pickler = cpickle.Pickler(file)
		pickler.dump(objects)

def unpickle(name):
	with open(name,'rb') as f:
		pickler=cpickle.Unpickler(f)
		objects=pickler.load()
	return objects

def writing(desired,strings,output):
	o=open('EDirect.sh','w') 
	o.write('#!/bin/bash'+'\n')
	for string in strings:
		o.write('\n'+'echo "#'+string+'" >> '+output)
		o.write('\n'+'esearch -db pubmed -query "'+string+'" -sort Relevance |efetch -format medline | grep -wf pattern >> '+output)
	o.write('\n'+'echo "DONE"')
	o.close()

def main():
	args=argsparse()
	islands=unpickle(args.pickle)
	desired=['ACCESSION','ORGANISM','START','END','SEQUENCE','INSERTION','REFERENCE','DETECTION']
	Query=[]
	for island in islands:
		organism=island.get_ajusted(desired)[1]
		print(organism)
		ID=island.ID.split('_')
		print(ID)
		key=[ID[0],int(float(ID[1])-args.pad),int(float(ID[2])+args.pad)]
		print(key)
		Query.append(key)

	#writing(desired,strings,args.output)

if __name__ == "__main__":
	timestamp("STARTING")
	main()
	timestamp("DONE")