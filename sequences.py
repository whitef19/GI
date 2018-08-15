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
		self.ID=line[col.index('ACCESSION')]+'-'+str(line[col.index('START')])+'-'+str(line[col.index('END')])
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
	parser.add_argument('-i',           action='store',     dest='input', help='input file', required=True)
	parser.add_argument('-sh',          action='store',     dest='sh',    help='output sh script (default=sequences.sh)', default='sequences.sh')
	parser.add_argument('-o','--output',action='store',     dest='output',help='output prefix (default=sequence.)', default='sequence.')
	parser.add_argument('-c','--crop',  action='store_true',dest='crop',  help='create a cropped fasta file with interval of the island')
	parser.add_argument('-bp','--pad',  action='store',     dest='pad',   help='number of base pair to add',type=int,default=50)

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

def reading(file):
	query=[]
	with open(file,'r') as f:
		for line in f:
			if not line.startswith('#'):
				line=line.replace('\n','').split('\t')
				query.append([line[0],line[2],line[3]])
	return query
			
def writing(query):
	o=open('get_sequences.sh','w') 
	o.write( '#!/bin/bash'+'\n' )
	for island in query:
		if island[1]!='1':
			ID='-'.join(island)
			genome='genomes/sequence.'+island[0]+'.fasta'
			output='sequences/island.'+ID+'.fa'
			interval=str(int(island[1])-50)+'-'+str(int(island[2])+50)
			o.write('\n'+ 'if [ ! -f '+genome+' ] ; then') 
			o.write('\n'+'\t'+'echo "'+ID+'" >> no_fasta.out')
			o.write('\n'+'else')
			o.write('\n'+ '\t' +'header=` grep ">" '+genome+'`' )
			o.write('\n'+ '\t' +'echo "$header '+island[0]+' '+interval+'" > '+output)
			o.write('\n'+ '\t' +'grep -v ">" '+genome+" | sed ':a;N;$!ba;s/\\n//g'| cut -c "+interval+' >> '+output)
			o.write('\n'+'fi'+'\n')
	o.write('\n'+'echo "DONE"')
	o.close()

def main():
	args=argsparse()
	database=args.input
	query=reading(database)
	writing(query)

if __name__ == "__main__":
	timestamp("STARTING")
	main()
	timestamp("DONE")