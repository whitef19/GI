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
	parser.add_argument('-db',          action='store',     dest='pickle',help='pickle file of genomic islands', required=True)
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

def writing(Query,IDs,sh,output,crop):
	o=open(sh,'w') 
	o.write( '#!/bin/bash'+'\n' )
	for organism in set(IDs):
		acc=organism.split('-')[0]
		file=output+acc+'.fasta'
		o.write( '\n'+'\n'+'#'+acc ) # verify if genome file already exist
		o.write('\n'+ 'if [ ! -f '+file+' ] ; then esearch -db Nucleotide -query "('+acc+')"|efetch -format fasta > '+file+' ;fi') 
	o.write('\n'+'echo "DONE"')
	o.close()


	if crop: # creation of cropped fasta 
		c=open('cropped.'+sh,'w')
		c.write( '#!/bin/bash'+'\n' )
		for query in Query:
			file=output+query[0]+'.fasta'
			o.write('\n'+ 'if [ ! -f '+output+'crop.'+query[3]+'.fasta ] ; then ') 
			c.write('\n'+ '\t' +'header=` grep ">" '+file+'`' )
			c.write('\n'+ '\t' +'echo "$header '+query[0]+' '+query[1]+'-'+query[2]+'" >'+output+'crop.'+query[3]+'.fasta' )
			c.write('\n'+ '\t' +'grep -v ">" '+file+" | sed ':a;N;$!ba;s/\\n//g'| cut -c "+query[1]+'-'+query[2]+' >>'+output+'crop.'+query[3]+'.fasta ; fi' )
		c.write('\n'+'echo "DONE"')
		c.close()
			
def main():
	args=argsparse()
	islands=unpickle(args.pickle)
	desired=['ACCESSION','ORGANISM','START','END','SEQUENCE','INSERTION','REFERENCE','DETECTION']
	Query=[]
	IDs=[]
	for island in islands:
		ID=island.ID.split('-')
		key=[ID[0],str(int(ID[1])-args.pad),str(int(ID[2])+args.pad),island.ID]
		IDs.append(island.ID)
		Query.append(key)
	writing(Query,IDs,args.sh,args.output,args.crop)

if __name__ == "__main__":
	timestamp("STARTING")
	main()
	timestamp("DONE")