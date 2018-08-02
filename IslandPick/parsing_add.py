# -*- coding: utf-8 -*-
# python 3.5
import os
import sys
import csv
import datetime
import pandas as pd
import _pickle as cpickle
import argparse
import optparse

class Ilot:
	def __init__(self,line,col,desired):
		self.ID=line[col.index('ACCESSION')]+'_'+str(line[col.index('START')])+'_'+str(line[col.index('END')])
		self.positif=False
		self.line=line
		self.col=col
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
	parser.add_argument('-i', action='store', nargs='*', dest='data', metavar='file', help='xlsx files list to parse')
	parser.add_argument('-acc', action='store', dest='acc', metavar='file', help='list of accession number with test number')
	parser.add_argument('-add', action='store', dest='add', metavar='file', help='additional table with detection method')
	args=parser.parse_args()
	return args

def reading(table,desired):
	df=pd.read_excel(table)
	columns=[c.upper() for c in df.columns] # upper case columns name
	df.columns=columns	#replace columns name for the upper one
	col=[c for c in desired if c in columns]
	df=df[col]
	df=df.values.tolist()

	element=[]
	for line in df:
		line=Ilot(line,col,desired)
		element.append(line)
	return element

def adding(desired,data):
	df=pd.read_excel(data)
	columns=[c.upper() for c in df.columns] # upper case columns name
	df.columns=columns	#replace columns name for the upper one
	col=[c for c in desired if c in columns]
	df=df[col] # cut desired columns
	data_list=df.values.tolist()
	island=[Ilot(line,col,desired) for line in data_list]
	return island

def detection(desired,data,accession,islands):
	df=pd.read_excel(data,index_col=0,header=2,usecols=[0,1,2,9,11,13,15,17,19])
	df=df.filter(like='landp',axis=0) 
	df=df.reset_index()	
	columns=df.columns[3:]
	detection=df.values.tolist() 
	for line in detection:
		line[0]=accession[str(line[0].replace('islandpick_','').split('_')[0])]
		line=[line[0], int(line[1]), int(line[2]),':'.join((columns[p]+'='+str(pour)) for p,pour in enumerate(line[3:]))]
		#line=[line[0], int(line[1]), int(line[2]),'islanpick']
		ID=line[0]+'_'+str(int(line[1]))+'_'+str(int(line[2]))
		if ID in islands:
			islands[islands.index(ID)].col.append('DETECTION')
			islands[islands.index(ID)].line.append('islanpick')
		else:
			line=Ilot(line,['ACCESSION','START','END','DETECTION'],desired)
			islands.append(line)
	return islands

def writing(desired,islands):
	timestamp("WRITING")
	output=open('output.out','w') # output file
	#output_pos=open('output.pos.out','w') # output file
	output.write('#ID'+'\t'+'\t'.join(desired))
	#output_pos.write('#ID'+'\t'+'\t'.join(desired))
	count=1
	for line in islands:
		output.write('\n'+str(count)+'\t'+'\t'.join(line.get_ajusted()))
		#if line.positif:
		#	output_pos.write('\n'+str(count)+'\t'+'\t'.join(line.ajusted))
		count+=1

	output.close()
	#output_pos.close()

def main():
	desired=['ACCESSION','ORGANISM','START','END','SEQUENCE','INSERTION','REFERENCE','DETECTION']
	args=argsparse()
	islands=[ adding(desired,file) for file in args.data ]
	accession={line.replace('\n','').split('\t')[0]:line.replace('\n','').split('\t')[1] for line in open(args.acc,'r')}
	additional=args.add 
	detection(desired,args.add,accession,islands)

#	pickle('islandpick.pickle',GI)
#	writing(desired,GI)	


if __name__ == "__main__":
	timestamp("STARTING")
	main()
	timestamp("DONE")

"""
args=optsparse()

def optsparse():
	parser=optparse.OptionParser()
	parser.add_option('-i', action='store', nargs='*', dest='data', metavar='file', help='xlsx files list to parse')
	parser.add_option('-acc', action='store', dest='acc', metavar='file', help='list of accession number with test number')
	parser.add_option('-add', action='store', dest='add', metavar='file', help='additional table with detection method')
	args, _ = parser.parse_args()
	return args
"""