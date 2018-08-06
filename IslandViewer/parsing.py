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
		self.ID=line[col.index('ACCESSION')]+'-'+str(line[col.index('START')])+'-'+str(line[col.index('END')])
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
	parser.add_argument('-i', action='store', dest='data', metavar='file', help='db files list to parse',required=True)
	parser.add_argument('-f', action='store', dest='format', metavar='file', help='format of input file(txt,xlsx)',required=True)
#	parser.add_argument('-acc', action='store', dest='acc', metavar='file', help='list of accession number with test number')
#	parser.add_argument('-add', action='store', dest='add', metavar='file', help='additional table with detection method')
#	parser.add_argument('-o', action='store', dest='output', help='output file database (default=ouput.out',default='output.out')
	parser.add_argument('-bn', action='store', dest='pickle', help='output binary file (default=input.pickle)',default='input.pickle')
	args=parser.parse_args()
	return args

def excel(desired,data):
	df=pd.read_excel(data)
	columns=[c.upper() for c in df.columns] # upper case columns name
	df.columns=columns	#replace columns name for the upper one
	col=[c for c in desired if c in columns]
	df=df[col] # cut desired columns
	data_list=df.values.tolist()
	island=[Ilot(line,col,desired) for line in data_list]
	return island

def txt(desired,data):
	f=open(data,'r')
	islands=[]
	for line in f:
		if line.startswith('#'):
			col=[c.upper() for c in line.replace('\n','').replace('#','').split('\t')]
		else:
			line=line.replace('\n','').split('\t')
			islands.append(Ilot(line,col,desired))
	return islands

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
	add=[]
	df=pd.read_excel(data,index_col=0,header=2,usecols=[0,1,2,5,9,11,13,15,17,19])
	df=df.filter(like='landp',axis=0) 
	df=df.reset_index()	
	columns=[col.replace(' (%)','') for col in df.columns[3:]]
	columns[0]='islandpick'
	detection=df.values.tolist() 
	for line in detection:
		line[0]=accession[str(line[0].replace('islandpick_','').split('_')[0])]
		line=[line[0], int(line[1]), int(line[2]),';'.join((columns[p]+'='+str(pour)) for p,pour in enumerate(line[3:]))]
		ID=line[0]+'_'+str(int(line[1]))+'_'+str(int(line[2]))
		if ID in islands:
			islands[islands.index(ID)].col.append('DETECTION')
			islands[islands.index(ID)].line.append('islanpick')
		else:
			line=Ilot(line,['ACCESSION','START','END','DETECTION'],desired)
			add.append(line)
	return add

def writing(desired,islands,IDs,ouput):
	timestamp("WRITING")
	o=open(ouput,'w') # output file
	#o_pos=open('output.pos.out','w') # output file
	o.write('#ID'+'\t'+'\t'.join(desired))
	#o_pos.write('#ID'+'\t'+'\t'.join(desired))
	count=1
	for ID in IDs:
		o.write('\n'+str(count)+'\t'+'\t'.join(islands[islands.index(ID)].get_ajusted(desired)))
		#if line.positif:
		#	o_pos.write('\n'+str(count)+'\t'+'\t'.join(line.ajusted))
		count+=1

	o.close()
	#o_pos.close()
	#my_file = Path("/path/to/file")
	#if my_file.is_file():
def main():
	desired=['ACCESSION','ORGANISM','START','END','SEQUENCE','INSERTION','REFERENCE','DETECTION']
	args=argsparse()

	if args.format =='xlsx':
		islands=excel(desired,args.data)
	if args.format =='txt' :
		islands=txt(desired,args.data)
	pickle(args.pickle,islands)

#	accession={line.replace('\n','').split('\t')[0]:line.replace('\n','').split('\t')[1] for line in open(args.acc,'r')}
#	additional=args.add 
#	islands+=detection(desired,args.add,accession,islands)
#	IDs=set([line.ID for line in islands])
#	writing(desired,islands,IDs,args.output)
	pickle(args.pickle,islands)
	
if __name__ == "__main__":
	timestamp("STARTING")
	main()
	timestamp("DONE")