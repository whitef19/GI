#python64/3.5.2
#Usage : python3 parsing table.xlsx 

import os
import sys
import csv
import datetime
import pandas as pd
#import numpy as np
#from math import log,exp
#import _pickle as cpickle

class Ilot:
	def __init__(self,liste,col,desired):
		self.col=col
		self.acc=liste[col.index('ACCESSION')]
		self.start=str(liste[col.index('START')])
		self.stop=str(liste[col.index('END')])
		self.ID=self.acc+'_'+self.start
		self.positif=False
		self.line=liste
		self.ref=''
		self.detection=''

		if 'REFERENCE' in col:
			self.ref=liste[col.index('REFERENCE')]
			self.positif=True 
		self.ajusted=[str(self.line[self.col.index(c)]) if c in self.col else '' for c in desired]

		
	def get_ajusted(self,desired):
		self.ajusted=[str(self.line[self.col.index(c)]) if c in self.col else '' for c in desired]
		return self.ajusted

	def __eq__ (self, other):
		return self.ID==other.ID
	def __eq__ (self, ID):
		return self.ID==ID

def timestamp(name):
	print("{0} : Timestamp: {1:%Y-%m-%d %H:%M:%S}".format(name, datetime.datetime.now()))

def reading(table,desired):
	timestamp("READING")
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

def detection(additional,accession_nb,GI,desired):
	IDs={}
	for line in accession_nb:
		line=line.replace('\n','').split('\t')
		IDs[line[0]]=line[1]

	df=pd.read_excel(additional,index_col=0,header=2,usecols=[0,1,2,9,11,13,15,17,19])
	df=df.filter(like='landp',axis=0)
	df=df.reset_index()	
	detection=df.values.tolist()
	for line in detection:
		line[0]=IDs[str(line[0].replace('islandpick_','').split('_')[0])]
		#line=[line[0], int(line[1]), int(line[2]),':'.join(str(p) for p in line[3:])]
		line=[line[0], int(line[1]), int(line[2]),'islanpick']

		ID=line[0]+'_'+str(int(line[1]))
		if ID in GI:
			GI[GI.index(ID)].col.append('DETECTION')
			GI[GI.index(ID)].line.append('islanpick')
			ajusted=GI[GI.index(ID)].get_ajusted(desired)
		else:
			line=Ilot(line,['ACCESSION','START','END','DETECTION'],desired)
			GI.append(line)
	return GI

def writing(desired,element):
	timestamp("WRITING")

	output=open('output.out','w') # output file
	output_pos=open('output.pos.out','w') # output file
	output.write('#ID'+'\t'+'\t'.join(desired))
	output_pos.write('#ID'+'\t'+'\t'.join(desired))
	count=1
	for line in element:
		output.write('\n'+str(count)+'\t'+'\t'.join(line.ajusted))
		if line.positif:
			output_pos.write('\n'+str(count)+'\t'+'\t'.join(line.ajusted))
		count+=1

	output.close()
	output_pos.close()

def main():
	accession_nb=[line.replace('\n','') for line in open(sys.argv[1],'r')]
	additional=sys.argv[2] # input file
	tables=sys.argv[3:] # input file
	desired=['ACCESSION','ORGANISM','START','END','SEQUENCE','INSERTION','REFERENCE','DETECTION']
	

	GI=[]
	for file in tables:
		data=reading(file,desired)
		GI=GI+data
	
	GI=detection(additional,accession_nb,GI,desired)
	writing(desired,GI)


if __name__ == "__main__":
	timestamp("STARTING")
	main()
	timestamp("DONE")

	#GI=df=pd.read_excel(additional,index_col=0,usecols=[0,1,5,7,9,11,13,15,17,19],header=2,nrows=10)

	#print(IDs)

	#print(df.head())
	#acces=(df.at[4,0]).split(',')[0]
	#df[0]=(acces+'_'+str(df[1]))
	#print(str(df[1]))
	#print(df.head(20))
