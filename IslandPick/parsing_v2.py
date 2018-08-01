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
		
		self.acc=liste[col.index('ACCESSION')]
		self.start=str(liste[col.index('START')])
		self.stop=str(liste[col.index('END')])
		self.ID=self.acc+'_'+self.start
		self.positif=False
		self.line=liste
		self.ref=''

		if 'REFERENCE' in col:
			self.ref=liste[col.index('REFERENCE')]
			self.positif=True 

		self.ajusted=[str(liste[col.index(c)]) if c in col else '' for c in desired]

	def __eq__ (self, other):
		return self.ID==other.ID
	def __eq__ (self, pos):
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

def writing(desired,element):
	timestamp("WRITING")

	output=open('output.class.out','w') # output file
	output_pos=open('output.class.pos.out','w') # output file
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

	tables=sys.argv[1:] # input file
	desired=['ACCESSION','ORGANISM','START','END','SEQUENCE','INSERTION','REFERENCE']
	element=[]
	for file in tables:
		data=reading(file,desired)
		element=element+data
	writing(desired,element)

if __name__ == "__main__":
	timestamp("STARTING")
	main()
	timestamp("DONE")