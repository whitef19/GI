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
	return df,col

def writing(desired,data,columns):
	timestamp("WRITING")

	exist=False
	file=open('output.out','r')
	exist=True
	for line in file:
		line=line.replace('\n','').split('\t')
		count=line[0]
	if not exist:
		count=1

	positif=False
	if 'REFERENCE' in columns:
		positif=True
	
	output=open('output.out','a') # output file
	
	output_pos=open('output.pos.out','a') # output file
	if not exist:
		output.write('#ID'+'\t'+'\t'.join(desired))
	output_pos.write('#ID'+'\t'+'\t'.join(desired))
	for line in data:
		ajusted=[str(line[columns.index(c)]) if c in columns else '' for c in desired]
		#output.write('\n'+str(count)+'\t'+'\t'.join(ajusted))
		output.write('\n'+'\t'.join(ajusted))
		#count+=1
		if positif:
			output_pos.write('\n'+'\t'.join(line))

	output.close()
	output_pos.close()


def main():

	table=sys.argv[1] # input file
	desired=['ACCESSION','ORGANISM','START','END','SEQUENCE','INSERTION','REFERENCE']
	data,column=reading(table,desired)
	writing(desired,data,column)

if __name__ == "__main__":
	timestamp("STARTING")
	main()
	timestamp("DONE")



def xlsx():
	key = pd.DataFrame({'variants': snp})
	for idx,pop in enumerate(population) :
		colonne = pd.DataFrame({'variants': snp,(str(pop)): frequency[idx]})
		key = pd.merge(key, colonne, on='variants')
	key.to_csv('reports/report.'+str(ID)+'.csv',  sep ="\t")
	labels = ['#SNP','AC_AFR','AC_AMR','AC_ASJ','AC_EAS','AC_FIN','AC_NFE','AC_OTH','AC_SAS','AN_AFR','AN_AMR','AN_ASJ','AN_EAS','AN_FIN','AN_NFE','AN_OTH','AN_SAS']
	df = pd.DataFrame.from_records(Snps, columns=labels)
	cat = df.groupby(['#SNP']).sum()
	cat.to_csv(path_output+'gnomad.coding.cat.csv', sep ="\t")
	cat=cat.reset_index()
	cat=cat.values.tolist()
	return cat