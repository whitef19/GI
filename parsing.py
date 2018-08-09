# -*- coding: utf-8 -*-
# python 3.5
import os
import sys
import csv
import datetime
import argparse
import pandas as pd
import _pickle as cpickle
from bs4 import BeautifulSoup,SoupStrainer

class Ilot:
	def __init__(self,line,col,desired):
		self.ID=line[col.index('ACCESSION')]+'-'+str(line[col.index('START')])+'-'+str(line[col.index('END')])
		self.positif=False
		self.line=line
		self.col=col
		if 'REFERENCE' in col:
			self.positif=True 
		self.info=[]
			
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
	parser.add_argument('-db', action='store', dest='db', help='islandviewer,paidb,iceberg',required=True)
	parser.add_argument('-bn', action='store', dest='pickle', help='output binary file (default=input.pickle)',default='input.pickle')
	parser.add_argument('-acc', action='store', dest='acc', help='list of accession number and organism')
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

def IV4(desired,data,org_dict):
	col=['ACCESSION','ORGANISM','START','END','DETECTION']
	f=open('table.iv4.out','w')
	f.write('#'+'\t'.join(col))
	
	file=open(data,'r')
	islands=[]
	for line in file:
		if not line.startswith('#'):
			line=line.replace('\n','').split('\t')
			organism=org_dict[line[0]] if line[0] in org_dict else None
			line=[line[0],organism]+line[1:]
			print(line)
			f.write('\n'+'\t'.join(line))
			islands.append(Ilot(line,col,desired))
	f.close()
	#return islands

def PAIDB(desired,data,acc_dict):
	import re
	import urllib
	import requests as req

	col=['ACCESSION','ORGANISM','START','END','SEQUENCE','INSERTION','DETECTION']
	f=open('table.PAIDB.out','w')
	f.write('#'+'\t'.join(col))
	file=open(data,'r')
	table=SoupStrainer("table")
	soup=BeautifulSoup(file,'html.parser',parse_only=table).find_all(bordercolordark='white')
	islands=[]
	for species in soup:
		for island in species.findAll(valign='top'):
			cells=island.findAll('td') # list
			if len(cells) >0:
				link=cells[1].find('a').get('href')
				strain=cells[2].text.strip()
				site=cells[4].text.strip()
				locus=cells[5].text.strip().split('(')[0].replace(' ','')
				end=int(float(cells[5].text.strip().split('(')[1].split('kb')[0])*1000)

				organism=' '.join(strain.split(' ')[:2])
				if organism in acc_dict:
					accession=acc_dict[organism]

				resp=req.get('http://www.paidb.re.kr/'+link) # crawling island page
				file=resp.text
				soup=BeautifulSoup(file,'html.parser').select('table') #.findAll('b')
				reference=[]
				for i,row in enumerate(soup):
					if str(i)=='1':
						section=row.findAll('td')
						publication=section[0].text.strip().split('.')
						for pub in publication:
							pub=pub.split(' ')
							if 'PUBMED' in pub :
								reference.append('PMID:'+str(pub[pub.index('PUBMED')+1]))

				line=[accession,strain,'1',str(end),locus,site,'PAIDB-REI',','.join(reference)]

				#islands.append(Ilot(line,col,desired))


				f.write('\n'+'\t'.join(line))
	f.close()
#	return islands		

def ICEberg(desired,data,acc_dict):
	islands=[]
	col=['ACCESSION','ORGANISM','START','END','SEQUENCE','INSERTION','REFERENCE','DETECTION']
	columns=['accession','Organism','Genome coordinates','Nucleotide Sequence','Insertion site','detection']
	f=open('table.ICEberg.out','w')
	f.write('#'+'\t'.join(col))
	for page in range(1,467): # nombre de page 
		info={}
		reference=[]
		file=open(data+'page_'+str(page)+'.html','r')
		soup=BeautifulSoup(file,'html.parser').findAll('tr')
		for row in soup:
			cells=row.findAll('td')
			if len(cells)>0 :
				if cells[0].text.strip() in columns:
					info[cells[0].text.strip()]=cells[1].text.strip()
				elif cells[0].text.strip().startswith('This is a') :
					info['detection']=cells[0].text.strip()
				elif cells[0].text.strip().startswith('(') :
					reference.append(cells[0].text.strip().split('[')[1].replace(']',''))
		info['reference']=','.join(reference)
		organism=' '.join(info['Organism'].split(' ')[:2]) if 'Organism' in info else None
		if organism in acc_dict:
			info['accession']=acc_dict[organism]
		for column in columns :
			if column not in info:
				info[column]='- ..-'
		line=[info['accession'],info['Organism'],info['Genome coordinates'].split('..')[0],info['Genome coordinates'].split('..')[1],info['Nucleotide Sequence'],info['Insertion site'], info['reference'],info['detection']]
		f.write('\n'+'\t'.join(line))
		"""
		for column in columns :
			if column not in info:
				info[column]='- ..-'
		line=[info['Organism'],info['Genome coordinates'].split('..')[0],info['Genome coordinates'].split('..')[1],info['Nucleotide Sequence'],info['Insertion site'], info['reference'],info['detection']]
		islands.append(Ilot(line,col,desired))
		"""
	f.close()
#	return islands

def Islander(desired,data,org_dict):

	islands=[]
	values=[]
	tables=['`accessions`','`smpb`']

	file=open(data,'r')
	for line in file:
		line=line.replace('\n','').split(' ')
		if 'INSERT' in line:
			if line[2] in tables:
				values.append(line[4].replace('),',')|').replace('(','').replace('\'','').replace(')','').split('|'))
	values= [values[0],(values[1]+values[2]+values[3])]

	# section smpb
	col=['ACCESSION','START','END','SEQUENCE']
	for line in values[1]:
		line=line.split(',')
		info=[line[0],str(line[1]),str(line[2]),line[5]]
		ID=info[0]+'-'+info[1]+'-'+info[2]
		if ID not in islands :
			islands.append(Ilot(info,col,desired))

	# section accessions
	col=['ACCESSION','START','END']
	for line in values[0]:
		line=line.split(',')
		info=[line[1],line[2].split('-')[0],line[2].split('-')[1]]
		ID='-'.join(info)
		if ID not in islands :
			islands.append(Ilot(info,col,desired))

	#writing
	f=open('table.islander.out','w')
	f.write('#'+'\t'.join(['ACCESSION','ORGANISM','START','END','SEQUENCE']))
	for island in islands:
		organism=org_dict[island.line[0]] if island.line[0] in org_dict else ''
		f.write('\n'+island.line[0]+'\t'+organism+'\t'+'\t'.join(island.line[1:]))
	f.close()



def main():
	desired=['ACCESSION','ORGANISM','START','END','SEQUENCE','INSERTION','REFERENCE','DETECTION']
	args=argsparse()

	
	if args.db.lower() =='xlsx':
		islands=excel(desired,args.data)
	if args.db.lower() =='islandviewer' :
		organisms={}
		for nb in open(args.acc,'r'):
			organisms[nb.replace('\n','').split('\t')[0]]=nb.replace('\n','').split('\t')[1]
		islands=IV4(desired,args.data,organisms)
	if args.db.lower() =='paidb' :
		accession_nb={}
		for nb in open(args.acc,'r'):
			accession_nb[nb.replace('\n','').split('\t')[1]]=nb.replace('\n','').split('\t')[0]
		islands=PAIDB(desired,args.data,accession_nb)
	if args.db.lower() =='iceberg' :
		accession_nb={}
		for nb in open(args.acc,'r'):
			accession_nb[nb.replace('\n','').split('\t')[1]]=nb.replace('\n','').split('\t')[0]
		islands=ICEberg(desired,args.data,accession_nb)

	if args.db.lower() =='islander' :
		organisms={}
		for nb in open(args.acc,'r'):
			organisms[nb.replace('\n','').split('\t')[0]]=nb.replace('\n','').split('\t')[1]
		islands=Islander(desired,args.data,organisms)
	
#	pickle(args.pickle,islands)


	
if __name__ == "__main__":
	timestamp("STARTING")
	main()
	timestamp("DONE")

"""

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

"""
#	accession={line.replace('\n','').split('\t')[0]:line.replace('\n','').split('\t')[1] for line in open(args.acc,'r')}
#	additional=args.add 
#	islands+=detection(desired,args.add,accession,islands)
#	IDs=set([line.ID for line in islands])
#	writing(desired,islands,IDs,args.output)