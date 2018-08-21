# -*- coding: utf-8 -*-
# python 3.5
import os
import sys
import csv
import datetime
import argparse
import pandas as pd
from bs4 import BeautifulSoup,SoupStrainer

def timestamp(name):
	print("{0} : Timestamp: {1:%Y-%m-%d %H:%M:%S}".format(name, datetime.datetime.now()))

def writing(desired,data,columns,source):
	f=open('table.'+source+'.tsv','w')
	f.write('#'+'\t'.join(desired))
	for line in data:
		ajusted=[str(line[columns.index(c)]) if c in columns else '' for c in desired]
		f.write('\n'+'\t'.join(ajusted))
	f.close()

def argsparse():
	parser=argparse.ArgumentParser(description="The parsing.py tool should be used with files from genomic island databases\nsuch as Islander (.sql) or ICEberg (.html) to parse them and extract the\ninformation into a tsv output file.\n\nThe user must enter the name of the database from where it comes from.\nThe user may indicate the number of base pairs that must be added to the\nends of the genomic island sequence.\n\nFor examples and more information, find the README document on the tool's github page.\n(https://github.com/whitef19/GI) ",formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('-i', action='store', dest='data', metavar='input', help='input file from database',required=True)
	parser.add_argument('-db', action='store', dest='db' , metavar='source', help='islander,paidb,iceberg,iv or tsv',required=True)
	parser.add_argument('-bp', action='store', dest='pad', metavar='base pairs',help='number of base pair to add',type=int, default=50)
	args=parser.parse_args()
	return args

def sequences(islands,index,bp):
	coordinate=open('no_coordinate.out','a')
	for island in islands:
		if (str(island[2])!='1') and ((int(island[3])-int(island[2]))>0):
			ID=island[0]+'-'+island[2]+'-'+island[3]
			sequence='sequences/island.'+ID+'.fa'
			if not os.path.isfile(sequence):
				interval=str(int(island[2])-bp)+'-'+str(int(island[3])+bp)
				genome='genomes/sequence.'+island[0]+'.fasta'
				head=os.popen('grep ">" '+genome).read()
				os.system('echo "'+str(head.replace('\n',''))+' '+str(island[0])+' '+str(interval)+'" >> '+str(sequence))
				os.system('grep -v ">" '+str(genome)+" | sed ':a;N;$!ba;s/\\n//g'| cut -c "+str(interval)+' >> '+str(sequence))
			if island[index]=='':
				with open(sequence,'r') as f:
					for line in f:
						if not line.startswith('>') :
							island[index]=line.replace('\n','')
		else :
			coordinate.write('\n'+'\t'.join(island))
	coordinate.close()
	return islands

def form(desired,data,bp):
	islands=[]
	for line in open(data,'r'):
		if not line.startswith('#'):
			line=line.replace('\n','').split('\t')
			islands.append(line)
	islands=sequences(islands,desired.index('SEQUENCE'),bp)
	writing(desired,islands,desired,'iv4')

def IV4(desired,data,org_dict,bp):
	islands=[]
	columns=['ACCESSION','ORGANISM','START','END','DETECTION','SEQUENCE']
	file=open(data,'r')
	for line in file:
		if not line.startswith('#'):
			line=line.replace('\n','').split('\t')
			organism=org_dict[line[0]] if line[0] in org_dict else None
			line=[line[0],organism]+line[1:]+['']
			islands.append(line)
	islands=sequences(islands,columns.index('SEQUENCE'),bp)
	writing(desired,islands,columns,'iv4')

def PAIDB(desired,data,acc_dict,bp):
	import re
	import urllib
	import requests as req

	islands=[]
	file=open(data,'r')
	table=SoupStrainer("table")
	soup=BeautifulSoup(file,'html.parser',parse_only=table).find_all(bordercolordark='white')
	for species in soup:
		for island in species.findAll(valign='top'):
			cells=island.findAll('td') # list
			if len(cells) >0:
				link=cells[1].find('a').get('href')
				strain=cells[2].text.strip()
				insertion=cells[4].text.strip()
				locus=cells[5].text.strip().split('(')[0].replace(' ','')
				end=int(float(cells[5].text.strip().split('(')[1].split('kb')[0])*1000)

				organism=' '.join(strain.split(' '))
				if organism in acc_dict:
					accession=acc_dict[organism]
				else:
					accession=''
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
				line=[accession,strain,'1',str(end),insertion,'PAIDB',','.join(reference),locus]
				islands.append(line)
	#islands=sequences(islands,desired.index('SEQUENCE'),bp)
	writing(desired,islands,desired,'PAIDB')

def ICEberg(desired,data,acc_dict,bp):
	islands=[]
	query=[]
	for page in range(1,467): # nombre de page
		line=[ '' for i in desired]
		line[5]='ICEberg'
		reference=[]
		file=open(data+'page_'+str(page)+'.html','r')
		soup=BeautifulSoup(file,'html.parser').findAll('tr')
		for row in soup:
			cells=row.findAll('td')
			if len(cells)>0 :
				if cells[0].text.strip() == 'Organism':
					line[1]=cells[1].text.strip()
				elif cells[0].text.strip() == 'Genome coordinates':
					line[2]=cells[1].text.strip().split('..')[0]
					line[3]=cells[1].text.strip().split('..')[1]
				elif cells[0].text.strip() == 'Insertion site':
					line[4]=cells[1].text.strip()
				elif cells[0].text.strip() == 'Nucleotide Sequence':
					line[7]=cells[1].text.strip()
					query.append(line[7].split(';')[0])
				elif cells[0].text.strip().startswith('This is a') :
					reference.append(cells[0].text.strip())
				elif cells[0].text.strip().startswith('(') :
					reference.append(cells[0].text.strip().split('[')[1].replace(']','').replace('PudMed','PMID'))

		line[6]=','.join(reference)
		organism=' '.join(line[1].split(' ')[:2]) if line[1]!='' else None
		if organism in acc_dict:
			line[0]=acc_dict[organism]
		islands.append(line)
	islands=sequences(islands,desired.index('SEQUENCE'),bp)
	writing(desired,islands,desired,'ICEberg')

def Islander(desired,data,bp):

	islands=[]
	tables=['`island_sequence`','`islander`','`literature_islands`','`island`']
	values=[[] for i in tables]
	columns=['ACCESSION','ORGANISM','START','END','DETECTION','REFERENCE','SEQUENCE']
	file=open(data,'r')
	for line in file:
		section=line.replace('\n','').split(' ')
		if 'INSERT' in section:
			if section[2] in tables:
				values[tables.index(section[2])]+=line.split('VALUES ')[1].split('),(')

	# island_sequence
	sequences={}
	for line in values[0]:
		line=line.replace('\'','').replace('(','').split(',')
		sequences[(line[0])]=line[1]

	# literature_islands
	IDs=[line.replace('\'','').split(',')[5] for line in values[2]]

	# island
	for line in values[3]:
		line=line.replace('(','').replace('\'','').split(',')
		if line[0] in IDs:
			if line[0] in sequences:
				islands.append([line[23],' '.join(line[24].split('_')),line[7],line[8],'islander','PMID:14681358',sequences[line[0]]])
			else:
				islands.append([line[23],' '.join(line[24].split('_')),line[7],line[8],'islander','PMID:14681358',''])
		elif line[0] in sequences:
			islands.append([line[23],' '.join(line[24].split('_')),line[7],line[8],'islander','',sequences[line[0]]])
		else:
			islands.append([line[23],' '.join(line[24].split('_')),line[7],line[8],'islander','',''])

	# islander
	for line in values[1]:
		line=line.replace('(','').replace('\'','').split(',')
		for i,ind in enumerate(line):
			if ind.startswith('NC_'):
				names=line[i:]
				break
		if line[0] in IDs:
			if line[0] in sequences:
				islands.append([names[0],names[2],line[7],line[8],'islander','PMID:14681358',sequences[line[0]]])
			else:
				islands.append([names[0],names[2],line[7],line[8],'islander','PMID:14681358',''])
		elif line[0] in sequences:
			islands.append([names[0],names[2],line[7],line[8],'islander','',sequences[line[0]]])
		else:
			islands.append([names[0],names[2],line[7],line[8],'islander','',''])

	islands=sequences(islands,columns.index('SEQUENCE'),bp)
	writing(desired,islands,columns,'Islander')


def main():
	accesion='index.accession_number.txt'
	desired=['ACCESSION','ORGANISM','START','END','INSERTION','DETECTION','REFERENCE','SEQUENCE']
	args=argsparse()
	bp=args.pad

	if args.db.lower() =='iv' :
		organisms={}
		for nb in open(accesion,'r'):
			organisms[nb.replace('\n','').split('\t')[0]]=nb.replace('\n','').split('\t')[1]
		IV4(desired,args.data,organisms,bp)

	if args.db.lower() =='paidb' :
		accession_nb={}
		for nb in open(accesion,'r'):
			accession_nb[nb.replace('\n','').split('\t')[1]]=nb.replace('\n','').split('\t')[0]
		PAIDB(desired,args.data,accession_nb,bp)

	if args.db.lower() =='iceberg' :
		accession_nb={}
		for nb in open(accesion,'r'):
			accession_nb[nb.replace('\n','').split('\t')[1]]=nb.replace('\n','').split('\t')[0]
		ICEberg(desired,args.data,accession_nb,bp)

	if args.db.lower() =='islander' :
		Islander(desired,args.data,bp)

	if args.db.lower() == 'tsv' :
		form(desired,args.data,bp)
	#if args.db.lower() =='xlsx':
		#excel(desired,args.data)

if __name__ == "__main__":
	timestamp("STARTING")
	main()
	timestamp("DONE")

"""

def excel(desired,data):
	df=pd.read_excel(data)
	columns=[c.upper() for c in df.columns] # upper case columns name
	df.columns=columns	#replace columns name for the upper one
	col=[c for c in desired if c in columns]
	df=df[col] # cut desired columns
	data_list=df.values.tolist()
	island=[Ilot(line,col,desired) for line in data_list]
	return island


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

"""



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
