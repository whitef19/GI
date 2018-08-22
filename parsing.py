# -*- coding: utf-8 -*-
# python 3.5
import os
import sys
import csv
import logging as log
from Bio import SeqIO
import datetime
#import argparse
import click
import pandas as pd
from bs4 import BeautifulSoup,SoupStrainer
log.basicConfig(filename='parse.log',level=log.DEBUG,format='%(asctime)s %(message)s')

#concat: result = pd.concat([df1, df4], axis=1, sort=False)

def writing(desired,data,columns,source):
	# Essayer avec pandas to.csv

	#import json
	#with open('data.json') as f:
	#data = json.load(f)
	#with open('data.json', 'w') as outfile:
	#json.dump(data,outfile)

	f=open('table.'+source+'.tsv','w')
	f.write('#'+'\t'.join(desired))
	for line in data:
		ajusted=[str(line[columns.index(c)]) if c in columns else '' for c in desired]
		f.write('\n'+'\t'.join(ajusted))
	f.close()

def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    bases = list(seq) 
    bases = [complement[base] for base in bases] 
    return ''.join(bases)

def sequences(islands,basepairs):
	# Ajouter le brin reverse
	# Blast quand déjà séquence https://krother.gitbooks.io/biopython-tutorial/content/BLAST.html
	for island in islands:
		start = int(island[2])	
		end =  int(island[3])
		if start != 1:
			ID = '{}-{}-{}'.format(island[0], island[2], island[3])
			sequence_file = 'island_sequences/island.{}.fa'.format(ID)

			if not os.path.isfile(sequence_file):
				genome_file = 'genomes/sequence.{}.fasta'.format(island[0])
				genome = SeqIO.read(genome_file, "fasta")
				with open(sequence_file, "w") as out:
					if end-start > 0:
						SeqIO.write(genome[(start-basepairs):(end+basepairs)], out, "fasta")
					#else:
						#reverse = complement(genome[(end-basepairs):(start+basepairs)])
						#SeqIO.write([reverse], out, "fasta")
					"""
					cmd_header = 'echo "{} {} {}-{}" >> {}'.format(head, island[0], end-basepairs, start+basepairs, sequence_file)
					os.system(cmd_header)
					cmd_seq = 'grep -v ">" '+genome_file+" | sed ':a;N;$!ba;s/\\n//g'| cut -c {}-{}".format(end-basepairs, start+basepairs)
					sequence = os.popen(cmd_seq).read()
					reverse = complement(sequence)
					os.system('echo "{}" >> {}'.format(reverse, sequence_file))
					"""
			if island[8]=='NA':
				with open(sequence_file,'r') as f:
					for line in f:
						if not line.startswith('>') :
							island[8]=line.replace('\n','')
	return islands

def sequences_sed_cut(islands,bp):
	# Ajouter le brin reverse
	# Blast quand déjà séquence https://krother.gitbooks.io/biopython-tutorial/content/BLAST.html
	for island in islands:
		if (str(island[2])!='1'):
			ID=island[0]+'-'+island[2]+'-'+island[3]
			sequence='sequences_test/island.'+ID+'.fa'
			if not os.path.isfile(sequence):
				interval=str(int(island[2])-bp)+'-'+str(int(island[3])+bp)
				genome='genomes/sequence.'+island[0]+'.fasta'
				head=os.popen('grep ">" '+genome).read()
				os.system('echo "'+str(head.replace('\n',''))+' '+str(island[0])+' '+str(interval)+'" >> '+str(sequence))
				if ((int(island[3])-int(island[2]))>0):
					os.system('grep -v ">" '+str(genome)+" | sed ':a;N;$!ba;s/\\n//g'| cut -c "+str(interval)+' >> '+str(sequence))
				#else:
					#os.system('grep -v ">" '+str(genome)+" | sed ':a;N;$!ba;s/\\n//g'| cut -c "+(str(int(island[3])-bp)+'-'+str(int(island[2])+bp))+' >> '+str(sequence))
			if island[8]=='':
				with open(sequence,'r') as f:
					for line in f:
						if not line.startswith('>') :
							island[8]=line.replace('\n','')
	return islands


def form(desired,data,bp):
	islands=[]
	for line in open(data,'r'):
		if not line.startswith('#'):
			line=line.replace('\n','').split('\t')
			islands.append(line)
	islands=sequences(islands,desired.index('SEQUENCE'),bp)
	writing(desired,islands,desired,'iv4')

def island_viewer(columns_of_interest, file, accession_nb_index, basepairs):
	"""parsing of files from island viewer.
		file: txt file from the database.
	"""
	islands = []
	columns_in_db = [
				'ACCESSION',
				'ORGANISM',
				'START',
				'END',
				'DETECTION',
				]

	with open(file,'r') as f:
		for line in f:
			if not line.startswith('#'):
				island = line.replace('\n','').split('\t')
				organism = accession_nb_index[island[0]] if island[0] in accession_nb_index else 'NA'
				island.insert(1, organism)
				island = [island[columns_in_db.index(c)] if c in columns_in_db else 'NA' for c in columns_of_interest]
				islands.append(island)
	Islands = pd.DataFrame.from_records(islands, columns=columns_of_interest)

#	log.info('lecture terminé')
#	islands = sequences(islands, basepairs)
	log.info('seqIO')
	islands = sequences_sed_cut(islands, basepairs)
	log.info('terminé')
#	writing(desired,islands,columns_in_db,'iv4')

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


@click.command()
@click.option('--file', '-i', help='input file from database')
@click.option('--database', '-db', help='islander,paidb,iceberg,iv or tsv')
@click.option('--basepairs', '-bp', help='number of base pair to add', default=50, type=int)
def main(file, database, basepairs):
	log.info('INPUT : '+file)
	accession_nb_index = {nb.replace('\n','').split('\t')[0]:nb.replace('\n','').split('\t')[1] for nb in open('index.accession_number.txt','r')} 
	columns_of_interest = [
					'ACCESSION',
					'ORGANISM',
					'START',
					'END',
					'INSERTION',
					'DETECTION',
					'REFERENCE',
					'TYPE',
					'SEQUENCE',
					'OTHER'
					]

	if database.lower() == 'iv' :
		island_viewer(columns_of_interest, file, accession_nb_index, basepairs)

	if database.lower() =='paidb' :
		accession_nb={}
		for nb in open(accesion,'r'):
			accession_nb[nb.replace('\n','').split('\t')[1]]=nb.replace('\n','').split('\t')[0]
		PAIDB(desired,args.data,accession_nb,bp)

	if database.lower() =='iceberg' :
		accession_nb={}
		for nb in open(accesion,'r'):
			accession_nb[nb.replace('\n','').split('\t')[1]]=nb.replace('\n','').split('\t')[0]
		ICEberg(desired,args.data,accession_nb,bp)

	if database.lower() =='islander' :
		Islander(desired,args.data,bp)

	if database.lower() == 'tsv' :
		form(desired,args.data,bp)
	#if args.db.lower() =='xlsx':
		#excel(desired,args.data)

if __name__ == "__main__":
	log.info('STARTING')
	main()
	log.info('DONE')
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
"""