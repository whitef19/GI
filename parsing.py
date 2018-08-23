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
log.basicConfig(filename='log_parsing',level=log.DEBUG,format='%(asctime)s %(message)s')

#concat: result = pd.concat([df1, df4], axis=1, sort=False)

def writing(columns_of_interest, islands, source):
	Islands = pd.DataFrame.from_records(islands, columns=columns_of_interest)
	Islands.to_csv(source)


	# Essayer avec pandas to.csv

	#import json
	#with open('data.json') as f:
	#data = json.load(f)
	#with open('data.json', 'w') as outfile:
	#json.dump(data,outfile)
	"""
	f=open('table.'+source+'.tsv','w')
	f.write('#'+'\t'.join(desired))
	for line in data:
		ajusted=[str(line[columns.index(c)]) if c in columns else '' for c in desired]
		f.write('\n'+'\t'.join(ajusted))
	f.close()
	"""
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


def island_viewer(columns_of_interest, file, accession_nb_index, basepairs):
	""" parsing of files from island viewer.
		
		file: txt file from the database.
		# NC_009925.1 2268064 2280662 Islander
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
	
	islands = sequences(islands, basepairs)
	Islands = pd.DataFrame.from_records(islands, columns=columns_of_interest)
	Islands.to_csv('table.{}.csv'.format(file.split('/')[1]))

def pai_db(columns_of_interest, file, strain_index, basepairs):
	""" parsing of files from pai-db, pai and rei
		extract only the table sections 

		file: html file from the web page of the database.
		<td width=200><a href=view_pai.php?pn=Bar.G1.NC_010161_P1&m=p>Not named</a></td>
		<td width=150><I>Bartonella tribocorum</I> CIP 105476</td>
		<td width=200>Type IV secretion system</td>
		<td width=120 align=left><B>NC_010161_P1</B> (25.6kb, complete PAI in the sequenced genome)<BR></td>
	"""
	import re
	import urllib
	import requests as req

	islands = []
	with open (file, 'r') as f:
		subset = SoupStrainer("table")
		tables = BeautifulSoup(f, 'html.parser', parse_only=subset).find_all(bordercolordark='white')
		for species in tables:
			for island in species.findAll(valign='top'): # only table with data, not header
				cells = island.findAll('td') 
				if len(cells) >0:
					line = ['NA' for c in columns_of_interest]
					line[1] = cells[2].text.strip()
					line[4] = cells[4].text.strip()
					line[9] = cells[5].text.strip().split('(')[0].replace(' ','')
					line[0] = strain_index[line[1]] if line[1] in strain_index else 'NA'
					line[5] = 'paidb'
					link = cells[1].find('a').get('href')

					#references
					page = req.get('http://www.paidb.re.kr/{}'.format(link)) # crawling island page
					page_file = page.text
					page_tables = BeautifulSoup(page_file,'html.parser').select('table')
					references = []
					publications = page_tables[1].findAll('i')
					for ref in publications: 
						ref = ref.text.strip().split(' ')
						if 'PUBMED' in ref:
							references.append('PMID:{}'.format(ref[ref.index('PUBMED')+1]))
					line[6] = ','.join(references)
					islands.append(line)

	#islands = sequences(islands, basepairs)
	Islands = pd.DataFrame.from_records(islands, columns=columns_of_interest)
	Islands.to_csv('table.{}.csv'.format(file.split('/')[1]))

def iceberg(columns_of_interest, file, basepairs):
	""" parsing of html pages from ICEberg 
		file: directory of html files from the web page of the database.
	"""
	islands = []
	query = []
	keywords = {'Organism':'ORGANISM','Insertion site':'INSERTION','Nucleotide Sequence':'OTHER','Replicon':'ACCESSION','Genome coordinates':'START'}
	for page in range(1,468): # page number in directory
		with open('{}page_{}.html'.format(file, page), 'r') as f:
			data = BeautifulSoup(f, 'html.parser').select('table')
			if len(data) > 5: # empty page
				table = data[0].select('td')
				references = data[4].select('a')
				information = [keywords[column.text.strip()] if column.text.strip() in keywords else column.text.strip() for column in table ]
				information = information + ['REFERENCE', ','.join([ref.text.strip().replace('PudMed','PMID') for ref in references])]
				line = [information[information.index(column)+1] if (column in information) and (information[information.index(column)+1]) != '-' else 'NA' for column in columns_of_interest]
				line[5] = 'ICEberg'
				if line[2] != 'NA':
					line[3] = line[2].split('..')[1]
					line[2] = line[2].split('..')[0]
				if line[0] != 'NA':
					line[0] = line[0].split('[')[1].replace(']','')
					print(line[0])
				line[9] = ','.join([col for col in information if col.startswith('This')])
		islands.append(line)
		
	Islands = pd.DataFrame.from_records(islands, columns=columns_of_interest)
	Islands.to_csv('table.{}.tsv'.format(file.split('/')[1]), seq='\t')

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

	#islands = sequences(islands, basepairs)
	writing(desired,islands,columns,'Islander')


@click.command()
@click.option('--file', '-i', help='input file from database')
@click.option('--database', '-db', help='islander,paidb,iceberg,iv or tsv')
@click.option('--basepairs', '-bp', help='number of base pair to add', default=50, type=int)
def main(file, database, basepairs):
	log.info('INPUT : '+file)
	accession_nb_index = {nb.replace('\n','').split('\t')[0]:nb.replace('\n','').split('\t')[1] for nb in open('index.accession_number.txt','r')} 
	strain_index = {nb.replace('\n','').split('\t')[1]:nb.replace('\n','').split('\t')[0] for nb in open('index.accession_number.txt','r')} 
	columns_of_interest = [
					'ACCESSION', #0
					'ORGANISM',  #1
					'START',     #2
					'END',       #3
					'INSERTION', #4
					'DETECTION', #5
					'REFERENCE', #6
					'TYPE',      #7
					'SEQUENCE',  #8
					'OTHER'      #9
					]

	if database.lower() == 'iv' :
		island_viewer(columns_of_interest, file, accession_nb_index, basepairs)

	if database.lower() =='paidb' :
		pai_db(columns_of_interest, file, strain_index, basepairs)

	if database.lower() =='iceberg' :
		iceberg(columns_of_interest, file, basepairs)

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

def form(desired,data,bp):
	islands=[]
	for line in open(data,'r'):
		if not line.startswith('#'):
			line=line.replace('\n','').split('\t')
			islands.append(line)
	islands=sequences(islands,desired.index('SEQUENCE'),bp)
	writing(desired,islands,desired,'iv4')

"""