# -*- coding: utf-8 -*-
# python 3.5
import os
import sys
import csv
import click
import pandas as pd
import logging as log
from Bio import SeqIO
from bs4 import BeautifulSoup,SoupStrainer
log.basicConfig(filename='log_parsing',level=log.DEBUG,format='%(asctime)s %(message)s')

def sequences(islands,basepairs):
	validation = []
	for island in islands:
		ID = '{}-{}-{}'.format(island[0], island[2], island[3])
		sequence_file = 'island_sequences/island.{}.fa'.format(ID)
		if not os.path.isfile(sequence_file):
			start = int(island[2])	
			end =  int(island[3]) 
			genome_file = 'genomes/sequence.{}.fasta'.format(island[0])
			genome = SeqIO.read(genome_file, "fasta")
			with open(sequence_file, "w") as outfile:
				if end-start > 0:
					sequence = str(genome.seq[(start-basepairs):(end+basepairs)])
					outfile.write('>{} {} {}-{}'.format(genome.description, island[0], start-basepairs, end+basepairs))
					outfile.write('\n{}'.format(sequence))
				else:
					reverse = complement(genome.seq[(end-basepairs):(start+basepairs)])
					outfile.write('>{} {} {}-{}'.format(genome.description, island[0], start+basepairs, end-basepairs))
					outfile.write('\n{}'.format(reverse))
		
		with open(sequence_file,'r') as f:
			for line in f:
				if not line.startswith('>'):
					if island[9]=='NA':
						island[9] = line.replace('\n','')
					else :
						validation.append([ID, island[5], island[9], line.replace('\n','')]) # to compare sequences

	with open('sequences_validation.txt', 'a') as out:
		for seq in validation:
			out.write('\n'+'\t'.join(seq))
	return islands

def complement(seq):
    complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'} 
    bases = list(seq) 
    bases = [complement[base] for base in bases] 
    bases = bases[::-1] # reverse
    return ''.join(bases)

def excel(columns_of_interest,file):
	df = pd.read_excel(file, engine='xlrd')
	df.columns = [c.upper() for c in df.columns] #replace columns name for the upper one
	columns = [c for c in columns_of_interest if c in df.columns] # order
	data = (df[columns]).values.tolist() # cut desired columns
	islands ={'{}-{}-{}'.format(d[0], d[2], d[3]):d[1] for d in data }
	return islands

def island_viewer(columns_of_interest,file,positif_dataset,negative_dataset,accession_nb_index,basepairs):
	""" parsing of files from island viewer.
		
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
				island[7] = 'VP' if '{}-{}-{}'.format(island[0], island[2], island[3]) in positif_dataset else 'NA'
				island[7] = 'NEG' if '{}-{}-{}'.format(island[0], island[2], island[3]) in negative_dataset else 'NA'
				islands.append(island)
	islands = sequences(islands, basepairs)
	Islands = pd.DataFrame.from_records(islands, columns=columns_of_interest, index='ACCESSION')
	Islands.to_csv('table.{}.tsv'.format(file.split('/')[1]), sep='\t')

def pai_db(columns_of_interest,file,positif_dataset,negative_dataset,strain_index,basepairs):
	""" parsing of files from pai-db, pai and rei
		extract only the table sections 

		file: html file from the web page of the database.
		No coordinate in this database
	"""
	import re
	import urllib
	import requests as req
	islands = []
	with open (file, 'r') as f:
		subset = SoupStrainer("table")
		tables = BeautifulSoup(f, 'html.parser', parse_only=subset).find_all(bordercolordark='white')
		for species in tables:
			for row in species.findAll(valign='top'): # only table with data, not header
				cells = row.findAll('td') 
				if len(cells) >0:
					island = ['NA' for c in columns_of_interest]
					island[1] = cells[2].text.strip()
					island[4] = cells[4].text.strip()
					island[8] = cells[5].text.strip().split('(')[0].replace(' ','')
					island[0] = strain_index[island[1]] if island[1] in strain_index else 'NA'
					island[5] = 'paidb'
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
					island[6] = ','.join(references)
					island[7] = 'VP' if '{}-{}-{}'.format(island[0], island[2], island[3]) in positif_dataset else 'NA'
					island[7] = 'NEG' if '{}-{}-{}'.format(island[0], island[2], island[3]) in negative_dataset else 'NA'
					islands.append(island)

	#islands = sequences(islands, basepairs)
	Islands = pd.DataFrame.from_records(islands, columns=columns_of_interest, index='ACCESSION')
	Islands.to_csv('table.{}.tsv'.format(file.split('/')[1]), sep='\t')

def iceberg(columns_of_interest,path,positif_dataset,negative_dataset,strain_index,basepairs):
	""" parsing of html pages from ICEberg 
		file: directory of html files from the web page of ICEberg.
		Some island miss the coordinate in this database
	"""
	islands = []
	keywords = {'Organism':'ORGANISM','Insertion site':'INSERTION','Nucleotide Sequence':'OTHER','Replicon':'ACCESSION','Genome coordinates':'START'}
	for page in range(1,468): # page number in directory
		with open('{}page_{}.html'.format(path, page), 'r') as f:
			html_tables = BeautifulSoup(f, 'html.parser').select('table')
			if len(html_tables) > 5: # empty page
				table = html_tables[0].select('td')
				references = html_tables[4].select('a')
				information = [keywords[column.text.strip()] if column.text.strip() in keywords else column.text.strip() for column in table ]
				information = information + ['REFERENCE', ','.join([ref.text.strip().replace('PudMed','PMID') for ref in references])]
				island = [information[information.index(column)+1] if (column in information) and (information[information.index(column)+1]) != '-' else 'NA' for column in columns_of_interest]
				island[5] = 'ICEberg'
				if island[0] != 'NA':
					island[0] = island[0].split('[')[1].replace(']','')
				else:
					island[0] = strain_index[island[1]] if island[1] in strain_index else 'NA' 
				if island[2] != 'NA':
					island[3] = island[2].split('..')[1]
					island[2] = island[2].split('..')[0]
				island[8] = island[8]+','.join([col for col in information if col.startswith('This')])
				island[7] = 'VP' if '{}-{}-{}'.format(island[0], island[2], island[3]) in positif_dataset else 'NA'
				island[7] = 'NEG' if '{}-{}-{}'.format(island[0], island[2], island[3]) in negative_dataset else 'NA'
		islands.append(island)
	Islands = pd.DataFrame.from_records(islands, columns=columns_of_interest, index='ACCESSION')
	Islands.to_csv('table.iceberg.tsv', sep='\t')

def islander(columns_of_interest,file,positif_dataset,negative_dataset,basepairs):
	islands = []
	tables = ['`island_sequence`','`islander`','`literature_islands`','`island`']
	values = [[] for i in tables]
	with open(file, 'r') as f:
		for line in f:
			section = line.replace('\n','').split(' ')
			if 'INSERT' in section:
				if section[2] in tables:
					values[tables.index(section[2])] += line.split('VALUES ')[1].split('),(')
	# island_sequence
	sequences = {line.replace('\'','').replace('(','').split(',')[0] : line.replace('\'','').replace('(','').split(',')[1] for line in values[0]}
	# literature_islands
	reference = [line.replace('\'','').split(',')[5] for line in values[2]] # only one reference, list of island ID with the reference
	# island
	for line in values[3]:
		line = line.replace('(','').replace('\'','').split(',')
		columns_index = {'ACCESSION':line[23],'ORGANISM':' '.join(line[24].split('_')), 'START':line[7], 'END':line[8], 'DETECTION': 'islander'}
		island = [columns_index[col] if col in columns_index else 'NA' for col in columns_of_interest]
		island[6] = 'PMID:14681358' if island[0] in reference else 'NA'
		island[9] = sequences[island[0]] if island[0] in sequences else 'NA'
		island[7] = 'VP' if '{}-{}-{}'.format(island[0], island[2], island[3]) in positif_dataset else 'NA'
		island[7] = 'NEG' if '{}-{}-{}'.format(island[0], island[2], island[3]) in negative_dataset else 'NA'
		islands.append(island)
	# islander
	for line in values[1]:
		line = line.replace('(','').replace('\'','').split(',')
		strain = line[line.index(''.join([element for element in line if element.startswith('NC_')])):] # different lenght of list
		columns_index = {'ACCESSION':strain[0],'ORGANISM':strain[2], 'START':line[7], 'END':line[8], 'DETECTION': 'islander'}
		island = [columns_index[col] if col in columns_index else 'NA' for col in columns_of_interest]
		island[6] = 'PMID:14681358' if island[0] in reference else 'NA'
		island[9] = sequences[island[0]] if island[0] in sequences else 'NA'
		island[7] = 'VP' if '{}-{}-{}'.format(island[0], island[2], island[3]) in positif_dataset else 'NA'
		island[7] = 'NEG' if '{}-{}-{}'.format(island[0], island[2], island[3]) in negative_dataset else 'NA'
		islands.append(island)

	Islands = pd.DataFrame.from_records(islands, columns=columns_of_interest, index='ACCESSION')
	Islands.to_csv('table.{}.tsv'.format(file.split('/')[1]), sep='\t')

def uniq(columns_of_interest,file): #concat: result = pd.concat([df1, df4], axis=1, sort=False)
	"""
	files = (pd.concat([pd.read_csv(f, sep='\t', index_col='ACCESSION' ) for f in file ])).reset_index().values.tolist()
	Islands=[]
	IDs=[]
	sources=['ICEberg','islander','PAIDB','Dimob','islandviewer','Islander','Sigi','Islandpick']
	islands=[{} for i in sources]
	for line in files:
		if 'NA' not in line[:4]:
			source = line[5]
			ID = '{}-{}-{}'.format(line[0],line[2],line[3])
			islands[sources.index(source)][ID] = (line)
			IDs.append(ID)
		else :
			Islands.append(line)

	for ID in list(set(IDs)):
		table = [source[ID] for source in islands if ID in source]
		if len(table)>1:
			colums = [(list(set([row[i] for row in table]))) for i in range(8) ]
			line = [','.join(col) if idx!=7 else col[0] for idx,col in enumerate(colums)]
			Islands.append(line)
		else: 
			Islands.append(table[0])

	Islands = pd.DataFrame.from_records(Islands, columns=columns_of_interest)
	Islands.to_csv('database.tsv', sep='\t')


################################
"""
	
	ilot = []
	for db in file:
		with open (db,'r') as f:
			for line in f:
				if not line.startswith("ACCESSION"):
					ilot.append(line.replace('\n','').split('\t'))

	Islands=[]
	IDs=[]
	sources=['ICEberg','islander','PAIDB','Dimob','islandviewer','Islander','Sigi','Islandpick']
	islands=[{} for i in sources]
	for line in ilot:
		if ('' not in line[:4]) and (len(line)==8):
			source=line[5]
			ID=line[0]+'-'+line[2]+'-'+line[3]
			islands[sources.index(source)][ID]=(line)
			IDs.append(ID)
		else:
			Islands.append(line)
	for ID in list(set(IDs)):
		table=[source[ID] for source in islands if ID in source]
		if len(table)>1:
			colums=[(list(set([row[i] for row in table]))) for i in range(8) ]
			line=[','.join(col) if idx!=7 else col[0] for idx,col in enumerate(colums)]
			Islands.append(line)
		else: 
			Islands.append(table[0])

	Islands = (Islands)
	Islands = pd.DataFrame.from_records(Islands, columns=columns_of_interest)
	Islands.to_csv('database.tsv', sep='\t')


###################################
# OPTIONS
###################################
@click.command()
@click.option('--file', '-i', help='input file from database', multiple=True)
@click.option('--positif', '-pos', help='positif dataset')
@click.option('--negative', '-neg', help='negative dataset')
@click.option('--database', '-db', help='islander,paidb,iceberg,iv or tsv')
@click.option('--basepairs', '-bp', help='number of base pair to add', default=50, type=int)
@click.option('--cat', '-cat', help='concatenate all the csv', is_flag=True)
def main(file, positif, negative, database, basepairs, cat):
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
					'OTHER',      #8
					'SEQUENCE'  #9
					]

	positif_dataset = excel(columns_of_interest, positif)
	negative_dataset = excel(columns_of_interest, negative)

	if database:
		log.info('INPUT : '+file[0])

		if database.lower() == 'iv' :
			island_viewer(columns_of_interest, file[0], positif_dataset, negative_dataset, accession_nb_index, basepairs)

		if database.lower() =='paidb' :
			pai_db(columns_of_interest, file[0], positif_dataset, negative_dataset, strain_index, basepairs)

		if database.lower() =='iceberg' :
			iceberg(columns_of_interest, file[0], positif_dataset, negative_dataset, strain_index, basepairs)

		if database.lower() =='islander' :
			islander(columns_of_interest, file[0], positif_dataset, negative_dataset, basepairs)

	if cat :
		uniq(columns_of_interest, file)

	log.info('DONE')

if __name__ == "__main__":
	log.info('STARTING')
	main()
	


"""
	#if args.db.lower() =='xlsx':
		#excel(desired,args.data)
def excel(desired,data):
	df=pd.read_excel(data)
	columns=[c.upper() for c in df.columns] # upper case columns name
	df.columns=columns	#replace columns name for the upper one
	col=[c for c in desired if c in columns]
	df=df[col] # cut desired columns
	data_list=df.values.tolist()
	island=[Ilot(line,col,desired) for line in data_list]
	return island

	if database.lower() == 'tsv' :
		form(desired,args.data,bp)
def form(desired,data,bp):
	islands=[]
	for line in open(data,'r'):
		if not line.startswith('#'):
			line=line.replace('\n','').split('\t')
			islands.append(line)
	islands=sequences(islands,desired.index('SEQUENCE'),bp)
	writing(desired,islands,desired,'iv4')

"""