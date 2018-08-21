# -*- coding: utf-8 -*-
# python 3.5
import os
import sys
import argparse
import datetime

def timestamp(name):
	print("{0} : Timestamp: {1:%Y-%m-%d %H:%M:%S}".format(name, datetime.datetime.now()))

def argsparse():
	parser=argparse.ArgumentParser(description='Fetch genomic sequences from Nucleotide')
	parser.add_argument('-i',  action='store', dest='input', help='input file with accession number', required=True)
	parser.add_argument('-o',  action='store', dest='sh',    help='.sh file to fetch genomes', default='get_genomes.sh')
	args=parser.parse_args()
	return args
	
def writing(query,output):
	o=open(output,'w') 
	o.write( '#!/bin/bash'+'\n' )
	for island in query:
		genome='genomes_new/sequence.'+island+'.fasta'
		o.write('\n'+ 'if [ ! -f '+genome+' ] ; then') 
		o.write('\n'+'\t'+'esearch -db Nucleotide -query "'+island+'" |efetch -format fasta >> '+genome)
		o.write('\n'+'fi'+'\n')
	o.write('\n'+'echo "DONE"')
	o.close()

def main():
	args=argsparse()
	query=[nb.replace('\n','') for nb in open(args.input,'r')]
	writing(query,args.sh)

if __name__ == "__main__":
	timestamp("STARTING")
	main()
	timestamp("DONE")

"""
def sequences(file):
	Islands=[]
	errors=[]
	with open (file,'r') as f:
		for island in f:
			if not island.startswith("#"):
				island=island.replace('\n','').split('\t')
				ID=island[0]+'-'+island[2]+'-'+island[3]

				header=[]
				error_fasta=False
				if os.path.isfile('sequences/island.'+ID+'.fa'):
					with open('sequences/island.'+ID+'.fa','r') as f:
						for line in f:
							if line.startswith('>') :
								header.append(line.replace('\n',''))
								if len(header)!=1:
									errors.append(ID)
									error_fasta=True

							elif (not error_fasta) and (island[2]!='1'): 
								if island[7]=='':
									island[7]=line.replace('\n','')
									Islands.append(island)
							else: 
								Islands.append(island)
				else:
					Islands.append(island)
					if island[7]!='':
						o=open('sequences/island.'+ID+'.fa','w')
						o.write('>'+island[0]+' '+island[1]+' db:'+island[5]+' '+ID)
						o.write('\n'+island[7])
						o.close()

	return Islands,errors
"""