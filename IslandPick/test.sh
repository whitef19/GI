#!/bin/bash

if [ ! -f sequence.NC_004431.1.fasta ] ; then	
	esearch -db Nucleotide -query "(NC_004431.1)"|efetch -format fasta > sequence.NC_004431.1.fasta
	grep ">" sequence.NC_004431.1.fasta |wc -l ;fi

header=` grep ">" sequence.NC_004431.1.fasta `
fasta=` grep ">" sequence.NC_004431.1.fasta |sed 's/,/\t/g'|cut -f1 |sed 's/ /\t/g' | cut -f2- |sed 's/\t/ /g' `

if [ "$fasta" = "Escherichia coli CFT073" ] ; then
	echo "$header 258023-271043" >sequence.crop.NC_004431.1-258073-270993.fasta
	grep -v ">" sequence.NC_004431.1.fasta |sed ':a;N;$!ba;s/\n//g'|cut -c 258023-271043 >>sequence.crop.NC_004431.1-258073-270993.fasta
else
	echo "NC_004431.1-258073-270993	Escherichia coli CFT073" >> error.not_found ; fi