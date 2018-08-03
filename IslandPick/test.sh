#!/bin/bash

echo "##Escherichia coli CFT073" >> sequence.NC_004431.1-258073-270993.fasta
echo "#NC_004431.1-258073-270993" >> sequence.NC_004431.1-258073-270993.fasta

esearch -db Nucleotide -query "(NC_004431.1)"|efetch -format fasta >> sequence.NC_004431.1-258073-270993.fasta

organism=grep "##" sequence.NC_004431.1-258073-270993.fasta |sed 's/##//g'
echo "$organism"
fasta=grep ">" sequence.NC_004431.1-258073-270993.fasta|sed 's/,/\t/g'|cut -f1 |sed 's/ /\t/g' | cut -f2- |sed 's/\t/ /g'
echo "$fasta"
if [$organism -eq $fasta]
then
	grep -v ">" sequence.NC_004431.1-258073-270993.fasta | grep -v "#" |sed ':a;N;$!ba;s/\n//g'|cut -c 258073-270993 >sequence.crop.NC_004431.1-258073-270993.fasta
fi
echo "DONE"