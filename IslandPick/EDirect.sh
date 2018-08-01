#!/bin/bash

for gi in `cat $1`
do
	echo "$gi"
	#esearch -db pubmed -query "$gi" -sort Relevance |efetch -format medline | grep -wf pattern >>output/EDirect_references
done

echo "DONE"
