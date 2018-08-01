#!/bin/bash

esearch -db pubmed -query "Baba et al., 2002, The Lancet" -sort Relevance |efetch -format medline | grep -wf pattern >>output/EDirect_references 
esearch -db pubmed -query "Beres et al., 2002, PNAS" -sort Relevance |efetch -format medline | grep -wf pattern >>output/EDirect_references 
esearch -db pubmed -query "Boyd et al., 2008, BMC Microbiology" -sort Relevance |efetch -format medline | grep -wf pattern >>output/EDirect_references 
esearch -db pubmed -query "Hayashi et al., 2001, DNA Res" -sort Relevance |efetch -format medline | grep -wf pattern >>output/EDirect_references 
esearch -db pubmed -query "Hurley et al., 2006, BMC Genomics" -sort Relevance |efetch -format medline | grep -wf pattern >>output/EDirect_references 
esearch -db pubmed -query "Lloyd et al., 2007 J Bacteriol" -sort Relevance |efetch -format medline | grep -wf pattern >>output/EDirect_references 
esearch -db pubmed -query "Okuda et al., 1997, J Clin Microbiol" -sort Relevance |efetch -format medline | grep -wf pattern >>output/EDirect_references 
esearch -db pubmed -query "Parkhill et al., 2001, Nature" -sort Relevance |efetch -format medline | grep -wf pattern >>output/EDirect_references 
esearch -db pubmed -query "Reen et al., 2006 Nat Rev Microbiol" -sort Relevance |efetch -format medline | grep -wf pattern >>output/EDirect_references 

echo "DONE"