#!/bin/bash


for gi in `cat $1`
do
	echo "$gi"
	python scholar.py -A "$gi" -A "genomic island"
done

echo "DONE"