# Project GI

*Élaboration d'une base de donnée détaillée des îlots génomique à partir de bases de données existante tel que IslandViewer*
 - [ ] IslandViewer
 - [ ] PAI-DB
 - [ ] ICEberg
 - [ ] Islander

**Colonnes**
```
ID	ACC_NB	ORGANISME	START	STOP	SEQUENCE	SITE D'INSERTION*	REFERENCES DOI/PMID	EXP/DÉTECTION*
```

**IslandViewer**

[Article](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6022643/#sup1)
Inputs :
* Supp1: positive dataset 
* Supp2: negative dataset 
* Supp3: 80 GIs ranging in size from 3 to 133 kb
* Additional_file_1.xls
* Additional_file_2.xls
* Additional_file_3.xls
* Additional_file_8.xls

```python
# Extraction de l'information des tables supplémentaires 1-3 

python3 parsing.py -i input/SupTable1.xlsx -b islandpick_suptable1.pickle

# Création des scripts sh pour executer des recherches de génomes avec **EDirect** sur NCBI Nucleotide

python3 sequences.py -p islandpick_suptable1.pickle -o sequences/sequence. --crop -bp 50 

sh sequences.sh

sh cropped.sequences.sh

