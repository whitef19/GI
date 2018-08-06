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
* all_gis_islander_iv4.txt
* all_gis_islandpick_iv4.txt
* all_gis_islandpath_dimob_iv4.txt
* all_gis_sigi_hmm_iv4.txt
* all_gis_islandviewer_iv4.txt 
* Supp1: positive dataset 
* Supp2: negative dataset 
* Supp3: 80 GIs ranging in size from 3 to 133 kb
* Additional_file_1.xls
* Additional_file_2.xls
* Additional_file_3.xls
* Additional_file_8.xls

```bash

# Extraction de l'information des tables supplémentaires 1-3 
python3 parsing.py -i input/all_gis_islander_iv4.txt         -f txt -bn input/iv4_islander.pickle
python3 parsing.py -i input/all_gis_islandpick_iv4.txt       -f txt -bn input/iv4_islandpick.pickle
python3 parsing.py -i input/all_gis_islandpath_dimob_iv4.txt -f txt -bn input/iv4_islandpath.pickle
python3 parsing.py -i input/all_gis_sigi_hmm_iv4.txt         -f txt -bn input/iv4_sigi.pickle
python3 parsing.py -i input/all_gis_islandviewer_iv4.txt     -f txt -bn input/iv4.pickle
python3 parsing.py -i input/SupTable1.xlsx -f xlsx -bn input/iv4_supp1.pickle
python3 parsing.py -i input/SupTable2.xlsx -f xlsx -bn input/iv4_supp2.pickle
python3 parsing.py -i input/SupTable3.xlsx -f xlsx -bn input/iv4_supp3.pickle

# Ajout dans la base de données
python3 update.py -db database.pickle -i input/iv4_islander.pickle   -new
python3 update.py -db database.pickle -i input/iv4_islandpick.pickle
python3 update.py -db database.pickle -i input/iv4_supp1.pickle
python3 update.py -db database.pickle -i input/iv4_supp3.pickle

python3 update.py -db database.pickle -i input/iv4_islandpath.pickle 

python3 update.py -db database.pickle -i input/iv4_sigi.pickle
python3 update.py -db database.pickle -i input/iv4.pickle

# Création des scripts sh pour executer des recherches de génomes avec **EDirect** sur NCBI Nucleotide
python3 sequences.py -db database.pickle -o sequences/sequence.  

# Download des genomes
sh sequences.sh

# Sélection de l'interval +/- 50 pb
sh cropped.sequences.sh
```



