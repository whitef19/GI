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
### à refaire avec la liste d'organismes
# Extraction de l'information des tables supplémentaires 1-3 
python3 parsing.py -i input/all_gis_islander_iv4.txt         -f txt -bn input/iv4_islander.pickle
python3 parsing.py -i input/all_gis_islandpick_iv4.txt       -f txt -bn input/iv4_islandpick.pickle
python3 parsing.py -i input/all_gis_islandpath_dimob_iv4.txt -f txt -bn input/iv4_islandpath.pickle
python3 parsing.py -i input/all_gis_sigi_hmm_iv4.txt         -f txt -bn input/iv4_sigi.pickle
python3 parsing.py -i input/all_gis_islandviewer_iv4.txt     -f txt -bn input/iv4.pickle
python3 parsing.py -i input/SupTable1.xlsx -f xlsx -bn input/iv4_supp1.pickle
python3 parsing.py -i input/SupTable2.xlsx -f xlsx -bn input/iv4_supp2.pickle
python3 parsing.py -i input/SupTable3.xlsx -f xlsx -bn input/iv4_supp3.pickle



### à refaire avec l'information complète
# Ajout dans la base de données
python3 update.py -db database.pickle -i input/iv4_islander.pickle   -new
python3 update.py -db database.pickle -i input/iv4_islandpick.pickle
python3 update.py -db database.pickle -i input/iv4_supp1.pickle
python3 update.py -db database.pickle -i input/iv4_supp3.pickle
python3 update.py -db database.pickle -i input/iv4_islandpath.pickle 
python3 update.py -db database.pickle -i input/iv4_sigi.pickle
python3 update.py -db database.pickle -i input/iv4.pickle # (3 heures)

# Création des scripts sh pour executer des recherches de génomes avec **EDirect** sur NCBI Nucleotide
python3 sequences.py -db database.pickle -o sequences/get_sequence.  

# Download des genomes
sh get_sequences.sh

```




**PAI-DB**

Inputs :   
* PAIDB_REI.html : REsistance Islands    
* PAIDB_PAI.html : PAthogenicity Islands

```bash 
# Fonction PAIDB ajouté au script parsing.py
# Création de tables des informations "table.PAIDB-REI" , "table.PAIDB-PAI"

python3 IslandViewer/parsing.py -i PAIDB_PAI.html -db paidb -acc IslandViewer/accession-organisms.txt
> table.PAIDB-PAI.out
python3 IslandViewer/parsing.py -i PAIDB_REI.html -db paidb -acc IslandViewer/accession-organisms.txt
> table.PAIDB-REI.out
```



**ICEberg**
Inputs:
* page_x.html 467 pages html de la base de données
```bash
# Fonction ICEberg ajouté au script parsing
python3 IslandViewer/parsing.py -i ICEberg/ -db iceberg

```
