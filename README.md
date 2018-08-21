# Project GI

*Élaboration d'une base de donnée détaillée des îlots génomiques à partir de bases de données existantes.*
 - [x] IslandViewer
 - [x] PAI-DB
 - [x] ICEberg
 - [x] Islander

```bash
# Colonnes de la base de données
ID	ACC_NB	ORGANISME	START	STOP	SEQUENCE	INSERTION*	REFERENCES  DETECTION*
```
-------------------
# Sources de données

**IslandViewer**   
[ site web ](http://www.pathogenomics.sfu.ca/islandviewer/)   
[publication](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6022643/#sup1)  

>>
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
>>

```bash
# Extraction de l'information et présentation uniforme dans des tables tsv avec les colonnes d'intéret

python3 parsing.py -i input/all_gis_islander_iv4.txt -db iv -acc accessionNb.organisms.txt
python3 parsing.py -i input/all_gis_islandpick_iv4.txt -db iv -acc accessionNb.organisms.txt
python3 parsing.py -i input/all_gis_islandpath_dimob_iv4.txt -db iv -acc accessionNb.organisms.txt
python3 parsing.py -i input/all_gis_sigi_hmm_iv4.txt -db iv -acc accessionNb.organisms.txt
python3 parsing.py -i input/all_gis_islandviewer_iv4.txt -db iv -acc accessionNb.organisms.txt

# Création des scripts sh pour executer des recherches de génomes avec EDirect sur NCBI Nucleotide
python3 sequences.py -db database.pickle -o sequences/get_sequence.sh  

# Téléchargement des genomes
sh get_sequences.sh
```

**PAI-DB**   
[ site web ](http://www.paidb.re.kr/browse_pais.php?m=p)
>>
Inputs :   
* PAIDB_REI.html : REsistance Islands    
* PAIDB_PAI.html : PAthogenicity Islands
>>

```bash
# Fonction PAIDB ajouté au script parsing.py
# Création de tables
python3 parsing.py -i PAIDB_PAI.html -db paidb -acc accessionNb.organisms.txt
python3 parsing.py -i PAIDB_REI.html -db paidb -acc accessionNb.organisms.txt
```

**ICEberg**   
[site web ](http://db-mml.sjtu.edu.cn/ICEberg/)   
>>
Inputs:
* page_x.html 467 pages html de la base de données
>>

```bash
# Fonction ICEberg ajouté au script parsing.py
python3 parsing.py -i input/ICEberg/ -db iceberg -acc accessionNb.organisms.txt
```

**Islander**   
[site web](https://bioinformatics.sandia.gov/islander/about.html)   
>>
Inputs:
* island_sequence.06.03.2015.fasta
* islander.08.03.2015.sql
>>

```bash
# Fonction Islander ajouté au script parsing.py
python3 parsing.py -i input/islander.08.03.2015.sql -db islander
```

# Base de données


>>
Inputs
* table.iv4.islander.txt
* table.iv4.islandpath.txt
* table.iv4.islandpick.txt
* table.iv4.sigi_hmm.txt
* table.iv4.islandviewer.txt
* table.paidb-pai.txt
* table.paidb-rei.txt
* table.iceberg.txt
* table.islander.txt
>>

```bash

cat Tables/table.* |sort -k2,2 |grep -v # >mastertable.txt

python3 update.py -i mastertable.txt -o database2.txt
>database.txt

python3 sequences.txt -i database.txt

sh get_sequences.sh

python3 update.py -i database2.txt -o database3.txt -seq


```

grep ">" genomes/* >header
