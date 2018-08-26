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

# ICEberg 
[site web ](http://db-mml.sjtu.edu.cn/ICEberg/)   

**Inputs :**
* page_x.html 467 pages html de la base de données

```bash
# ICEberg
python3 parsing.py -i input/iceberg/ -pos input/SupTable1.xlsx -neg input/SupTable2.xlsx  -db iceberg -bp 50
```
*output : table.iceberg.tsv*
> Islands  :   468   
> Uniques  :   354    
> NA       :   258 (39 ont des souches inconnus ('Unknown strain'))    


# Islander 
[site web](https://bioinformatics.sandia.gov/islander/about.html)   

**Inputs:**
* island_sequence.06.03.2015.fasta
* islander.08.03.2015.sql

```bash
# Islander
python3 parsing.py -i input/islander.08.03.2015.sql -pos input/SupTable1.xlsx -neg input/SupTable2.xlsx  -db islander -bp 50
```
*output :table.islander.08.03.2015.sql.tsv*
> Islands  : 4 208   
> Uniques  : 4 208    
> NA       :     1   

# IslandViewer
[ site web ](http://www.pathogenomics.sfu.ca/islandviewer/)   
[publication](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6022643/#sup1)  

**Inputs :**
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
# IslandViewer 4 islander
python3 parsing.py -i input/all_gis_islander_iv4.txt -pos input/SupTable1.xlsx -neg input/SupTable2.xlsx  -db iv -bp 50
```
*output : table.all_gis_islander_iv4.txt.tsv*
> Islands  : 3 170   
> Uniques  : 3 170   
> NA       :     0    

```bash
# IslandViewer 4 islandpath dimob
python3 parsing.py -i input/all_gis_islandpath_dimob_iv4.txt -pos input/SupTable1.xlsx -neg input/SupTable2.xlsx  -db iv -bp 50
```
*output : table.all_gis_islandpath_dimob_iv4.txt.tsv*
> Islands  : 51 425    
> Uniques  : 51 425     
> NA       :      0   

```bash
# IslandViewer 4 islandpick
python3 parsing.py -i input/all_gis_islandpick_iv4.txt -pos input/SupTable1.xlsx -neg input/SupTable2.xlsx  -db iv -bp 50
```
*output :table.all_gis_islandpick_iv4.txt.tsv*
> Islands  : 43 469
> Uniques  : 43 469
> NA       :      0

```bash
# IslandViewer 4 islandviewer
python3 parsing.py -i input/all_gis_islandviewer_iv4.txt -pos input/SupTable1.xlsx -neg input/SupTable2.xlsx  -db iv -bp 50
```
*output :table.all_gis_islandviewer_iv4.txt.tsv*
> Islands  :110 914     
> Uniques  :110 914    
> NA       :      0   

```bash
# IslandViewer 4 sigi hmm
python3 parsing.py -i input/all_gis_sigi_hmm_iv4.txt -pos input/SupTable1.xlsx -neg input/SupTable2.xlsx  -db iv -bp 50
```
*output :table.all_gis_sigi_hmm_iv4.txt.tsv*
> Islands  : 63 917   
> Uniques  : 63 917    
> NA       :      0     

# PAIDB
[ site web ](http://www.paidb.re.kr/browse_pais.php?m=p)

**Inputs :**   
* PAIDB_REI.html : REsistance Islands    
* PAIDB_PAI.html : PAthogenicity Islands

```bash
# PAIDB PAI
python3 parsing.py -i input/PAIDB_PAI.html -pos input/SupTable1.xlsx -neg input/SupTable2.xlsx  -db paidb -bp 50
```
*output :table.PAIDB_PAI.html.tsv*
> Islands  :   986    
> Uniques  :   986    
> NA       :   904   

```bash
# PAIDB REI
python3 parsing.py -i input/PAIDB_REI.html -pos input/SupTable1.xlsx -neg input/SupTable2.xlsx  -db paidb -bp 50
```
*output :table.PAIDB_REI.html.tsv*
> Islands  :   108   
> Uniques  :   108   
> NA       :    95   


# Analyse 

Temps de calcul : voir log_parsing_1.0

```bash
 
# Nombre d'island 
wc -l table.all_gis_islander_iv4.txt.tsv
# Nombre d'island uniques
sort -k3 table.all_gis_islander_iv4.txt.tsv | sort -k1 |uniq | wc -l
# Island avec numéro d'accession et positions
awk '$1=="NA" || $3=="NA" || $4=="NA"' table.all_gis_islander_iv4.txt.tsv | wc -l
mv table.all_gis_islander_iv4.txt.tsv Tables/.
```


# Concatenation des sources 
```bash
# L'option "cat" permet donner 9 arguments et lance la fonction concatenate  
python3 parsing.py -cat Tables/*
```
*output: database_1.0.tsv*
> Islands  : 181 699   
> NA       :   9 221   

