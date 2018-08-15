# Project GI

*Élaboration d'une base de donnée détaillée des îlots génomique à partir de bases de données existantes.*
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
# Extraction de l'information et présentation uniforme dans des tables
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


# Arranger les génomes accessionNb_fasta_dup

sed -n '311259,352653p' genomes/sequence.NC_006840.2.fasta > genomes_fix/sequence.NC_006840.2.fasta

## Fait
     16 NC_002971.3 - - # genomes manquants, fichier contient mauvais génomes
    368 NC_000962.3	ok 2393945,2456967
  19821 NC_006840.2	sed -n '311259,352653p' genomes/sequence.NC_006840.2.fasta > genomes_fix/sequence.NC_006840.2.fasta
   1610 NC_003197.2	sed -n '817599,886992p' genomes/sequence.NC_003197.2.fasta > genomes_fix/sequence.NC_003197.2.fasta
    585 NC_003143.1	sed -n '65339,131821p' genomes/sequence.NC_003143.1.fasta > genomes_fix/sequence.NC_003143.1.fasta
    490 NC_009613.3	sed -n '109210,150073p' genomes/sequence.NC_009613.3.fasta > genomes_fix/sequence.NC_009613.3.fasta
    379 NC_002162.1	sed -n '28680,39419p' genomes/sequence.NC_002162.1.fasta > genomes_fix/sequence.NC_002162.1.fasta
    354 NC_006510.1	sed -n '47599,98239p' genomes/sequence.NC_006510.1.fasta > genomes_fix/sequence.NC_006510.1.fasta
    241 NC_015312.1	sed -n '202673,304053p' genomes/sequence.NC_015312.1.fasta > genomes_fix/sequence.NC_015312.1.fasta
    175 NC_009089.1	sed -n '123189,184479p' genomes/sequence.NC_009089.1.fasta > genomes_fix/sequence.NC_009089.1.fasta
    157 NC_000964.3	sed -n '56048,116271p' genomes/sequence.NC_000964.3.fasta > genomes_fix/sequence.NC_000964.3.fasta
    132 NC_007434.1	sed -n '1,58949p' genomes/sequence.NC_007434.1.fasta > genomes_fix/sequence.NC_007434.1.fasta
    129 NC_013446.2	sed -n '82854,159621p' genomes/sequence.NC_013446.2.fasta > genomes_fix/sequence.NC_013446.2.fasta
     83 NC_009656.1	sed -n '1,94121p' genomes/sequence.NC_009656.1.fasta > genomes_fix/sequence.NC_009656.1.fasta
     60 NC_002488.3	sed -n '37262,75538p' genomes/sequence.NC_002488.3.fasta > genomes_fix/sequence.NC_002488.3.fasta
     44 NZ_CP012275.1	sed -n '123070,153776p' genomes/sequence.NZ_CP012275.1.fasta > genomes_fix/sequence.NZ_CP012275.1.fasta
     45 NZ_CP014873.1	sed -n '154408,192587p' genomes/sequence.NZ_CP014873.1.fasta > genomes_fix/sequence.NZ_CP014873.1.fasta
     43 NC_007086.1	sed -n '1,73554p' genomes/sequence.NC_007086.1.fasta > genomes_fix/sequence.NC_007086.1.fasta
     21 NC_010999.1	sed -n '88626,132615p' genomes/sequence.NC_010999.1.fasta > genomes_fix/sequence.NC_010999.1.fasta
     21 NC_004668.1	sed -n '161532,207504p' genomes/sequence.NC_004668.1.fasta > genomes_fix/sequence.NC_004668.1.fasta
     14 NC_011770.1	sed -n '1,94312p' genomes/sequence.NC_011770.1.fasta > genomes_fix/sequence.NC_011770.1.fasta
     12 NC_010612.1	sed -n '1,94813p' genomes/sequence.NC_010612.1.fasta > genomes_fix/sequence.NC_010612.1.fasta
      9 NC_002516.2	sed -n '1,89493p' genomes/sequence.NC_002516.2.fasta > genomes_fix/sequence.NC_002516.2.fasta
      7 NC_002620.2	sed -n '45988,61316p' genomes/sequence.NC_002620.2.fasta > genomes_fix/sequence.NC_002620.2.fasta


## À faire 
      8 NC_008563.1
      9 NC_008531.1
      9 NC_010117.1
      7 NC_010572.1
      7 NC_017196.2
      7 NC_017638.1
      7 NC_021214.1
      7 NC_012590.1
      7 NC_007296.1
      7 NC_007793.1
      6 NC_007613.1
      5 NC_002944.2
      5 NC_002978.6
      5 NC_007384.1
      5 NC_009725.1
      5 NC_011528.1
      5 NC_012560.1
      5 NC_016856.1
      5 NC_018101.1
      5 NC_022569.1
      4 NC_009648.1
      3 NC_007432.1
      3 NC_007492.2
      3 NC_007530.2
      3 NC_007626.1
      3 NC_003919.1
      3 NC_003997.3
      3 NC_004461.1
      3 NC_004603.1
      3 NC_004605.1
      3 NC_005126.1
      3 NC_005957.1
      3 NC_001318.1
      3 NC_008095.1
      3 NC_008319.1
      3 NC_008570.1
      3 NC_009614.1
      3 NC_009632.1
      3 NC_009882.1
      3 NC_010103.1
      3 NC_010104.1
      3 NC_010184.1
      3 NC_010397.1
      3 NC_010473.1
      3 NC_010501.1
      3 NC_010610.1
      3 NC_011149.1
      3 NC_011978.1
      3 NC_012892.2
      3 NC_012925.1
      3 NC_012926.1
      3 NC_012973.1
      3 NC_013209.1
      3 NC_014171.1
      3 NC_014329.1
      3 NC_014639.1
      3 NC_015725.1
      3 NC_015957.1
      3 NC_016810.1
      3 NC_017250.1
      3 NC_017251.1
      3 NC_017331.1
      3 NC_017522.1
      3 NC_017549.1
      3 NC_017646.1
      3 NC_017660.1
      3 NC_017730.3
      3 NC_018221.1
      3 NC_018525.1
      3 NC_019556.1
      3 NC_020819.1
      3 NC_021251.1
      3 NC_021725.1
      3 NC_022198.1
      3 NZ_CP008816.1
      3 NZ_CP009709.1
      2 NC_006368.1
      2 NC_006814.3
      2 NC_004741.1
      2 NC_006932.1
      2 NC_002952.2
      2 NC_003028.3
      2 NC_003098.1
      2 NC_007946.1
      2 NC_008258.1
      2 NC_012207.1
      2 NC_015758.1
      2 NC_016445.1
      2 NC_016446.1
      2 NC_016804.1
      2 NC_016937.1
      2 NC_017625.1
# Modifier parsing.py pour ajouter les séquences directement avant l'écriture, 'SEQUENCE' doit être ajouter aux liste 'columns'
# Ajouter les fasta manquants
# Refaire les tables avec le bon accessionNB.organisms.txt et en ajoutant tous de suite les séquences.
python3 parsing.py -i input/all_gis_islander_iv4.txt -db iv -acc accessionNb.organisms.txt
python3 parsing.py -i input/all_gis_islandpick_iv4.txt -db iv -acc accessionNb.organisms.txt
python3 parsing.py -i input/all_gis_islandpath_dimob_iv4.txt -db iv -acc accessionNb.organisms.txt
python3 parsing.py -i input/all_gis_sigi_hmm_iv4.txt -db iv -acc accessionNb.organisms.txt
python3 parsing.py -i input/all_gis_islandviewer_iv4.txt -db iv -acc accessionNb.organisms.txt
python3 parsing.py -i PAIDB_PAI.html -db paidb -acc accessionNb.organisms.txt
python3 parsing.py -i PAIDB_REI.html -db paidb -acc accessionNb.organisms.txt
python3 parsing.py -i input/ICEberg/ -db iceberg -acc accessionNb.organisms.txt

Ensuite on peut rouler update.py pour combiner les sources en une seule table. 

# Références