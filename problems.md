# Problemes rencontres

# ISLANDVIEWER
Documentation Pandas [https://pandas.pydata.org/pandas-docs/stable/generated/pandas.read_excel.html]   
Lecture des fichiers binaires SupTable1.xlsx SupTable2.xlsx avec pandas.read_csv.

# PAIDB 
[BeautifulSoup](https://www.crummy.com/software/BeautifulSoup/bs4/doc/#searching-the-tree)   
[Table HTML](https://gist.github.com/phillipsm/0ed98b2585f0ada5a769)

*Contient plusieurs tables, dont l'entête qui ne nous interesse pas.*   
> Pour se débarasser de l'entête : .find_all(bordercolordark='white')   
> Pour avoir seulement les données sans le nom des colonnes : .findAll(valign='top') 

# ICEBERG

*Quand on sort toute les objets tables des fichiers HTML, on peut juste sortir jusqu'à la ligne "organism".*
> C'est causé par une erreur dans le script HTML, à la fin de la ligne il y a deux balises pour fermer "a" </a></a>
```bash
# script python qui génère un script bash pour retirer la balise de trop
python3 iceberg_fixed.py
sh iceberg_fixed.sh
```


# ISLANDER

*plusieurs section, pas égal, pas la même information ni les mêmes îlots.*
> choix des sections les plus pertinantes pour couvrir les colonnes d'intérêt.

- [ ] accession_ids
- [ ] accessions
  * id
  * acc_id
  * start_stop
  * accession
  * tmrna_id
- [ ] bioprojects
- [ ] false_positives
- [ ] high_phylogeny
- [x] island
  * Island
  * Bacteria_nickname
  * tRNA
  * Portion
  * Int_Subfam
  * GC_Content
  * IR_Version
  * Island_Genome_L
  * Island_Genome_R
  * Orientation 
  * tRNA_L
  * tRNA_R
  * Dupli_L
  * Dupli_R
  * Frag_L
  * Frag_R
  * Damage_L
  * Damage_R
  * Extend
  * Mismatch
  * Index_Point
  * Comment
  * Markup
  * nc_number
  * ftp_folder
- [ ] island_annotation
- [ ] island_markup
- [ ] island_search
- [x] island_sequence
  * island
  * sequence
- [x] islander
  * island
  * bact_nickname
  * tRNA
  * portion
  * int_family
  * gc_content
  * ir_version
  * island_genome_l
  * island_genome_r
  * orientation
  * tRNA_L
  * tRNA_R
  * Dupli_L
  * Dupli_R
  * Frag_L
  * Frag_R
  * Mismatch
  * IndexPoint
  * nc
  * lineage
  * strain
  * island_index
- [ ] islander_supplemental
- [x] literature_islands
  * litterature_name
  * citation
  * accession
  * left 
  * right
  * island
- [ ] phast_overlap
- [ ] phylogeny
- [ ] potential_fragments
- [ ] potential_integrases
- [ ] project_acc
- [ ] smpb
- [ ] smpb_nicknames
- [ ] tandem
- [ ] tmrna
- [ ] tmrna_lineage
- [ ] tmrna_nicknames
- [ ] trna
- [ ] trna_aa
- [ ] unique_tmrna
- [ ] usage_monitor


# Fetch references**
* Google Scholar [https://github.com/ckreibich/scholar.py ]
```bash
 "accession_nb" --csv
 (en moyenne 10 par # d'accession)

 "accession_nb island" --csv

 "reference" --csv
```

>Ne fonctionne pas, aucun résultat ou trop vague.

* eDirect
[Documentations](https://www.ncbi.nlm.nih.gov/books/NBK179288/)
```bash
#Installation de EDirect 
cd ~
/bin/bash
perl -MNet::FTP -e \
  '$ftp = new Net::FTP("ftp.ncbi.nlm.nih.gov", Passive => 1);
   $ftp->login; $ftp->binary;
   $ftp->get("/entrez/entrezdirect/edirect.tar.gz");'
gunzip -c edirect.tar.gz | tar xf -
rm edirect.tar.gz
builtin exit
export PATH=${PATH}:$HOME/edirect >& /dev/null || setenv PATH "${PATH}:$HOME/edirect"
./edirect/setup.sh
echo "export PATH=\${PATH}:/home/local/USHERBROOKE/whif3801/edirect" >> $HOME/.bashrc
```

Creation d'un script bash pour executer des recherches avec **EDirect** avec des querys sélectionnés. 
```bash
# Exemple de recherche avec eDirect
esearch -db pubmed -query "Reen et al., 2006 Nat Rev Microbiol" -sort Relevance |efetch -format medline| grep -wf pattern

# script (pas à jour)
python3 references.py -p islandpick.pickle -q 'organism' 'start' 'end' -o 'output/references.organims_position_island' -key 'island'

```
> Beaucoup de résultats, peut⁻être pas assez précis pour identifier le bon article


# Fetch sequence

Creation d'un script bash pour executer des recherches avec **EDirect** avec des querys sélectionnés. 
```bash
# Exemple de recherche avec eDirect
esearch -db Nucleotide -query "NC_000919.1" |efetch -format fasta > genome/sequence.NC_000919.1.fasta

# script
python3 sequences.py -i list_of_accession_number 
sh get_genomes.sh

# Validation des génomes téléchargés
grep ">" genomes/* >validation
sed 's/.fasta:>/\t/g' validation | sed 's%genomes/sequence.%%g'| sed 's/ /\t/g'|cut -f1,2| awk '$1!=$2' |wc -l

```


*Plusieurs séquences avec mauvais # d'accession dans les fichiers fasta de certain génome*
> Imprimer les lignes contenants ">" et leur # de ligne   
> Identifier le bon génome manuellement    
> Sélectionner les lignes désirées   
```bash
grep -n genomes/sequence.NC_006840.2.fasta
sed -n '311259,352653p' genomes/sequence.NC_006840.2.fasta > genomes_fix/sequence.NC_006840.2.fasta
mv genomes_fix/* genomes/.
```

*Dans les fichiers fasta, le génome est coupé à chaque changement de ligne*
> 2 options trouvés
```python
## Option 1 : sed & cut
# en ligne de commande
grep -v ">" genome.fasta |sed ':a;N;$!ba;s/\n//g' |cut -c 271673-289384 >genome_island.fasta

# ou dans un script python
head = os.system(grep ">" genome_file).read()
sequence = os.system(grep -v ">" genome_file |sed ':a;N;$!ba;s/\n//g' |cut -c 271673-289384).read()

## Option 2 : SeqIO
from Bio import SeqIO
genome = SeqIO.read(genome_file, "fasta")
sequence = genome.seq[ (271673-basepairs):(289384+basepairs) ]
```
> SeqIO beaucoup plus rapide ! (6 min vs 17 min pour sed & cut ) 

# Index des numéros d'accession

```bash
# Obtenir les noms des organismes liés aux numéros d'accession dans les fasta.
grep ">" genomes/* > header
sed 's/.fasta:>/\t/g' header |cut -f2|awk '{print $1"\t"$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' > index.acc_num.txt

```


# à faire et problèmes restants
* déplacer les lignes manquant de l'information dans autre fichier
* modifier les numéros d'accession avec ou sans .1, .2 ... (solution temporaire dans le script parsing.py ligne 27-30)
* ajouter les génomes manquants (python3 sequences.py -i genomes_to_add.txt && sh get_sequences.sh )
* effectuer des blast pour comparer les séquences des bases de données avec celle ajouté des fasta (sequences_validation.txt)
* ajouter les références
