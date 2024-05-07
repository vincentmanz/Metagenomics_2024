
# Genome-resolved metagenomics with anvi'o

`anvi'o` is an analysis and visualization platform for omics data.  
You can read more from their [webpage](http://merenlab.org/software/anvio/).

### Rename the scaffolds and select those >5000nt.

`anvi'o` wants sequence IDs in your FASTA file as simple as possible.  
Therefore we need to reformat the headers to remove spaces and non-numeric characters.  
Also contigs shorter than 5000 bp will be removed.

```bash
mamba activate anvio-8

mkdir ANVIO
mkdir ANVIO/01_ASSEMBLIES_5000nt

for i in {518..547}
do
    anvi-script-reformat-fasta ASSEMBLIES/SRR15276"$i"-RAW.fa/final.contigs.fa \
                              --min-len 5000 \
                              --simplify-names \
                              --prefix SRR15276"$i" \
                              -r ANVIO/01_ASSEMBLIES_5000nt/REPORT \
                              -o ANVIO/01_ASSEMBLIES_5000nt/SRR15276"$i"_5000nt.fa
done
```

### Generate CONTIGS.db
The contigs database (`Sample03_5000nt_CONTIGS.db`) contains information on contig length, open reading frames (searched with `Prodigal`) and kmer composition.  
See the [anvi'o webpage](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#creating-an-anvio-contigs-database) for more information.  

```bash
mkdir ANVIO/02_CONTIGS_DB
for i in {518..547}
do
    anvi-gen-contigs-database --contigs-fasta  ANVIO/01_ASSEMBLIES_5000nt/SRR15276"$i"_5000nt.fa \
                          --output-db-path ANVIO/02_CONTIGS_DB/SRR15276"$i"_5000nt_CONTIGS.db \
                          -n SRR15276"$i" \
                          --num-threads 60
done

```

### Run HMMs to identify single-copy core genes for Bacteria, Archaea and Eukarya, plus rRNAs



First annotate the SCGs.
```bash

for i in {518..547}
do
    anvi-run-hmms -c ANVIO/02_CONTIGS_DB/SRR15276"$i"_5000nt_CONTIGS.db \
                  --num-threads 60
done
```


### Run the SCG taxonomy 

#### setup

Please first run diamond --version to make sure you have the right version. If not, install it first. If you are in a conda environment, you can simply run conda install diamond=0.9.14.
To setup your SCG taxonomy databases, run the following program, which will not take longer than a minute:


```bash
anvi-setup-scg-taxonomy
```



#### Populating contigs db with SCG taxonomy

This is something you will do once for every contigs database you wish to work with.




```bash
for i in {518..547}
do
    anvi-run-scg-taxonomy -c ANVIO/02_CONTIGS_DB/SRR15276"$i"_5000nt_CONTIGS.db --num-threads 20 --num-parallel-processes 3
done
```


#### Estimating taxonomy in the terminal

We are working with a contigs database that represents a metagenome, anvi’o will complaining that there is too much redundancy of SCGs for this to be a single genome. You can avoid that by using the flag --metagenome-mode. 

```bash
for i in {518..547}
do
    anvi-estimate-scg-taxonomy -c ANVIO/02_CONTIGS_DB/SRR15276"$i"_5000nt_CONTIGS.db \
                           --metagenome-mode
                           --num-threads 60 
done
```

First let’s take a look at the output, and then discuss what is happening behind the scenes:

<!---
Contigs DB ...................................: ANVIO/02_CONTIGS_DB/SRR15276544_5000nt_CONTIGS.db                                                                                                                                                                 
Metagenome mode ..............................: True
SCG for metagenome ...........................: None
                                                                                                                                                                                                                                                                  
* A total of 48 single-copy core genes with taxonomic affiliations were
  successfully initialized from the contigs database 🎉 Following shows the
  frequency of these SCGs: Ribosomal_S3_C (3), Ribosomal_S8 (3), Ribosomal_S20p
  (3), Ribosomal_L6 (3), Ribosomal_L16 (3), Ribosomal_L22 (3), ribosomal_L24
  (3), Ribosomal_L27A (3), Ribosomal_S2 (2), Ribosomal_S6 (2), Ribosomal_S7 (2),
  Ribosomal_S11 (2), Ribosomal_L2 (2), Ribosomal_L3 (2), Ribosomal_L4 (2),
  Ribosomal_L9_C (2), Ribosomal_L20 (2), Ribosomal_L21p (2), Ribosomal_S9 (1),
  Ribosomal_L1 (1), Ribosomal_L13 (1), Ribosomal_L17 (1).

WARNING
===============================================
Anvi'o automatically set 'Ribosomal_S3_C' to be THE single-copy core gene to
survey your metagenome for its taxonomic composition. If you are not happy with
that, you could change it with the parameter `--scg-name-for-metagenome-mode`.

                                                                                                                                                                                                                                                                  
Taxa in metagenome "SRR15276544"
===============================================
+---------------------------------+--------------------+--------------------------------------------------------------------------------------------------------------------+
|                                 |   percent_identity | taxonomy                                                                                                           |
+=================================+====================+====================================================================================================================+
| SRR15276544_Ribosomal_S3_C_2202 |               98.6 | Bacteria / Actinomycetota / Actinomycetia / Mycobacteriales / Mycobacteriaceae / Rhodococcus / Rhodococcus rhodnii |
+---------------------------------+--------------------+--------------------------------------------------------------------------------------------------------------------+
| SRR15276544_Ribosomal_S3_C_3182 |                 99 | Bacteria / Bacillota / Bacilli / Lactobacillales / Enterococcaceae / Enterococcus / Enterococcus faecalis          |
+---------------------------------+--------------------+--------------------------------------------------------------------------------------------------------------------+

-->

One of the first messages you see in that output is this, where anvi’o reports the frequencies of SCGs:

<!---
* A total of 48 single-copy core genes with taxonomic affiliations were
  successfully initialized from the contigs database 🎉 Following shows the
  frequency of these SCGs: Ribosomal_S3_C (3), Ribosomal_S8 (3), Ribosomal_S20p
  (3), Ribosomal_L6 (3), Ribosomal_L16 (3), Ribosomal_L22 (3), ribosomal_L24
  (3), Ribosomal_L27A (3), Ribosomal_S2 (2), Ribosomal_S6 (2), Ribosomal_S7 (2),
  Ribosomal_S11 (2), Ribosomal_L2 (2), Ribosomal_L3 (2), Ribosomal_L4 (2),
  Ribosomal_L9_C (2), Ribosomal_L20 (2), Ribosomal_L21p (2), Ribosomal_S9 (1),
  Ribosomal_L1 (1), Ribosomal_L13 (1), Ribosomal_L17 (1).
-->

Because there are 3 Ribosomal_S3 genes, anvi’o automatically used this SCG to estimate taxonomy. 

For instance, Ribosomal_S3 is a gene that encodes for a ribosomal protein, and it is a single-copy core gene in anvi’o’s database. This gene is present in all genomes, and it is a good gene to use to estimate the taxonomy of a genome.

Ribosomal_S3 occurs 9 times instead of 10, we can ask anvi’o to use that one instead the following way:

```bash
anvi-estimate-scg-taxonomy -c ANVIO/02_CONTIGS_DB/SRR15276544_5000nt_CONTIGS.db \
                           --metagenome-mode \
                           --scg-name Ribosomal_S3
```

<!---
Contigs DB ...................................: ANVIO/02_CONTIGS_DB/SRR15276544_5000nt_CONTIGS.db                                                                                                                                                                 
Metagenome mode ..............................: True
SCG for metagenome ...........................: None
                                                                                                                                                                                                                                                                  
* A total of 48 single-copy core genes with taxonomic affiliations were
  successfully initialized from the contigs database 🎉 Following shows the
  frequency of these SCGs: Ribosomal_S3_C (3), Ribosomal_S8 (3), Ribosomal_S20p
  (3), Ribosomal_L6 (3), Ribosomal_L16 (3), Ribosomal_L22 (3), ribosomal_L24
  (3), Ribosomal_L27A (3), Ribosomal_S2 (2), Ribosomal_S6 (2), Ribosomal_S7 (2),
  Ribosomal_S11 (2), Ribosomal_L2 (2), Ribosomal_L3 (2), Ribosomal_L4 (2),
  Ribosomal_L9_C (2), Ribosomal_L20 (2), Ribosomal_L21p (2), Ribosomal_S9 (1),
  Ribosomal_L1 (1), Ribosomal_L13 (1), Ribosomal_L17 (1).

WARNING
===============================================
Anvi'o automatically set 'Ribosomal_S3_C' to be THE single-copy core gene to
survey your metagenome for its taxonomic composition. If you are not happy with
that, you could change it with the parameter `--scg-name-for-metagenome-mode`.

                                                                                                                                                                                                                                                                  
Taxa in metagenome "SRR15276544"
===============================================
+---------------------------------+--------------------+--------------------------------------------------------------------------------------------------------------------+
|                                 |   percent_identity | taxonomy                                                                                                           |
+=================================+====================+====================================================================================================================+
| SRR15276544_Ribosomal_S3_C_2202 |               98.6 | Bacteria / Actinomycetota / Actinomycetia / Mycobacteriales / Mycobacteriaceae / Rhodococcus / Rhodococcus rhodnii |
+---------------------------------+--------------------+--------------------------------------------------------------------------------------------------------------------+
| SRR15276544_Ribosomal_S3_C_3182 |                 99 | Bacteria / Bacillota / Bacilli / Lactobacillales / Enterococcaceae / Enterococcus / Enterococcus faecalis          |
+---------------------------------+--------------------+--------------------------------------------------------------------------------------------------------------------+
| SRR15276544_Ribosomal_S3_C_7443 |                 99 | Bacteria / Pseudomonadota / Gammaproteobacteria / Enterobacterales / Enterobacteriaceae /  /                       |
+---------------------------------+--------------------+--------------------------------------------------------------------------------------------------------------------+
-->

These may not be conclusive taxonomic insights, but it’s still very useful to have a rapid idea about what you have. After all what is a conclusive taxonomic insight anyway?

