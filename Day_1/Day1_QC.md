
# Downloading the TRYPANO metagenomes

We are going to use the data from the following article: 

Eberhard, F.E., Klimpel, S., Guarneri, A.A. et al. Exposure to Trypanosoma parasites induces changes in the microbiome of the Chagas disease vector Rhodnius prolixus. Microbiome 10, 45 (2022). https://doi.org/10.1186/s40168-022-01240-z

***Abstract***
***Background***
***The causative agent of Chagas disease, Trypanosoma cruzi, and its nonpathogenic relative, Trypanosoma rangeli, are transmitted by haematophagous triatomines and undergo a crucial ontogenetic phase in the insect’s intestine. In the process, the parasites interfere with the host immune system as well as the microbiome present in the digestive tract potentially establishing an environment advantageous for development. However, the coherent interactions between host, pathogen and microbiota have not yet been elucidated in detail. We applied a metagenome shotgun sequencing approach to study the alterations in the microbiota of Rhodnius prolixus, a major vector of Chagas disease, after exposure to T. cruzi and T. rangeli focusing also on the functional capacities present in the intestinal microbiome of the insect.***
***Results***
***The intestinal microbiota of R. prolixus was dominated by the bacterial orders Enterobacterales, Corynebacteriales, Lactobacillales, Clostridiales and Chlamydiales, whereas the latter conceivably originated from the blood used for pathogen exposure. The anterior and posterior midgut samples of the exposed insects showed a reduced overall number of organisms compared to the control group. However, we also found enriched bacterial groups after exposure to T. cruzi as well as T rangeli. While the relative abundance of Enterobacterales and Corynebacteriales decreased considerably, the Lactobacillales, mainly composed of the genus Enterococcus, developed as the most abundant taxonomic group. This applies in particular to vectors challenged with T. rangeli and at early timepoints after exposure to vectors challenged with T. cruzi. Furthermore, we were able to reconstruct four metagenome-assembled genomes from the intestinal samples and elucidate their unique metabolic functionalities within the triatomine microbiome, including the genome of a recently described insect symbiont, Candidatus Symbiopectobacterium, and the secondary metabolites producing bacteria Kocuria spp.***
***Conclusions***
***Our results facilitate a deeper understanding of the processes that take place in the intestinal tract of triatomine vectors during colonisation by trypanosomal parasites and highlight the influential aspects of pathogen-microbiota interactions. In particular, the mostly unexplored metabolic capacities of the insect vector’s microbiome are clearer, underlining its role in the transmission of Chagas disease.***


![communities](http://url/to/img.png)
![pangenomes](http://url/to/img.png)


The following link contains FTP URLs for each raw data file for 30 samples:

> https://www.ncbi.nlm.nih.gov/Traces/study/?WebEnv=MCID_661792592760e53abf6a1827&query_key=2

To download each of the raw sequencing data file from the NCBI servers


```bash
for i in {482..511}; do prefetch SRX11581$i; done
```

And the NCBI SRA Toolkit to convert the SRA files to FASTQ files:

>https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit


```bash
NUM_CORES=64

# Define the function to be executed in each folder
execute_commands() {
    cd "$1" || return 1  # Change to the folder or return if unsuccessful
    for sra_file in *.sra; do
        fastq-dump --split-3 "$sra_file"
        pigz *.fastq &
    done
    wait  # Wait for all background jobs to finish
    cd ..
}

# Export the function for use by parallel
export -f execute_commands

# List directories in RAWDATA and execute the function in parallel
ls -d RAWDATA/*/ | parallel -j $NUM_CORES execute_commands {}
```


# Defining metagenomic sets, setting sample names, and linking those with the raw data

The dataset contains 5 timepoints T0, T1, T2,T3, T7 (numbers are weeks) and for each sample the anterior midgut (AM) and posterior midgut (PM) were sequenced individualy. 
The Triatomines were infected with 2 different Trypanosoma species, **T. cruzi** (cruz) and **T. rangeli** (rang) and a control group (Cont).

The samples were named as follows:
T0ContAM
T0ContPM
T0cruzAM
T0cruzPM
T0rangAM
T0rangPM
...

```bash
wget http://github.com/vincentmanz/integrative_metgenomics/data/metadatas.txt
```

The contents of the Trypanosoma set file should look like this:

```bash
$ head metadatas.txt
Sample	SRA identifier	Type	Time	Number-of-days	Midgut
T0ContAM	SRR15276521	Control	T0	0	AM
T0ContPM	SRR15276520	Control	T0	0	PM
T0cruzAM	SRR15276539	T._cruzi	T0	0	AM
T0cruzPM	SRR15276528	T._cruzi	T0	0	PM
T0rangAM	SRR15276527	T._rangeli	T0	0	AM

```

# Quality control and trimming of the raw reads

### Running fastQC


```bash
USER=your_username
mamba activate QC_env
```
> QC_env is a conda environment with fastqc, multiqc and cutadapt installed.

Run `fastQC` to the files stored in the RAWDATA folder. What does the `-o` and `-t` flags refer to?

Running the QC step on all sequence files would take too long, so they are already done and you can just copy them from the course files ''
Make sure you're on your own folder before copying.


```bash
mkdir RAWDATA/FASTQC_RAW

# for i in {518..547}; do
#     fastqc RAWDATA/SRR15276"$i"/SRR15276"$i"_*.fastq.gz -o RAWDATA/FASTQC_RAW/ -t 4
# done

NUM_CORES=15

# Define the function to be executed in parallel
run_fastqc() {
    fastqc RAWDATA/SRR15276"$1"/SRR15276"$1"_*.fastq.gz -o RAWDATA/FASTQC_RAW/ -t 5
}

# Export the function for use by parallel
export -f run_fastqc

# Run fastqc commands in parallel for the specified range of SRR numbers
seq 518 547 | parallel -j $NUM_CORES run_fastqc
```
Running the QC step on all sequence files would take too long, so they are already done and you can just copy them from the course files ''
Make sure you're on your own folder before copying.

Then combine the reports in FASTQC folder with multiQC:
MultiQC will run in the same virtual environment QC_env. 

```bash
multiqc RAWDATA/FASTQC_RAW/* -o RAWDATA/MULTIQC_RAW --interactive
```

To leave the interactive node, type `exit`.  

Open the resulting HTML file.   
Have a look at the QC report with your favourite browser.  

After inspecting the output, it should be clear that we need to do some trimming.  
__What kind of trimming do you think should be done?__



### Running Cutadapt
For trimming we have an array script that runs `Cutadapt` for each file in the `RAWDATA` folder.  
Go to your working directory and copy the `CUTADAPT.sh` script from `/scratch/project_2001499/COURSE_FILES/SBATCH_SCRIPTS`.  
Check the script for example with the command `less`.  
The adapter sequences that you want to trim are located after `-a` and `-A`.  
What is the difference with `-a` and `-A`?  
And what is specified with option `-p` or `-o`?
And how about `-m` and `-j`?  
You can find the answers from Cutadapt [manual](http://cutadapt.readthedocs.io).

Before running the script, we need to create the directory where the trimmed data will be written:

```bash
mkdir TRIMMEDDATA
```

Then we need to submit our jos to the SLURM system.  
Make sure to submit it from your own folder.  
More about CSC batch jobs here: https://docs.csc.fi/computing/running/creating-job-scripts-puhti/.  

```bash
bash HELPER/CUTADAPT.sh
```

### Running fastQC on the trimmed reads
Go to the folder containing the trimmed reads (`TRIMMED`) and view the `Cutadapt` log. Can you answer:

* How many read pairs we had originally?
* How many reads contained adapters?
* How many read pairs were removed because they were too short?
* How many base calls were quality-trimmed?
* Overall, what is the percentage of base pairs that were kept?

Then make a new folder (`FASTQC`) for the QC files of the trimmed data and run fastQC and multiQC again as you did before trimming:
Again the QC part would take too long, so we have created the files for you to copy and run only the multiQC part.

```bash
mkdir TRIMMEDDATA/FASTQC_TRIMMED
#fastqc TRIMMEDDATA/*.fastq.gz -o TRIMMEDDATA/FASTQC_TRIMMED/ -t 64


NUM_CORES=15

# Define the function to be executed in parallel
run_fastqc() {
    fastqc TRIMMEDDATA/SRR15276"$1".*.fastq.gz -o TRIMMEDDATA/FASTQC_TRIMMED/ -t 5
}

# Export the function for use by parallel
export -f run_fastqc

# Run fastqc commands in parallel for the specified range of SRR numbers
seq 518 547 | parallel -j $NUM_CORES run_fastqc

```

```bash
multiqc TRIMMEDDATA/FASTQC_TRIMMED/* -o TRIMMEDDATA/MULTIQC_TRIMMED --interactive --title "Trimmed data"
# To comapre the raw and trimmed data
multiqc RAWDATA/FASTQC_RAW/* TRIMMEDDATA/FASTQC_TRIMMED/* -o TRIMMEDDATA/MULTIQC_TRIMMED --interactive --title "Raw vs Trimmed data"
# deactivate the environment
mamba deactivate
```

Look how well the trimming went.


# Removing host reads

### Downloading the host genome

We will use the genome of the kissing bug, Triatoma infestans GCA_011037195.1 [Genome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_011037195.1/).

To downoload the **datasets** package, you need to have an account in the [NCBI website](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/).

```bash
mkdir HOST_MAPPING

./HOST_MAPPING/datasets download genome accession GCA_011037195.1 --include gff3,rna,cds,protein,genome,seq-report --filename Triatoma_infestans.fasta.zip
unzip HOST_MAPPING/Triatoma_infestans.fasta.zip
```

### Mapping the reads to the host genome
note: switch to [Sambamba](https://lomereiter.github.io/sambamba/) 

```bash
NUM_CORES=60

bowtie2-build --threads $NUM_CORES HOST_MAPPING/ncbi_dataset/data/GCA_011037195.1/GCA_011037195.1_UVM_Tinf_1.0_genomic.fna HOST_MAPPING/Triatoma_infestans

for i in {518..547}
do
    echo "Mapping: SRR15276"$i""
    # 1. map the reads to the host genome
    bowtie2 --threads $NUM_CORES \
            -x HOST_MAPPING/Triatoma_infestans \
            -1 TRIMMEDDATA/SRR15276"$i".R1.fastq.gz \
            -2 TRIMMEDDATA/SRR15276"$i".R2.fastq.gz \
            -S HOST_MAPPING/SRR15276"$i".sam
    # 2. covert the resulting SAM file to a BAM file:
    samtools view -@ $NUM_CORES -f 4 -bS  HOST_MAPPING/SRR15276"$i".sam > HOST_MAPPING/SRR15276"$i"_RAW.bam

    # 3. sort the BAM file:
    samtools sort -@ $NUM_CORES HOST_MAPPING/SRR15276"$i"_RAW.bam -o HOST_MAPPING/SRR15276"$i"_sorted.bam

    # 4. index the BAM file:
    samtools index -@ $NUM_CORES HOST_MAPPING/SRR15276"$i"_sorted.bam

    # 5. Convert unmapped reads to FastQ
    bamToFastq -i HOST_MAPPING/SRR15276"$i"_sorted.bam -fq HOST_MAPPING/SRR15276"$i"_unmapped_R1.fastq -fq2 HOST_MAPPING/SRR15276"$i"_unmapped_R2.fastq

    # Optional: Remove intermediate files to save disk space
    #rm HOST_MAPPING/SRR15276"$i".sam HOST_MAPPING/SRR15276"$i"_RAW.bam HOST_MAPPING/SRR15276"$i"_sorted.bam
done
```

# Assembling metagenomic data

We used the program MEGAHIT v1.0.3 (Li et al., 2014) to perform 30 metagenomic each experiment. 

We used MEGAHIT the following way to individually co-assemble each metagenomic set:

```bash
mamba activate MEGAHIT_env
mkdir ASSEMBLIES

for i in {518..547}
do
    megahit -1 TRIMMEDDATA/SRR15276"$i".R1.fastq.gz \
            -2 TRIMMEDDATA/SRR15276"$i".R2.fastq.gz \
            --min-contig-len 1000 \
            --k-min 27 \
            --k-max 127 \
            --k-step 10 \
            --memory 0.9 \
            -o ASSEMBLIES/SRR15276"$i"-RAW.fa \
            --num-cpu-threads 64
done
```


### Assembly quality statistics

Let's take a look at the assemblies in a bit more detail with [MetaQUAST](http://bioinf.spbau.ru/metaquast).

```bash
metaquast.py -o ASSEMBLIES/QUAST_REPORTS -t 64 ASSEMBLIES/SRR15276*/*.fa  
python3 ~/Downloads/quast-5.2.0/metaquast.py -f -b -t 60 ASSEMBLIES/SRR15276*.fa/final.contigs.fa   -o ASSEMBLIES/QUAST_REPORTS 
mamba deactivate
```

Copy the folder called METAQUAST_FAST to your computer.
You can view the results (report.html) in your favorite browser.

Questions about the assembly QC:

Which assembly has the longest contig when also long reads assemblies are included?
Which assembly had the most contigs?
Were the long read assemblies different from the corresponding short read assemblies?
If yes, in what way?




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



# Mapping the reads back to the assembly

Next thing to do is mapping all the reads back to the assembly. We use the renamed >5000 nt contigs and do it sample-wise, so each sample is mapped separately using the trimmed R1 & R2 reads.


```bash

# Make a directory to house the mapping results
mkdir -p MAPPING

# use a for loop to map the recipient gut metagenomes from PRE and POST FMT metagenomes
# against our MAG reference

for i in {518..547}
do
    # 1. perform read recruitment with bowtie2 to get a SAM file:
    echo -e "Mapping: SRR15276"$i""

    bowtie2-build ANVIO/01_ASSEMBLIES_5000nt/SRR15276"$i"_5000nt.fa MAPPING/SRR15276"$i"_5000nt

    bowtie2 --threads $NUM_CORES \
            -x MAPPING/SRR15276"$i"_5000nt \
            -1 TRIMMEDDATA/SRR15276"$i".R1.fastq.gz \
            -2 TRIMMEDDATA/SRR15276"$i".R2.fastq.gz \
            -S MAPPING/SRR15276"$i".sam

    # 2. covert the resulting SAM file to a BAM file:
    samtools view -F 4 -bS  MAPPING/SRR15276"$i".sam > MAPPING/SRR15276"$i"_RAW.bam

    # 3. sort the BAM file:
    samtools sort MAPPING/SRR15276"$i"_RAW.bam -o MAPPING/SRR15276"$i"_sorted.bam

    # 4. index the BAM file:
    samtools index MAPPING/SRR15276"$i"_sorted.bam

done
```



### Profile the mapping results

Previously, in "Rename the scaffolds and select those >5000nt", we have selected contigs > 5000bp only, to be consistant, we will profile the mapping results with the same minimum size of contigs. 
We will build profiles for each sample that was mapped against the assembly. The mapping output from each sample is the $SAMPLE.bam file.

```bash
mkdir PROFILE
for i in {518..547}
do    anvi-profile -c ANVIO/02_CONTIGS_DB/SRR15276"$i"_5000nt_CONTIGS.db \
                 -i MAPPING/SRR15276"$i"_sorted.bam \
                 --num-threads $NUM_CORES \
                 --min-contig-length 5000 \
                 -o PROFILE/SRR15276"$i".db
done
```

### Merging the profiles

When the profiling is done, you can merge them with this command. Remember to re-attach to you screen and run the command in there.

```bash

mkdir ANVIO/03_MERGED_PROFILES
for i in {518..547}
do    
    anvi-merge PROFILE/SRR15276"$i".db/PROFILE.db \
                --output-dir ANVIO/03_MERGED_PROFILES \
                --contigs-db ANVIO/02_CONTIGS_DB/SRR15276"$i"_5000nt_CONTIGS.db \
                --force-overwrite \
                --enforce-hierarchical-clustering  &> ANVIO/03_MERGED_PROFILES/SRR15276"$i"_merge.log
done

anvi-merge PROFILE/SRR15276*.db/PROFILE.db \
                --output-dir ANVIO/03_MERGED_PROFILES \
                --contigs-db ANVIO/02_CONTIGS_DB/SRR15276*_5000nt_CONTIGS.db \

for i in {518..547}
do   
anvi-estimate-scg-taxonomy -p  PROFILE/SRR15276"$i".db/PROFILE.db \
                --contigs-db ANVIO/02_CONTIGS_DB/SRR15276"$i"_5000nt_CONTIGS.db \
                --metagenome-mode \
                --compute-scg-coverages
done

``` 


# Visualizing the metagenomic data

### Interactive interface introduction

To start the interactive interface, you need to run the following command:

> anvi-interactive -c ANVIO/02_CONTIGS_DB/SRR15276518_5000nt_CONTIGS.db  -p PROFILE/SRR15276518.db/PROFILE.db 


anvi-estimate-scg-taxonomy --metagenomes METADATA/metagenomes.tsv \
      --output-file-prefix METAGENOME-EXAMPLE






# read based analysis with host reference genome

### Diamond blastx


The data base is the nr db from NCBI with the Trypanosoma and Triatome genomes added.

GCA_002087225.1	Trypanosoma theileri
GCA_000691245.1	Trypanosoma grayi
GCA_003719475.1	Trypanosoma rangeli
GCA_003719485.1	Trypanosoma conorhini
GCA_001457755.2	Trypanosoma equiperdum
GCA_003543875.1	Trypanosoma brucei equiperdum
GCA_015033655.1	Trypanosoma cruzi
GCA_015033625.1	Trypanosoma cruzi
GCA_013358655.1	Trypanosoma cruzi
GCA_011037195.1	Triatoma infestans

https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_011037195.1/

```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip


cat  triatome.fasta Trypanosoma.fasta >> nr
mv nr  tt.fasta
diamond makedb \
    --threads 60 \
    --taxonnodes taxdmp/nodes.dmp \
    --taxonnames taxdmp/names.dmp \
    --in tt.fasta \
    --db tt_nr_diamond

tar --use-compress-program="pigz  --best --recursive " -cf Triatome_genomes.tar.gz Triatome\ genomes/
```



```bash
conda activate read_based_env
mkdir READBASED
mkdir READBASED/RESAMPLED
mkdir READBASED/MEGAN
mkdir READBASED/METAXA

for i in {518..547}; do
# resample 2 million reads
  seqtk sample -s100 TRIMMEDDATA/SRR15276"$i".R1.fastq.gz 2000000 > READBASED/RESAMPLED/SRR15276"$i".R1.fastq
  seqtk sample -s100 TRIMMEDDATA/SRR15276"$i".R2.fastq.gz 2000000 >> READBASED/RESAMPLED/SRR15276"$i".R1.fastq
# run diamond blastx
  diamond blastx --query READBASED/RESAMPLED/SRR15276"$i".R1.fastq \
                 --out MEGAN/$SAMPLE.blastx.txt \
                 --db ~/Share/DBs/nr \
                 --outfmt 0 \
                 --threads 4 > MEGAN/$SAMPLE.diamond.log.txt
# convert blastx to rma6
  ~/Share/megan/tools/blast2rma --in MEGAN/$SAMPLE.blastx.txt \
                                --out MEGAN/$SAMPLE.rma6 \
                                --mapDB ~/Share/DBs/megan-map-Jan2021.db \
                                --format BlastText \
                                --threads 4 > MEGAN/$SAMPLE.megan.log.txt
# import rma6 to metaxa2
  metaxa2 -1 RESAMPLED/$SAMPLE.R1.fastq \
          -2 RESAMPLED/$SAMPLE.R2.fastq \
          -o METAXA/$SAMPLE \
          --align none \
          --graphical F \
          --cpu 4 \
          --plus > METAXA/$SAMPLE.metaxa.log.txt

  metaxa2_ttt -i METAXA/$SAMPLE.taxonomy.txt \
              -o METAXA/$SAMPLE >> METAXA/$SAMPLE.metaxa.log.txt
done

metaxa2_dc -o METAXA/metaxa_genus.txt METAXA/*level_6.txt
metaxa2_dc -o METAXA/metaxa_species.txt METAXA/*level_7.txt
```



# With nf-core MAGs

we are using a pipeline called [MAGs](https://nf-co.re/mag/2.5.4) to assemble and bin the metagenomic data.

[Installation](https://nf-co.re/docs/usage/installation)

data/Trypanosoma_exposure/HOST_MAPPING/ncbi_dataset/data/GCA_011037195.1/GCA_011037195.1_UVM_Tinf_1.0_genomic.fna \
      --host_removal_verysensitive True \

```bash
nf-core download mags
```

```bash
nextflow run ~/nf-core/nf-core-mag_2.5.4/2_5_4/ \
      -resume \
      -profile singularity \
      --input '/home/vincent/Documents/project/Metagenomics_2024/data/Trypanosoma_exposure/TRIMMEDDATA/SRR152765*.R{1,2}.fastq.gz' \
      --outdir /home/vincent/Documents/project/Metagenomics_2024/data/Trypanosoma_exposure/MAGs \
      --multiqc_title "MAGs Trypanosoma exposure" \
      # reproducibility \
      --megahit_fix_cpu_1 FALSE \ 
      #  Quality control for short reads options \
      --reads_minlength 130 \
      --host_genome GCA_011037195_Triatoma_infestans \
      --host_fasta /home/vincent/Documents/project/Metagenomics_2024/ \
      # Taxonomic profiling options \
      --bbnorm /media/vincent/Data/DB/centrifuge/p_compressed_2018_4_15.tar.gz \
      --kraken2_db /media/vincent/Data/DB/kraken/k2_viral_20240112 \
      --cat_db_generate TRUE \
      --save_cat_db TRUE \
      --skip_gtdbtk TRUE \ 
      # Assembly options \
      --coassemble_group FALSE \ 
      --skip_spades TRUE \
      --skip_spadeshybrid TRUE \
      --skip_megahit FALSE  \
      --skip_quast FALSE \
      # gene prediction and annotation options \
      --skip_prodigal FALSE \
      --skip_prokka FALSE \
      --skip_metaeus TRUE \
      #  Virus identification options \
      --run_virus_identification TRUE \
      # Binning options \
      --bowtie2_mode="--very-sensitive" \ 
      # Bin quality check options
      --skip_binqc FLASE \
      --busco_auto_lineage_prok TRUE \
      --save_busco_db TRUE \
      --busco_clean TRUE \
      --save_checkm_data TRUE \
      --refine_bins_dastool TRUE  \
      --run_gunc FALSE \
      --ancient_dna FALSE








