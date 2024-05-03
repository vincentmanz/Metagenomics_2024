# Read based analyses

The first question to answer in a metagenomics study is who is there. It is to identify the members of the microbial community. There are two approaches that could be used to achieve this goal: the first is to utilize the reads and the other is to assemble the metagenomes before using homology search to the database. We will look at the read based approach in this section. There are three main algorithms to classify the reads to taxa: the first is to do homology search (for example, using blast: MEGAN) of the reads against huge reference databases, the second is k-mer based and the third is marker gene based. In this step, we are going to focus on the k-mer based classification usually more accurate and faster to compute. 

<p align="center">
      <img src="https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_1/pictures/benchmark_calssifier.jpg"  width="400">
      <img src="https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_1/pictures/benchmark_calssifier_abundance.jpg"  width="400">
</p>


[Reference](https://doi.org/10.1016/j.cell.2019.07.010)

![benchmark_calssifier_abundance](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_1/pictures/benchmark.png)

[Reference](https://doi.org/10.1186/s40793-024-00561-w)


### Taxonomy profiling using Kraken2 & Bracken abundance estimation

K-mer based methods count the k-mer frequency of the reads, and compare it to a model trained with sequences from known genomes. There are a few tools belong to this group of method: Kraken, Centrifuge, Kraken2, Clark, Kaiju. For the Read based analyses we will use [Kraken2](https://doi.org/10.1186/s13059-019-1891-0). 

#### Download the database for Kraken2 and Bracken
In order to run Kraken2, one has to build corresponding database first, the command to build the standard Kraken2 database is kraken2-build –standard –threads 24 –db kraken.db. This will download NCBI taxonomic information, as well as the complete genomes in RefSeq for the bacterial, archaeal, and viral domains, along with the human genome and a collection of known vectors (UniVec_Core). The build process is the most time-consuming, so we are not going to perform it in this workshop. We will link to prebuild databases are [here](https://benlangmead.github.io/aws-indexes/k2).

k2_standard_08gb = Standard with DB capped at 8 GB: Refeq archaea, bacteria, viral, plasmid, human1, UniVec_Core

```bash
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240112.tar.gz

```


#### Kraken2 classification 

```bash
mkdir READBASED

for i in {518..547}
do
    echo "Mapping: SRR15276"$i""
    kraken2 --db READBASED/k2_standard_20240112/ --threads 60  --output READBASED/SRR15276"$i".kraken.out --report READBASED/SRR15276"$i".kraken_report.out --paired TRIMMEDDATA/SRR15276"$i".R1.fastq.gz TRIMMEDDATA/SRR15276"$i".R2.fastq.gz
done
```


### Bracken abundance estimation

What Kraken2 has produced is the classification of each read to a taxonomic rank. This result needs to be further refined to generate species (genus, phylum) level abundance table for downstream statistical analysis. This process needs to be done properly. Kraken2 only classifies the reads to the Lowest Common Ancestor (LCA) because there are many genomes present in the database have large fractions of their sequences identical to other genomes. This leads to the result that for well-populated clades with low genome diversity, Kraken only reports species-level assignments for reads from unique regions, while for many species the majority of reads might be classified at a higher level of the taxonomy and the number of reads classified directly to a species may be far lower than the actual number present. Therefore, Kraken’s raw read assignments cannot be directly translated into species- or strain-level abundance estimates. Bracken has been designed to perform sophisticated probabilistically re-distribution of reads to estimate the abundance.

[Bracken](https://github.com/jenniferlu717/Bracken) takes the output from Kraken and estimate the abundance at user specified level: species, genus, or phylum.

Run the script for 3 levels one with species, Genus and Families.
    # sed -i '/Viruses/,$d'  READBASED/SRR15276"$i".kraken_report.out  # remove viruses / bug in the DB formatting


```bash
mkdir READBASED
mkdir READBASED/BRACKEN

for i in {518..547}
do
    sed -i '/Viruses/,$d'  READBASED/SRR15276"$i".kraken_report.out  # remove viruses / bug in the DB formatting

    bracken -d READBASED/k2_standard_20240112/ -t 1000 -i READBASED/SRR15276"$i".kraken_report.out -o READBASED/BRACKEN/SRR15276"$i"_genus.bracken -l G
    bracken -d READBASED/k2_standard_20240112/ -t 1000 -i READBASED/SRR15276"$i".kraken_report.out -o READBASED/BRACKEN/SRR15276"$i"_family.bracken -l F
    bracken -d READBASED/k2_standard_20240112/ -t 1000 -i READBASED/SRR15276"$i".kraken_report.out -o READBASED/BRACKEN/SRR15276"$i"_phylum.bracken -l P
done
```

This step runs very fast, a few seconds. It generates two files for each sample in its corresponding subdirectory inside 03-Kraken: samplename_report_species.txt and samplename.kraken_report_bracken.out. Please take a look at both files to understand what they contain.


**Q: look at the report do you see un expected results?**

<details>
<summary>
HINT
</summary>

> there is human associated reads, whish is more likely can come from contamination. 

</details>  

We need to filter the filnal report and remove possible contamination. 


```bash 
for i in {518..547}
do
    echo "Processing: SRR15276"$i""
    python3 HELPER/filter_bracken.out.py \
    --input READBASED/BRACKEN/SRR15276"$i"_genus.bracken \
    -o READBASED/BRACKEN/SRR15276"$i"_genus_filtered.bracken \
    --exclude  9604
        python3 HELPER/filter_bracken.out.py \
    --input READBASED/BRACKEN/SRR15276"$i"_family.bracken \
    -o READBASED/BRACKEN/SRR15276"$i"_family_filtered.bracken \
    --exclude  9605
        python3 HELPER/filter_bracken.out.py \
    --input READBASED/BRACKEN/SRR15276"$i"_phylum.bracken \
    -o READBASED/BRACKEN/SRR15276"$i"_phylum_filtered.bracken \
    --exclude  7711
done
```

Finally, we are going to combine the abundance estimation for each sample into an abundance table.

```bash 
python3 HELPER/combine_bracken_outputs.py --files READBASED/BRACKEN/SRR152765*_genus_filtered.bracken -o READBASED/BRACKEN/bracken_merged_genus.csv
python3 HELPER/combine_bracken_outputs.py --files READBASED/BRACKEN/SRR152765*_family_filtered.bracken -o READBASED/BRACKEN/bracken_merged_family.csv
python3 HELPER/combine_bracken_outputs.py --files READBASED/BRACKEN/SRR152765*_phylum_filtered.bracken -o READBASED/BRACKEN/bracken_merged_phylum.csv
```

This produces 2 files in the same directory where the input files are:

- *merged_rel_abund.csv*: contains table for all samples with bracken relative abundances and taxonimic assignments
- *merged_read_numbers.csv*: contains table for all samples with bracken read counts and taxonimic assignments

