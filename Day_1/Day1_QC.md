
# Data and QC

We are going to use the data from the following article: 

Eberhard, F.E., Klimpel, S., Guarneri, A.A. et al. Exposure to Trypanosoma parasites induces changes in the microbiome of the Chagas disease vector Rhodnius prolixus. Microbiome 10, 45 (2022). https://doi.org/10.1186/s40168-022-01240-z

***Abstract***

***Background***

***The causative agent of Chagas disease, Trypanosoma cruzi, and its nonpathogenic relative, Trypanosoma rangeli, are transmitted by haematophagous triatomines and undergo a crucial ontogenetic phase in the insect’s intestine. In the process, the parasites interfere with the host immune system as well as the microbiome present in the digestive tract potentially establishing an environment advantageous for development. However, the coherent interactions between host, pathogen and microbiota have not yet been elucidated in detail. We applied a metagenome shotgun sequencing approach to study the alterations in the microbiota of Rhodnius prolixus, a major vector of Chagas disease, after exposure to T. cruzi and T. rangeli focusing also on the functional capacities present in the intestinal microbiome of the insect.***

***Results***
***The intestinal microbiota of R. prolixus was dominated by the bacterial orders Enterobacterales, Corynebacteriales, Lactobacillales, Clostridiales and Chlamydiales, whereas the latter conceivably originated from the blood used for pathogen exposure. The anterior and posterior midgut samples of the exposed insects showed a reduced overall number of organisms compared to the control group. However, we also found enriched bacterial groups after exposure to T. cruzi as well as T rangeli. While the relative abundance of Enterobacterales and Corynebacteriales decreased considerably, the Lactobacillales, mainly composed of the genus Enterococcus, developed as the most abundant taxonomic group. This applies in particular to vectors challenged with T. rangeli and at early timepoints after exposure to vectors challenged with T. cruzi. Furthermore, we were able to reconstruct four metagenome-assembled genomes from the intestinal samples and elucidate their unique metabolic functionalities within the triatomine microbiome, including the genome of a recently described insect symbiont, Candidatus Symbiopectobacterium, and the secondary metabolites producing bacteria Kocuria spp.***

***Conclusions***
***Our results facilitate a deeper understanding of the processes that take place in the intestinal tract of triatomine vectors during colonisation by trypanosomal parasites and highlight the influential aspects of pathogen-microbiota interactions. In particular, the mostly unexplored metabolic capacities of the insect vector’s microbiome are clearer, underlining its role in the transmission of Chagas disease.***



![communities](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_1/pictures/communities.png)

![pangenomes](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_1/pictures/pangenomes.png)


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

mamba activate meta_quito_2024
Sample,SRA identifier,Type,Time,Number-of-days,Gut,Reads
T1rangPM,SRR15276518,T._rangeli,T1,1,PM,7.70124
T1rangAM,SRR15276519,T._rangeli,T1,1,AM,42.786669
T0ContPM,SRR15276520,Control,T0,0,PM,12.945862
T0ContAM,SRR15276521,Control,T0,0,AM,33.104564
T1cruzPM,SRR15276522,T._cruzi,T1,1,PM,9.564389
T1cruzAM,SRR15276523,T._cruzi,T1,1,AM,18.341077
T1ContAM,SRR15276524,Control,T1,1,AM,19.498048
T1ContPM,SRR15276525,Control,T1,1,PM,9.575527
T0rangPM,SRR15276526,T._rangeli,T0,0,PM,8.858336
```

# Quality control and trimming of the raw reads

### Running fastQC


```bash
USER=your_username
mamba activate meta_quito_2024
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

**Q: What kind of trimming do you think should be done?**



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

![multiqc1](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_1/pictures/fastqc_sequence_counts_plot.png)

![multiqc2](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_1/pictures/fastqc_per_base_sequence_quality_plot.png)

# Removing host reads


For samples that are suppected to have very high levels of host material, it is worth considering adding a step in the library prep stage to remove host material. 
Even with host depletion at the library prep stage, it is important to remove leftover host reads bioinformatically.

A [publication](https://doi.org/10.1099/mgen.0.000393) has done a very thorough study using different software/approaches to detect/remove human reads in microbial datasets. If we are to pick one method, Bowtie2 performs the best overall to remove human reads from the microbial datasets. So, we have chosen to use bowtie2 to remove the host DNA from our dataset.


F-measure is the measure of predictive performance.

![host1](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_1/pictures/rmhost.gif)

![host2](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_1/pictures/rmhost2.gif)

![host3](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_1/pictures/rmhost3.gif)

F-measure is a measure of predictive performance
Based on the findings from the above paper, we are going to use Bowtie2 to remove host DNA. The basic idea is to map/align the preprocessed reads to the bovine reference genome. Then we extract the reads that do not align for downstream analysis.



### Downloading the host genome

We will use the genome of the kissing bug, Triatoma infestans GCA_011037195.1 [Genome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_011037195.1/).

To downoload the **datasets** package, you need to have an account in the [NCBI website](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/).

```bash
mkdir HOST_MAPPING

datasets download genome accession GCA_011037195.1 --include gff3,rna,cds,protein,genome,seq-report --filename Triatoma_infestans.fasta.zip
unzip HOST_MAPPING/Triatoma_infestans.fasta.zip
```

### Mapping the reads to the host genome

note: alternative software (faster) [Sambamba](https://lomereiter.github.io/sambamba/) 


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
    rm HOST_MAPPING/SRR15276"$i".sam HOST_MAPPING/SRR15276"$i"_RAW.bam HOST_MAPPING/SRR15276"$i"_sorted.bam
done
```



