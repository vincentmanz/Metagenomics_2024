# Assembling metagenomic data


![binning](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/binning.png)


### Assemblies with MEGAHIT

We used [MEGAHIT](https://github.com/voutcn/megahit) the following way to individually co-assemble each metagenomic set:

```bash
mamba activate ASSEMBLIES_env
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

### Spades

Alternatively [SPADES](https://github.com/ablab/spades) has specific pipelines for metagenomes. 

```bash
mamba activate MEGAHIT_env
mkdir ASSEMBLIES

for i in {518..547l.}
do
    spades   \
        --meta \
        --threads 60 \
        --memory 60 \ 
        --careful \
        --cov-cutoff auto \
        -1 TRIMMEDDATA/SRR15276"$i".R1.fastq.gz \
        -2 TRIMMEDDATA/SRR15276"$i".R2.fastq.gz \
        -k 55,77 \
        -o ASSMBLIES/ 
done
```

### Assembly quality statistics

Let's take a look at the assemblies in a bit more detail with [MetaQUAST](http://bioinf.spbau.ru/metaquast).

```bash
HELPER/metaquast.py -o ASSEMBLIES/QUAST_REPORTS -t 64 ASSEMBLIES/SRR15276*/*.fa  
python3 metaquast.py -f -b -t 60 ASSEMBLIES/SRR15276*.fa/final.contigs.fa   -o ASSEMBLIES/QUAST_REPORTS 
mamba deactivate
```

You can view the results (report.html) in your favorite browser.

Questions about the assembly QC:

> Which assembly has the longest contig when also long reads assemblies are included?

> Which assembly had the most contigs?

> Were the long read assemblies different from the corresponding short read assemblies?

> If yes, in what way?


# Back Mapping: mapping the reads back to the assembly

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



### METABAT2


```bash
mkdir TEST
bowtie2-build --threads 64 ASSEMBLIES/SRR15276518-RAW.fa/final.contigs.fa TEST/SRR15276518_1000nt

bowtie2 --threads 64 -x TEST/SRR15276518_1000nt -1 TRIMMEDDATA/SRR15276518.R1.fastq.gz -2 TRIMMEDDATA/SRR15276518.R2.fastq.gz -S TEST/SRR15276518.sam
samtools view --threads 64 -F 4 -bS TEST/SRR15276518.sam > TEST/SRR15276518_RAW.bam
samtools sort  --threads 64 TEST/SRR15276518_RAW.bam -o TEST/SRR15276518_sorted.bam


runMetaBat.sh ASSEMBLIES/SRR15276518-RAW.fa/final.contigs.fa TEST/SRR15276518_sorted.bam
metabat2 -i ASSEMBLIES/SRR15276518-RAW.fa/final.contigs.fa -a final.contigs.fa.depth.txt -o TEST/SRR15276518_depth_matrix.csv -v

# add the DB in the path: https://github.com/Ecogenomics/CheckM/wiki/Installation#required-reference-data
checkm lineage_wf -t 60 -x fa final.contigs.fa.metabat-bins-20240506_120635/ TEST/


```