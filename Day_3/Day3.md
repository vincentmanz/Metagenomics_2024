


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

We are using a pipeline called [MAGs](https://nf-co.re/mag/2.5.4) to assemble and bin the metagenomic data.

Please find an example here: 

[Installation](https://nf-co.re/docs/usage/installation)


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