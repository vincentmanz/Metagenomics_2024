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