---
title: "Processing data and stats"
author: "Vincent Manzanilla, INTERTRYP"
output: github_document
---

# Day 2

| Time      | Activity                      | Slides                                | Hands-on                                                |
|-----------|-------------------------------|---------------------------------------|---------------------------------------------------------|
| Morning   | Read Base analysis            |                                       |       [Link here](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/Day2_READBASED_stats.md)                         | 




# Metagenomic read based profiles

In this tutorial we will explore and analyse the results of kraken2 read profiling data using the re-estimatd abundances of bracken. These exercises are meant to show how to conceptually approach your data analysis but there are many more and different ways to explore your data. The most important thing to keep in mind is that you have to understand your own data and analyses. One way to achieve this is to perform visual explorations that help you to judge whether the data are appropriate for your question.

Now let’s start the fun!

In R studio.

Load libraries:

**Data Manipulation and Wrangling**


```r
# dplyr: A grammar of data manipulation, providing a consistent set of verbs to solve common data manipulation challenges
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
# tibble: Provides a modern reimagining of data frames, making them more user-friendly
library(tibble)

# tidyverse: A collection of R packages designed for data science, all sharing an underlying design philosophy, grammar, and data structures
library(tidyverse)
```

```
## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
## ✔ forcats   1.0.0     ✔ readr     2.1.5
## ✔ ggplot2   3.5.0     ✔ stringr   1.5.1
## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
## ✔ purrr     1.0.2
```

```
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
```

**Microbiome Analysis**


```r
# microbiome: Provides tools for microbiome data analysis and visualization
library(microbiome)
```

```
## Loading required package: phyloseq
```

```
## 
## microbiome R package (microbiome.github.com)
##     
## 
## 
##  Copyright (C) 2011-2022 Leo Lahti, 
##     Sudarshan Shetty et al. <microbiome.github.io>
```

```
## 
## Attaching package: 'microbiome'
```

```
## The following object is masked from 'package:ggplot2':
## 
##     alpha
```

```
## The following object is masked from 'package:base':
## 
##     transform
```

```r
# microbiomeutilities: Extends functionalities of the microbiome package with additional utilities
library(microbiomeutilities)
```

```
## Warning: replacing previous import 'ggplot2::alpha' by 'microbiome::alpha' when
## loading 'microbiomeutilities'
```

```
## 
## Attaching package: 'microbiomeutilities'
```

```
## The following object is masked from 'package:microbiome':
## 
##     add_refseq
```

```r
# vegan: Functions for ecological analysis, including ordination methods, diversity analysis and other functions for community ecologists
library(vegan)
```

```
## Loading required package: permute
```

```
## Loading required package: lattice
```

```
## This is vegan 2.6-4
```

```
## 
## Attaching package: 'vegan'
```

```
## The following object is masked from 'package:microbiome':
## 
##     diversity
```

```r
# compositions: Provides methods for analyzing compositional data, which is common in fields such as microbiome research
library(compositions)
```

```
## Welcome to compositions, a package for compositional data analysis.
## Find an intro with "? compositions"
```

```
## 
## Attaching package: 'compositions'
```

```
## The following objects are masked from 'package:stats':
## 
##     anova, cor, cov, dist, var
```

```
## The following object is masked from 'package:graphics':
## 
##     segments
```

```
## The following objects are masked from 'package:base':
## 
##     %*%, norm, scale, scale.default
```

```r
# mia: Microbiome Analysis, provides functionalities for microbiome data analysis and visualization
library(mia)
```

```
## Loading required package: SummarizedExperiment
```

```
## Loading required package: MatrixGenerics
```

```
## Loading required package: matrixStats
```

```
## 
## Attaching package: 'matrixStats'
```

```
## The following object is masked from 'package:dplyr':
## 
##     count
```

```
## 
## Attaching package: 'MatrixGenerics'
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
##     colWeightedMeans, colWeightedMedians, colWeightedSds,
##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
##     rowWeightedSds, rowWeightedVars
```

```
## Loading required package: GenomicRanges
```

```
## Loading required package: stats4
```

```
## Loading required package: BiocGenerics
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:compositions':
## 
##     normalize, var
```

```
## The following objects are masked from 'package:lubridate':
## 
##     intersect, setdiff, union
```

```
## The following objects are masked from 'package:dplyr':
## 
##     combine, intersect, setdiff, union
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unsplit, which.max, which.min
```

```
## Loading required package: S4Vectors
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following objects are masked from 'package:lubridate':
## 
##     second, second<-
```

```
## The following object is masked from 'package:tidyr':
## 
##     expand
```

```
## The following objects are masked from 'package:dplyr':
## 
##     first, rename
```

```
## The following object is masked from 'package:utils':
## 
##     findMatches
```

```
## The following objects are masked from 'package:base':
## 
##     expand.grid, I, unname
```

```
## Loading required package: IRanges
```

```
## 
## Attaching package: 'IRanges'
```

```
## The following object is masked from 'package:microbiome':
## 
##     coverage
```

```
## The following object is masked from 'package:phyloseq':
## 
##     distance
```

```
## The following object is masked from 'package:lubridate':
## 
##     %within%
```

```
## The following object is masked from 'package:purrr':
## 
##     reduce
```

```
## The following objects are masked from 'package:dplyr':
## 
##     collapse, desc, slice
```

```
## Loading required package: GenomeInfoDb
```

```
## Loading required package: Biobase
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## 
## Attaching package: 'Biobase'
```

```
## The following object is masked from 'package:MatrixGenerics':
## 
##     rowMedians
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     anyMissing, rowMedians
```

```
## The following object is masked from 'package:phyloseq':
## 
##     sampleNames
```

```
## Loading required package: SingleCellExperiment
```

```
## Loading required package: TreeSummarizedExperiment
```

```
## Loading required package: Biostrings
```

```
## Loading required package: XVector
```

```
## 
## Attaching package: 'XVector'
```

```
## The following object is masked from 'package:purrr':
## 
##     compact
```

```
## 
## Attaching package: 'Biostrings'
```

```
## The following object is masked from 'package:base':
## 
##     strsplit
```

```
## Loading required package: MultiAssayExperiment
```

```
## 
## Attaching package: 'mia'
```

```
## The following objects are masked from 'package:dplyr':
## 
##     full_join, inner_join, left_join, right_join
```

```r
# miaViz: Visualization utilities for microbiome analysis, part of the mia package ecosystem
library(miaViz)
```

```
## Loading required package: ggraph
```

```r
# phyloseq: Microbiome Analysis, provides functionalities for microbiome data analysis and visualization
library(phyloseq)
```

**Data Visualization**


```r
# ggplot2: A system for creating elegant and versatile data visualizations based on the grammar of graphics
library(ggplot2)

# hrbrthemes: Contains additional themes, theme components, and utilities for ggplot2
library(hrbrthemes)

# pheatmap: Pretty heatmaps, provides more control over the heatmap visualization
library(pheatmap)

# RColorBrewer: Provides color palettes for visualizing data, particularly useful in ggplot2 visualizations
library(RColorBrewer)

# ggrepel: Extends ggplot2 by adding better text label placement to avoid overlaps
library(ggrepel)

# patchwork: Makes it easier to combine multiple ggplot2 plots into one overall plot layout
library(patchwork)
```

## 1. load the data

First we need to load our data. Usually the biggest bottleneck between raw data and analyses is to get the data in the right shape for your purpose. Often this requires a little bit of data mingling. On this road - google is your best friend to master the R universe :)

Let’s first load the relative abundance table of the bracken results.

For this part we are using the [phyloseq](https://joey711.github.io/phyloseq/) package. The phyloseq package is a tool to import, store, analyze, and graphically display complex phylogenetic sequencing data that has already been clustered into Operational Taxonomic Units (OTUs), especially when there is associated sample data, phylogenetic tree, and/or taxonomic assignment of the OTUs. 

First st your working directory: 


```r
setwd("~/Documents/project/Metagenomics_2024/Metagenomics_2024/Day_2/")
```

And load your two files:


```r
merged_metagenomes <- phyloseq::import_biom("../DATA//merge_species.biom") 
meta <- read.csv(file = "../DATA/tryp_metadata.csv", sep = ",")
```

#### phyloseq-ize Data

Any data already in an R session can be annoated/coerced to be recognized by phyloseq’s functions and methods. This is important, because there are lots of ways you might receive data related to a microbiome project, and not all of these will come from a popular server or workflow that is already supported by a phyloseq import function. 

![](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/phyloseq.png)

- **otu_table** - Works on any numeric matrix. You must also specify if the species are rows or columns
- **tax_table** - Works on any character matrix. The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
- **sample_data** - Works on any data.frame. The rownames must match the sample names in the otu_table if you plan to combine them as a phyloseq-object

 Let's look at the phlyseq object. 


```r
merged_metagenomes
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 2251 taxa and 30 samples ]
## tax_table()   Taxonomy Table:    [ 2251 taxa by 7 taxonomic ranks ]
```

Summarize the data.  


```r
# Summarize the phyloseq object 'merged_metagenomes'
#microbiome::summarize_phyloseq(merged_metagenomes)

# Display the first few rows of the OTU (Operational Taxonomic Unit) table
head(otu_table(merged_metagenomes))
```

```
## OTU Table:          [6 taxa and 30 samples]
##                      taxa are rows
##         SRR15276518.kraken_report_bracken_species
## 1351                                       225235
## 2920934                                      3837
## 1352                                          522
## 417368                                         32
## 118060                                         25
## 167972                                         22
##         SRR15276519.kraken_report_bracken_species
## 1351                                       910112
## 2920934                                     16203
## 1352                                          399
## 417368                                        125
## 118060                                        101
## 167972                                         98
##         SRR15276520.kraken_report_bracken_species
## 1351                                       475076
## 2920934                                      9116
## 1352                                          795
## 417368                                         63
## 118060                                         62
## 167972                                         56
##         SRR15276521.kraken_report_bracken_species
## 1351                                         4024
## 2920934                                        67
## 1352                                          816
## 417368                                          0
## 118060                                          0
## 167972                                          0
##         SRR15276522.kraken_report_bracken_species
## 1351                                      2205946
## 2920934                                     40456
## 1352                                          822
## 417368                                        287
## 118060                                        318
## 167972                                        293
##         SRR15276523.kraken_report_bracken_species
## 1351                                      1860911
## 2920934                                     30836
## 1352                                          575
## 417368                                        229
## 118060                                        246
## 167972                                        203
##         SRR15276524.kraken_report_bracken_species
## 1351                                       793332
## 2920934                                     14336
## 1352                                          418
## 417368                                         98
## 118060                                         85
## 167972                                         85
##         SRR15276525.kraken_report_bracken_species
## 1351                                       207171
## 2920934                                      3854
## 1352                                          826
## 417368                                         29
## 118060                                         23
## 167972                                         21
##         SRR15276526.kraken_report_bracken_species
## 1351                                        81067
## 2920934                                      1571
## 1352                                          731
## 417368                                          0
## 118060                                          0
## 167972                                          0
##         SRR15276527.kraken_report_bracken_species
## 1351                                         7105
## 2920934                                       111
## 1352                                            0
## 417368                                          0
## 118060                                          0
## 167972                                          0
##         SRR15276528.kraken_report_bracken_species
## 1351                                        66193
## 2920934                                      1216
## 1352                                          784
## 417368                                          0
## 118060                                          0
## 167972                                          0
##         SRR15276529.kraken_report_bracken_species
## 1351                                      2364001
## 2920934                                     41302
## 1352                                         1471
## 417368                                        346
## 118060                                        324
## 167972                                        333
##         SRR15276530.kraken_report_bracken_species
## 1351                                      1314916
## 2920934                                     24422
## 1352                                          993
## 417368                                        175
## 118060                                        157
## 167972                                        203
##         SRR15276531.kraken_report_bracken_species
## 1351                                      1010148
## 2920934                                     18448
## 1352                                          341
## 417368                                        133
## 118060                                        118
## 167972                                        144
##         SRR15276532.kraken_report_bracken_species
## 1351                                      2343002
## 2920934                                     40276
## 1352                                          679
## 417368                                        295
## 118060                                        334
## 167972                                        330
##         SRR15276533.kraken_report_bracken_species
## 1351                                       718800
## 2920934                                     12881
## 1352                                          585
## 417368                                         89
## 118060                                        125
## 167972                                         87
##         SRR15276534.kraken_report_bracken_species
## 1351                                       902748
## 2920934                                     16141
## 1352                                          963
## 417368                                        132
## 118060                                        132
## 167972                                        100
##         SRR15276535.kraken_report_bracken_species
## 1351                                      3662690
## 2920934                                     67130
## 1352                                         1743
## 417368                                        481
## 118060                                        477
## 167972                                        431
##         SRR15276536.kraken_report_bracken_species
## 1351                                      3425155
## 2920934                                     62794
## 1352                                         2033
## 417368                                        451
## 118060                                        519
## 167972                                        444
##         SRR15276537.kraken_report_bracken_species
## 1351                                      1666832
## 2920934                                     30963
## 1352                                          987
## 417368                                        215
## 118060                                        223
## 167972                                        262
##         SRR15276538.kraken_report_bracken_species
## 1351                                      3626804
## 2920934                                     66022
## 1352                                         1534
## 417368                                        544
## 118060                                        516
## 167972                                        455
##         SRR15276539.kraken_report_bracken_species
## 1351                                         2086
## 2920934                                         0
## 1352                                          938
## 417368                                          0
## 118060                                          0
## 167972                                          0
##         SRR15276540.kraken_report_bracken_species
## 1351                                       255083
## 2920934                                      4507
## 1352                                          445
## 417368                                         46
## 118060                                         33
## 167972                                         31
##         SRR15276541.kraken_report_bracken_species
## 1351                                      1139637
## 2920934                                     20546
## 1352                                          434
## 417368                                        136
## 118060                                        164
## 167972                                        125
##         SRR15276542.kraken_report_bracken_species
## 1351                                      1036494
## 2920934                                     19109
## 1352                                          600
## 417368                                        176
## 118060                                        117
## 167972                                        165
##         SRR15276543.kraken_report_bracken_species
## 1351                                      1838154
## 2920934                                     34128
## 1352                                         1056
## 417368                                        220
## 118060                                        255
## 167972                                        219
##         SRR15276544.kraken_report_bracken_species
## 1351                                      1641627
## 2920934                                     29517
## 1352                                         1194
## 417368                                        191
## 118060                                        203
## 167972                                        204
##         SRR15276545.kraken_report_bracken_species
## 1351                                      2230928
## 2920934                                     39380
## 1352                                         1195
## 417368                                        281
## 118060                                        328
## 167972                                        281
##         SRR15276546.kraken_report_bracken_species
## 1351                                       363473
## 2920934                                      6730
## 1352                                          406
## 417368                                         50
## 118060                                         39
## 167972                                         49
##         SRR15276547.kraken_report_bracken_species
## 1351                                      1198563
## 2920934                                     21856
## 1352                                          603
## 417368                                        138
## 118060                                        159
## 167972                                        119
```

```r
# Display the first few rows of the taxonomy table
head(tax_table(merged_metagenomes))
```

```
## Taxonomy Table:     [6 taxa by 7 taxonomic ranks]:
##         Rank1         Rank2          Rank3        Rank4               
## 1351    "k__Bacteria" "p__Bacillota" "c__Bacilli" "o__Lactobacillales"
## 2920934 "k__Bacteria" "p__Bacillota" "c__Bacilli" "o__Lactobacillales"
## 1352    "k__Bacteria" "p__Bacillota" "c__Bacilli" "o__Lactobacillales"
## 417368  "k__Bacteria" "p__Bacillota" "c__Bacilli" "o__Lactobacillales"
## 118060  "k__Bacteria" "p__Bacillota" "c__Bacilli" "o__Lactobacillales"
## 167972  "k__Bacteria" "p__Bacillota" "c__Bacilli" "o__Lactobacillales"
##         Rank5                Rank6             Rank7                           
## 1351    "f__Enterococcaceae" "g__Enterococcus" "s__faecalis"                   
## 2920934 "f__Enterococcaceae" "g__Enterococcus" "s__sp. LX10"                   
## 1352    "f__Enterococcaceae" "g__Enterococcus" "s__faecium"                    
## 417368  "f__Enterococcaceae" "g__Enterococcus" "s__thailandicus"               
## 118060  "f__Enterococcaceae" "g__Enterococcus" "s__rotai"                      
## 167972  "f__Enterococcaceae" "g__Enterococcus" "s__uncultured Enterococcus sp."
```

```r
# Display the first few rows of the sample data associated with 'merged_metagenomes'
#head(sample_data(merged_metagenomes))

# Get the sample variables of the phyloseq object 'merged_metagenomes'
#sample_variables(merged_metagenomes)
```



We need to add the metadata to the phyloseq object. 

## 2. format the data


```r
# Sort the 'meta' data frame by the 'SRA.identifier' column
meta <- meta %>% arrange(row_number(SRA.identifier))

# Associate the sorted metadata to the phyloseq object as sample data
merged_metagenomes@sam_data <- sample_data(meta)

# Extract the sample names from the 'meta' data frame
column_name <- meta %>% pull(Sample)

# Associate the extracted sample names to the phyloseq object
sample_names(merged_metagenomes) <- column_name

# Remove the unnecessary 'k_' prefix in the taxonomy data
merged_metagenomes@tax_table@.Data <- substring(merged_metagenomes@tax_table@.Data, 4)

# Rename the columns of the taxonomy table to represent taxonomic ranks
colnames(merged_metagenomes@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
```

## 3. Basic stats

Before we start anything, let’s just check out or data a little bit (sanity check). Never go blind into your analyses.


```r
head(psmelt(merged_metagenomes))
```

```
## Warning in psmelt(merged_metagenomes): The sample variables: 
## Sample
##  have been renamed to: 
## sample_Sample
## to avoid conflicts with special phyloseq plot attribute names.
```

```
##        OTU   Sample Abundance sample_Sample SRA.identifier       Type Time
## 5029  1351 T3rangAM   3662690      T3rangAM    SRR15276535 T._rangeli   T3
## 5039  1351 T3cruzAM   3626804      T3cruzAM    SRR15276538   T._cruzi   T3
## 5016  1351 T3rangPM   3425155      T3rangPM    SRR15276536 T._rangeli   T3
## 18354 2035 T0rangAM   2602231      T0rangAM    SRR15276527 T._rangeli   T0
## 5036  1351 T7rangAM   2364001      T7rangAM    SRR15276529 T._rangeli   T7
## 5021  1351 T7cruzAM   2343002      T7cruzAM    SRR15276532   T._cruzi   T7
##       Number.of.days Gut     Reads  Kingdom         Phylum         Class
## 5029               3  AM 21.104708 Bacteria      Bacillota       Bacilli
## 5039               3  AM 11.656514 Bacteria      Bacillota       Bacilli
## 5016               3  PM  8.744824 Bacteria      Bacillota       Bacilli
## 18354              0  AM 35.555314 Bacteria Actinomycetota Actinomycetes
## 5036               7  AM 12.138032 Bacteria      Bacillota       Bacilli
## 5021               7  AM 13.686443 Bacteria      Bacillota       Bacilli
##                 Order            Family          Genus        Species
## 5029  Lactobacillales   Enterococcaceae   Enterococcus       faecalis
## 5039  Lactobacillales   Enterococcaceae   Enterococcus       faecalis
## 5016  Lactobacillales   Enterococcaceae   Enterococcus       faecalis
## 18354   Micrococcales Microbacteriaceae Curtobacterium flaccumfaciens
## 5036  Lactobacillales   Enterococcaceae   Enterococcus       faecalis
## 5021  Lactobacillales   Enterococcaceae   Enterococcus       faecalis
```

</details>  

**Q: How many species and samples are detected in our dataset?**

<details>
<summary>
HINT
</summary>

> ntaxa(merged_metagenomes), nsamples.

</details>  

**Q: Look at the taxonomy  is ther some problem with the taxonomy?**

<details>
<summary>
HINT
</summary>

>  We have contaminant in the data. Which one? (head(psmelt(merged_metagenomes)))

</details>  

**Q: what are the taxa with the more abundant reads?**

<details>
<summary>
HINT
</summary>

>  head(otu_table(merged_metagenomes))

</details> 

Extraction of the sample's names. 


```r
sample_names(merged_metagenomes)
```

```
##  [1] "T1rangPM" "T1rangAM" "T0ContPM" "T0ContAM" "T1cruzPM" "T1cruzAM"
##  [7] "T1ContAM" "T1ContPM" "T0rangPM" "T0rangAM" "T0cruzPM" "T7rangAM"
## [13] "T7rangPM" "T7cruzPM" "T7cruzAM" "T7ContPM" "T7ContAM" "T3rangAM"
## [19] "T3rangPM" "T3cruzPM" "T3cruzAM" "T0cruzAM" "T3ContPM" "T3ContAM"
## [25] "T2rangPM" "T2rangAM" "T2cruzPM" "T2cruzAM" "T2ContPM" "T2ContAM"
```

### Aggregation

Microbial species can be called at multiple taxonomic resolutions. We can easily agglomerate the data based on taxonomic ranks. Here, we agglomerate the data at Family level.


```r
# Aggregate rare taxa at the family level for the phyloseq object 'merged_metagenomes'
merged_metagenomes_family <- aggregate_rare(merged_metagenomes, level = "Genus", detection = 0/100, prevalence = 0/100)

# Display the dimensionality of the abundances of 'merged_metagenomes_family'
dim(abundances(merged_metagenomes_family))
```

```
## [1] 788  30
```

**Q: How many sample and tax do we have now?**

<details>
<summary>
HINT
</summary>

> There are XX taxa and 30 samples, meaning that there are 29 different Phylum level taxonomic groups. Looking at the rowData after agglomeration shows all Enterococcaceae are combined together, and all lower rank information is lost.


```r
head(tax_table(merged_metagenomes_family))
```

```
## Taxonomy Table:     [6 taxa by 2 taxonomic ranks]:
##                        Genus                    unique                  
## Abyssalbus             "Abyssalbus"             "Abyssalbus"            
## Achromobacter          "Achromobacter"          "Achromobacter"         
## Acidipropionibacterium "Acidipropionibacterium" "Acidipropionibacterium"
## Acidothermus           "Acidothermus"           "Acidothermus"          
## Acidovorax             "Acidovorax"             "Acidovorax"            
## Acinetobacter          "Acinetobacter"          "Acinetobacter"
```

</details> 


## 4. QC and Pre-process data

Now that we know a little bit about our data we can start the pre-processing. 


###  Library size / read count

Let us check for distribution of number of sequences retained from the Kraken/Bracken approach.


```r
plot_read_distribution(merged_metagenomes, "Type", plot.type = "density")
```
![](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/density_reads.png)

*You can try to plot with diffrent metadata.*


Plotting the read count per sample


```r
df <- psmelt(merged_metagenomes)  %>%  group_by(Sample, Time) %>%  
  summarise(sum_reads = sum(Abundance)) %>% arrange(sum_reads) 

ggplot(df) +
  geom_bar(aes(reorder(Sample, -sum_reads), sum_reads, fill=Time),
           col="red", alpha = .2, stat="identity") 
```


![](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/reads_count.png)


### Prevalence - Detection 


Prevalence quantifies the frequency of samples where certain microbes were detected (above a given detection threshold). The prevalence can be given as sample size (N) or percentage (unit interval).

The population prevalence (frequency) at a 1% relative abundance threshold (detection = 1/100 and count = false), can look like this.


```r
head(prevalence(merged_metagenomes, detection=1/100, sort=TRUE, count=FALSE))
```

```
##  85471   5855   3218   3369 210225   3469 
##      1      1      1      1      1      1
```

The function arguments detection and count can also be used to access, how many samples do pass a threshold for raw counts. Here, the population prevalence (frequency) at the absolute abundance threshold (count=true) at read count 1 (detection = 1) is accessed.


```r
# Calculate the prevalence of taxa in the phyloseq object 'merged_metagenomes'
# Detection threshold is set to 1, taxa are sorted by prevalence, and count of taxa is returned
head(prevalence(merged_metagenomes, detection = 1, sort = TRUE, count = TRUE))

# Plot the prevalence of taxa at the Phylum level with a detection threshold of 10000
plot_taxa_prevalence(merged_metagenomes, level = "Phylum", detection = 10000)

# Alternative plot: Create a scatter plot of the coefficient of variation (CV) of taxa in 'merged_metagenomes'
# The CV is a measure of relative variability, calculated as the standard deviation divided by the mean
p1 <- plot_taxa_cv(merged_metagenomes, plot.type = "scatter")

# Apply a log10 scale to the x-axis of the scatter plot
p1 + scale_x_log10()
```

Each point corresponds to a different or unique taxon. The y-axis represents the fraction of samples, these taxa are present. The low prevalence suggests there is a low overlap across samples. 

![](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/prevalence.png)

**Q: Which taxa is present in all samples and with a high abundance?**

<details>
<summary>
HINT
</summary>

> Homo sapiens

</details> 


#### Contaminant sequences
Samples might be contaminated with exogenous sequences. We have observed 1 contaminants Homo sapiens.



```r
#check 
head(psmelt(subset_taxa(merged_metagenomes, Genus == "Trypanosoma")))
```

```
## Warning in psmelt(subset_taxa(merged_metagenomes, Genus == "Trypanosoma")): The sample variables: 
## Sample
##  have been renamed to: 
## sample_Sample
## to avoid conflicts with special phyloseq plot attribute names.
```

```
##     OTU   Sample Abundance sample_Sample SRA.identifier       Type Time
## 11 5691 T1rangAM        69      T1rangAM    SRR15276519 T._rangeli   T1
## 24 5691 T3rangPM        51      T3rangPM    SRR15276536 T._rangeli   T3
## 1  5691 T0ContAM        47      T0ContAM    SRR15276521    Control   T0
## 17 5691 T2rangAM        40      T2rangAM    SRR15276543 T._rangeli   T2
## 5  5691 T0rangAM        38      T0rangAM    SRR15276527 T._rangeli   T0
## 23 5691 T3rangAM        32      T3rangAM    SRR15276535 T._rangeli   T3
##    Number.of.days Gut     Reads   Kingdom     Phylum         Class
## 11              1  AM 42.786669 Eukaryota Euglenozoa Kinetoplastea
## 24              3  PM  8.744824 Eukaryota Euglenozoa Kinetoplastea
## 1               0  AM 33.104564 Eukaryota Euglenozoa Kinetoplastea
## 17              2  AM 25.297229 Eukaryota Euglenozoa Kinetoplastea
## 5               0  AM 35.555314 Eukaryota Euglenozoa Kinetoplastea
## 23              3  AM 21.104708 Eukaryota Euglenozoa Kinetoplastea
##              Order           Family       Genus Species
## 11 Trypanosomatida Trypanosomatidae Trypanosoma  brucei
## 24 Trypanosomatida Trypanosomatidae Trypanosoma  brucei
## 1  Trypanosomatida Trypanosomatidae Trypanosoma  brucei
## 17 Trypanosomatida Trypanosomatidae Trypanosoma  brucei
## 5  Trypanosomatida Trypanosomatidae Trypanosoma  brucei
## 23 Trypanosomatida Trypanosomatidae Trypanosoma  brucei
```

```r
head(psmelt(subset_taxa(merged_metagenomes, Genus == "Homo")))
```

```
## Warning in psmelt(subset_taxa(merged_metagenomes, Genus == "Homo")): The sample variables: 
## Sample
##  have been renamed to: 
## sample_Sample
## to avoid conflicts with special phyloseq plot attribute names.
```

```
##     OTU   Sample Abundance sample_Sample SRA.identifier       Type Time
## 11 9606 T1rangAM    769548      T1rangAM    SRR15276519 T._rangeli   T1
## 5  9606 T0rangAM    631939      T0rangAM    SRR15276527 T._rangeli   T0
## 1  9606 T0ContAM    579258      T0ContAM    SRR15276521    Control   T0
## 17 9606 T2rangAM    440572      T2rangAM    SRR15276543 T._rangeli   T2
## 7  9606 T1ContAM    326725      T1ContAM    SRR15276524    Control   T1
## 3  9606 T0cruzAM    303760      T0cruzAM    SRR15276539   T._cruzi   T0
##    Number.of.days Gut    Reads   Kingdom   Phylum    Class    Order    Family
## 11              1  AM 42.78667 Eukaryota Chordata Mammalia Primates Hominidae
## 5               0  AM 35.55531 Eukaryota Chordata Mammalia Primates Hominidae
## 1               0  AM 33.10456 Eukaryota Chordata Mammalia Primates Hominidae
## 17              2  AM 25.29723 Eukaryota Chordata Mammalia Primates Hominidae
## 7               1  AM 19.49805 Eukaryota Chordata Mammalia Primates Hominidae
## 3               0  AM 19.60417 Eukaryota Chordata Mammalia Primates Hominidae
##    Genus Species
## 11  Homo sapiens
## 5   Homo sapiens
## 1   Homo sapiens
## 17  Homo sapiens
## 7   Homo sapiens
## 3   Homo sapiens
```

```r
#Keep only the kingdom of interest
merged_metagenomes <- subset_taxa(merged_metagenomes, Kingdom %in% c("Archaea", "Bacteria", "Fungi", "Viruses"))
```

Now recheck that your data are clean before continuing the analysis. 

**Q: Which taxa is present in all samples and with a high abundance?**

<details>
<summary>
HINT
</summary>

> head(tax_table(merged_metagenomes)) or plot_taxa_prevalence(merged_metagenomes, level="Phylum", detection = 10000)

</details> 

### Run rarefaction curves

Before normalization by sub-sampling, let’s have a look at rarefaction curves, evaluate your sequencing effort and make decisions. 
Rarefaction is a statistical technique employed to evaluate species richness based on sampling results. The process involved subsampling reads without replacement to a defined sequencing depth, thereby creating a standardized library size across samples. Any sample with a total read count less than the defined sequencing depth used to rarefy will be discarded. What we want to see in those curves is that it reaches a plateau meaning no new OTUs are discovered (going up on the Y axis) as sequencing is deeper (the X axis). If this happens we can be pretty confident that we sequenced all the OTUs that were present in the sample.


```r
install.packages("remotes")
remotes::install_github("gauravsk/ranacapa")
```


```r
ranacapa::ggrare(pseq, step = 10, color = "Sample", se = FALSE) +
    geom_vline(xintercept = min(sample_sums(pseq)), color = "gray60")
```

![](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/rare.png)
 ***Q: What do you observe?***

### Data Normalization

Data transformations are common in (microbial) ecology (Legendre 2001) and used to improve compatibility with assumptions related to specific statistical methods, mitigate biases, enhance the comparability of samples or features, or to obtain more interpretable values.

Examples include the logarithmic transformation, calculation of relative abundances (percentages), and compositionality-aware transformations such as the centered log-ratio transformation (clr).

Let us summarize some commonly used transformations in microbiome data science; further details and benchmarkings available in the references:

- *relabundance* relative transformation; also known as total sum scaling (TSS) and compositional transformation. This converts counts into percentages (at the scale [0, 1]) that sum up to. Much of the currently available taxonomic abundance data from high-throughput assays (16S, metagenomic sequencing) is compositional by nature, even if the data is provided as counts (Gloor et al. 2017).

- *clr* Centered log ratio transformation (Aitchison 1986) is used to reduce data skewness and compositionality bias in relative abundances, while bringing the data to the logarithmic scale. This transformation is frequently applied in microbial ecology (Gloor et al. 2017). However, this transformation only applies to positive values. Usual solution is to add pseudocount, which adds another type of bias in the data. The robust clr transformation (‘rclr’) aims to circumvent the need to add a pseudocount. While the resulting values from these transformations are difficult interpret directly, this transformation may enhance comparability of relative differences between samples. It is part of a broader Aitchison family of transformations; the additive log ratio transformation (`alr’) is also available. The robust clr (“rclr”) is similar to regular clr (see above) but allows data with zeroes and avoids the need to add pseudocount Martino et al. (2019). See library *compositions*. 

- *pa* presence/absence transformation ignores abundances and only indicates whether the given feature is detected above the given threshold (default: 0). This simple transformation is relatively widely used in ecological research. It has shown good performance in microbiome-based classification performance (Giliberti et al. 2022, Karwowska2024).

- *z* Z transformation scales data to zero mean and unit variance; this us used to bring features (or samples) to more comparable levels in terms of mean and scale of the values. This can enhance visualization and interpretation of the data

- *log*, *log2*, *log10* Logarithmic transformations; used e.g. to reduce data skewness; with compositional data the clr (or rclr) transformation is often preferred.

- *hellinger* Hellinger transformation equals to the square root of relative abundances. This ecological transformation can be useful if we are interested in changes in relative abundances.

- *rank* Rank transformation replaces each value by its rank. Also see ‘rrank’ (relative rank transformation). This has use for instance in non-parametric statistics.

- Other available transformations include Chi square (‘chi.square’), Frequency transformation (‘frequency’), and Make margin sum of squares equal to one (‘normalize’)


The data contains read counts. We can convert these into relative abundances and other formats. Compare abundance of a given taxonomic group using the example data before and after the compositionality transformation (with a cross-plot, for instance). You can also compare the results to CLR-transformed data (see e.g. Gloor et al. 2017)

Have a look at the function *?microbiome::transform*.



```r
microbiome::transform(merged_metagenomes, transform = "compositional")
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 2025 taxa and 30 samples ]
## sample_data() Sample Data:       [ 30 samples by 7 sample variables ]
## tax_table()   Taxonomy Table:    [ 2025 taxa by 7 taxonomic ranks ]
```

In order to assess the effect of the different transformations, you can use the *clr* or the *log10* for the  downstream analysis. 

## 5. Microbiome composition

Microbial abundances are typically ‘compositional’ (relative) in the current microbiome profiling data sets. This is due to technical aspects of the data generation process (see e.g. Gloor et al., 2017).

The next example calculates relative abundances as these are usually easier to interpret than plain counts. For some statistical models we need to transform the data into other formats as explained in above link (and as we will see later).

*detection*: Detection threshold for absence/presence (percentage reads).
*prevalence*: Prevalence threshold (in [0, 1]) (presence accross samples).


```r
pseq <- aggregate_rare(merged_metagenomes, level = "Genus", detection = 0.05/100, prevalence = 20/100)
pseq <- microbiome::transform(pseq, transform = "compositional")
```


**Q: what is the meaning and effect of compositional?**

<details>
<summary>
HINT
</summary>

>  abundances(pseq)

</details>  


#### Visualization

**Bar plot** are useful to have a broad look at the data. 


```r
# Create a plot of the taxonomic composition at the Family level
p <- plot_composition(pseq,
                      taxonomic.level = "Species", # Specify the taxonomic level for the plot
                      sample.sort = "Sample",     # Sort samples by sample identifier
                      x.label = "Sample",         # Label for the x-axis
                      group_by = "Gut") +          # Group samples by the 'Time' variable

  # Use a color palette from the Color Brewer for filling the Family groups
  scale_fill_brewer("Species", palette = "magma-viridis") +

  # Arrange the legend items in a single column
  guides(fill = guide_legend(ncol = 1)) +

  # Convert the y-axis to represent relative abundance as a percentage
  scale_y_percent() +

  # Add labels and a title to the plot
  labs(x = "Samples", y = "Relative abundance (%)",
       title = "Relative abundance data") + 

  # Apply a clean theme with a horizontal grid
  theme_ipsum(grid = "Y") +

  # Customize the theme for x-axis text and legend text
  theme(axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis text 90 degrees for readability
        legend.text = element_text(face = "italic"),    # Italicize the legend text
        legend.position = "none")

# Print the plot
print(p)
```

![](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/barplot_family.png)

**Q: Does the profiles look similar between samples? Can you spot any trends?**

<details>
<summary>
HINT
</summary>

> A: We see that the same one or two OTU dominates all samples but there is some variability between samples. We can observe some trend on the Entrococcus proption over time.

</details>  


*If you have time you can also visualize the other taxonomic levels (e.g. species) with the same approach. Try to come up with the code yourself.*


```r
# Aggregate rare taxa at the Genus level for the phyloseq object 'merged_metagenomes'
# Taxa are included if they have a detection threshold of 0.05% and a prevalence of 20%
pseq <- aggregate_rare(merged_metagenomes, level = "Phylum", detection = 0.05/100, prevalence = 20/100)

# Transform the data to compositional (relative abundance) format
pseq <- microbiome::transform(pseq, transform = "compositional")



n <- dim(tax_table(pseq))[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))



# Create a plot of the taxonomic composition at the Genus level
p <- plot_composition(pseq,
                      taxonomic.level = "Phylum", # Specify the taxonomic level for the plot
                      sample.sort = "Sample",    # Sort samples by sample identifier
                      x.label = "Sample",        # Label for the x-axis
                      otu.sort = "abundance",    # Sort OTUs by their abundance
                      group_by = "Type") +       # Group samples by the 'Type' variable

  # Use a color palette from the Color Brewer for filling the Genus groups
  scale_fill_manual(values = col_vector) +

  # Arrange the legend items in a single column
  guides(fill = guide_legend(ncol = 1)) +

  # Convert the y-axis to represent relative abundance as a percentage
  scale_y_percent() +

  # Add labels and a title to the plot
  labs(x = "Samples", y = "Relative abundance (%)",
       title = "Relative abundance data") + 

  # Apply a clean theme with a horizontal grid
  theme_ipsum(grid = "Y") +

  # Customize the theme for x-axis text and legend text
  theme(axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis text 90 degrees for readability
        legend.text = element_text(face = "italic"))       # Italicize the legend text

# Print the plot
print(p)
```

![](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/barplot_genus.png)


## 6 Diversity

Species diversity, in its simplest definition, is the number of species in a particular area and their relative abundance (evenness). Once we know the taxonomic composition of our metagenomes, we can do diversity analyses. Here we will discuss the two most used diversity metrics, α diversity (within one metagenome) and β (across metagenomes).

- *α* Diversity: Can be represented only as richness (, i.e., the number of different species in an environment), or it can be measured considering the abundance of the species in the environment as well (i.e., the number of individuals of each species inside the environment). To measure α-diversity, we use indexes such as Shannon’s, Simpson’s, Chao1, etc.

![](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/diversity1.png)
*Alpha diversity is calculated according to fish diversity in a pond. Here, alpha diversity is represented in its simplest way: Richness.*

In the next example, we will look at the α and the β components of the diversity of a dataset of fishes in three lakes. The most simple way to calculate the β-diversity is to calculate the distinct species between two lakes (sites). Let us take as an example the diversity between Lake A and Lake B. The number of species in Lake A is 3. To this quantity, we will subtract the number of these species that are shared with the Lake B: 2. So the number of unique species in Lake A compared to Lake B is (3-2) = 1. To this number, we will sum the result of the same operations but now take Lake B as our reference site. In the end, the β diversity between Lake A and Lake B is (3-2) + (3-2) = 2. This process can be repeated, taking each pair of lakes as the focused sites.

![](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/diversity2.png)
*Alpha and beta diversity indexes of fishes in a pond.*

- *β* diversity mesures how different two or more communities are, either in their composition (richness) or in the abundance of the organisms that compose it (abundance).

- *Bray-Curtis dissimilarity*: The difference in richness and abundance across environments (samples). Weight on abundance. Measures the differences from 0 (equal communities) to 1 (different communities)
- *Jaccard distance*: Based on the presence/absence of species (diversity). It goes from 0 (same species in the community) to 1 (no species in common)
- *UniFrac*: Measures the phylogenetic distance; how alike the trees in each community are. There are two types, without weights (diversity) and with weights (diversity and abundance)
There are different ways to plot and show the results of such analysis. Among others, PCA, PCoA, or NMDS analysis are widely used.

![](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/diversity3.png)

**Q: Which of the options below is true for the alpha diversity in lakes A, B, and beta diversity between lakes A and B, respectively?**

- 4, 3, 1
- 4, 3, 5
- 9, 7, 16


<details>
<summary>
HINT
</summary>

> 4, 3, 5 Alpha diversity in this case, is the sum of different species. Lake A has 4 different species and lake B has 3 different species. Beta diversity refers to the difference between lake A and lake B. If we use the formula in Figure 2 we can see that to calculate beta diversity, we have to detect the number of species and the number of shared species in both lakes. There is only one shared species, so we have to subtract the number of shared species to the total species and sum the result. In this case, in lake A, we have 4 different species and one shared species with lake B (4-1)=3, and in lake B we have three species and one shared species with lake A (3-1)=2. If we add 3+2, the result is 5.

</details>  


### Alpha diversity

A comprehensive list of global indicators of the ecosystem state can be obtained as follows. This includes various measures of richness, evenness, diversity, dominance, and rarity with default parameters. See the individual functions for more options regarding parameter tuning.


```r
# Aggregate rare taxa at the Species level for the phyloseq object 'merged_metagenomes'
# Taxa are included if they have a detection threshold of 0.05% and a prevalence of 20%
pseq <- aggregate_rare(merged_metagenomes, level = "Species", detection = 0.05/100, prevalence = 20/100)

# Calculate the Shannon diversity index for the aggregated phyloseq object 'pseq'
# The Shannon index measures the diversity within a sample, considering both richness and evenness
tab <- microbiome::alpha(pseq, index = 'shannon')
```

```
## Observed richness
```

```
## Other forms of richness
```

```
## Diversity
```

```
## Evenness
```

```
## Dominance
```

```
## Rarity
```

```r
# Display the first few rows of the calculated Shannon diversity index table
head(tab)
```

```
##          diversity_shannon
## T1rangPM         1.0597104
## T1rangAM         0.8341934
## T0ContPM         3.3703144
## T0ContAM         0.2441067
## T1cruzPM         0.3210316
## T1cruzAM         0.7733807
```
This returns observed richness with given detection threshold(s).


```r
# Calculate the richness of the phyloseq object 'pseq' with a detection threshold of 1000
# Richness is a measure of the number of different taxa present in a sample
tab <- richness(pseq, detection = 1000)

# Display the first few rows of the calculated richness table
head(tab)
```

```
##          1000 chao1
## T1rangPM    4   192
## T1rangAM   11   248
## T0ContPM   93   867
## T0ContAM   12   167
## T1cruzPM    5   312
## T1cruzAM   10   216
```


```r
p.shannon <- boxplot_alpha(pseq, 
                           index = "shannon",
                           x_var = "Type",
                          fill.colors = c(Control="cyan4", T._cruzi="deeppink4",T._rangeli="darkorange1" ))
p.shannon
```

![](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/shanon.png)

Alternative: 


```r
# get the metadata out as seprate object
hmp.meta <- meta(pseq)

# Add the rownames as a new colum for easy integration later.
hmp.meta$sam_name <- (hmp.meta$Type)

hmp.div <- microbiome::alpha(pseq)
# Add the rownames to diversity table
hmp.div$Sample <- rownames(hmp.div)

# merge these two data frames into one
div.df <- merge(hmp.div,hmp.meta, by = "Sample")

# check the tables
colnames(div.df)

div.df2 <- div.df[, c("sam_name", "diversity_inverse_simpson", "diversity_gini_simpson", "diversity_shannon", "diversity_fisher", "chao1")]

# the names are not pretty. we can replace them

colnames(div.df2) <- c("Type", "Inverse Simpson", "Gini-Simpson", "Shannon", "Fisher", "Chao1")

# check
colnames(div.df2)
div_df_melt <- reshape2::melt(div.df2)
## Using Location as id variables

head(div_df_melt)

library(ggpubr)
# Now use this data frame to plot 
p <- ggboxplot(div_df_melt, x = "Type", y = "value",
               fill = "Type", 
               palette = "jco", 
               legend= "right",
               facet.by = "variable", 
               scales = "free")

p <- p + rotate_x_text() 
# we will remove the x axis lables

p <- p + rremove("x.text")
p
```

![](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/indices.png)

Each of these metrics can give an insight into the distribution of the OTUs inside our samples. For example, the Chao1 diversity index gives more weight to singletons and doubletons observed in our samples, while Shannon is an entropy index remarking the impossibility of taking two reads out of the metagenome “bag” and that these two will belong to the same OTU.

**Q: What do you observe?**

####  Normality test: Check the Normal or not normal distribution  


Check normality of data: Shapiro Test & QQ-plots. 
Shapiro: H0 is «data follow normal distribution», H1 is «data do not follow normal distribution». 

![](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/shap1.png)

Use the custom function indices_normality() (defined in HELPER/indices_normality.R) plots the results of Shapiro test (Theoritical Quantiles) as well as Q-Qplots.


```r
tab <- microbiome::alpha(pseq)


tab |>
  dplyr::select(observed,
                diversity_gini_simpson,
                diversity_shannon,
                chao1) |>
  indices_normality(nrow = 2, ncol = 2)
```

![](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/shap.png)



***Q: What are your conclusions?***




### Investigate the top factors

Show coefficients for the top taxa separating the groups
 

```r
pseq <- aggregate_rare(merged_metagenomes, level = "Genus", detection = 0.1/100, prevalence = 50/100)
pseq <- microbiome::transform(pseq, transform = "compositional")
head(abundances(pseq))

p <- plot_landscape(pseq, method = "NMDS", distance = "bray", col = "Time", size = 3)
```
![association](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/pca.png)


### Estimating associations with an external variable

Next to visualizing whether any variable is associated with differences between samples, we can also quantify the strength of the association between community composition (beta diversity) and external factors.

Permutational Analysis of Variance (PERMANOVA; (2001)) is a widely used non-parametric multivariate method that aims to estimate the actual statistical significance of differences in the observed community composition between two groups of samples. This method takes as input the abundance table, which measure of distance you want to base the test on and a formula that tells the model how you think the variables are associated with each other.

PERMANOVA tests the hypothesis that the centroids and dispersion of the community are equivalent between the compared groups. A p-value smaller than the significance threshold indicates that the groups have a different community composition. This method is implemented with the adonis2 function from the vegan package.


```r
otu <- abundances(pseq)
meta <- meta(pseq)

permanova <- adonis(t(otu) ~ Time,
                    by = "margin",
                    data = meta, 
                    permutations = 9999, 
                    method = "euclidean")
```

```
## 'adonis' will be deprecated: use 'adonis2' instead
```

P-value: 


```r
print(as.data.frame(permanova$aov.tab))
```

```
##           Df    SumsOfSqs      MeanSqs  F.Model        R2 Pr(>F)
## Time       4 2.011394e+13 5.028484e+12 4.444902 0.4156094 0.0015
## Residuals 25 2.828231e+13 1.131292e+12       NA 0.5843906     NA
## Total     29 4.839624e+13           NA       NA 1.0000000     NA
```

The time variable is significantly associated with microbiota composition (p-value is below 0.05).

Let us visualize the model coefficients for species that exhibit the largest differences between the groups. This gives some insights into how the groups tend to differ from each other in terms of community composition.



```r
coef <- coefficients(permanova)["Time1",]
top.coef <- coef[rev(order(abs(as.numeric(coef))))[1:20]]
par(mar = c(3, 14, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")
```

![association](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/associations.png)

**Enterococcus** (Bacillota), which plays a crucial role in metabolic adaptability against pathogenic or plant toxins and anti-herbivore defense, was found to be one of the predominant gut microorganism of lepidopteran insects, including B. mori, Helicoverpa zea, and Porthetria dispar (Paniagua Voirol et al., 2018; Zhang et al., 2022). 

**Symbiopectobacterium** (Enterobacteriaceae) has recently been described for the first time as an intracellular bacterial symbiont, responsible for to the biosynthesis of vitamins and cofactors.  This bacteria may be boosting the parsite fitness, for example, by aiding in evading the Triatome immune response, or providing a novel function, such as supplementing nutrition or metabolism

**Rhodococcus** (Nocardiaceae) in the triatomine gut are believed to play an important role in the metabolism of the vector, such as by participating in the synthesis of group B vitamins or by being digested by the bugs directly to provide missing nutrients (Sassera et al., 2013). Moreover, the most attractive aspect is the host-symbiont relationship between triatomines and Rhodococcus; since Rhodococcus bacteria can be easily cultured and genetically modified to harm the pathogen in vector gut, they are probably suitable tools for the control of trypanosomiasis (Sassera et al., 2013). as the blood is poor in B vitamins compared to what is generally required for insect development. The blood is poor in B vitamins compared to what is generally required for insect development, Kissing bugs, Rhodnius prolixus, notably require *Rhodococcus* bacteria for nymph development, but the addition of B vitamins in the diet can rescue nymph development in the absence of Rhodococcus  (Serrato-Salas and Gendrin 2023).

**Wolbachia** (Ehrlichiaceae)  The obligate intracellular bacteria Wolbachia spp. are common in a wide range of insects, including sand flies, bed bugs, fleas and mosquitoes, and can cause reproduction alterations such as feminization, male killing and cytoplasmic incompatibility ([Landmann 20219](https://doi.org/10.1128/microbiolspec.BAI-0018-2019.)). In triatomines, Wolbachia has been solely reported for the genus Rhodnius, where it occurs in the intestine, salivary glands and gonads.

**Curtobacterium** (Microbacteriaceae) *C. flaccumfaciens* is the only species of Curtobacterium associated with plant pathogenesis (Young et al., 1996), the presence of *C. flaccumfaciens* in the rhizosphere induced a systematic resistance in cucumber plants to pathogens.


**Exercises** 

Community-level comparisons: Use PERMANOVA to investigate whether the community composition differs between two groups of individuals (e.g. times, or some other grouping of your choice). You can also include covariates such as type, gut, and see how this affects the results?


### Beta diversity

Beta diversity quantifies the dissimilarity between communities (multiple samples), as opposed to alpha diversity which focuses on variation within a community (one sample). In microbiome research, commonly used metrics of beta diversity include the Bray-Curtis index (for compositional data), Jaccard index (for presence/absence data, ignoring abundance information), Aitchison distance (Euclidean distance for clr transformed abundances, aiming to avoid the compositionality bias), and the Unifrac distance (that takes into account the phylogenetic tree information). Notably, only some of these measures are actual distances, as this is a mathematical concept whose definition is not satisfied by certain ecological measure, such as the Bray-Curtis index. Therefore, the terms dissimilarity and beta diversity are preferred.

| Method description             | Assay type          | Beta diversity metric |
|--------------------------------|---------------------|-----------------------|
| Quantitative profiling         | Absolute counts     | Bray-Curtis           |
| Relative profiling             | Relative abundances | Bray-Curtis           |
| Aitchison distance             | Absolute counts     | Aitchison             |
| Aitchison distance             | clr                 | Euclidean             |
| Robust Aitchison distance      | rclr                | Euclidean             |
| Presence/Absence similarity    | Relative abundances | Jaccard               |
| Presence/Absence similarity    | Absolute counts     | Jaccard               |
| Phylogenetic distance          | Rarefied counts     | Unifrac               |

In practice, beta diversity is usually represented as a *dist* object, a triangular matrix where the distance between each pair of samples is encoded by a specific cell. This distance matrix can then undergo ordination, which is an important ecological tool to reduce the dimensionality of data for a more efficient analysis and visualization. Ordination techniques aim to capture as much essential information from the data as possible and turn it into a lower dimensional representation. Dimension reduction is bound to lose information but commonly used ordination techniques can preserve relevant information of sample similarities in an optimal way, which is defined in different ways by different methods.

Based on the type of algorithm, ordination methods in microbiome research can be generally divided in two categories: unsupervised and supervised ordination. The former includes Principal Coordinate Analysis (PCoA), Principal Component Analysis (PCA) and Uniform Manifold Approximation and Projection for Dimension Reduction (UMAP), whereas the latter is mainly represented by distance-based Redundancy Analysis (dbRDA). We will first discuss unsupervised ordination methods and then proceed to supervised ones.


#### Unsupervised ordination

Unsupervised ordination methods variation in the data without additional information on covariates or other supervision of the model. Among the different approaches, Multi-Dimensional Scaling (MDS) and non-metric MDS (NMDS) can be regarded as the standard. They are jointly referred to as PCoA. 

A typical comparison of community compositions starts with a visual representation of the groups by a 2D ordination. Then we estimate relative abundances and MDS ordination based on Bray-Curtis index between the groups, and visualize the results.


```r
# To ensure reproducibility we can fix the seed here. This will ensure you always get the same result each time you run your data.
set.seed(34521)

pseq <- aggregate_rare(merged_metagenomes, level = "Genus", detection = 0.1/100, prevalence = 50/100)
pseq <- microbiome::transform(pseq, transform = "compositional")
# abundances(pseq)
pseq  <- abundances(pseq) %>% as.data.frame() 
pseq <- pseq <- pseq[order(rowSums(pseq), decreasing = TRUE), ] %>% t()

# Calculate distance matrix
species_frac_filtered_dist <- vegdist(pseq, method = "bray")

# Perform NMDS on distance matrix
nmds_spec <- metaMDS(species_frac_filtered_dist, distance = "bray",k = 2)
```

```
## Run 0 stress 0.03362996 
## Run 1 stress 0.0336299 
## ... New best solution
## ... Procrustes: rmse 0.0001257671  max resid 0.0002233539 
## ... Similar to previous best
## Run 2 stress 0.0339897 
## ... Procrustes: rmse 0.01002435  max resid 0.04134134 
## Run 3 stress 0.03362999 
## ... Procrustes: rmse 5.825382e-05  max resid 0.0001048328 
## ... Similar to previous best
## Run 4 stress 0.06288395 
## Run 5 stress 0.04144533 
## Run 6 stress 0.04452223 
## Run 7 stress 0.05896259 
## Run 8 stress 0.03398968 
## ... Procrustes: rmse 0.009997987  max resid 0.04149452 
## Run 9 stress 0.05778166 
## Run 10 stress 0.05752581 
## Run 11 stress 0.03362988 
## ... New best solution
## ... Procrustes: rmse 3.457753e-05  max resid 6.164393e-05 
## ... Similar to previous best
## Run 12 stress 0.03398975 
## ... Procrustes: rmse 0.01002627  max resid 0.04135569 
## Run 13 stress 0.03398965 
## ... Procrustes: rmse 0.01000719  max resid 0.04146811 
## Run 14 stress 0.05573411 
## Run 15 stress 0.05535289 
## Run 16 stress 0.05817008 
## Run 17 stress 0.05810237 
## Run 18 stress 0.0336299 
## ... Procrustes: rmse 2.72224e-05  max resid 4.87934e-05 
## ... Similar to previous best
## Run 19 stress 0.0614943 
## Run 20 stress 0.06068536 
## *** Best solution repeated 2 times
```

Check the output. 


```r
# Check the output
nmds_spec
```

```
## 
## Call:
## metaMDS(comm = species_frac_filtered_dist, distance = "bray",      k = 2) 
## 
## global Multidimensional Scaling using monoMDS
## 
## Data:     species_frac_filtered_dist 
## Distance: bray 
## 
## Dimensions: 2 
## Stress:     0.03362988 
## Stress type 1, weak ties
## Best solution was repeated 2 times in 20 tries
## The best solution was from try 11 (random start)
## Scaling: centring, PC rotation, halfchange scaling 
## Species: scores missing
```
Here you see a kind of summary of the analysis. For example, you can see that you used 2 dimensions and the stress was approx. 0.03. In general if a stress is above 0.2 then the clustering is not reliably representing the data and should be interpreted with caution. But here the stress is below 0.2, so we are okay.

Now let’s look at the ordination. To plot the data with ggplot, we need to extract the coordinaties of each point from nmds_spec$points.


```r
# Extract and reshape the data to plot ordination as ggplot  and add the metadata
# Convert the NMDS points to a data frame
nmds_spec_gg <- as.data.frame(nmds_spec$points)

# Add the "Sample" column based on row names
nmds_spec_gg <- nmds_spec_gg %>% rownames_to_column("Sample")

# Merge the NMDS data with the metadata by the "Sample" column
merged_data <- dplyr::left_join(meta, nmds_spec_gg, by = "Sample")
```
Then we can create the plot easily and color according to the metadata. We are choosing timepoint and mocktreat for the coloring respectively. But feel free to explore other parameters.


```r
ggplot(merged_data, aes(x = MDS1, y = MDS2)) +
    geom_point(aes(color = Time), size = 3, alpha = 0.5) +
    geom_text_repel(aes(label = Type), size = 3, color = "black") +
    ggtitle("NMDS colored according to Time") +
    theme_minimal()
```
![nmds_time](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/nmds_time.png)


```r
ggplot(merged_data, aes(x = MDS1, y = MDS2)) +
    geom_point(aes(color = Reads), size = 3, alpha = 0.5) +
    geom_text_repel(aes(label = Time), size = 3, color = "black") +
    scale_color_continuous(name = "Reads") +  # Add a continuous color scale
    ggtitle("NMDS colored according to Reads numbers") +
    theme_minimal()
```

![nmds_type](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/nmds_reads.png)

A few combinations of beta diversity metrics and assay types are typically used. For instance, Bray-Curtis dissimilarity and Euclidean distance are often applied to the relative abundance and the clr assays, respectively. Besides beta diversity metric and assay type, the PCoA algorithm is also a variable that should be considered. Below, we show how the choice of these three factors can affect the resulting lower-dimensional data.


```r
# Run NMDS on relabundance assay with Bray-Curtis distances
pseq <- aggregate_rare(merged_metagenomes, level = "Genus", detection = 0.1/100, prevalence = 50/100)
pseq <- microbiome::transform(pseq, transform = "compositional")
# abundances(pseq)
pseq  <- abundances(pseq) %>% as.data.frame() 
pseq <- pseq <- pseq[order(rowSums(pseq), decreasing = TRUE), ] %>% t()
# Calculate distance matrix
species_frac_filtered_dist_bray <- vegdist(pseq, method = "bray")
# Perform NMDS on distance matrix
nmds_spec_comp_bray <- metaMDS(species_frac_filtered_dist_bray,distance = "bray",k = 2)
```

```
## Run 0 stress 0.03362996 
## Run 1 stress 0.03362994 
## ... New best solution
## ... Procrustes: rmse 2.517756e-05  max resid 7.99797e-05 
## ... Similar to previous best
## Run 2 stress 0.03398969 
## ... Procrustes: rmse 0.009980697  max resid 0.04169512 
## Run 3 stress 0.04150819 
## Run 4 stress 0.04150804 
## Run 5 stress 0.05924 
## Run 6 stress 0.04144544 
## Run 7 stress 0.0446174 
## Run 8 stress 0.04144544 
## Run 9 stress 0.03362999 
## ... Procrustes: rmse 2.613272e-05  max resid 4.637768e-05 
## ... Similar to previous best
## Run 10 stress 0.05831549 
## Run 11 stress 0.03362988 
## ... New best solution
## ... Procrustes: rmse 8.202762e-05  max resid 0.0001496796 
## ... Similar to previous best
## Run 12 stress 0.05817021 
## Run 13 stress 0.04144542 
## Run 14 stress 0.05748917 
## Run 15 stress 0.03398972 
## ... Procrustes: rmse 0.0100232  max resid 0.04137259 
## Run 16 stress 0.05810249 
## Run 17 stress 0.04144545 
## Run 18 stress 0.05778144 
## Run 19 stress 0.04144542 
## Run 20 stress 0.05817676 
## *** Best solution repeated 1 times
```

```r
# Run NMDS on compositional assay with Euclidean distances
# Calculate distance matrix
species_frac_filtered_dist_euclidean <- vegdist(pseq, method = "euclidean")
# Perform NMDS on distance matrix
nmds_spec_comp_euclidean <- metaMDS(species_frac_filtered_dist_euclidean,distance = "euclidean",k = 2)
```

```
## Run 0 stress 0.007739837 
## Run 1 stress 0.02284928 
## Run 2 stress 0.02292606 
## Run 3 stress 0.007740427 
## ... Procrustes: rmse 0.002822558  max resid 0.01294398 
## Run 4 stress 0.01336195 
## Run 5 stress 0.007740287 
## ... Procrustes: rmse 0.002826063  max resid 0.01292524 
## Run 6 stress 0.00819932 
## ... Procrustes: rmse 0.009693796  max resid 0.03835989 
## Run 7 stress 0.008282441 
## Run 8 stress 0.007740289 
## ... Procrustes: rmse 0.002826031  max resid 0.01292548 
## Run 9 stress 0.007740324 
## ... Procrustes: rmse 0.002824945  max resid 0.01293063 
## Run 10 stress 0.02293236 
## Run 11 stress 0.008199333 
## ... Procrustes: rmse 0.009694176  max resid 0.03834877 
## Run 12 stress 0.008273883 
## Run 13 stress 0.022853 
## Run 14 stress 0.01336205 
## Run 15 stress 0.02101797 
## Run 16 stress 0.01336186 
## Run 17 stress 0.01327813 
## Run 18 stress 0.008273854 
## Run 19 stress 0.00774005 
## ... Procrustes: rmse 7.269943e-05  max resid 0.0001752729 
## ... Similar to previous best
## Run 20 stress 0.01336166 
## *** Best solution repeated 1 times
```

```r
# Run NMDS on compositional assay with Aitchison distances
# Calculate distance matrix
species_frac_filtered_dist_aitchison <- vegdist(pseq, method = "robust.aitchison")
# Perform NMDS on distance matrix
nmds_spec_comp_aitchison <- metaMDS(species_frac_filtered_dist_aitchison,distance = "robust.aitchison",k = 2)
```

```
## Run 0 stress 0.1172405 
## Run 1 stress 0.1465789 
## Run 2 stress 0.1513655 
## Run 3 stress 0.1172406 
## ... Procrustes: rmse 0.0001538783  max resid 0.000660751 
## ... Similar to previous best
## Run 4 stress 0.1131166 
## ... New best solution
## ... Procrustes: rmse 0.02751694  max resid 0.1211573 
## Run 5 stress 0.1172406 
## Run 6 stress 0.1131165 
## ... New best solution
## ... Procrustes: rmse 3.099105e-05  max resid 0.0001323215 
## ... Similar to previous best
## Run 7 stress 0.1644752 
## Run 8 stress 0.1131167 
## ... Procrustes: rmse 4.010913e-05  max resid 0.0001276984 
## ... Similar to previous best
## Run 9 stress 0.1172405 
## Run 10 stress 0.1923506 
## Run 11 stress 0.1550482 
## Run 12 stress 0.113127 
## ... Procrustes: rmse 0.004699189  max resid 0.01863022 
## Run 13 stress 0.1172407 
## Run 14 stress 0.1131166 
## ... Procrustes: rmse 1.209206e-05  max resid 3.604119e-05 
## ... Similar to previous best
## Run 15 stress 0.1131166 
## ... Procrustes: rmse 7.759492e-05  max resid 0.0002872779 
## ... Similar to previous best
## Run 16 stress 0.1552812 
## Run 17 stress 0.1131165 
## ... Procrustes: rmse 1.068509e-05  max resid 4.409835e-05 
## ... Similar to previous best
## Run 18 stress 0.1898461 
## Run 19 stress 0.1131166 
## ... Procrustes: rmse 3.188427e-05  max resid 8.271452e-05 
## ... Similar to previous best
## Run 20 stress 0.1172406 
## *** Best solution repeated 6 times
```

```r
# Run NMDS on clr assay with Euclidean distances
pseq <- aggregate_rare(merged_metagenomes, level = "Genus", detection = 0.1/100, prevalence = 50/100)
pseq <- microbiome::transform(pseq, transform = "clr")
pseq  <- abundances(pseq) %>% as.data.frame() 
pseq <- pseq <- pseq[order(rowSums(pseq), decreasing = TRUE), ] %>% t()
# Calculate distance matrix
species_frac_filtered_dist_euclidean <- vegdist(pseq, method = "euclidean")
# Perform NMDS on distance matrix
nmds_spec_clr_euclidean <- metaMDS(species_frac_filtered_dist_euclidean,distance = "euclidean",k = 2)
```

```
## Run 0 stress 0.09713366 
## Run 1 stress 0.1160723 
## Run 2 stress 0.1097908 
## Run 3 stress 0.1352095 
## Run 4 stress 0.1038031 
## Run 5 stress 0.09713361 
## ... New best solution
## ... Procrustes: rmse 8.804707e-05  max resid 0.0003260719 
## ... Similar to previous best
## Run 6 stress 0.09713361 
## ... New best solution
## ... Procrustes: rmse 7.101807e-06  max resid 2.043321e-05 
## ... Similar to previous best
## Run 7 stress 0.1012603 
## Run 8 stress 0.09713362 
## ... Procrustes: rmse 3.842243e-05  max resid 0.0001411418 
## ... Similar to previous best
## Run 9 stress 0.1012603 
## Run 10 stress 0.1012603 
## Run 11 stress 0.09896957 
## Run 12 stress 0.09896957 
## Run 13 stress 0.1366336 
## Run 14 stress 0.09896957 
## Run 15 stress 0.09896957 
## Run 16 stress 0.1198454 
## Run 17 stress 0.1029259 
## Run 18 stress 0.09896957 
## Run 19 stress 0.1306306 
## Run 20 stress 0.09896957 
## *** Best solution repeated 2 times
```

```r
# List of NMDS objects and their corresponding titles
nmds_list <- list(
  list(data = nmds_spec_comp_bray, title = "Comp Bray"),
  list(data = nmds_spec_comp_euclidean, title = "Comp Euclidean"),
  list(data = nmds_spec_comp_aitchison, title = "Comp Aitchison"),
  list(data = nmds_spec_clr_euclidean, title = "Clr Euclidean")
)

# Initialize an empty list to store plots
plot_list <- list()

# Loop through each NMDS object and generate the corresponding plot
for (nmds_item in nmds_list) {
  nmds_spec_gg <- as.data.frame(nmds_item$data$points) %>%
    rownames_to_column("Sample") %>%
    dplyr::left_join(meta, by = "Sample")
  
  plot <- ggplot(nmds_spec_gg, aes(x = MDS1, y = MDS2)) +
    geom_point(aes(color = Time), size = 3, alpha = 0.5) +
    geom_text_repel(aes(label = Type), size = 3, color = "black") +
    ggtitle(paste("NMDS colored according to Time -", nmds_item$title)) +
    theme_minimal()
  
  # Add the plot to the plot list
  plot_list[[nmds_item$title]] <- plot
}

# Combine all plots using patchwork
combined_plot <- wrap_plots(plot_list) +
  plot_layout(guides = "collect")

# Print the combined plot
#print(combined_plot)
```

![nmds_type](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/nmds_all.png)


#### Unsupervised ordination

dbRDA is a supervised counterpart of PCoA. It maximize the variance with respect to the covariates provided by the user. This can be used to quantify associations between each covariate and community composition (beta diversity). The table below summarizes the relations between the supervised and unsupervised ordination methods.

| Supervised ordination  | Unsupervised ordination          |
|------------------------|----------------------------------|
| Euclidean distance     | RDA                              | PCA                         |
| Non-Euclidean distance | dbRDA                            | PCoA/MDS, NMDS, UMAP        |


In summary, the dbRDA is the more general method that allows a wider variety dissimilarity, or beta diversity, indices. 

Let's use the package mia for this analysis, fist we transform our phyloseq object into a mia compatible object. The colData lists the covariates such as Gut and Type...


```r
pseq <-  mia::makeTreeSummarizedExperimentFromPhyloseq(merged_metagenomes)

# Apply relative transform
pseq <- mia::transformAssay(pseq,
                       method = "relabundance")
```

dbRDA can be perfomed with the runRDA function. In addition to the arguments previously defined for unsupervised ordination, this function takes a formula to control for variables and an action to treat missing values. Along with metadata, which is the main outcome,we can treat observations missing values (not the case here). 


```r
pseq <- mia::runRDA(pseq,
               assay.type = "relabundance",
               formula = assay ~ Gut + Time + Type + Reads,
               distance = "bray",
               na.action = na.exclude)
```
The importance of each variable on the similarity between samples can be assessed from the results of PERMANOVA, automatically provided by the runRDA function. We see that Time explain more than 42% of the variance and it is also significant.


```r
# Store results of PERMANOVA test
rda_info <- attr(reducedDim(pseq, "RDA"), "significance")
rda_info$permanova
```

```
##          Df   SumOfSqs          F Pr(>F) Total variance Explained variance
## Model     8 2.98947546 13.0670703  0.001       3.590021         0.83271806
## Gut       1 0.03886857  1.3591639  0.230       3.590021         0.01082684
## Time      4 1.51258308 13.2230752  0.001       3.590021         0.42132985
## Type      2 0.05388041  0.9420504  0.410       3.590021         0.01500838
## Reads     1 0.05231320  1.8292984  0.156       3.590021         0.01457184
## Residual 21 0.60054571         NA     NA       3.590021         0.16728194
```

|          | Df | SumOfSqs |         F |  Pr(>F) | Total variance | Explained variance |
|----------|----|----------|-----------|---------|----------------|-------------------|
| Model    |  8 | 2.989475 | 13.067070 |   0.001 |       3.590021 |         0.832718  |
| Gut      |  1 | 0.038869 |  1.359164 |   0.237 |       3.590021 |         0.010827  |
| Time     |  4 | 1.512583 | 13.223075 |   0.001 |       3.590021 |         0.421330  |
| Type     |  2 | 0.053880 |  0.942050 |   0.405 |       3.590021 |         0.015008  |
| Reads    |  1 | 0.052313 |  1.829298 |   0.163 |       3.590021 |         0.014572  |
| Residual | 21 | 0.600546 |     NA    |     NA  |       3.590021 |         0.167282  |


Next, we proceed to visualize the weight and significance of each variable on the similarity between samples with an RDA plot, which can be generated with the plotRDA function from the miaViz package.


```r
# Load packages for plotting function
library(miaViz)

# Generate RDA plot coloured by clinical status
plotRDA(pseq, "RDA", colour_by = "Time")
```

![nmds_type](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/dbrda_plot.png)


### Confounding effects
Confounders can be defined as variables that are related to and affect the apparent dynamics between the response and the main independent variable. They are common in experimental studies. Generally, they can be classified into 3 groups:

- Biological confounders, such as age and sex

- Technical confounders produced during sample collection, processing and analysis

- Confounders resulting from experimental models, such as batch effects and sample history

Controlling for confounders is an important practice to reach an unbiased conclusion. To perform causal inference, it is crucial that the method is able to include confounders in the model. This is not possible with statistical tests of general use, such as the Wilcoxon test.
