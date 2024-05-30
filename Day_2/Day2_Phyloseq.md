Processing data and stats
================
Vincent Manzanilla, INTERTRYP

# Day 2

| Time    | Activity           | Slides | Hands-on                                                                                              |
|---------|--------------------|--------|-------------------------------------------------------------------------------------------------------|
| Morning | Read Base analysis |        | [Link here](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/Day2_READBASED_stats.md) |

# Metagenomic read based profiles

In this tutorial we will explore and analyse the results of kraken2 read
profiling data using the re-estimatd abundances of bracken. These
exercises are meant to show how to conceptually approach your data
analysis but there are many more and different ways to explore your
data. The most important thing to keep in mind is that you have to
understand your own data and analyses. One way to achieve this is to
perform visual explorations that help you to judge whether the data are
appropriate for your question.

Now let’s start the fun!

In R studio.

Load libraries:

**Data Manipulation and Wrangling**

``` r
# dplyr: A grammar of data manipulation, providing a consistent set of verbs to solve common data manipulation challenges
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
# tibble: Provides a modern reimagining of data frames, making them more user-friendly
library(tibble)

# tidyverse: A collection of R packages designed for data science, all sharing an underlying design philosophy, grammar, and data structures
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ forcats   1.0.0     ✔ readr     2.1.5
    ## ✔ ggplot2   3.5.0     ✔ stringr   1.5.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

**Microbiome Analysis**

``` r
# microbiome: Provides tools for microbiome data analysis and visualization
library(microbiome)
```

    ## Loading required package: phyloseq

    ## 
    ## microbiome R package (microbiome.github.com)
    ##     
    ## 
    ## 
    ##  Copyright (C) 2011-2022 Leo Lahti, 
    ##     Sudarshan Shetty et al. <microbiome.github.io>

    ## 
    ## Attaching package: 'microbiome'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     alpha

    ## The following object is masked from 'package:base':
    ## 
    ##     transform

``` r
# microbiomeutilities: Extends functionalities of the microbiome package with additional utilities
library(microbiomeutilities)
```

    ## Warning: replacing previous import 'ggplot2::alpha' by 'microbiome::alpha' when
    ## loading 'microbiomeutilities'

    ## 
    ## Attaching package: 'microbiomeutilities'

    ## The following object is masked from 'package:microbiome':
    ## 
    ##     add_refseq

``` r
# vegan: Functions for ecological analysis, including ordination methods, diversity analysis and other functions for community ecologists
library(vegan)
```

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.6-4

    ## 
    ## Attaching package: 'vegan'

    ## The following object is masked from 'package:microbiome':
    ## 
    ##     diversity

``` r
# compositions: Provides methods for analyzing compositional data, which is common in fields such as microbiome research
library(compositions)
```

    ## Welcome to compositions, a package for compositional data analysis.
    ## Find an intro with "? compositions"

    ## 
    ## Attaching package: 'compositions'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     anova, cor, cov, dist, var

    ## The following object is masked from 'package:graphics':
    ## 
    ##     segments

    ## The following objects are masked from 'package:base':
    ## 
    ##     %*%, norm, scale, scale.default

``` r
# mia: Microbiome Analysis, provides functionalities for microbiome data analysis and visualization
library(mia)
```

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count

    ## 
    ## Attaching package: 'MatrixGenerics'

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

    ## Loading required package: GenomicRanges

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:compositions':
    ## 
    ##     normalize, var

    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     intersect, setdiff, union

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     second, second<-

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:microbiome':
    ## 
    ##     coverage

    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     distance

    ## The following object is masked from 'package:lubridate':
    ## 
    ##     %within%

    ## The following object is masked from 'package:purrr':
    ## 
    ##     reduce

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## Loading required package: GenomeInfoDb

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     sampleNames

    ## Loading required package: SingleCellExperiment

    ## Loading required package: TreeSummarizedExperiment

    ## Loading required package: Biostrings

    ## Loading required package: XVector

    ## 
    ## Attaching package: 'XVector'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     compact

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    ## Loading required package: MultiAssayExperiment

    ## 
    ## Attaching package: 'mia'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     full_join, inner_join, left_join, right_join

``` r
# miaViz: Visualization utilities for microbiome analysis, part of the mia package ecosystem
library(miaViz)
```

    ## Loading required package: ggraph

``` r
# phyloseq: Microbiome Analysis, provides functionalities for microbiome data analysis and visualization
library(phyloseq)
```

**Data Visualization**

``` r
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

First we need to load our data. Usually the biggest bottleneck between
raw data and analyses is to get the data in the right shape for your
purpose. Often this requires a little bit of data mingling. On this
road - google is your best friend to master the R universe :)

Let’s first load the relative abundance table of the bracken results.

For this part we are using the
[phyloseq](https://joey711.github.io/phyloseq/) package. The phyloseq
package is a tool to import, store, analyze, and graphically display
complex phylogenetic sequencing data that has already been clustered
into Operational Taxonomic Units (OTUs), especially when there is
associated sample data, phylogenetic tree, and/or taxonomic assignment
of the OTUs.

First st your working directory:

``` r
setwd("~/Documents/project/Metagenomics_2024/Metagenomics_2024/Day_2/")
```

And load your two files:

``` r
merged_metagenomes <- phyloseq::import_biom("../DATA//merge_species.biom") 
meta <- read.csv(file = "../DATA/tryp_metadata.csv", sep = ",")
```

#### phyloseq-ize Data

Any data already in an R session can be annoated/coerced to be
recognized by phyloseq’s functions and methods. This is important,
because there are lots of ways you might receive data related to a
microbiome project, and not all of these will come from a popular server
or workflow that is already supported by a phyloseq import function.

![](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/phyloseq.png)

- **otu_table** - Works on any numeric matrix. You must also specify if
  the species are rows or columns
- **tax_table** - Works on any character matrix. The rownames must match
  the OTU names (taxa_names) of the otu_table if you plan to combine it
  with a phyloseq-object.
- **sample_data** - Works on any data.frame. The rownames must match the
  sample names in the otu_table if you plan to combine them as a
  phyloseq-object

Let’s look at the phlyseq object.

``` r
merged_metagenomes
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 2251 taxa and 30 samples ]
    ## tax_table()   Taxonomy Table:    [ 2251 taxa by 7 taxonomic ranks ]

Summarize the data.

``` r
# Summarize the phyloseq object 'merged_metagenomes'
#microbiome::summarize_phyloseq(merged_metagenomes)

# Display the first few rows of the OTU (Operational Taxonomic Unit) table
head(otu_table(merged_metagenomes))
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

``` r
# Display the first few rows of the taxonomy table
head(tax_table(merged_metagenomes))
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

``` r
# Display the first few rows of the sample data associated with 'merged_metagenomes'
#head(sample_data(merged_metagenomes))

# Get the sample variables of the phyloseq object 'merged_metagenomes'
#sample_variables(merged_metagenomes)
```

We need to add the metadata to the phyloseq object.

## 2. format the data

``` r
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

Before we start anything, let’s just check out or data a little bit
(sanity check). Never go blind into your analyses.

``` r
head(psmelt(merged_metagenomes))
```

    ## Warning in psmelt(merged_metagenomes): The sample variables: 
    ## Sample
    ##  have been renamed to: 
    ## sample_Sample
    ## to avoid conflicts with special phyloseq plot attribute names.

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

</details>

**Q: How many species and samples are detected in our dataset?**

<details>
<summary>
HINT
</summary>

> ntaxa(merged_metagenomes), nsamples.

</details>

**Q: Look at the taxonomy is ther some problem with the taxonomy?**

<details>
<summary>
HINT
</summary>

> We have contaminant in the data. Which one?
> (head(psmelt(merged_metagenomes)))

</details>

**Q: what are the taxa with the more abundant reads?**

<details>
<summary>
HINT
</summary>

> head(otu_table(merged_metagenomes))

</details>

Extraction of the sample’s names.

``` r
sample_names(merged_metagenomes)
```

    ##  [1] "T1rangPM" "T1rangAM" "T0ContPM" "T0ContAM" "T1cruzPM" "T1cruzAM"
    ##  [7] "T1ContAM" "T1ContPM" "T0rangPM" "T0rangAM" "T0cruzPM" "T7rangAM"
    ## [13] "T7rangPM" "T7cruzPM" "T7cruzAM" "T7ContPM" "T7ContAM" "T3rangAM"
    ## [19] "T3rangPM" "T3cruzPM" "T3cruzAM" "T0cruzAM" "T3ContPM" "T3ContAM"
    ## [25] "T2rangPM" "T2rangAM" "T2cruzPM" "T2cruzAM" "T2ContPM" "T2ContAM"

### Aggregation

Microbial species can be called at multiple taxonomic resolutions. We
can easily agglomerate the data based on taxonomic ranks. Here, we
agglomerate the data at Family level.

``` r
# Aggregate rare taxa at the family level for the phyloseq object 'merged_metagenomes'
merged_metagenomes_family <- aggregate_rare(merged_metagenomes, level = "Genus", detection = 0/100, prevalence = 0/100)

# Display the dimensionality of the abundances of 'merged_metagenomes_family'
dim(abundances(merged_metagenomes_family))
```

    ## [1] 788  30

**Q: How many sample and tax do we have now?**

<details>
<summary>
HINT
</summary>

> There are XX taxa and 30 samples, meaning that there are 29 different
> Phylum level taxonomic groups. Looking at the rowData after
> agglomeration shows all Enterococcaceae are combined together, and all
> lower rank information is lost.

``` r
head(tax_table(merged_metagenomes_family))
```

    ## Taxonomy Table:     [6 taxa by 2 taxonomic ranks]:
    ##                        Genus                    unique                  
    ## Abyssalbus             "Abyssalbus"             "Abyssalbus"            
    ## Achromobacter          "Achromobacter"          "Achromobacter"         
    ## Acidipropionibacterium "Acidipropionibacterium" "Acidipropionibacterium"
    ## Acidothermus           "Acidothermus"           "Acidothermus"          
    ## Acidovorax             "Acidovorax"             "Acidovorax"            
    ## Acinetobacter          "Acinetobacter"          "Acinetobacter"

</details>

## 4. QC and Pre-process data

Now that we know a little bit about our data we can start the
pre-processing.

### Library size / read count

Let us check for distribution of number of sequences retained from the
Kraken/Bracken approach.

``` r
plot_read_distribution(merged_metagenomes, "Type", plot.type = "density")
```

![](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/density_reads.png)

*You can try to plot with diffrent metadata.*

Plotting the read count per sample

``` r
df <- psmelt(merged_metagenomes)  %>%  group_by(Sample, Time) %>%  
  summarise(sum_reads = sum(Abundance)) %>% arrange(sum_reads) 

ggplot(df) +
  geom_bar(aes(reorder(Sample, -sum_reads), sum_reads, fill=Time),
           col="red", alpha = .2, stat="identity") 
```

![](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/reads_count.png)

### Prevalence - Detection

Prevalence quantifies the frequency of samples where certain microbes
were detected (above a given detection threshold). The prevalence can be
given as sample size (N) or percentage (unit interval).

The population prevalence (frequency) at a 1% relative abundance
threshold (detection = 1/100 and count = false), can look like this.

``` r
head(prevalence(merged_metagenomes, detection=1/100, sort=TRUE, count=FALSE))
```

    ##  85471   5855   3218   3369 210225   3469 
    ##      1      1      1      1      1      1

The function arguments detection and count can also be used to access,
how many samples do pass a threshold for raw counts. Here, the
population prevalence (frequency) at the absolute abundance threshold
(count=true) at read count 1 (detection = 1) is accessed.

``` r
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

Each point corresponds to a different or unique taxon. The y-axis
represents the fraction of samples, these taxa are present. The low
prevalence suggests there is a low overlap across samples.

![](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/prevalence.png)

**Q: Which taxa is present in all samples and with a high abundance?**

<details>
<summary>
HINT
</summary>

> Homo sapiens

</details>

#### Contaminant sequences

Samples might be contaminated with exogenous sequences. We have observed
1 contaminants Homo sapiens.

``` r
#check 
head(psmelt(subset_taxa(merged_metagenomes, Genus == "Trypanosoma")))
```

    ## Warning in psmelt(subset_taxa(merged_metagenomes, Genus == "Trypanosoma")): The sample variables: 
    ## Sample
    ##  have been renamed to: 
    ## sample_Sample
    ## to avoid conflicts with special phyloseq plot attribute names.

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

``` r
head(psmelt(subset_taxa(merged_metagenomes, Genus == "Homo")))
```

    ## Warning in psmelt(subset_taxa(merged_metagenomes, Genus == "Homo")): The sample variables: 
    ## Sample
    ##  have been renamed to: 
    ## sample_Sample
    ## to avoid conflicts with special phyloseq plot attribute names.

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

``` r
#Keep only the kingdom of interest
merged_metagenomes <- subset_taxa(merged_metagenomes, Kingdom %in% c("Archaea", "Bacteria", "Fungi", "Viruses"))
```

Now recheck that your data are clean before continuing the analysis.

**Q: Which taxa is present in all samples and with a high abundance?**

<details>
<summary>
HINT
</summary>

> head(tax_table(merged_metagenomes)) or
> plot_taxa_prevalence(merged_metagenomes, level=“Phylum”, detection =
> 10000)

</details>

### Run rarefaction curves

Before normalization by sub-sampling, let’s have a look at rarefaction
curves, evaluate your sequencing effort and make decisions. Rarefaction
is a statistical technique employed to evaluate species richness based
on sampling results. The process involved subsampling reads without
replacement to a defined sequencing depth, thereby creating a
standardized library size across samples. Any sample with a total read
count less than the defined sequencing depth used to rarefy will be
discarded. What we want to see in those curves is that it reaches a
plateau meaning no new OTUs are discovered (going up on the Y axis) as
sequencing is deeper (the X axis). If this happens we can be pretty
confident that we sequenced all the OTUs that were present in the
sample.

``` r
install.packages("remotes")
remotes::install_github("gauravsk/ranacapa")
```

``` r
ranacapa::ggrare(pseq, step = 10, color = "Sample", se = FALSE) +
    geom_vline(xintercept = min(sample_sums(pseq)), color = "gray60")
```

![](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/rare.png)
***Q: What do you observe?***

### Data Normalization

Data transformations are common in (microbial) ecology (Legendre 2001)
and used to improve compatibility with assumptions related to specific
statistical methods, mitigate biases, enhance the comparability of
samples or features, or to obtain more interpretable values.

Examples include the logarithmic transformation, calculation of relative
abundances (percentages), and compositionality-aware transformations
such as the centered log-ratio transformation (clr).

Let us summarize some commonly used transformations in microbiome data
science; further details and benchmarkings available in the references:

- *relabundance* relative transformation; also known as total sum
  scaling (TSS) and compositional transformation. This converts counts
  into percentages (at the scale \[0, 1\]) that sum up to. Much of the
  currently available taxonomic abundance data from high-throughput
  assays (16S, metagenomic sequencing) is compositional by nature, even
  if the data is provided as counts (Gloor et al. 2017).

- *clr* Centered log ratio transformation (Aitchison 1986) is used to
  reduce data skewness and compositionality bias in relative abundances,
  while bringing the data to the logarithmic scale. This transformation
  is frequently applied in microbial ecology (Gloor et al. 2017).
  However, this transformation only applies to positive values. Usual
  solution is to add pseudocount, which adds another type of bias in the
  data. The robust clr transformation (‘rclr’) aims to circumvent the
  need to add a pseudocount. While the resulting values from these
  transformations are difficult interpret directly, this transformation
  may enhance comparability of relative differences between samples. It
  is part of a broader Aitchison family of transformations; the additive
  log ratio transformation (\`alr’) is also available. The robust clr
  (“rclr”) is similar to regular clr (see above) but allows data with
  zeroes and avoids the need to add pseudocount Martino et al. (2019).
  See library *compositions*.

- *pa* presence/absence transformation ignores abundances and only
  indicates whether the given feature is detected above the given
  threshold (default: 0). This simple transformation is relatively
  widely used in ecological research. It has shown good performance in
  microbiome-based classification performance (Giliberti et al. 2022,
  Karwowska2024).

- *z* Z transformation scales data to zero mean and unit variance; this
  us used to bring features (or samples) to more comparable levels in
  terms of mean and scale of the values. This can enhance visualization
  and interpretation of the data

- *log*, *log2*, *log10* Logarithmic transformations; used e.g. to
  reduce data skewness; with compositional data the clr (or rclr)
  transformation is often preferred.

- *hellinger* Hellinger transformation equals to the square root of
  relative abundances. This ecological transformation can be useful if
  we are interested in changes in relative abundances.

- *rank* Rank transformation replaces each value by its rank. Also see
  ‘rrank’ (relative rank transformation). This has use for instance in
  non-parametric statistics.

- Other available transformations include Chi square (‘chi.square’),
  Frequency transformation (‘frequency’), and Make margin sum of squares
  equal to one (‘normalize’)

The data contains read counts. We can convert these into relative
abundances and other formats. Compare abundance of a given taxonomic
group using the example data before and after the compositionality
transformation (with a cross-plot, for instance). You can also compare
the results to CLR-transformed data (see e.g. Gloor et al. 2017)

Have a look at the function *?microbiome::transform*.

``` r
microbiome::transform(merged_metagenomes, transform = "compositional")
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 2025 taxa and 30 samples ]
    ## sample_data() Sample Data:       [ 30 samples by 7 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 2025 taxa by 7 taxonomic ranks ]

In order to assess the effect of the different transformations, you can
use the *clr* or the *log10* for the downstream analysis.

## 5. Microbiome composition

Microbial abundances are typically ‘compositional’ (relative) in the
current microbiome profiling data sets. This is due to technical aspects
of the data generation process (see e.g. Gloor et al., 2017).

The next example calculates relative abundances as these are usually
easier to interpret than plain counts. For some statistical models we
need to transform the data into other formats as explained in above link
(and as we will see later).

*detection*: Detection threshold for absence/presence (percentage
reads). *prevalence*: Prevalence threshold (in \[0, 1\]) (presence
accross samples).

``` r
pseq <- aggregate_rare(merged_metagenomes, level = "Genus", detection = 0.05/100, prevalence = 20/100)
pseq <- microbiome::transform(pseq, transform = "compositional")
```

**Q: what is the meaning and effect of compositional?**

<details>
<summary>
HINT
</summary>

> abundances(pseq)

</details>

#### Visualization

**Bar plot** are useful to have a broad look at the data.

``` r
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

**Q: Does the profiles look similar between samples? Can you spot any
trends?**

<details>
<summary>
HINT
</summary>

> A: We see that the same one or two OTU dominates all samples but there
> is some variability between samples. We can observe some trend on the
> Entrococcus proption over time.

</details>

*If you have time you can also visualize the other taxonomic levels
(e.g. species) with the same approach. Try to come up with the code
yourself.*

``` r
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
