# Day 2

| Time      | Activity                      | Slides                                | Hands-on                                                |
|-----------|-------------------------------|---------------------------------------|---------------------------------------------------------|
| Morning   | Read Base analysis            |                                       |       [Link here](Readbased.md)                            | 
| Afternoon | Assembly                      | [Link here](linkhere.pdf)                | [Link here](Day2/Day2_assembly.md) |




# Metagenomic read based profiles

In this tutorial we will explore and analyse the results of kraken2 read profiling data using the re-estimatd abundances of bracken. These exercises are meant to show how to conceptually approach your data analysis but there are many more and different ways to explore your data. The most important thing to keep in mind is that you have to understand your own data and analyses. One way to achieve this is to perform visual explorations that help you to judge whether the data are appropriate for your question.

Now let’s start the fun!



### Load libraries

```R
library(tidyverse)
library(vegan)
library(coin)
library(pheatmap)
```


##### 1. load the data

First we need to load our data. Usually the biggest bottleneck between raw data and analyses is to get the data in the right shape for your purpose. Often this requires a little bit of data mingling. On this road - google is your best friend to master the R universe :)

Let’s first load the relative abundance table of the bracken results.

```R
merge_rel_abund <- read_csv(file = "READBASED/merged_rel_abund.csv") 

meta <- read_csv(file = "DATA/tryp_metadata.csv")
```


#### 2. Basic stats

Before we start anything, let’s just check out or data a little bit. Never go blind into your analyses.

Q1: How many samples do we have in the matadata and in out data tibble


Q2: How many families were detected in our dataset?


#### 3. Pre-process data

Now that we know a little bit about our data we can start selecting what we want to look at. Let’s extract a separate tibble for species, family and phylum level. At each tibble we remove all the taxa that are not present in any of the samples.

Let’s start with species:

```R
merge_rel_abund_spe <- subset(brel, taxlevel == "S") %>% 
  select(-c("taxid","taxlevel")) %>%
  filter(rowSums(select_if(., is.numeric)) != 0)
head(brel_spec)
```

Can you come up with the command for family (F) and phylum (P) level? Tip: make sure you save it to a different variable. Otherwise you are overwriting your species data.



#### 4. Plot relative abundances
Let’s plot and explore microbial composition and relative abundances in our samples by visualizing them.

As a first overview of coarse differences we can create a stacked bar plot of phyla. We need to mingle our data structure a bit to make it compliant with ggplot and add the metadata to get an extra level of information.

```r
# data mingling
merge_rel_abund_gg1 <- merge_rel_abund_phy %>% gather(key="Sample",value="rel_abun",-taxa)
# now join with metadata by column 'Sample'. We are using left join in case the metadata file contains additional samples not included in our dataset
merge_rel_abund_gg1 <- left_join(merge_rel_abund_gg1,meta,by="Sample") 
# check how it looks
head(brel_phy_gg1) 
```


```r
# bar plot phyla relative abundances
ggplot(merge_rel_abund_gg1, aes(x=Sample, y=rel_abun)) +
  geom_bar(stat="identity",position="stack", aes(fill=taxa)) + # chose bar plot
  theme(axis.text.x = element_text(angle=45, hjust=1)) + # put x-axis label at 45 degree angle
  facet_grid(. ~ mocktreat, scales="free_x",space = "free_x") # produce two panels according to metatadata category 'mocktreat'
```
Q3: Do the phyla profiles look similar between samples? Can you spot any trends?


When we want to look at species composition, the communities are often very complex and looking at abundance profiles is not as informative. One way to explore species-level community composition is to filter you data so you focus only on the most abundant taxa across samples. What cutoff and approach you use here of course depends on your question and data.

Let’s create a heatmap of the 20 most abundant spexies in our samples.

```r
# sort species by abundance across samples and select top 20
merge_rel_abund_spe_gg1<-brel_spec[order(rowSums(select_if(merge_rel_abund_spe, is.numeric)),decreasing=T),] %>%
  head(20) %>%
  column_to_rownames("taxa")

# shape the metadata
meta_h <- subset(meta, Sample %in% colnames(merge_rel_abund_spe_gg1)) %>% 
  column_to_rownames("Sample") # shoft the 'Sample' column to rownames

# plot the heatmap
pheatmap::pheatmap(brel_spec_gg1,
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   annotation_col = meta_h[,c(1,3,4)],
                   annotation_names_col=TRUE)
```

#### 5. Beta diversity


Often we want to know whether the microbiomes are different between conditions or groups. One way to explore this is to look at the beta-diversity in an ordination. There are different distances and approaches that can be done and explored. We will perform an NMDS on bray curtis dissimilarities of the species profiles.

```r
# To ensure reproducibility we can fix the seed here. This will ensure you always get the same result each time you run your data.
set.seed(34521)

# Data mingling
merge_rel_abund_gg2 <- merge_rel_abund_spe %>% 
  column_to_rownames("taxa") %>% 
  t() # transpose

# Calculate distance matrix
merge_rel_abund_gg2_dist <- vegdist(merge_rel_abund_gg2, method = "bray")

# Perform NMDS on distance matrix
nmds_spec <- metaMDS(merge_rel_abund_gg2_dist,distance = "bray",k = 2)
```
#### 6. Differential abundance
