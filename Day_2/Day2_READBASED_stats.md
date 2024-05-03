# Day 2

| Time      | Activity                      | Slides                                | Hands-on                                                |
|-----------|-------------------------------|---------------------------------------|---------------------------------------------------------|
| Morning   | Read Base analysis            |                                       |       [Link here](Readbased.md)                            | 
| Afternoon | Assembly                      | [Link here](linkhere.pdf)                | [Link here](Day2/Day2_assembly.md) |




# Metagenomic read based profiles

In this tutorial we will explore and analyse the results of kraken2 read profiling data using the re-estimatd abundances of bracken. These exercises are meant to show how to conceptually approach your data analysis but there are many more and different ways to explore your data. The most important thing to keep in mind is that you have to understand your own data and analyses. One way to achieve this is to perform visual explorations that help you to judge whether the data are appropriate for your question.

Now let’s start the fun!

In R studio.

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

bracken_merged <- read.csv(file = "bracken_merged.csv", sep = "\t") 

# Select columns ending with ".bracken_num" for bracken_merged_reads
bracken_merged_reads <- bracken_merged[grep("\\.bracken_num$", names(bracken_merged))]

# Select columns ending with ".bracken_frac" for bracken_merged_frac
bracken_merged_frac <- bracken_merged[grep("\\.bracken_frac$", names(bracken_merged))]

# Remove the suffixes from the column names
colnames(bracken_merged_reads) <- sub("\\.bracken_num$", "", colnames(bracken_merged_reads))
colnames(bracken_merged_frac) <- sub("\\.bracken_frac$", "", colnames(bracken_merged_frac))

# Add name, taxonomy_id, and taxonomy_lvl to both data frames
bracken_merged_reads <- cbind(bracken_merged[, c("name", "taxonomy_id", "taxonomy_lvl")], bracken_merged_reads)
bracken_merged_frac <- cbind(bracken_merged[, c("name", "taxonomy_id", "taxonomy_lvl")], bracken_merged_frac)

meta <- read.csv(file = "../../Metagenomics_2024/DATA/tryp_metadata.csv", sep = ",")
```


#### 2. Basic stats

Before we start anything, let’s just check out or data a little bit. Never go blind into your analyses.

**Q1: How many samples do we have in the metadata and in out data tibble**

<details>
<summary>
HINT
</summary>

> nrow(meta)

</details>  

**Q2: How many families were detected in our dataset?**

<details>
<summary>
HINT
</summary>

> nrow(bracken_merged)

</details>  

#### 3. Pre-process data

Now that we know a little bit about our data we can start selecting what we want to look at. Let’s extract a separate tibble for species, family and phylum level. At each tibble we remove all the taxa that are not present in any of the samples.

Let’s start with species:

```R
merge_rel_abund_spe <- subset(bracken_merged_reads, taxonomy_lvl == "S") %>% 
  select(-c("taxonomy_id","taxonomy_lvl")) %>%
  filter(rowSums(select_if(., is.numeric)) != 0)
head(merge_rel_abund_spe)
```

Can you come up with the command for family (F) and phylum (P) level? 

**Tip**: make sure you save it to a different variable. Otherwise you are overwriting your species data.

<details>
<summary>
HINT
</summary>

```r
merge_rel_abund_family <- subset(bracken_merged_reads, taxonomy_lvl == "F") %>% 
  select(-c("taxonomy_id","taxonomy_lvl")) %>%
  filter(rowSums(select_if(., is.numeric)) != 0)
head(merge_rel_abund_spe)

merge_rel_abund_genus <- subset(bracken_merged_reads, taxonomy_lvl == "G") %>% 
  select(-c("taxonomy_id","taxonomy_lvl")) %>%
  filter(rowSums(select_if(., is.numeric)) != 0)
head(merge_rel_abund_spe)
```

</details>  




#### 4. Plot relative abundances
Let’s plot and explore microbial composition and relative abundances in our samples by visualizing them.

As a first overview of coarse differences we can create a stacked bar plot of phyla. We need to mingle our data structure a bit to make it compliant with ggplot and add the metadata to get an extra level of information.

```r
# data mingling
merge_rel_abund_spe1 <- merge_rel_abund_spe %>% gather(key="SRA.identifier",value="rel_abun",-name)
# now join with metadata by column 'Sample'. We are using left join in case the metadata file contains additional samples not included in our dataset
merge_rel_abund_spe1 <- left_join(merge_rel_abund_spe1,meta,by="SRA.identifier") 
# check how it looks
head(merge_rel_abund_spe1) 
```


```r
ggplot(merge_rel_abund_spe1, aes(x=SRA.identifier, y=rel_abun)) +
  geom_bar(stat="identity", position="stack", aes(fill=name)) + # chose bar plot
  theme(axis.text.x = element_text(angle=45, hjust=1)) + # put x-axis label at 45 degree angle
  facet_grid(. ~ Time, scales="free_x",space = "free_x") # produce two panels according to metatadata category 'mocktreat'
```

**Q3: Do the phyla profiles look similar between samples? Can you spot any trends?**

<details>
<summary>
HINT
</summary>

> A3: We see that the same two Phyla dominate all samples but there is some variability between samples. From a first look there does not seem to be a major difference between the two groups A and B at this taxonomic level.

</details>  


*If you have time you can also visualize the other taxonomic levels (e.g. species) with the same approach. Try to come up with the code yourself. Hint: Omit legend using legend.position.*


<details>
<summary>
HINT
</summary>

```r
# bar plot phyla relative abundances
brel_spec_gg1 <- brel_spec %>% gather(key="Sample",value="rel_abun",-taxa)
brel_spec_gg1 <- left_join(brel_spec_gg1,meta,by="Sample")

ggplot(brel_spec_gg1, aes(x=Sample, y=rel_abun)) +
  geom_bar(stat="identity",position="stack", aes(fill=taxa)) +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "none") +
  facet_grid(. ~ mocktreat, scales="free_x",space = "free_x")
```

</details>  


**Q4: Does the profile change according to taxonomic level? Is the stacked bar plot helpful in all scenarios?**


<details>
<summary>
HINT
</summary>

> A4: When too many taxa are present, such as at species level, it becomes difficult to distinguish the colors. As you might realize, when there are too many taxa it becomes very difficult to spot anything in the stacked bar plot. 

</details>  


Another way to visualize the relative abundances is by creating a bubble plot. Let’s do that for the family composition.


```r
# bubble plot
brel_fam_gg1 <- brel_fam %>% gather(key="Sample",value="rel_abun",-taxa)
brel_fam_gg1 <- left_join(brel_fam_gg1,meta,by="Sample")

ggplot(brel_fam_gg1, aes(x=Sample, y=taxa)) +
  geom_point(aes(size=rel_abun, color=timepoint), alpha=0.7) + # this time we use points
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  scale_size_continuous(limits = c(0.00001,max(brel_phy_gg1$rel_abun))) + # sets minimum above '0' 
  facet_grid(. ~ mocktreat, scales="free_x",space = "free_x")
```


**Q5: Did you notice that here we added an extra level of information? Can you spot what it is?**

<details>
<summary>
HINT
</summary>

> A5: Now we produced panels according to mocktreat and colored according to time point. Like this you can see that we only have two time points in one of the groups. This plot allows us therefore to combine multiple metadata layers.

</details>  

If you have time you can play with the different information levels and see what this can tell you about your data.

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

**Q6: what can you learn from the heatmap. Are there any informative clusters?**

<details>
<summary>
HINT
</summary>

A6: We can see that 4 species dominate the communities in most samples. Adding the metadata we can also see that data do not cluster strongly according to group or time point, but there is some degree of structuring in mocktreat.

</details>  

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

Check the output. 

```r
# Check the output
nmds_spec
```
Here you see a kind of summary of the analysis. For example, you can see that you used 2 dimensions and the stress was approx. 0.15. In general if a stress is above 0.2 then the clustering is not reliably representing the data and should be interpreted with caution. But here the stress is below 0.2, so we are okay.

Now let’s look at the ordination. To plot the data with ggplot, we need to extract the coordinaties of each point from nmds_spec$points.

```r
# Extract and reshape the data to plot ordination as ggplot  and add the metadata
nmds_spec_gg<-as.data.frame(nmds_spec$points) %>%
  rownames_to_column("Sample") %>%
  left_join(meta, by="Sample")
```
Then we can greate the plot easily and color according to the metadata. We are choosing timepoint and mocktreat for the coloring respectively. But feel free to explore other parameters.

```r
# Let's plot and color according to time point
ggplot(nmds_spec_gg, aes(x=MDS1,y=MDS2)) +
  geom_point(aes(color=timepoint), size=3, alpha=0.5) +
  ggtitle("NMDS colored according to timepoint")
```



```r
# Now let's plot according to mocktreat
ggplot(nmds_spec_gg, aes(x=MDS1,y=MDS2)) +
  geom_point(aes(color=mocktreat), size=3, alpha=0.5)+
  ggtitle("NMDS colored according to mocktreat")
```
#### 6. Differential abundance









