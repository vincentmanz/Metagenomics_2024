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
library(ggplot2)

```


##### 1. load the data

First we need to load our data. Usually the biggest bottleneck between raw data and analyses is to get the data in the right shape for your purpose. Often this requires a little bit of data mingling. On this road - google is your best friend to master the R universe :)

Let’s first load the relative abundance table of the bracken results.

```R
bracken_species <- read.csv(file = "READBASED/BRACKEN/bracken_merged_species.csv", sep = "\t") 
bracken_genus <- read.csv(file = "READBASED/BRACKEN/bracken_merged_genus.csv", sep = "\t") 
bracken_family <- read.csv(file = "READBASED/BRACKEN/bracken_merged_family.csv", sep = "\t") 
bracken_phylum <- read.csv(file = "READBASED/BRACKEN/bracken_merged_phylum.csv", sep = "\t") 

meta <- read.csv(file = "../../Metagenomics_2024/DATA/tryp_metadata.csv", sep = ",")
```


#### 2. Basic stats

Before we start anything, let’s just check out or data a little bit (sanity check). Never go blind into your analyses.

**Q: How many samples do we have in the metadata and in out data tibble**

<details>
<summary>
HINT
</summary>

> nrow(meta)

</details>  

**Q: How many families were detected in our dataset?**

<details>
<summary>
HINT
</summary>

> nrow(bracken_merged)

</details>  

**Q: Look at the data frame is ther some problem with the taxonomy**

<details>
<summary>
HINT
</summary>

>  The family level has been recoded automatically form F for family to FALSE. 

</details>  



#### 3. Pre-process data

Now that we know a little bit about our data we can start selecting what we want to look at. Let’s extract a separate tibble for species, family and phylum level. At each tibble we remove all the taxa that are not present in any of the samples.

Let’s start with species:

```r


# Change "FALSE" to "F" in the column "taxonomy_lvl"
bracken_family$taxonomy_lvl <- ifelse(bracken_family$taxonomy_lvl == "FALSE", "F", bracken_family$taxonomy_lvl)

# Select columns ending with ".bracken_num" for bracken_merged_reads
species_reads <- bracken_species[grep("\\.bracken_num$", names(bracken_species))]
genus_reads <- bracken_genus[grep("\\.bracken_num$", names(bracken_genus))]
family_reads <- bracken_family[grep("\\.bracken_num$", names(bracken_family))]
phylum_reads <- bracken_phylum[grep("\\.bracken_num$", names(bracken_phylum))]

# Select columns ending with ".bracken_frac" for bracken_merged_frac
species_frac <- bracken_species[grep("\\.bracken_frac$", names(bracken_species))]
genus_frac <- bracken_genus[grep("\\.bracken_frac$", names(bracken_genus))]
family_frac <- bracken_family[grep("\\.bracken_frac$", names(bracken_family))]
phylum_frac <- bracken_phylum[grep("\\.bracken_frac$", names(bracken_phylum))]

# Remove the suffixes from the column names
colnames(species_reads) <- sub("\\_species_filtered.bracken_num$", "", colnames(species_reads))
colnames(species_frac) <- sub("\\_species_filtered.bracken_frac$", "", colnames(species_frac))
colnames(genus_reads) <- sub("\\_genus_filtered.bracken_num$", "", colnames(genus_reads))
colnames(genus_frac) <- sub("\\_genus_filtered.bracken_frac$", "", colnames(genus_frac))
colnames(family_reads) <- sub("\\_family_filtered.bracken_num$", "", colnames(family_reads))
colnames(family_frac) <- sub("\\_family_filtered.bracken_frac$", "", colnames(family_frac))
colnames(phylum_reads) <- sub("\\_phylum_filtered.bracken_num$", "", colnames(phylum_reads))
colnames(phylum_frac) <- sub("\\_phylum_filtered.bracken_frac$", "", colnames(phylum_frac))


# Add name, taxonomy_id, and taxonomy_lvl to both data frames
species_reads <- cbind(bracken_species[, c("name", "taxonomy_id", "taxonomy_lvl")], species_reads)
species_frac <- cbind(bracken_species[, c("name", "taxonomy_id", "taxonomy_lvl")], species_frac)
genus_reads <- cbind(bracken_genus[, c("name", "taxonomy_id", "taxonomy_lvl")], genus_reads)
genus_frac <- cbind(bracken_genus[, c("name", "taxonomy_id", "taxonomy_lvl")], genus_frac)
family_reads <- cbind(bracken_family[, c("name", "taxonomy_id", "taxonomy_lvl")], family_reads)
family_frac <- cbind(bracken_family[, c("name", "taxonomy_id", "taxonomy_lvl")], family_frac)
phylum_reads <- cbind(bracken_phylum[, c("name", "taxonomy_id", "taxonomy_lvl")], phylum_reads)
phylum_frac <- cbind(bracken_phylum[, c("name", "taxonomy_id", "taxonomy_lvl")], phylum_frac)
head(phylum_frac)

# filter data
species_frac_filtered <- subset(species_frac, taxonomy_lvl == "S") %>% 
  select(-c("taxonomy_id","taxonomy_lvl")) %>%
  filter(rowSums(select_if(., is.numeric)) >= 0.001)
genus_frac_filtered <- subset(genus_frac, taxonomy_lvl == "G") %>% 
  select(-c("taxonomy_id","taxonomy_lvl")) %>%
  filter(rowSums(select_if(., is.numeric)) >= 0.001)
head(genus_frac_filtered)
family_frac_filtered <- subset(family_frac, taxonomy_lvl == "F") %>% 
  select(-c("taxonomy_id","taxonomy_lvl")) %>%
  filter(rowSums(select_if(., is.numeric)) >= 0.001)
head(family_frac_filtered)
phylum_frac_filtered <- subset(phylum_frac, taxonomy_lvl == "P") %>% 
  select(-c("taxonomy_id","taxonomy_lvl")) %>%
  filter(rowSums(select_if(., is.numeric)) >= 0.001)
head(phylum_frac_filtered)
```

**Q Can you come up with the command for family (S) and phylum (P) level?**

**Tip**: make sure you save it to a different variable. Otherwise you are overwriting your species data.

<details>
<summary>
HINT
</summary>

> If not done, run bracken on diffrent levels, load the merged file into rstudio. 

</details>  




#### 4. Plot relative abundances
Let’s plot and explore microbial composition and relative abundances in our samples by visualizing them.

As a first overview of coarse differences we can create a stacked bar plot of phyla. We need to mingle our data structure a bit to make it compliant with ggplot and add the metadata to get an extra level of information.

```r

# data mingling
species_level <- species_frac_filtered %>% gather(key="SRA.identifier",value="rel_abun",-name)
# now join with metadata by column 'Sample'. We are using left join in case the metadata file contains additional samples not included in our dataset
species_level <- left_join(species_level,meta,by="SRA.identifier") 
# check how it looks
head(species_level) 
genus_level <- genus_frac_filtered %>% gather(key="SRA.identifier",value="rel_abun",-name)
# now join with metadata by column 'Sample'. We are using left join in case the metadata file contains additional samples not included in our dataset
genus_level <- left_join(genus_level,meta,by="SRA.identifier") 
# check how it looks
head(genus_level) 
family_level <- family_frac_filtered %>% gather(key="SRA.identifier",value="rel_abun",-name)
# now join with metadata by column 'Sample'. We are using left join in case the metadata file contains additional samples not included in our dataset
family_level <- left_join(family_level,meta,by="SRA.identifier") 
# check how it looks
head(family_level) 
phylum_level <- phylum_frac_filtered %>% gather(key="SRA.identifier",value="rel_abun",-name)
# now join with metadata by column 'Sample'. We are using left join in case the metadata file contains additional samples not included in our dataset
phylum_level <- left_join(phylum_level,meta,by="SRA.identifier") 
# check how it looks
head(phylum_level) 


```


```r
ggplot(phylum_level, aes(x=Sample, y=rel_abun)) +
  geom_bar(stat="identity", position="stack", aes(fill=name)) + # chose bar plot
  theme(axis.text.x = element_text(angle=45, hjust=1)) + # put x-axis label at 45 degree angle
  facet_grid(. ~ Time, scales="free_x",space = "free_x") # produce two panels according to metatadata category 'Time' 
  
```

![barplot](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/barplot_phylum.png)

**Q: Do the phyla profiles look similar between samples? Can you spot any trends?**

<details>
<summary>
HINT
</summary>

> A: We see that the same two Phyla dominate all samples but there is some variability between samples. From a first look there does not seem to be a major difference between the two groups A and B at this taxonomic level.

</details>  


*If you have time you can also visualize the other taxonomic levels (e.g. species) with the same approach. Try to come up with the code yourself. Hint: Omit legend using legend.position (guides(fill = FALSE)).*


<details>
<summary>
HINT
</summary>

```r
ggplot(genus_level, aes(x=Sample, y=rel_abun)) +
  geom_bar(stat="identity", position="stack", aes(fill=name)) + # chose bar plot
  theme(axis.text.x = element_text(angle=45, hjust=1)) + # put x-axis label at 45 degree angle
  facet_grid(. ~ Time, scales="free_x",space = "free_x") +# produce two panels according to metatadata category 'Time' 
  guides(fill = FALSE)
ggplot(family_level, aes(x=Sample, y=rel_abun)) +
  geom_bar(stat="identity", position="stack", aes(fill=name)) + # chose bar plot
  theme(axis.text.x = element_text(angle=45, hjust=1)) + # put x-axis label at 45 degree angle
  facet_grid(. ~ Time, scales="free_x",space = "free_x") + # produce two panels according to metatadata category 'Time' 
  guides(fill = FALSE)
```

</details>  

![barplot](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/barplot_family.png)


**Q: Does the profile change according to taxonomic level? Is the stacked bar plot helpful in all scenarios?**


<details>
<summary>
HINT
</summary>

> A: When too many taxa are present, such as at species level, it becomes difficult to distinguish the colors. As you might realize, when there are too many taxa it becomes very difficult to spot anything in the stacked bar plot. 

</details>  


Another way to visualize the relative abundances is by creating a bubble plot. Let’s do that for the family composition.


```r
# bubble plot
ggplot(phylum_level, aes(x=Sample, y=name)) +
  geom_point(aes(size=rel_abun, color=Type), alpha=0.7) + # this time we use points
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  scale_size_continuous(limits = c(0.00001,max(phylum_level$rel_abun))) + # sets minimum above '0'
  facet_grid(. ~ Time, scales="free_x",space = "free_x")

```
![bubble](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/bubble_plot.png)


**Q: Did you notice that here we added an extra level of information? Can you spot what it is?**

<details>
<summary>
HINT
</summary>

> A: Now we produced panels according to Time and colored according to Type. This plot allows us to combine multiple metadata layers. Can you see the change for Enterococcus? 

</details>  

If you have time you can play with the different information levels and see what this can tell you about your data.

When we want to look at species composition, the communities are often very complex and looking at abundance profiles is not as informative. One way to explore species-level community composition is to filter you data so you focus only on the most abundant taxa across samples. What cutoff and approach you use here of course depends on your question and data.

Let’s create a heatmap of the 20 most abundant spexies in our samples.


```r
# sort species by abundance across samples and select top 20
species_level_abundance <- species_frac %>% select(-taxonomy_id, -taxonomy_lvl) %>% head(20) %>% as_tibble() %>%  column_to_rownames( "name")
genus_level_abundance <- genus_frac %>% select(-taxonomy_id, -taxonomy_lvl) %>% head(20) %>% as_tibble() %>%  column_to_rownames( "name")
phylum_level_abundance <- phylum_frac %>% select(-taxonomy_id, -taxonomy_lvl) %>% head(20) %>% as_tibble() %>%  column_to_rownames( "name")


# shape the metadata
meta_s <- subset(meta, SRA.identifier %in% colnames(species_level_abundance)) %>%  as_tibble() %>% 
  column_to_rownames("Sample") # shoft the 'Sample' column to rownames
meta_h <- subset(meta, SRA.identifier %in% colnames(genus_level_abundance)) %>%  as_tibble() %>% 
  column_to_rownames("Sample") # shoft the 'Sample' column to rownames
meta_p <- subset(meta, SRA.identifier %in% colnames(phylum_level_abundance)) %>%  as_tibble() %>% 
  column_to_rownames("Sample") # shoft the 'Sample' column to rownames

column_name <- meta_h %>% rownames_to_column() %>%  arrange(row_number(SRA.identifier)) %>% pull(rowname)
colnames(species_level_abundance) <- column_name
colnames(genus_level_abundance) <- column_name
colnames(phylum_level_abundance) <- column_name


# plot the heatmap
pheatmap::pheatmap(genus_level_abundance,
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   annotation_col = meta_h[,c(2,3,4,5)],
                   annotation_names_col=TRUE)

pheatmap::pheatmap(phylum_level_abundance,
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   annotation_col = meta_h[,c(2,3,4,5)],
                   annotation_names_col=TRUE)
pheatmap::pheatmap(species_level_abundance,
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   annotation_col = meta_h[,c(2,3,4,5)],
                   annotation_names_col=TRUE)
```
![bubble](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/heatmap_phy.png)


**Q: what can you learn from the heatmap. Are there any informative clusters?**

<details>
<summary>
HINT
</summary>

> A: We can see that 3 phylum dominate the communities in most samples. Adding the metadata we can also see that data do not cluster strongly according to group or time point, but there is some degree of structuring.

</details>  




#### 5. Beta diversity


Often we want to know whether the microbiomes are different between conditions or groups. One way to explore this is to look at the beta-diversity in an ordination. There are different distances and approaches that can be done and explored. We will perform an NMDS on bray curtis dissimilarities of the species profiles.


```r
# To ensure reproducibility we can fix the seed here. This will ensure you always get the same result each time you run your data.
set.seed(34521)

# Data mingling
species_frac_filtered_f <- species_frac_filtered %>% 
  column_to_rownames("name") %>% 
  t() # transpose

# Calculate distance matrix
species_frac_filtered_dist <- vegdist(species_frac_filtered_f, method = "bray")

# Perform NMDS on distance matrix
nmds_spec <- metaMDS(species_frac_filtered_dist,distance = "bray",k = 2)
```

Check the output. 

```r
# Check the output
nmds_spec
```
Here you see a kind of summary of the analysis. For example, you can see that you used 2 dimensions and the stress was approx. 0.05. In general if a stress is above 0.2 then the clustering is not reliably representing the data and should be interpreted with caution. But here the stress is below 0.2, so we are okay.

Now let’s look at the ordination. To plot the data with ggplot, we need to extract the coordinaties of each point from nmds_spec$points.

```r
# Extract and reshape the data to plot ordination as ggplot  and add the metadata
nmds_spec_gg<-as.data.frame(nmds_spec$points) %>%
  rownames_to_column("Sample") %>%
  left_join(meta, by="Sample")
```
Then we can create the plot easily and color according to the metadata. We are choosing timepoint and mocktreat for the coloring respectively. But feel free to explore other parameters.

```r
# Let's plot and color according to time point
ggplot(nmds_spec_gg, aes(x=MDS1,y=MDS2)) +
  geom_point(aes(color=Time), size=3, alpha=0.5) +
  ggtitle("NMDS colored according to Time")
```

![nmds_time](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/nmds_time.png)


```r
# Let's plot and color according to type
ggplot(nmds_spec_gg, aes(x=MDS1,y=MDS2)) +
  geom_point(aes(color=Type), size=3, alpha=0.5) +
  ggtitle("NMDS colored according to Type")
```
![nmds_type](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/nmds_type.png)



#### 6. Differential abundance









