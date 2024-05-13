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

```r
library(microbiome)
library(dplyr)
library(hrbrthemes)
library(tibble)
library(tidyverse)
library(vegan)

library(mia)
library(tidyverse)
library(ggplot2)
library(knitr)
library(phyloseq)
library(gcookbook)
```


##### 1. load the data

First we need to load our data. Usually the biggest bottleneck between raw data and analyses is to get the data in the right shape for your purpose. Often this requires a little bit of data mingling. On this road - google is your best friend to master the R universe :)

Let’s first load the relative abundance table of the bracken results.

For this part we are using the [phyloseq](https://joey711.github.io/phyloseq/) package. The phyloseq package is a tool to import, store, analyze, and graphically display complex phylogenetic sequencing data that has already been clustered into Operational Taxonomic Units (OTUs), especially when there is associated sample data, phylogenetic tree, and/or taxonomic assignment of the OTUs. 

```R
merged_metagenomes <- import_biom("READBASED/merge_species.biom") 
meta <- read.csv(file = "DATA/tryp_metadata.csv", sep = ",")
```

#### phyloseq-ize Data

Any data already in an R session can be annoated/coerced to be recognized by phyloseq’s functions and methods. This is important, because there are lots of ways you might receive data related to a microbiome project, and not all of these will come from a popular server or workflow that is already supported by a phyloseq import function. 

![phyloseq](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/phyloseq.png)

- **otu_table** - Works on any numeric matrix. You must also specify if the species are rows or columns
- **tax_table** - Works on any character matrix. The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
- **sample_data** - Works on any data.frame. The rownames must match the sample names in the otu_table if you plan to combine them as a phyloseq-object

 Let's look at the phlyseq object. 

```r
merged_metagenomes
```

We need to add the metadata to the phyloseq object. 

#### 2. format the data

```r
meta <- meta  %>%  arrange(row_number(SRA.identifier)) # sort the data frame
merged_metagenomes@sam_data <- sample_data(meta) # associate the metadata to the to the phyloseq object
column_name <-  meta %>% pull(Sample) # extract the sample names
sample_names(merged_metagenomes) <- column_name # associate the sample names to the phyloseq object

merged_metagenomes@tax_table@.Data <- substring(merged_metagenomes@tax_table@.Data, 4) # remove the unnecessary 'k_' in the taxonomy.
colnames(merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") # change the rank name 
```


#### 3. Basic stats

Before we start anything, let’s just check out or data a little bit (sanity check). Never go blind into your analyses.
 



Summurize the data.  
```r
summarize_phyloseq(merged_metagenomes)

head(otu_table(merged_metagenomes))
head(tax_table(merged_metagenomes))
head(sample_data(merged_metagenomes))

sample_variables(merged_metagenomes)


```


```r
head(psmelt(merged_metagenomes))
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

Etraction of the sample's names. 

```r
sample_names(merged_metagenomes)
```

#### Aggregation

Microbial species can be called at multiple taxonomic resolutions. We can easily agglomerate the data based on taxonomic ranks. Here, we agglomerate the data at Family level.

```r

merged_metagenomes_family <- aggregate_rare(merged_metagenomes, level = "Family", detection = 0/100, prevalence = 0/100) # no filter except family
# Show dimensionality
dim(abundances(merged_metagenomes_family))

```

**Q: How many sample and tax do we have now?**

<details>
<summary>
HINT
</summary>

> Now there are 53 taxa and 30 samples, meaning that there are 53 different Phylum level taxonomic groups. Looking at the rowData after agglomeration shows all Enterococcaceae are combined together, and all lower rank information is lost.

```r
knitr::kable(head(tax_table(merged_metagenomes_family))) %>% 
  kableExtra::kable_styling("striped", 
                            latex_options="scale_down") %>% 
  kableExtra::scroll_box(width = "100%")
```

</details> 



### 4. QC and Pre-process data

Now that we know a little bit about our data we can start the pre-processing. 


####  Library size / read count


#### Prevalence - Detection 


Prevalence quantifies the frequency of samples where certain microbes were detected (above a given detection threshold). The prevalence can be given as sample size (N) or percentage (unit interval).

The population prevalence (frequency) at a 1% relative abundance threshold (detection = 1/100 and count = false), can look like this.

```r
prevalence(merged_metagenomes, detection=1/100, sort=TRUE, count=FALSE)
```

The function arguments detection and count can also be used to access, how many samples do pass a threshold for raw counts. Here, the population prevalence (frequency) at the absolute abundance threshold (count=true) at read count 1 (detection = 1) is accessed.

```r
prevalence(merged_metagenomes, detection=1, sort=TRUE, count=true)

plot_taxa_prevalence(merged_metagenomes, level="Phylum", detection = 10000)

```

Each point corresponds to a different or unique taxon. The y-axis represents the fraction of samples, these taxa are present. The low prevalence suggests there is a low overlap across samples. 

![prevalence](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/prevalence.png)

**Q: Which taxa is present in all samples and with a high abundance?**

<details>
<summary>
HINT
</summary>

> Homo sapiens

</details> 


##### Contaminant sequences
Samples might be contaminated with exogenous sequences. We have observed 1 contaminants Homo sapiens.


```r
# remove contaminants
contaminants <- c("9606") 
allTaxa = taxa_names(merged_metagenomes)
allTaxa <- allTaxa[!(allTaxa %in% contaminants)]
merged_metagenomes = prune_taxa(allTaxa, merged_metagenomes)
head(tax_table(merged_metagenomes))
```

Now recheck that your data are clean before continuing the analysis. 



### 5. Microbiome composition

Microbial abundances are typically ‘compositional’ (relative) in the current microbiome profiling data sets. This is due to technical aspects of the data generation process (see e.g. Gloor et al., 2017).

The next example calculates relative abundances as these are usually easier to interpret than plain counts. For some statistical models we need to transform the data into other formats as explained in above link (and as we will see later).

*detection*: Detection threshold for absence/presence (percentage reads).
*prevalence*: Prevalence threshold (in [0, 1]) (presence accross samples).

```r
pseq <- aggregate_rare(merged_metagenomes, level = "Family", detection = 0.1/100, prevalence = 50/100)

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
p <- plot_composition(pseq,
                      taxonomic.level = "Family",
                      sample.sort = "Sample",
                      x.label = "Sample") +
  scale_fill_brewer("Family", palette = "Paired") +
  guides(fill = guide_legend(ncol = 1)) +
  scale_y_percent() +
  labs(x = "Samples", y = "Relative abundance (%)",
       title = "Relative abundance data") + 
  theme_ipsum(grid="Y") +
  theme(axis.text.x = element_text(angle=90, hjust=1),
        legend.text = element_text(face = "italic"))
print(p)  
```

![barplot](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/barplot_family.png)

**Q: Does the profiles look similar between samples? Can you spot any trends?**

<details>
<summary>
HINT
</summary>

> A: We see that the same one or two OTU dominates all samples but there is some variability between samples. We can observe some trend on the Entrococcus proption over time.

</details>  


*If you have time you can also visualize the other taxonomic levels (e.g. species) with the same approach. Try to come up with the code yourself.*


**Density plot** shows the overall abundance distribution for a given taxonomic group. Let us check the relative abundance of Enterococcus across the sample collection. The density plot is a smoothened version of a standard histogram.

```r
pseq <- microbiome::transform(merged_metagenomes, "compositional")
Enterococcus_abun = prune_taxa("1351", pseq)
Enterococcus_abun_df <- t(abundances(Enterococcus_abun)) 
Enterococcus_abun_df <- rownames_to_column(as.data.frame(Enterococcus_abun_df), var = "Sample") %>% arrange(Sample)
meta_df <- meta %>% arrange(Sample)
Enterococcus_abun_df_meta <- merge(Enterococcus_abun_df, meta_df, by = "Sample", all = TRUE) %>% select("Time", "Morning.Afternoon", "Type", "1351") %>% filter(Morning.Afternoon == "AM") %>% select(-c("Morning.Afternoon", "Time")) 
Enterococcus_abun_df_meta$"1351" <- Enterococcus_abun_df_meta$"1351"*100
colnames(Enterococcus_abun_df_meta) <- c("Type", "Enterococcus")

Enterococcus_abund_plot <- ggplot(Enterococcus_abun_df_meta, aes(x = as.numeric(Enterococcus), color = Type, fill = Type)) + 
  geom_density(alpha = 0.3) + 
  labs(x = "Relative abundance", title = "Enterococcus") +
  theme_classic()
  
```

![density_plot](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/density_Entero.png)




###  Exercises Normalization (optional) 

**Transformations** The data contains read counts. We can convert these into relative abundances and other formats. Compare abundance of a given taxonomic group using the example data before and after the compositionality transformation (with a cross-plot, for instance). You can also compare the results to CLR-transformed data (see e.g. Gloor et al. 2017)




# Alpha diversity

```r
pseq <- aggregate_rare(merged_metagenomes, level = "Species", detection = 0.1/100, prevalence = 50/100)

tab <-microbiome::alpha(pseq, index = "all")
kable(head(tab))
```
```r
tab <- richness(pseq)
kable(head(tab))
```

```r
p.shannon <- boxplot_alpha(pseq, 
                           index = "shannon",
                           x_var = "Type",
                          fill.colors = c(Control="cyan4", T._cruzi="deeppink4",T._rangeli="darkorange1" ))
p.shannon
```

![shanon](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/shanon.png)

### Investigate the top factors

Show coefficients for the top taxa separating the groups

```r
pseq <- aggregate_rare(merged_metagenomes, level = "Genus", detection = 0.1/100, prevalence = 50/100)
pseq <- microbiome::transform(pseq, transform = "compositional")
abundances(pseq)

p <- plot_landscape(pseq, method = "NMDS", distance = "bray", col = "Time", size = 3)
print(p)
```
![association](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/pcq.png)


### Estimating associations with an external variable
Next to visualizing whether any variable is associated with differences between samples, we can also quantify the strength of the association between community composition (beta diversity) and external factors.

The standard way to do this is to perform a so-called permutational multivariate analysis of variance (PERMANOVA). This method takes as input the abundance table, which measure of distance you want to base the test on and a formula that tells the model how you think the variables are associated with each other.


```r
otu <- abundances(pseq)
meta <- meta(pseq)

permanova <- adonis(t(otu) ~ Time,
                    data = meta, permutations=9999, method = "euclidean")
```

# P-value
```r
print(as.data.frame(permanova$aov.tab))
```

The time variable is significantly associated with microbiota composition (p-value is below 0.05).

We can, however, visualize those taxa whose abundances drive the differences between cohorts. We first need to extract the model coefficients of taxa:

```r
coef <- coefficients(permanova)["Time1",]
top.coef <- coef[rev(order(abs(coef)))[1:20]]
par(mar = c(3, 14, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")
```

![association](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/associations.png)

**Enterococcus** (Bacillota), which plays a crucial role in metabolic adaptability against pathogenic or plant toxins and anti-herbivore defense, was found to be one of the predominant gut microorganism of lepidopteran insects, including B. mori, Helicoverpa zea, and Porthetria dispar (Paniagua Voirol et al., 2018; Zhang et al., 2022). 

**Symbiopectobacterium** (Enterobacteriaceae) has recently been described for the first time as an intracellular bacterial symbiont, responsible for to the biosynthesis of vitamins and cofactors.  This bacteria may be boosting the parsite fitness, for example, by aiding in evading the Triatome immune response, or providing a novel function, such as supplementing nutrition or metabolism

**Rhodococcus** (Nocardiaceae) in the triatomine gut are believed to play an important role in the metabolism of the vector, such as by participating in the synthesis of group B vitamins or by being digested by the bugs directly to provide missing nutrients (Sassera et al., 2013). Moreover, the most attractive aspect is the host-symbiont relationship between triatomines and Rhodococcus; since Rhodococcus bacteria can be easily cultured and genetically modified to harm the pathogen in vector gut, they are probably suitable tools for the control of trypanosomiasis (Sassera et al., 2013). as the blood is poor in B vitamins compared to what is generally required for insect development. The blood is poor in B vitamins compared to what is generally required for insect development, Kissing bugs, Rhodnius prolixus, notably require *Rhodococcus* bacteria for nymph development, but the addition of B vitamins in the diet can rescue nymph development in the absence of Rhodococcus  (Serrato-Salas and Gendrin 2023).

**Wolbachia** (Ehrlichiaceae)  The obligate intracellular bacteria Wolbachia spp. are common in a wide range of insects, including sand flies, bed bugs, fleas and mosquitoes, and can cause reproduction alterations such as feminization, male killing and cytoplasmic incompatibility ([Landmann 20219](https://doi.org/10.1128/microbiolspec.BAI-0018-2019.)). In triatomines, Wolbachia has been solely reported for the genus Rhodnius, where it occurs in the intestine, salivary glands and gonads.

**Curtobacterium** (Microbacteriaceae) 



Exercises
Community-level comparisons: Use PERMANOVA to investigate whether the community composition differs between two groups of individuals (e.g. times, or some other grouping of your choice). You can also include covariates such as age and gender, and see how this affects the results.



























































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





