# Day 2

| Time      | Activity                      | Slides                                | Hands-on                                                |
|-----------|-------------------------------|---------------------------------------|---------------------------------------------------------|
| Morning   | Read Base analysis            |                                       |       [Link here](Readbased.md)                         | 




# Metagenomic read based profiles

In this tutorial we will explore and analyse the results of kraken2 read profiling data using the re-estimatd abundances of bracken. These exercises are meant to show how to conceptually approach your data analysis but there are many more and different ways to explore your data. The most important thing to keep in mind is that you have to understand your own data and analyses. One way to achieve this is to perform visual explorations that help you to judge whether the data are appropriate for your question.

Now let’s start the fun!

In R studio.

Load libraries:

```r
library(microbiome)
library(microbiomeutilities)
library(dplyr)
library(hrbrthemes)
library(tibble)
library(tidyverse)
library(vegan)
library(compositions)
library(pheatmap)
library(ggplot2)
library(mia)
library(tidyverse)
library(RColorBrewer)
```


## 1. load the data

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

Summarize the data.  

```r
# Summarize the phyloseq object 'merged_metagenomes'
summarize_phyloseq(merged_metagenomes)

# Display the first few rows of the OTU (Operational Taxonomic Unit) table
head(otu_table(merged_metagenomes))

# Display the first few rows of the taxonomy table
head(tax_table(merged_metagenomes))

# Display the first few rows of the sample data associated with 'merged_metagenomes'
head(sample_data(merged_metagenomes))

# Get the sample variables of the phyloseq object 'merged_metagenomes'
sample_variables(merged_metagenomes)
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

### Aggregation

Microbial species can be called at multiple taxonomic resolutions. We can easily agglomerate the data based on taxonomic ranks. Here, we agglomerate the data at Family level.

```r
# Aggregate rare taxa at the family level for the phyloseq object 'merged_metagenomes'
merged_metagenomes_family <- aggregate_rare(merged_metagenomes, level = "Family", detection = 0/100, prevalence = 0/100)

# Display the dimensionality of the abundances of 'merged_metagenomes_family'
dim(abundances(merged_metagenomes_family))
```

**Q: How many sample and tax do we have now?**

<details>
<summary>
HINT
</summary>

> There are 53 taxa and 28 samples, meaning that there are 29 different Phylum level taxonomic groups. Looking at the rowData after agglomeration shows all Enterococcaceae are combined together, and all lower rank information is lost.

```r
head(tax_table(merged_metagenomes_family))
```

</details> 


## 4. QC and Pre-process data

Now that we know a little bit about our data we can start the pre-processing. 


###  Library size / read count

Let us check for distribution of number of sequences retained from the Kraken/Bracken approach.



```r
plot_read_distribution(merged_metagenomes, "Type", plot.type = "density")

```
![density](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/density_reads.png)

*You can try to plot with diffrent metadata.*


Plotting the read count per sample

```r
df <- psmelt(merged_metagenomes)  %>%  group_by(Sample, Time) %>%  
  summarise(sum_reads = sum(Abundance)) %>% arrange(sum_reads) 

ggplot(df) +
  geom_bar(aes(reorder(Sample, -sum_reads), sum_reads, fill=Time),
           col="red", alpha = .2, stat="identity") 

```


![density](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/reads_count.png)


### Prevalence - Detection 


Prevalence quantifies the frequency of samples where certain microbes were detected (above a given detection threshold). The prevalence can be given as sample size (N) or percentage (unit interval).

The population prevalence (frequency) at a 1% relative abundance threshold (detection = 1/100 and count = false), can look like this.

```r
prevalence(merged_metagenomes, detection=1/100, sort=TRUE, count=FALSE)
```

The function arguments detection and count can also be used to access, how many samples do pass a threshold for raw counts. Here, the population prevalence (frequency) at the absolute abundance threshold (count=true) at read count 1 (detection = 1) is accessed.

```r
# Calculate the prevalence of taxa in the phyloseq object 'merged_metagenomes'
# Detection threshold is set to 1, taxa are sorted by prevalence, and count of taxa is returned
prevalence(merged_metagenomes, detection = 1, sort = TRUE, count = TRUE)

# Plot the prevalence of taxa at the Phylum level with a detection threshold of 10000
plot_taxa_prevalence(merged_metagenomes, level = "Phylum", detection = 10000)

# Alternative plot: Create a scatter plot of the coefficient of variation (CV) of taxa in 'merged_metagenomes'
# The CV is a measure of relative variability, calculated as the standard deviation divided by the mean
p1 <- plot_taxa_cv(merged_metagenomes, plot.type = "scatter")

# Apply a log10 scale to the x-axis of the scatter plot
p1 + scale_x_log10()
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


#### Contaminant sequences
Samples might be contaminated with exogenous sequences. We have observed 1 contaminants Homo sapiens.


```r
#check 
psmelt(subset_taxa(merged_metagenomes, Genus == "Trypanosoma"))
psmelt(subset_taxa(merged_metagenomes, Genus == "Homo"))

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

In order to assess the effect of the different transformations, you can use the *clr* or the *log10* for the  downstream analysis. 

## 5. Microbiome composition

Microbial abundances are typically ‘compositional’ (relative) in the current microbiome profiling data sets. This is due to technical aspects of the data generation process (see e.g. Gloor et al., 2017).

The next example calculates relative abundances as these are usually easier to interpret than plain counts. For some statistical models we need to transform the data into other formats as explained in above link (and as we will see later).

*detection*: Detection threshold for absence/presence (percentage reads).
*prevalence*: Prevalence threshold (in [0, 1]) (presence accross samples).

```r
pseq <- aggregate_rare(merged_metagenomes, level = "Family", detection = 0.05/100, prevalence = 20/100)
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

![barplot](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/barplot_family.png)

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

![barplot](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/barplot_genus.png)

## 6 Diversity


Species diversity, in its simplest definition, is the number of species in a particular area and their relative abundance (evenness). Once we know the taxonomic composition of our metagenomes, we can do diversity analyses. Here we will discuss the two most used diversity metrics, α diversity (within one metagenome) and β (across metagenomes).

- *α* Diversity: Can be represented only as richness (, i.e., the number of different species in an environment), or it can be measured considering the abundance of the species in the environment as well (i.e., the number of individuals of each species inside the environment). To measure α-diversity, we use indexes such as Shannon’s, Simpson’s, Chao1, etc.

![density_plot](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/diversity1.png)
*Alpha diversity is calculated according to fish diversity in a pond. Here, alpha diversity is represented in its simplest way: Richness.*

In the next example, we will look at the α and the β components of the diversity of a dataset of fishes in three lakes. The most simple way to calculate the β-diversity is to calculate the distinct species between two lakes (sites). Let us take as an example the diversity between Lake A and Lake B. The number of species in Lake A is 3. To this quantity, we will subtract the number of these species that are shared with the Lake B: 2. So the number of unique species in Lake A compared to Lake B is (3-2) = 1. To this number, we will sum the result of the same operations but now take Lake B as our reference site. In the end, the β diversity between Lake A and Lake B is (3-2) + (3-2) = 2. This process can be repeated, taking each pair of lakes as the focused sites.

![density_plot](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/diversity2.png)
*Alpha and beta diversity indexes of fishes in a pond.*

- *β* diversity mesures how different two or more communities are, either in their composition (richness) or in the abundance of the organisms that compose it (abundance).

- *Bray-Curtis dissimilarity*: The difference in richness and abundance across environments (samples). Weight on abundance. Measures the differences from 0 (equal communities) to 1 (different communities)
- *Jaccard distance*: Based on the presence/absence of species (diversity). It goes from 0 (same species in the community) to 1 (no species in common)
- *UniFrac*: Measures the phylogenetic distance; how alike the trees in each community are. There are two types, without weights (diversity) and with weights (diversity and abundance)
There are different ways to plot and show the results of such analysis. Among others, PCA, PCoA, or NMDS analysis are widely used.

![density_plot](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/diversity3.png)

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

# Display the first few rows of the calculated Shannon diversity index table
head(tab)
```
This returns observed richness with given detection threshold(s).

```r
# Calculate the richness of the phyloseq object 'pseq' with a detection threshold of 1000
# Richness is a measure of the number of different taxa present in a sample
tab <- richness(pseq, detection = 1000)

# Display the first few rows of the calculated richness table
head(tab)
```

```r
p.shannon <- boxplot_alpha(pseq, 
                           index = "shannon",
                           x_var = "Type",
                          fill.colors = c(Control="cyan4", T._cruzi="deeppink4",T._rangeli="darkorange1" ))
p.shannon
```

![shanon](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/shanon.png)

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
![indices](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/indices.png)

Each of these metrics can give an insight into the distribution of the OTUs inside our samples. For example, the Chao1 diversity index gives more weight to singletons and doubletons observed in our samples, while Shannon is an entropy index remarking the impossibility of taking two reads out of the metagenome “bag” and that these two will belong to the same OTU.

**Q: What do you observe?**


### Investigate the top factors

Show coefficients for the top taxa separating the groups

```r
pseq <- aggregate_rare(merged_metagenomes, level = "Genus", detection = 0.1/100, prevalence = 50/100)
pseq <- microbiome::transform(pseq, transform = "compositional")
abundances(pseq)

p <- plot_landscape(pseq, method = "NMDS", distance = "bray", col = "Time", size = 3)
print(p)
```
![association](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_2/img/pca.png)


### Estimating associations with an external variable

Next to visualizing whether any variable is associated with differences between samples, we can also quantify the strength of the association between community composition (beta diversity) and external factors.

The standard way to do this is to perform a so-called permutational multivariate analysis of variance (PERMANOVA). This method takes as input the abundance table, which measure of distance you want to base the test on and a formula that tells the model how you think the variables are associated with each other.

PERMANOVA tests the hypothesis that the centroids and dispersion of the community are equivalent between the compared groups. A p-value smaller than the significance threshold indicates that the groups have a different community composition. This method is implemented with the adonis2 function from the vegan package.

By default, the argument by is set to "terms", in which the order of variables in the formula matters. In this case, each variable is analyzed sequentially, and the result is different when more than 1 variable is introduced and their order differs. Therefore, it is recommended to set by = "margin", which specifies that the marginal effect of each variable is analyzed individually. You can view a comparison between the two designs in chapter ?sec-compare-permanova.



```r
otu <- abundances(pseq)
meta <- meta(pseq)

permanova <- adonis(t(otu) ~ Time,
                    data = meta, permutations=9999, method = "euclidean")
```

P-value: 

```r
print(as.data.frame(permanova$aov.tab))
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

Community-level comparisons: Use PERMANOVA to investigate whether the community composition differs between two groups of individuals (e.g. times, or some other grouping of your choice). You can also include covariates such as type, and see how this affects the results?




### Beta diversity


Often we want to know whether the microbiomes are different between conditions or groups. One way to explore this is to look at the beta-diversity in an ordination. There are different distances and approaches that can be done and explored. We will perform an NMDS on bray curtis dissimilarities of the species profiles.


```r
# To ensure reproducibility we can fix the seed here. This will ensure you always get the same result each time you run your data.
set.seed(34521)

https://telatin.github.io/microbiome-bioinformatics/data/kraken-r/2021-03-33-ExploreMGprofiles_solutions.html
species_frac_filtered <- psmelt(pseq) %>% head(20) %>% select("OTU", "Abundance", "Type", "Reads", "Time") %>%  as.data.frame()
rownames(species_frac_filtered) <- NULL
# Data mingling
species_frac_filtered_f <- species_frac_filtered %>% 
  column_to_rownames("OTU") %>% 
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










### Confounding effects
Confounders can be defined as variables that are related to and affect the apparent dynamics between the response and the main independent variable. They are common in experimental studies. Generally, they can be classified into 3 groups:

- Biological confounders, such as age and sex

- Technical confounders produced during sample collection, processing and analysis

- Confounders resulting from experimental models, such as batch effects and sample history

Controlling for confounders is an important practice to reach an unbiased conclusion. To perform causal inference, it is crucial that the method is able to include confounders in the model. This is not possible with statistical tests of general use, such as the Wilcoxon test. In contrast, methods that target Differential Abundance Analysis (DAA), such as those described in this chapter, allow controlling for confounders. In the following examples, we will perform DAA with a main independent variable and a few confounders.

