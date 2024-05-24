# Mofa - on Trypanosoma data


Separate the data in 3 groups, control, *T. cruzi* and *T. rangeli*.
Load the data. 
```{r}
library(microbiome)
library(dplyr)
library(ggpubr)
library(MOFA2)
library(cowplot)
library(pheatmap)
library(grid)
library(gridExtra)
```

Load, formating and exclude contaminants from the data. 

```r
# Import the merged metagenomes data from a BIOM file
merged_metagenomes <- import_biom("DATA/merge_species.biom")

# Read the metadata from a CSV file
meta <- read.csv(file = "DATA/tryp_metadata.csv", sep = ",")

# Sort the metadata data frame by the SRA.identifier column
meta <- meta %>% arrange(row_number(SRA.identifier))

# Associate the sorted metadata to the phyloseq object as sample data
merged_metagenomes@sam_data <- sample_data(meta)

# Extract the sample names from the metadata
column_name <- meta %>% pull(Sample)

# Assign the extracted sample names to the phyloseq object
sample_names(merged_metagenomes) <- column_name

# Remove the unnecessary 'k_' prefix in the taxonomy data
merged_metagenomes@tax_table@.Data <- substring(merged_metagenomes@tax_table@.Data, 4)

# Rename the columns of the taxonomy table to represent taxonomic ranks
colnames(merged_metagenomes@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Keep the kingdom of interest
merged_metagenomes <- subset_taxa(merged_metagenomes, Kingdom %in% c("Archaea", "Bacteria", "Fungi", "Viruses"))
```



# Formating the data

We concider that we have 3 differents block of data, viruses, fungi and bacteria. We need to separate the data.  


```r
# Subset the merged metagenomes to only include morning samples (AM)
AM_metagenomes <- phyloseq::subset_samples(merged_metagenomes, Gut == "AM")

# Further subset the morning metagenomes by sample type
block_control <- phyloseq::subset_samples(AM_metagenomes, Type == "Control")
block_cruzi <- phyloseq::subset_samples(AM_metagenomes, Type == "T._cruzi")
block_rangeli <- phyloseq::subset_samples(AM_metagenomes, Type == "T._rangeli")

# Aggregate rare taxa at the genus level for each subset, using specific detection and prevalence thresholds
block_control_aggre <- aggregate_rare(block_control, level = "Genus", detection = 0.1 / 100, prevalence = 50 / 100)
block_cruzi_aggre <- aggregate_rare(block_cruzi, level = "Genus", detection = 0.1 / 100, prevalence = 50 / 100)
block_rangeli_aggre <- aggregate_rare(block_rangeli, level = "Genus", detection = 0.1 / 100, prevalence = 50 / 100)
```
## Normalization

You can use different normalisation and plot the density of the data to check the distribution of the data. 

```r
# Transform the aggregated data using centered log-ratio (clr) transformation
block_control_aggre <- microbiome::transform(block_control_aggre, transform = "clr") 
block_cruzi_aggre <- microbiome::transform(block_cruzi_aggre, transform = "clr")
block_rangeli_aggre <- microbiome::transform(block_rangeli_aggre, transform = "clr")
```


```r
# Melt the transformed data into long format and select relevant columns
block_control_aggre_df <- psmelt(block_control_aggre) %>% select(Time, OTU, Abundance) %>%  mutate(view = "Control")
block_cruzi_aggre_df <- psmelt(block_cruzi_aggre)  %>% select(Time, OTU, Abundance) %>%  mutate(view = "T. cruzi")
block_rangeli_aggre_df <- psmelt(block_rangeli_aggre) %>% select(Time, OTU, Abundance) %>%  mutate(view = "T. rangeli")

# Rename the columns for consistency
colnames(block_control_aggre_df) <-  c("sample", "feature", "value", "view")
colnames(block_cruzi_aggre_df) <- c("sample", "feature", "value", "view")
colnames(block_rangeli_aggre_df) <- c("sample", "feature", "value", "view")

head(block_control_aggre_df)
head(block_cruzi_aggre_df)
head(block_rangeli_aggre_df)


merged_df <- bind_rows(block_cruzi_aggre_df, block_rangeli_aggre_df, block_control_aggre_df)
head(merged_df)

```

The data display real-valued distributions that are appropiately modelled using the gaussian likelihood. The normalization is well done. 


```r
ggdensity(rbind(block_control_aggre_df, block_cruzi_aggre_df,block_rangeli_aggre_df), x="value", fill="gray70") +
  facet_wrap(~view, nrow=1, scales="free")
```

CLR nomalization
![prevalence](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_3/img/density_clr.png)

Z nomalization
![prevalence](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_3/img/density_Z.png)


# Integration in Mofa

## Create MOFA object

```r
# Create a MOFA (Multi-Omics Factor Analysis) object from the list of matrices
mofa <- create_mofa(data = matrix_list)
```

Visualise data structure, sanity check. 

```r
# Plot an overview of the data in the MOFA object
plot_data_overview(mofa)
```
![data_overview_mofa)](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_3/img/data_overview_mofa.png)


In that case we do not have missing data, if it was the case, you will have grey line on the graph. It helps to see if the data are correctly aligned. 

## Prepare MOFA object

```r
# Get the default model options for the MOFA object
model_opts <- get_default_model_options(mofa)

# Set the number of factors to be inferred in the model to 5
model_opts$num_factors <- 5

# Prepare the MOFA object with the specified model options
mofa <- prepare_mofa(mofa, model_options = model_opts)

# Run the MOFA model using the Basilisk environment (for reproducibility and isolation of dependencies)
mofa <- run_mofa(mofa, use_basilisk = TRUE)
```

Mofa will complain about the sample size. We have only 5 samples, we should have at list 20.

# Downstream analysis
## Variance decomposition


```r
plot_variance_explained(mofa, plot_total = T)[[2]]
```

![variance_explained1)](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_3/img/variance_explained1.png)

Tha variance is well represented accross the data blocks. 

```r
plot_variance_explained(mofa, max_r2=5)
```

![variance_explained2)](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_3/img/variance_explained2.png)

We have only 2 factors but they capture coordinated variability across all data blocks. 


# Visualise the MOFA Factors



```r
time.colors <- c(
  "T0" = "#66C2A5", 
  "T1" = "#8DA0CB",
  "T2" = "#E78AC3",
  "T3" = "#FC8D62",
  "T7" = "#D8EC82"
)
plot_factors(mofa, 
             factors = c(1,2), 
             color_by = c("T0", "T1", "T2", "T3", "T7"), 
             #shape_by="samples",
             dot_size = 4
) + scale_fill_manual(values=time.colors)
```

![factors)](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_3/img/factors.png)


We can see that Factor 1 is well explaining the variation through time. 



## Plot Weights 

Helper function. 

```r
plot_weights_fn <- function(mofa, factor=1, view=1, nfeatures=10) {
  p1 <- plot_weights(mofa, 
                     factors = factor, 
                     view = view,
                     nfeatures = nfeatures,
                     text_size = 4
  )
  
  p2 <- plot_top_weights(mofa, 
                         factors = factor, 
                         view = view,
                         nfeatures = nfeatures
  )
  
  p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow=1)
  return(p)
}
```

```r
# Generate the plots
plot_control <- plot_weights_fn(mofa, factor = 1, view = "Control", nfeatures = 8)
plot_cruzi <- plot_weights_fn(mofa, factor = 1, view = "T.cruzi", nfeatures = 8)
plot_rangeli <- plot_weights_fn(mofa, factor = 1, view = "T.rangeli", nfeatures = 8)

# Combine the plots in a single line
combined_plot <- plot_grid(plot_control, plot_cruzi, plot_rangeli, nrow = 3)

# Print the combined plot
print(combined_plot)
```

![Weights)](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_3/img/Weights.png)


## Plot feature values

Heatmaps of bacterial concentration, no denoising:


```r
# Generate the heatmap plots
plot_control <- plot_data_heatmap(mofa, factor = 1, view = "Control", features = 20, denoise = TRUE, cluster_rows = TRUE, cluster_cols = FALSE, show_colnames = TRUE, show_rownames = TRUE, scale = "row")
plot_cruzi <- plot_data_heatmap(mofa, factor = 1, view = "T.cruzi", features = 20, denoise = TRUE, cluster_rows = TRUE, cluster_cols = FALSE, show_colnames = TRUE, show_rownames = TRUE, scale = "row")
plot_rangeli <- plot_data_heatmap(mofa, factor = 1, view = "T.rangeli", features = 20, denoise = TRUE, cluster_rows = TRUE, cluster_cols = FALSE, show_colnames = TRUE, show_rownames = TRUE, scale = "row")

# Arrange the heatmaps in a single column
grid.arrange(plot_control[[4]], plot_cruzi[[4]], plot_rangeli[[4]], ncol = 1)
```

![heatmap)](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_3/img/heatmap.png)


**Enterococcus** (Bacillota), which plays a crucial role in metabolic adaptability against pathogenic or plant toxins and anti-herbivore defense, was found to be one of the predominant gut microorganism of lepidopteran insects, including B. mori, Helicoverpa zea, and Porthetria dispar (Paniagua Voirol et al., 2018; Zhang et al., 2022). 

**Symbiopectobacterium** (Enterobacteriaceae) has recently been described for the first time as an intracellular bacterial symbiont, responsible for to the biosynthesis of vitamins and cofactors.  This bacteria may be boosting the parsite fitness, for example, by aiding in evading the Triatome immune response, or providing a novel function, such as supplementing nutrition or metabolism

**Rhodococcus** (Nocardiaceae) in the triatomine gut are believed to play an important role in the metabolism of the vector, such as by participating in the synthesis of group B vitamins or by being digested by the bugs directly to provide missing nutrients (Sassera et al., 2013). Moreover, the most attractive aspect is the host-symbiont relationship between triatomines and Rhodococcus; since Rhodococcus bacteria can be easily cultured and genetically modified to harm the pathogen in vector gut, they are probably suitable tools for the control of trypanosomiasis (Sassera et al., 2013). as the blood is poor in B vitamins compared to what is generally required for insect development. The blood is poor in B vitamins compared to what is generally required for insect development, Kissing bugs, Rhodnius prolixus, notably require *Rhodococcus* bacteria for nymph development, but the addition of B vitamins in the diet can rescue nymph development in the absence of Rhodococcus  (Serrato-Salas and Gendrin 2023).

**Wolbachia** (Ehrlichiaceae)  The obligate intracellular bacteria Wolbachia spp. are common in a wide range of insects, including sand flies, bed bugs, fleas and mosquitoes, and can cause reproduction alterations such as feminization, male killing and cytoplasmic incompatibility ([Landmann 20219](https://doi.org/10.1128/microbiolspec.BAI-0018-2019.)). In triatomines, Wolbachia has been solely reported for the genus Rhodnius, where it occurs in the intestine, salivary glands and gonads.

**Curtobacterium** (Microbacteriaceae) *C. flaccumfaciens* is the only species of Curtobacterium associated with plant pathogenesis (Young et al., 1996), the presence of *C. flaccumfaciens* in the rhizosphere induced a systematic resistance in cucumber plants to pathogens.





# Mofa Integrative by kingdom


```r
library(data.table)
library(purrr)
library(ggplot2)
library(ggpubr)
library(MOFA2)
```

```r

setwd("~/Documents/project/Metagenomics_2024/Metagenomics_2024")

merged_metagenomes <- import_biom("DATA/merge_species.biom") 
meta <- read.csv(file = "DATA/tryp_metadata.csv", sep = ",")

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

merged_metagenomes <- subset_taxa(merged_metagenomes, Kingdom %in% c("Archaea", "Bacteria", "Fungi", "Viruses"))
```




```r
# Subset the merged metagenomes to only include morning samples (AM)
#AM_metagenomes <- phyloseq::subset_samples(merged_metagenomes, Gut == "AM")
AM_metagenomes <- merged_metagenomes
# Further subset the morning metagenomes by sample type
block_bacteria <- phyloseq::subset_taxa(AM_metagenomes, Kingdom == "Bacteria")
#block_fungi <- phyloseq::subset_taxa(AM_metagenomes, Kingdom == "Fungi")
block_virus <- phyloseq::subset_taxa(AM_metagenomes, Kingdom == "Viruses")

# Aggregate rare taxa at the genus level for each subset, using specific detection and prevalence thresholds
block_bacteria <- aggregate_rare(block_bacteria, level = "Genus", detection = 0.1 / 100, prevalence = 50 / 100)
#block_fungi <- aggregate_rare(block_fungi, level = "Genus", detection = 0.1 / 100, prevalence = 50 / 100)
block_virus <- aggregate_rare(block_virus, level = "Genus", detection = 0.1 / 100, prevalence = 50 / 100)
```



```r
# Further subset the morning metagenomes by sample type
block_bacteria <- phyloseq::subset_taxa(merged_metagenomes, Kingdom == "Bacteria")
#block_fungi <- phyloseq::subset_taxa(merged_metagenomes, Kingdom == "Fungi")
block_virus <- phyloseq::subset_taxa(merged_metagenomes, Kingdom == "Viruses")

# Aggregate rare taxa at the genus level for each subset, using specific detection and prevalence thresholds
block_bacteria <- aggregate_rare(block_bacteria, level = "Genus", detection = 0.05 / 100, prevalence = 20 / 100)
#block_fungi <- aggregate_rare(block_fungi, level = "Genus", detection = 0.1 / 100, prevalence = 50 / 100)
block_virus <- aggregate_rare(block_virus, level = "Genus", detection = 0.05 / 100, prevalence = 20 / 100)

# Transform the aggregated data using centered log-ratio (clr) transformation
block_bacteria <- microbiome::transform(block_bacteria, transform = "clr")
block_virus <- microbiome::transform(block_virus, transform = "clr")

# Melt the transformed data into long format and select relevant columns
block_bacteria_df <- psmelt(block_bacteria) %>% select(Sample, OTU, Abundance) %>%  mutate(view = "Bacteria")
block_virus_df <- psmelt(block_virus) %>% select(Sample, OTU, Abundance) %>% mutate(view = "Virus")

# Rename the columns for consistency
colnames(block_bacteria_df) <- c("sample", "feature", "value", "view")
colnames(block_virus_df) <- c("sample", "feature", "value", "view")

#head(block_bacteria_df)
#head(block_virus_df)

merged_df <- bind_rows(block_bacteria_df, block_virus_df)
```


# Density plot abundance

```r
ggdensity(rbind(block_control_aggre_df, block_cruzi_aggre_df,block_rangeli_aggre_df), x="value", fill="gray70") +
  facet_wrap(~view, nrow=1, scales="free")
```

```r
# Create a MOFA (Multi-Omics Factor Analysis) object from the list of matrices
mofa <- create_mofa(data = merged_df)
```


## Prepare MOFA object

```r
# Get the default model options for the MOFA object
model_opts <- get_default_model_options(mofa)

# Set the number of factors to be inferred in the model to 5
model_opts$num_factors <- 10

# Prepare the MOFA object with the specified model options
mofa <- prepare_mofa(mofa, model_options = model_opts)

# Run the MOFA model using the Basilisk environment (for reproducibility and isolation of dependencies)
mofa <- run_mofa(mofa, use_basilisk = TRUE)
```
# Add sample metadata to the model

```r
colnames(meta) <- c("sample", "SRA.identifier", "Type", "Time", "Number.of.days", "Gut", "Reads")
samples_metadata(mofa) <- meta
```

Mofa will complain about the sample size. We have only 5 samples, we should have at list 20.

# Downstream analysis

```r
plot_variance_explained(mofa, plot_total = T)[[2]]


plot_variance_explained(mofa, max_r2=15)
```

```r
gut.colors <- c(
  "AM" = "#66C2A5", 
  "PM" = "#8DA0CB"
  )
time.colors <- c(
  "T0" = "#66C2A5", 
  "T1" = "#8DA0CB",
  "T2" = "#E78AC3",
  "T3" = "#FC8D62",
  "T7" = "#D8EC82"
  )
type.colors <- c(
  "Control" = "#66C2A5",
  "T._rangeli" = "#FC8D62",
  "T._cruzi" = "#D8EC82"
  )
```


```r
plot_data_heatmap(mofa, 
                  factor = 1, 
                  view = "Bacteria", 
                  features = 20,
                  denoise = TRUE,
                  cluster_rows = T, cluster_cols = F,
                  show_colnames =T, show_rownames = T,
                  annotation_samples = c("Gut", "Time", "Type"),  
                  annotation_colors = list("Time"=time.colors), 
                  annotation_legend = TRUE,
                  scale = "row"
)
```

```r
p <- plot_factors(mofa, 
                  factors = c(1,2), 
                  color_by = "Time", 
                  dot_size = 4
) + scale_fill_manual(values=time.colors)
p +
  # geom_density_2d(aes_string(color="color_by")) +
  stat_ellipse(aes(color=color_by), geom = "polygon", alpha=0.25) +
  scale_color_manual(values=time.colors)
```



```r
plot_factor(mofa,
            factor = 1,
            color_by = "Time", 
            dot_size = 5,
            scale = TRUE, 
            legend = TRUE, 
            dodge = TRUE,
) +
  scale_color_manual(values=time.colors) + 
  scale_fill_manual(values=time.colors)
```

![factor)](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_3/img/factor_value.png)

```r
plot_data_scatter(mofa, 
                  factor = 1, 
                  view = "Bacteria", 
                  features = 6,
                  dot_size = 3,
                  color_by = "Time",
                  legend = T
)
```

![factor)](https://github.com/vincentmanz/Metagenomics_2024/blob/main/Day_3/img/correlation.png)
