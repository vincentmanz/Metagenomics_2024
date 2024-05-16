# Mofa - on Trypanomsoma data


Separate the dat in 3 groups, control, *T. cruzi* and *T. rangeli*.
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
merged_metagenomes <- import_biom("READBASED/merge_species.biom")

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

# Define a list of contaminants to be removed (e.g., human DNA with taxonomic ID '9606')
contaminants <- c("9606")

# Get the names of all taxa in the merged metagenomes object
allTaxa = taxa_names(merged_metagenomes)

# Remove the contaminants from the list of all taxa
allTaxa <- allTaxa[!(allTaxa %in% contaminants)]

# Prune the phyloseq object to exclude the contaminants
merged_metagenomes = prune_taxa(allTaxa, merged_metagenomes)
```



# Fromating the data

We concider that we have 3 differents experiments, Control, T. cruzi and T. rangeli. We will need to separate these experiments into 3 blocks of data. 


```r
# Subset the merged metagenomes to only include morning samples (AM)
AM_metagenomes <- phyloseq::subset_samples(merged_metagenomes, Morning.Afternoon == "AM")

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
block_control_aggre_df <- psmelt(block_control_aggre) %>% select(Time, OTU, Type, Abundance)
block_cruzi_aggre_df <- psmelt(block_cruzi_aggre) %>% select(Time, OTU, Type, Abundance)
block_rangeli_aggre_df <- psmelt(block_rangeli_aggre) %>% select(Time, OTU, Type, Abundance)

# Rename the columns for consistency
colnames(block_control_aggre_df) <- c("sample", "feature", "view", "value")
colnames(block_cruzi_aggre_df) <- c("sample", "feature", "view", "value")
colnames(block_rangeli_aggre_df) <- c("sample", "feature", "view", "value")
```

The data display real-valued distributions that are appropiately modelled using the gaussian likelihood. The normalization is well done. 


```r
ggdensity(rbind(block_control_aggre_df, block_cruzi_aggre_df,block_rangeli_aggre_df), x="value", fill="gray70") +
  facet_wrap(~view, nrow=1, scales="free")
```

CLR nomalization
![prevalence](Day_3/img/density_clr.png)

Z nomalization
![prevalence](Day_3/img/density_Z.png)


Back to the data transformation.

```r
# Pivot the data from long to wide format, remove the 'view' column, and set the 'feature' column as row names
block_control_aggre_df <- pivot_wider(data = block_control_aggre_df, names_from = sample, values_from = value) %>% select(-view) %>% tibble::column_to_rownames(var = "feature")
block_cruzi_aggre_df <- pivot_wider(data = block_cruzi_aggre_df, names_from = sample, values_from = value) %>% select(-view) %>% tibble::column_to_rownames(var = "feature")
block_rangeli_aggre_df <- pivot_wider(data = block_rangeli_aggre_df, names_from = sample, values_from = value) %>% select(-view) %>% tibble::column_to_rownames(var = "feature")

# Get the row names for each data frame
row_control <- row.names(block_control_aggre_df)
row_cruzi <- row.names(block_cruzi_aggre_df)
row_rangeli <- row.names(block_rangeli_aggre_df)

# Convert the data frames to numeric matrices
block_control_aggre_df_n <- sapply(block_control_aggre_df, as.numeric)
block_cruzi_aggre_df_n <- sapply(block_cruzi_aggre_df, as.numeric)
block_rangeli_aggre_df_n <- sapply(block_rangeli_aggre_df, as.numeric)

# Reassign the row names to the numeric matrices
row.names(block_control_aggre_df_n) <- row_control
row.names(block_cruzi_aggre_df_n) <- row_cruzi
row.names(block_rangeli_aggre_df_n) <- row_rangeli

# Order the columns of each numeric matrix alphabetically by column names
block_control_aggre_df_n <- block_control_aggre_df_n[, order(colnames(block_control_aggre_df_n))]
block_cruzi_aggre_df_n <- block_cruzi_aggre_df_n[, order(colnames(block_cruzi_aggre_df_n))]
block_rangeli_aggre_df_n <- block_rangeli_aggre_df_n[, order(colnames(block_rangeli_aggre_df_n))]

# Combine the numeric matrices into a list and name the elements accordingly
matrix_list <- list(block_control_aggre_df_n, block_cruzi_aggre_df_n, block_rangeli_aggre_df_n)
names(matrix_list) <- c("Control", "T.cruzi", "T.rangeli")
```


# Integration in Mofa

## Create MOFA object

```r
# Create a MOFA (Multi-Omics Factor Analysis) object from the list of matrices
mofa <- create_mofa_from_matrix(data = matrix_list)
```

Visualise data structure, sanity check. 

```r
# Plot an overview of the data in the MOFA object
plot_data_overview(mofa)
```
![data_overview_mofa)](Day_3/img/data_overview_mofa.png)


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

![variance_explained1)](Day_3/img/variance_explained1.png)

Tha variance is well represented accross the data blocks. 

```r
plot_variance_explained(mofa, max_r2=5)
```

![variance_explained2)](Day_3/img/variance_explained2.png)

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

![factors)](Day_3/img/factors.png)


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
combined_plot <- plot_grid(plot_control, plot_cruzi, plot_rangeli, ncol = 3)

# Print the combined plot
print(combined_plot)
```

![Weights)](Day_3/img/Weights.png)


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

![heatmap)](Day_3/img/heatmap.png)
```

