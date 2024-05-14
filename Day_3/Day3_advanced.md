# Tittle



#https://microbiome.github.io/OMA/docs/devel/pages/30_differential_abundance.html
#https://microbiome.github.io/OMA/docs/devel/pages/60_network_learning.html

Separate the dat in 3 groups, control, *T. cruzi* and *T. rangeli*.
Load the data. 
```{r}
merged_metagenomes <- import_biom("READBASED/merge_species.biom") 
meta <- read.csv(file = "DATA/tryp_metadata.csv", sep = ",")


meta <- meta  %>%  arrange(row_number(SRA.identifier)) # sort the data frame
merged_metagenomes@sam_data <- sample_data(meta) # associate the metadata to the to the phyloseq object
column_name <-  meta %>% pull(Sample) # extract the sample names
sample_names(merged_metagenomes) <- column_name # associate the sample names to the phyloseq object

merged_metagenomes@tax_table@.Data <- substring(merged_metagenomes@tax_table@.Data, 4) # remove the unnecessary 'k_' in the taxonomy.
colnames(merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") # change the rank name 

# remove contaminants
contaminants <- c("9606") 
allTaxa = taxa_names(merged_metagenomes)
allTaxa <- allTaxa[!(allTaxa %in% contaminants)]
merged_metagenomes = prune_taxa(allTaxa, merged_metagenomes)
```


```{r}
t_cruzi <- subset_samples(merged_metagenomes, Type== "T._cruzi")
t_rangeli <- subset_samples(merged_metagenomes, Type== "T._rangeli")
control <- subset_samples(merged_metagenomes, Type== "Control")



```



# MOFA


```{r}
library(dplyr)

pseq <- aggregate_rare(merged_metagenomes, level = "Genus", detection = 0.1/100, prevalence = 20/100)

# Transforming microbiome data with rclr
pseq <- microbiome::transform(pseq, transform = "clr")

mofa_df <- (psmelt(pseq))

# Number of features per Type

unique(mofa_df$Type)

mofa_df %>%
  group_by(Type) %>%
  summarize(unique_count = n_distinct(Genus))

```


```{r}
library(ggpubr)

ggdensity(dt, x="value", fill="gray70") +
  facet_wrap(~view, nrow=1, scales="free")
```

![prevalence](Day_3/img/density.png)


```{r}
library(MOFA2)


mofa_df <- mofa_df%>%select("Sample", "OTU", "Abundance", "Type", "Time", "Morning.Afternoon")


colnames(mofa_df) <- c( "sample", "feature", "value", "view", "Time", "Morning.Afternoon")
mofa <- create_mofa_from_df(mofa_df,  extract_metadata = TRUE)
mofa
```


library(MOFA2)

pserq_AM <- phyloseq::subset_samples(merged_metagenomes, Morning.Afternoon == "AM")
pseq <- aggregate_rare(pserq_AM, level = "Genus", detection = 0.1/100, prevalence = 20/100)

# Transforming microbiome data with rclr
pseq <- microbiome::transform(pseq, transform = "clr")

mofa_df <- (psmelt(pseq))
head(mofa_df)
# Create a new column 'Time_Type' by concatenating 'Time' and 'Type' columns
mofa_df$Time_Type <- paste(mofa_df$Time, mofa_df$Morning.Afternoon, sep = "_")

mofa_df <- mofa_df%>%select("Time", "OTU", "Abundance", "Type")



colnames(mofa_df) <- c( "sample", "feature", "value", "view")
head(mofa_df)
mofa <- create_mofa_from_df(mofa_df)
mofa
plot_data_overview(mofa)


model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- 5

mofa <- prepare_mofa(mofa, model_options = model_opts)

mofa <- run_mofa(mofa, use_basilisk = TRUE)