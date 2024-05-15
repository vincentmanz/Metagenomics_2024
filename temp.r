
library(microbiome)
library(dplyr)
library(MOFA2)

merged_metagenomes <- import_biom("READBASED/merge_species.biom") 
meta <- read.csv(file = "DATA/tryp_metadata.csv", sep = ",")


meta <- meta  %>%  arrange(row_number(SRA.identifier)) # sort the data frame
merged_metagenomes@sam_data <- sample_data(meta) # associate the metadata to the to the phyloseq object
column_name <-  meta %>% pull(Sample) # extract the sample names
sample_names(merged_metagenomes) <- column_name # associate the sample names to the phyloseq object

merged_metagenomes@tax_table@.Data <- substring(merged_metagenomes@tax_table@.Data, 4) # remove the unnecessary 'k_' in the taxonomy.
colnames(merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") # change the rank name 

contaminants <- c("9606") 
allTaxa = taxa_names(merged_metagenomes)
allTaxa <- allTaxa[!(allTaxa %in% contaminants)]
merged_metagenomes = prune_taxa(allTaxa, merged_metagenomes)

AM_metagenomes <- phyloseq::subset_samples(merged_metagenomes, Morning.Afternoon=="AM")

block_control <- phyloseq::subset_samples(AM_metagenomes, Type=="Control")
block_cruzi <- phyloseq::subset_samples(AM_metagenomes, Type=="T._cruzi")
block_rangeli <- phyloseq::subset_samples(AM_metagenomes, Type=="T._rangeli")


block_control_aggre <- aggregate_rare(block_control, level = "Genus", detection = 0.1/100, prevalence = 50/100)
block_cruzi_aggre <- aggregate_rare(block_cruzi, level = "Genus", detection = 0.1/100, prevalence = 50/100)
block_rangeli_aggre <- aggregate_rare(block_rangeli, level = "Genus", detection = 0.1/100, prevalence = 50/100)


block_control_aggre <- microbiome::transform(block_control_aggre, transform = "clr")
block_cruzi_aggre <-  microbiome::transform(block_cruzi_aggre, transform = "clr")
block_rangeli_aggre <- microbiome::transform(block_rangeli_aggre, transform = "clr")


block_control_aggre_df <- psmelt(block_control_aggre) %>% select(Time, OTU, Type, Abundance ) 
block_cruzi_aggre_df <- psmelt(block_cruzi_aggre)  %>% select(Time, OTU, Type, Abundance )
block_rangeli_aggre_df <- psmelt(block_rangeli_aggre)  %>% select(Time, OTU, Type, Abundance )

colnames(block_control_aggre_df) <- c("sample", "feature", "view", "value")
colnames(block_cruzi_aggre_df) <- c("sample", "feature", "view", "value")
colnames(block_rangeli_aggre_df) <- c("sample", "feature", "view", "value")

block_control_aggre_df  <- pivot_wider(data = block_control_aggre_df, names_from = sample, values_from = value) %>% select(-view) %>% tibble::column_to_rownames(var = "feature") 
block_cruzi_aggre_df <- pivot_wider(data = block_cruzi_aggre_df, names_from = sample, values_from = value) %>% select(-view) %>% tibble::column_to_rownames(var = "feature")
block_rangeli_aggre_df <- pivot_wider(data = block_rangeli_aggre_df, names_from = sample, values_from = value) %>% select(-view) %>% tibble::column_to_rownames(var = "feature")

row_control <- row.names(block_control_aggre_df)
row_cruzi <- row.names(block_cruzi_aggre_df)
row_rangeli <- row.names(block_rangeli_aggre_df)



block_control_aggre_df_n <- sapply(block_control_aggre_df,  as.numeric)
block_cruzi_aggre_df_n <- sapply(block_cruzi_aggre_df,  as.numeric)
block_rangeli_aggre_df_n <- sapply(block_rangeli_aggre_df,  as.numeric)

row.names(block_control_aggre_df_n) <- row_control
row.names(block_cruzi_aggre_df_n) <- row_cruzi
row.names(block_rangeli_aggre_df_n) <- row_rangeli

block_control_aggre_df_n <- block_control_aggre_df_n[, order(colnames(block_control_aggre_df_n))]
block_cruzi_aggre_df_n <- block_cruzi_aggre_df_n[, order(colnames(block_cruzi_aggre_df_n))]
block_rangeli_aggre_df_n <- block_rangeli_aggre_df_n[, order(colnames(block_rangeli_aggre_df_n))]


matrix_list <- list(block_control_aggre_df_n, block_cruzi_aggre_df_n, block_rangeli_aggre_df_n)
names(matrix_list) <- c("Control", "T.cruzi", "T.rangeli")




mofa <- create_mofa_from_matrix(data = matrix_list)

plot_data_overview(mofa)
model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- 5

mofa <- prepare_mofa(mofa, model_options = model_opts)
mofa <- run_mofa(mofa, use_basilisk = T)






plot_variance_explained(mofa, plot_total = T)[[2]]

plot_variance_explained(mofa, max_r2=5)

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

plot_weights_fn(mofa, factor=1, view="Control", nfeatures=8)
plot_weights_fn(mofa, factor=1, view="T.cruzi", nfeatures=8)
plot_weights_fn(mofa, factor=1, view="T.rangeli", nfeatures=8)

plot_data_heatmap(mofa, 
                  factor = 2, 
                  view = "Control", 
                  features = 20,
                  denoise = T,
                  cluster_rows = T, cluster_cols = F,
                  show_colnames = T, show_rownames = T,
                  scale = "row"
)

plot_data_heatmap(mofa, 
                  factor = 1, 
                  view = "T.cruzi", 
                  features = 20,
                  denoise = T,
                  cluster_rows = T, cluster_cols = F,
                  show_colnames = T, show_rownames = T,
                  scale = "row"
)

plot_data_heatmap(mofa, 
                  factor = 1, 
                  view = "T.rangeli", 
                  features = 20,
                  denoise = T,
                  cluster_rows = T, cluster_cols = F,
                  show_colnames = T, show_rownames = T,
                  scale = "row"
)



get_weights(mofa)


NEXT DO THE SAME THING WITH NORMAL MICROBIAL STATS:COMPARE THE 3 SET FOR THE SAME TIME. 
