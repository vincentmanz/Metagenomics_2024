library(microbiome)
library(dplyr)
library(hrbrthemes)
library(tibble)
library(tidyverse)
library(vegan)




merged_metagenomes <- import_biom("/home/vincent/Documents/project/Metagenomics_2024/data/Trypanosoma_exposure/READBASED/merge_species.biom") 
meta <- read.csv(file = "DATA/tryp_metadata.csv", sep = ",")


meta <- meta  %>%  arrange(row_number(SRA.identifier)) # sort the data frame
merged_metagenomes@sam_data <- sample_data(meta) # associate the metadata to the to the phyloseq object
column_name <-  meta %>% pull(Sample) # extract the sample names
sample_names(merged_metagenomes) <- column_name # associate the sample names to the phyloseq object

merged_metagenomes@tax_table@.Data <- substring(merged_metagenomes@tax_table@.Data, 4) # remove the unnecessary 'k_' in the taxonomy.
colnames(merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") # change the rank name 


merged_metagenomes <- subset_taxa(merged_metagenomes, Kingdom %in% c("Archaea", "Bacteria", "Fungi", "Viruses"))
 

pseq <- aggregate_rare(merged_metagenomes, level = "Phylum", detection = 0.1/100, prevalence = 50/100)

pseq <- microbiome::transform(pseq, transform = "compositional")
dim(tax_table(pseq))[1]

# Summarize the phyloseq object 'merged_metagenomes'
summarize_phyloseq(pseq)

# Display the first few rows of the OTU (Operational Taxonomic Unit) table
head(otu_table(pseq))

# Display the first few rows of the taxonomy table
head(tax_table(merged_metagenomes))

# Display the first few rows of the sample data associated with 'merged_metagenomes'
head(sample_data(merged_metagenomes))

# Get the sample variables of the phyloseq object 'merged_metagenomes'
sample_variables(merged_metagenomes)

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
head(psmelt(merged_metagenomes), n=100)

a <- psmelt(merged_metagenomes)

d<-psmelt(subset_taxa(merged_metagenomes, Genus == "Homo"))




non_eukaryota_taxa <- taxa_names(merged_metagenomes)[!tax_table(merged_metagenomes)[, "Domain"] %in% "Eukaryota"]
b <- subset_taxa(merged_metagenomes, Kingdom %in% c("Archaea", "Bacteria", "Fungi", "Viruses"))
tax_table(b)
# remove contaminants
contaminants <- c("9606") 
allTaxa = taxa_names(merged_metagenomes)
allTaxa <- allTaxa[!(allTaxa %in% contaminants)]
merged_metagenomes = prune_taxa(allTaxa, merged_metagenomes)


merged_metagenomes_family <- aggregate_rare(merged_metagenomes, level = "Family", detection = 0/100, prevalence = 0/100)

dim(psmelt(merged_metagenomes))

head(otu_table(merged_metagenomes))

plot_read_distribution(merged_metagenomes, "Type", plot.type = "density")



df <- psmelt(merged_metagenomes)  %>%  group_by(Sample, Time) %>%  
  summarise(sum_reads = sum(Abundance)) %>% arrange(sum_reads) 

ggplot(df) +
  geom_bar(aes(reorder(Sample, -sum_reads), sum_reads, fill=Time),
           col="red", alpha = .2, stat="identity") 





prevalence(merged_metagenomes, detection=1/100, sort=TRUE, count=FALSE)


prevalence(merged_metagenomes, detection=1, sort=TRUE, count=T)

plot_taxa_prevalence(prevalence(merged_metagenomes, detection=1, sort=TRUE, count=T))

plot_taxa_prevalence(merged_metagenomes, level="Family", detection = 10000)

# alternative 
p1 <- plot_taxa_cv(merged_metagenomes, plot.type = "scatter")
p1 + scale_x_log10(), level="Phylum", detection = 10000)

# alternative 
p1 <- plot_taxa_cv(merged_metagenomes, plot.type = "scatter")
p1 + scale_x_log10()

otu_table(pseq)


alpha(pseq, index = "diversity_fisher", zeroes = F)




























































subset_samples(merged_metagenomes, Type== "Control")






pseq <- aggregate_rare(merged_metagenomes, level = "Genus", detection = 0.1/100, prevalence = 50/100)

otu <- abundances(pseq)
meta <- meta(pseq)

permanova <- adonis(t(otu) ~ Time,
                    data = meta, permutations=9999, method = "euclidean")

print(as.data.frame(permanova$aov.tab))


coef <- coefficients(permanova)["Time1",]
top.coef <- coef[rev(order(abs(coef)))[1:20]]
par(mar = c(3, 14, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")



pseq <- microbiome::transform(pseq, transform = "compositional")
head(psmelt(pseq))

p <- plot_composition(pseq,
                      taxonomic.level = "Family",
                      sample.sort = "Sample",
                      x.label = "Sample",
                      group_by = "Type") +
  scale_fill_brewer("Family", palette = "Paired") +
  guides(fill = guide_legend(ncol = 1)) +
  scale_y_percent() +
  labs(x = "Samples", y = "Relative abundance (%)",
       title = "Relative abundance data") + 
  theme_ipsum(grid="Y") +
  theme(axis.text.x = element_text(angle=90, hjust=1),
        legend.text = element_text(face = "italic"))
print(p)  


dim(abundances(merged_metagenomes))
head(otu_table(merged_metagenomes) )


df <- psmelt(merged_metagenomes)  %>%  group_by(Sample, Time) %>%  
  summarise(sum_reads = sum(Abundance)) %>% arrange(sum_reads) 

ggplot(df) +
  geom_bar(aes(reorder(Sample, -sum_reads), sum_reads, fill=Time),
           col="red", alpha = .2, stat="identity") 
pseq <- microbiome::transform(merged_metagenomes, transform = "compositional")
otu_table(pseq)

library("FactoMineR")
library("factoextra")
?microbiome::transform

res.pca <- PCA(df, graph = FALSE)

get_eigenvalue(res.pca)











pseq <- aggregate_rare(merged_metagenomes, level = "Genus", detection = 0.1/100, prevalence = 50/100)
pseq <- microbiome::transform(pseq, transform = "clr")

heatmap_df <- psmelt(pseq)

tse_phylum_subset <- tse_phylum_subset[top_taxa, ]

# Gets the assay table
mat <- assay(tse_phylum_subset, "clr_z")



heatmap_df <- psmelt(pseq)

genus_level_abundance <- heatmap_df  %>% as_tibble() %>%  column_to_rownames( "Sample")

pheatmap::pheatmap(heatmap_df,
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   annotation_col = meta_h[,c(1,6,7)],
                   annotation_names_col=TRUE)


















pseq <- aggregate_rare(merged_metagenomes, level = "Genus", detection = 0.1/100, prevalence = 50/100)

pseq <- microbiome::transform(pseq, transform = "clr")


psmelt(merged_metagenomes)

library(mia)
ancombc2_out <- ancombc2(merged_metagenomes,
                         assay.type = "counts",
                         fix_formula = "Abundance + Type + Time + Morning.Afternoon",
                         p_adj_method = "fdr",
                         lib_cut = 0,
                         group = "Type", 
                         struc_zero = TRUE, 
                         neg_lb = TRUE,
                         alpha = 0.05,
                         # multi-group comparison is deactivated automatically
                         global = TRUE)















ntaxa(merged_metagenomes)
#
get_taxa_unique(merged_metagenomes, "Family")



merged_metagenomes_family <- aggregate_rare(merged_metagenomes, level = "Family", detection = 0/100, prevalence = 0/100) # no filter except family


plot_richness(merged_metagenomes) 


plot_taxa_prevalence(merged_metagenomes, level="Phylum", detection = 1000)

plot_read_distribution


plot_read_distribution(merged_metagenomes, "Type", plot.type = "density")



p1 <- plot_taxa_cv(merged_metagenomes, plot.type = "scatter")

p1 + scale_x_log10()




data.table(tax_table(merged_metagenomes),
           ASVabundance = taxa_sums(ps.ng.tax),
           ASV = taxa_names(ps.ng.tax))



bracken_species <- read.csv(file = "READBASED/BRACKEN/bracken_merged_species.csv", sep = "\t") 
bracken_genus <- read.csv(file = "READBASED/BRACKEN/bracken_merged_genus.csv", sep = "\t") 
bracken_family <- read.csv(file = "READBASED/BRACKEN/bracken_merged_family.csv", sep = "\t") 
bracken_phylum <- read.csv(file = "READBASED/BRACKEN/bracken_merged_phylum.csv", sep = "\t") 

meta <- read.csv(file = "../../Metagenomics_2024/DATA/tryp_metadata.csv", sep = ",")



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






ggplot(phylum_level, aes(x=Sample, y=rel_abun)) +
  geom_bar(stat="identity", position="stack", aes(fill=name)) + # chose bar plot
  theme(axis.text.x = element_text(angle=45, hjust=1)) + # put x-axis label at 45 degree angle
  facet_grid(. ~ Time, scales="free_x",space = "free_x") # produce two panels according to metatadata category 'Time' 
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

ggplot(species_level, aes(x=Sample, y=rel_abun)) +
         geom_bar(stat="identity", position="stack", aes(fill=name)) + # chose bar plot
         theme(axis.text.x = element_text(angle=45, hjust=1)) + # put x-axis label at 45 degree angle
         facet_grid(. ~ Time, scales="free_x",space = "free_x") + # produce two panels according to metatadata category 'Time' 
         guides(fill = FALSE)


  ggplot(species_level, aes(x=Sample, y=name)) +
  geom_point(aes(size=rel_abun, color=Type), alpha=0.7) + # this time we use points
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  scale_size_continuous(limits = c(0.00001,max(species_level$rel_abun))) + # sets minimum above '0'
  facet_grid(. ~ Time, scales="free_x",space = "free_x")



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

pheatmap::pheatmap(species_level_abundance,
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   annotation_col = meta_h[,c(2,3,4,5)],
                   annotation_names_col=TRUE)






# To ensure reproducibility we can fix the seed here. This will ensure you always get the same result each time you run your data.
set.seed(34521)

# Data mingling
species_frac_filtered_f <- species_frac_filtered %>% 
  column_to_rownames("name")
colnames(species_frac_filtered_f) <- column_name
species_frac_filtered_f <- species_frac_filtered_f %>% t() # transpose

# Calculate distance matrix
species_frac_filtered_dist <- vegdist(species_frac_filtered_f, method = "bray")

# Perform NMDS on distance matrix
nmds_spec <- metaMDS(species_frac_filtered_dist,distance = "bray",k = 2, try = 50)

# Extract and reshape the data to plot ordination as ggplot  and add the metadata
nmds_spec_gg<-as.data.frame(nmds_spec$points) %>%
  rownames_to_column("Sample") %>%
  left_join(meta, by="Sample")
# Let's plot and color according to time point
ggplot(nmds_spec_gg, aes(x=MDS1,y=MDS2)) +
  geom_point(aes(color=Type), size=3, alpha=0.5) +
  ggtitle("NMDS colored according to Type")


set.seed(34521) # set seed for reproducibility

# PERMANOVA
adonis2(species_frac_filtered_dist ~ Time, data = as.data.frame(nmds_spec_gg))























# reformat table taxa: https://carpentries-lab.github.io/metagenomics-analysis/07-phyloseq/index.html


merged_metagenomes <- import_biom("READBASED/merge_species.biom")
meta <- read.csv(file = "../../Metagenomics_2024/DATA/tryp_metadata.csv", sep = ",")


#format the data
meta <- meta  %>%  arrange(row_number(SRA.identifier))
merged_metagenomes@sam_data <- sample_data(meta)
column_name <-  meta %>% pull(Sample)
sample_names(merged_metagenomes) <- column_name


merged_metagenomes@tax_table@.Data <- substring(merged_metagenomes@tax_table@.Data, 4)
colnames(merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

head(otu_table(merged_metagenomes))
head(tax_table(merged_metagenomes))
head(sample_data(merged_metagenomes))
sample_variables(merged_metagenomes)

# remove contaminants
contaminants <- c("9606")
allTaxa = taxa_names(merged_metagenomes)
allTaxa <- allTaxa[!(allTaxa %in% contaminants)]
merged_metagenomes = prune_taxa(allTaxa, merged_metagenomes)
head(tax_table(merged_metagenomes))


# Community typing with Dirichlet Multinomial Mixtures

library(microbiome)
library(DirichletMultinomial)
library(reshape2)
library(magrittr)
library(dplyr)

pseq <- aggregate_rare(merged_metagenomes, level = "Genus", detection = 0.1/100, prevalence = 50/100)

# Pick the OTU count matrix
# and convert it into samples x taxa format
dat <- abundances(pseq)
count <- as.matrix(t(dat)) 
# Remove rows with zero sums
count <- count[rowSums(count) != 0, ]
fit <- lapply(1:10, dmn, count = count, verbose=TRUE)



lplc <- base::sapply(fit, DirichletMultinomial::laplace) # AIC / BIC / Laplace
aic  <- base::sapply(fit, DirichletMultinomial::AIC) # AIC / BIC / Laplace
bic  <- base::sapply(fit, DirichletMultinomial::BIC) # AIC / BIC / Laplace
plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
lines(aic, type="b", lty = 2)
lines(bic, type="b", lty = 3)


best <- fit[[which.min(unlist(lplc))]]
mixturewt(best)
ass <- apply(mixture(best), 1, which.max)



for (k in seq(ncol(fitted(best)))) {
  d <- melt(fitted(best))
  colnames(d) <- c("OTU", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange OTUs by assignment strength
    arrange(value) %>%
    mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.8))     
  
  p <- ggplot(d, aes(x = OTU, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top drivers: community type", k))
  print(p)
}




# Composition barplots
Same with compositional (relative) abundances; for each sample (left), or averafged by group (right).

# Try another theme
# from https://github.com/hrbrmstr/hrbrthemes
library(hrbrthemes)
library(gcookbook)
library(tidyverse)
#theme_set(theme_bw(21))
# Limit the analysis on core taxa and specific sample group
pseq <- aggregate_rare(merged_metagenomes, level = "Species", detection = 0.1/100, prevalence = 50/100)

pseq <- microbiome::transform(pseq, transform = "compositional")

abundances(pseq)

p <- plot_composition(pseq,
                      taxonomic.level = "Genus",
                      sample.sort = "Sample",
                      x.label = "Sample") +
  scale_fill_brewer("Genera", palette = "Paired") +
  guides(fill = guide_legend(ncol = 1)) +
  scale_y_percent() +
  labs(x = "Samples", y = "Relative abundance (%)",
       title = "Relative abundance data") + 
  theme_ipsum(grid="Y") +
  theme(axis.text.x = element_text(angle=90, hjust=1),
        legend.text = element_text(face = "italic"))
print(p)  

# Alpha diversity

https://microbiome.github.io/tutorials/Alphadiversity.html
library(microbiome)
library(knitr)

pseq <- aggregate_rare(merged_metagenomes, level = "Genus", detection = 0.1/100, prevalence = 50/100)

tab <-microbiome::alpha(pseq, index = "all")
kable(head(tab))


tab <- richness(pseq)
kable(head(tab))


p.shannon <- boxplot_alpha(pseq, 
                           index = "shannon",
                           x_var = "Type")

p.shannon <- p.shannon + theme_minimal() + 
  labs(x="Sex", y="Shannon diversity") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16))
p.shannon


Investigate the top factors
Show coefficients for the top taxa separating the groups

pseq <- aggregate_rare(merged_metagenomes, level = "Genus", detection = 0.1/100, prevalence = 50/100)
pseq <- microbiome::transform(pseq, transform = "CLR")
abundances(pseq)

p <- plot_landscape(pseq, method = "NMDS", distance = "bray", col = "Time", size = 3)
print(p)


8.3 Estimating associations with an external variable
Next to visualizing whether any variable is associated with differences between samples, we can also quantify the strength of the association between community composition (beta diversity) and external factors.

The standard way to do this is to perform a so-called permutational multivariate analysis of variance (PERMANOVA). This method takes as input the abundance table, which measure of distance you want to base the test on and a formula that tells the model how you think the variables are associated with each other.
library(vegan)
otu <- abundances(pseq)
meta <- meta(pseq)

permanova <- adonis(t(otu) ~ Time,
                    data = meta, permutations=9999, method = "euclidean")

# P-value

print(as.data.frame(permanova$aov.tab))


The cohort variable is not significantly associated with microbiota composition (p-value is over 0.05).

We can, however, visualize those taxa whose abundances drive the differences between cohorts. We first need to extract the model coefficients of taxa:

  coef <- coefficients(permanova)["Time1",]
top.coef <- coef[rev(order(abs(coef)))[1:20]]
par(mar = c(3, 14, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")

#For example, Enterococcus (Firmicutes), which plays a crucial role in metabolic adaptability against pathogenic or plant toxins and anti-herbivore defense, was found to be one of the predominant gut microorganism of lepidopteran insects, including B. mori, Helicoverpa zea, and Porthetria dispar (Paniagua Voirol et al., 2018; Zhang et al., 2022). 

#Rhodococcus in the triatomine gut are believed to play an important role in the metabolism of the vector, such as by participating in the synthesis of group B vitamins or by being digested by the bugs directly to provide missing nutrients (Sassera et al., 2013). Moreover, the most attractive aspect is the host-symbiont relationship between triatomines and Rhodococcus; since Rhodococcus bacteria can be easily cultured and genetically modified to harm the pathogen in vector gut, they are probably suitable tools for the control of trypanosomiasis (Sassera et al., 2013). A
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7667259/








i=nmds_spec_comp_aitchison
# Convert the NMDS points to a data frame
nmds_spec_gg <- as.data.frame(i$points)
# Add the "Sample" column based on row names
nmds_spec_gg <- nmds_spec_gg %>% rownames_to_column("Sample")
# Merge the NMDS data with the metadata by the "Sample" column
nmds_spec_gg <- left_join(meta, nmds_spec_gg, by = "Sample")
# Merge the metadata and NMDS data frames
merged_data <- dplyr::left_join(as.data.frame(meta), as.data.frame(nmds_spec_gg), by = "Sample")

data("enterotype", package = "mia")
tse2 <- enterotype


# Apply relative transform
tse2 <- transformAssay(tse2,
                       method = "relabundance")



AM_metagenomes <- phyloseq::subset_samples(merged_metagenomes, Gut == "PM")


# Run NMDS on relabundance assay with Bray-Curtis distances
pseq <- aggregate_rare(AM_metagenomes, level = "Genus", detection = 0.1/100, prevalence = 50/100)
pseq <- microbiome::transform(pseq, transform = "compositional")
# Perform RDA
tse2 <- mia::runRDA(pseq,
               assay.type = "relabundance",
               formula = assay ~ ClinicalStatus + Gender + Age,
               distance = "bray",
               na.action = na.exclude)

# Store results of PERMANOVA test
rda_info <- attr(reducedDim(tse2, "RDA"), "significance")






































merged_metagenomes <- import_biom("/home/vincent/Documents/project/Metagenomics_2024/data/Trypanosoma_exposure/READBASED/merge_species.biom") 
meta <- read.csv(file = "DATA/tryp_metadata.csv", sep = ",")


meta <- meta  %>%  arrange(row_number(SRA.identifier)) # sort the data frame
merged_metagenomes@sam_data <- sample_data(meta) # associate the metadata to the to the phyloseq object
column_name <-  meta %>% pull(Sample) # extract the sample names
sample_names(merged_metagenomes) <- column_name # associate the sample names to the phyloseq object

merged_metagenomes@tax_table@.Data <- substring(merged_metagenomes@tax_table@.Data, 4) # remove the unnecessary 'k_' in the taxonomy.
colnames(merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") # change the rank name 


merged_metagenomes <- subset_taxa(merged_metagenomes, Kingdom %in% c("Archaea", "Bacteria", "Fungi", "Viruses"))




pseq <- aggregate_rare(merged_metagenomes, level = "Genus", detection = 0.1/100, prevalence = 50/100)
pseq <- microbiome::transform(pseq, transform = "compositional")

otu <- otu_table(pseq) %>% t()











head(t(otu))
head(meta)

dbrda_pseq <- vegan::dbrda(otu ~ factor(meta$Type) + factor(meta$Time) + factor(meta$Gut) + factor(meta$Reads), data = otu, distance = "bray")




pseq <-  mia::makeTreeSummarizedExperimentFromPhyloseq(merged_metagenomes)

# Apply relative transform
pseq <- mia::transformAssay(pseq,
                       method = "relabundance")

pseq <- mia::runRDA(pseq,
               assay.type = "relabundance",
               formula = assay ~ Gut + Time + Type + Reads,
               distance = "bray",
               na.action = na.exclude)


# Store results of PERMANOVA test
rda_info <- attr(reducedDim(pseq, "RDA"), "significance")
rda_info$permanova

# Load packages for plotting function
library(miaViz)

# Generate RDA plot coloured by clinical status
plotRDA(pseq, "RDA", colour_by = "Time")


