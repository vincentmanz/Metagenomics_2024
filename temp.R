
phyloseq::plot_richness(pseq, x = "Type", color = "Type",
                        measures = c("Observed", "Shannon", "Simpson","Chao1")) +  ggplot2::geom_boxplot()

microbiome::alpha(pseq)

tab <- microbiome::alpha(pseq)

metadata <- data.frame(sample_data(pseq))


tab |>
  dplyr::select(observed,
                diversity_gini_simpson,
                diversity_shannon,
                chao1) |>
  indices_normality(nrow = 2, ncol = 2)


stats::kruskal.test(tab$diversity_shannon ~ Time, data = metadata)


aov_observed <- stats::aov(tab$diversity_shannon ~ Time, metadata)

signif_pairgroups <- stats::TukeyHSD(aov_observed, method = "bh")
signif_pairgroups

diversity_shannon <- tab$diversity_shannon
pairwise_test <- ggpubr::compare_means(diversity_shannon ~ Type,
                                       metadata,
                                       method = "wilcox.test")

#Boxplot as previously seen
graph_shan <- ggplot(metadata, aes(x = Type, y = tab$diversity_gini_simpson)) + 
  geom_boxplot(alpha=0.6
               #fill = c("#00AFBB", "#E7B800"),
               #color = c("#00AFBB", "#E7B800")
               ) +
  geom_jitter(aes(colour = Time),
              position = position_jitter(0.02) ,
              cex=2.2)+
  stat_summary(fun = mean, geom = "point",
               shape = 17, size = 3,
               color = "white")

#Add p-value on graph
graph_shan + ggpubr::stat_pvalue_manual(
  pairwise_test,
  y.position = 3.5,
  label = "p.adj = {p.adj}",
  color = "blue",
  linetype = 1,
  tip.length = 0.01
)
https://microbiome.github.io/OMA/docs/devel/pages/14_alpha_diversity.html

https://joey711.github.io/phyloseq/plot_network-examples