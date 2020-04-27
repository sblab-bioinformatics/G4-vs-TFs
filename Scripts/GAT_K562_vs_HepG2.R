# Compare GAT analysis results in K562 vs HepG2



# ======= Libraries and working directory
library(data.table)
library(ggplot2)
library(ggrepel)



setwd('ENCODE/')


# ==== Script parameters
savetables <- 1
plotdata <- 1

# ====== Load GAT results for K562 and HepG2
HepG2_GAT <- fread(file = "HepG2/GAT_analysis/GAT_HepG2_rep1-3_100Peaks_Gat-Analysis_trim.csv", header = TRUE)
K562_GAT <- fread(file = "GAT rerun 2019_new/GAT_analysis/GAT_K562_BG4_Rerun_100Peaks_Jan2019_Gat-Analysis_trim.csv", header = TRUE)

# load known G4 interactors
Known_G4_proteins <- fread(file = "K562_Rerun_Jan2019/GAT_analysis/G4IPD_Oct2019.csv", header = F)

# Strip any extension "eGFP", "3xFLAG" etc
# remove information on tag
K562_GAT$Experiment.target <- gsub("eGFP-", "", K562_GAT$Experiment.target)
K562_GAT$Experiment.target <- gsub("3xFLAG-", "", K562_GAT$Experiment.target)

HepG2_GAT$Experiment.target <- gsub("eGFP-", "", HepG2_GAT$Experiment.target)
HepG2_GAT$Experiment.target <- gsub("3xFLAG-", "", HepG2_GAT$Experiment.target)


# Focus on GAT_DHS 
K562_GAT <- K562_GAT[, c("Experiment.target", "Experiment.accession", "DHS_fold")]
colnames(K562_GAT)[3] <- "K562_DHS" 

HepG2_GAT <- HepG2_GAT[, c("Experiment.target","Experiment.accession", "DHS_fold")]  
colnames(HepG2_GAT)[3] <- "HepG2_DHS"



# There are different approaches when multiple experiments exist for the same target e.g. SP1
# Comparing overlap bed files for a particular TF often showed substantial variability even within a particular cell line, suggesting low quality ChIP-seq datasets that may skew the analysis
# Keeping the highest enrichment per target should give the cleanest and most representative


K562_GAT_max <- aggregate(K562_GAT, by = list(K562_GAT$Experiment.target), FUN = max)[, c("Experiment.target", "K562_DHS")]
HepG2_GAT_max <- aggregate(HepG2_GAT, by = list(HepG2_GAT$Experiment.target), FUN = max)[, c("Experiment.target", "HepG2_DHS")]

# Perform full outer join. Results in table containing all different combinations for one particular target. 
Merged_K562_HepG2_max  <- merge(K562_GAT_max, HepG2_GAT_max, by = "Experiment.target", all = TRUE) 

# Remove all NAs to focus on targets that are present in both 
Merged_K562_HepG2_max <- Merged_K562_HepG2_max[rowSums(is.na(Merged_K562_HepG2_max)) == F, ]


write.csv(Merged_K562_HepG2_max, file =paste("HepG2/GAT_analysis/GAT_HepG2_Merged_K562_HepG2_max.csv",  sep=""))



temp <- Merged_K562_HepG2_max
temp$G4IPD <- ifelse(temp$Experiment.target %in% Known_G4_proteins$V1, "#279615", "black")

gg_max <- ggplot(Merged_K562_HepG2_max, aes(x=K562_DHS, y=HepG2_DHS)) + 
  geom_point(colour=temp$G4IPD) + 
  theme_minimal()+
  # labs(title="TF enrichment at endogenous G4s - K562 vs HepG2", y="enrichment in HepG2 \n randomization in DHS", x="enrichment in K562 \n randomization in DHS") + 
  theme(axis.title.y = element_text(size = rel(1), vjust =2, angle = 90)) + 
  theme(axis.title.x = element_text(size = rel(1), angle = 0)) + 
  theme(plot.title = element_text(face = "bold", size = rel(1.2))) + 
  coord_fixed() +
  scale_y_continuous(breaks=seq(0,11,2), limits=c(0,8)) +
  scale_x_continuous(breaks=seq(0,11,2), limits=c(0,11)) +
  #coord_fixed(ratio = 1) +
  labs(color = "") +
  theme(legend.position="bottom") +
  geom_label_repel(aes(label=ifelse(
    (Merged_K562_HepG2_max$K562_DHS > 8.5 & Merged_K562_HepG2_max$HepG2_DHS > 3  ), as.character(Experiment.target),'')),
    box.padding = 0.4, point.padding = 0.4, size=rel(2.5), segment.color = 'grey80') +
  geom_label_repel(aes(label=ifelse(
    (Merged_K562_HepG2_max$K562_DHS > 4 & Merged_K562_HepG2_max$HepG2_DHS > 5  ), as.character(Experiment.target),'')),
    box.padding = 0.4, point.padding = 0.4, size=rel(2.5), segment.color = 'grey80') +
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA))


ggsave("HepG2/GAT_analysis/figures/K562_vs_HepG2_max_10x10.pdf", width = 10/2.54, height = 10/2.54, useDingbats=F)



# Chekc Spearman correlation
corMat_GAT_spearman <- cor(Merged_K562_HepG2[, c(3,5)], method= "spearman")
corMat_GAT_spearman_max <- cor(Merged_K562_HepG2_max[, 2:3], method= "spearman")

corMat_GAT_pearson <- cor(Merged_K562_HepG2[, c(3,5)], method= "pearson")
corMat_GAT_pearson_max <- cor(Merged_K562_HepG2_max[, 2:3], method= "pearson")

corMat_GAT_spearman_max <- cor.test(Merged_K562_HepG2_max$K562_DHS, Merged_K562_HepG2_max$HepG2_DHS, method= "spearman")
corMat_GAT_spearman_max$p.value
