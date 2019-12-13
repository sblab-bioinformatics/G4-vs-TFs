
# ======= Libraries and working directory
library(data.table)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(ggrepel)


setwd('ENCODE_Analysis/')


# ==== Script parameters
plotdata <- 1
savetables <- 1


# ====== Parameters for Filtering the data ====
qValue_cutoff <-  0.05
PEAK_NUMBER_Cutoff <- 500


Explicit_exclusion <- c("ENCFF782GWS") # ("HNRNPF", "RBM15", "RBM17", "RBM34") Datasets have been flagged in original paper: Cell 2019 (10.1016/j.cell.2019.06.001)


# ====== Load meta data
ENCODE_bed_meta <- read.table(file = "ENCODE_K562_Jan2019_BED_Meta.tsv", sep = '\t', header = TRUE)
ENCODE_bed_meta$Experiment.target <- gsub("-human", "", ENCODE_bed_meta$Experiment.target)
ENCODE_bed_meta <- ENCODE_bed_meta[rowSums(is.na(ENCODE_bed_meta))  != ncol(ENCODE_bed_meta), ] # Remove empty rows


#read in peak numbers
Peak_Number <- read.table(file="GAT_opOQ_noBG4_TSSK562_Jan2019_NumberOfPeaks.txt", sep= '\t', header = FALSE)
colnames(Peak_Number)[1] <- "ENCODE_ID"
colnames(Peak_Number)[2] <- "Peaks"
Peak_Number$ENCODE_ID <- sapply(as.vector(Peak_Number$ENCODE_ID), function(x) {gsub("\\..*","",x)}) # strip extension
ENCODE_bed_meta$Peaks <- Peak_Number$Peaks[match(ENCODE_bed_meta$File.accession, Peak_Number$ENCODE_ID)]


# ================================================================================================================================================================
# ==== Load result tables for four GAT runs using different workspaces
# ================================================================================================================================================================

# ====== open OQS without BG4 signal upstream of TSS
# This are the old sites 5kb upstream + 5UTR
#GAT_opOQ_noBG4_TSS <- read.table(file = "GAT_opOQ_noBG4_TSSGAT.opOQs_noBG4_atTSS.K562.DHS.Rerun_2019.DHS.tsv", sep = '\t', header = TRUE) 

# Use 1kb upstream + 5'UTR
GAT_opOQ_noBG4_TSS <- read.table(file = "GAT_opOQ_noBG4_TSSGAT.opOQs_noBG4_1kbupstreamTSS.K562.DHS.Rerun_2019.DHS.tsv", sep = '\t', header = TRUE) 
GAT_opOQ_noBG4_TSS <- GAT_opOQ_noBG4_TSS[, 2:11]
colnames(GAT_opOQ_noBG4_TSS)[1] <- "ENCODE_ID"
GAT_opOQ_noBG4_TSS$ENCODE_ID <- gsub(".cut.bed.gz", "", GAT_opOQ_noBG4_TSS$ENCODE_ID)
###Sort by rank and include and extra column 'rank' to see how high the marker scored
GAT_opOQ_noBG4_TSS <- GAT_opOQ_noBG4_TSS[order(GAT_opOQ_noBG4_TSS$fold, decreasing = TRUE), ]
GAT_opOQ_noBG4_TSS$rank <- seq.int(nrow(GAT_opOQ_noBG4_TSS))


# ====== Original shuffling of BG4 peaks in DHS
Gat_DHS <- read.table(file = "GAT_opOQ_noBG4_TSSGAT_K562_BG4_Rerun_Jan2019.DHS.tsv", sep = '\t', header = TRUE)
Gat_DHS <- Gat_DHS[, 2:11]
colnames(Gat_DHS)[1] <- "ENCODE_ID"
Gat_DHS$ENCODE_ID <- gsub(".cut.bed.gz", "", Gat_DHS$ENCODE_ID)
###Sort by rank and include and extra column 'rank' to see how high the marker scored
Gat_DHS <- Gat_DHS[order(Gat_DHS$fold, decreasing = TRUE), ]
Gat_DHS$rank <- seq.int(nrow(Gat_DHS))





# ================================================================================================================================================================
# === Generate a Merged data sheet comprising ENCODE meta data and GAT analysis for the different runs ====
# ================================================================================================================================================================

Merged_all <- ENCODE_bed_meta


GAT_opOQ_noBG4_TSS$Peaks <- Peak_Number$Peaks[match(GAT_opOQ_noBG4_TSS$ENCODE_ID, Peak_Number$ENCODE_ID)]


# wihtout BG4 signal at TSS 5UTR
Merged_all$noBG4atTSS_observed <- GAT_opOQ_noBG4_TSS$observed[match(Merged_all$File.accession, GAT_opOQ_noBG4_TSS$ENCODE_ID)]
Merged_all$noBG4atTSS_expected <- GAT_opOQ_noBG4_TSS$expected[match(Merged_all$File.accession, GAT_opOQ_noBG4_TSS$ENCODE_ID)]
Merged_all$noBG4atTSS_CI95low <- GAT_opOQ_noBG4_TSS$CI95low[match(Merged_all$File.accession, GAT_opOQ_noBG4_TSS$ENCODE_ID)]
Merged_all$noBG4atTSS_CI95high <- GAT_opOQ_noBG4_TSS$CI95high[match(Merged_all$File.accession, GAT_opOQ_noBG4_TSS$ENCODE_ID)]
Merged_all$noBG4atTSS_stddev <- GAT_opOQ_noBG4_TSS$stddev[match(Merged_all$File.accession, GAT_opOQ_noBG4_TSS$ENCODE_ID)]
Merged_all$noBG4atTSS_fold <- GAT_opOQ_noBG4_TSS$fold[match(Merged_all$File.accession, GAT_opOQ_noBG4_TSS$ENCODE_ID)]
Merged_all$noBG4atTSS_l2fold <- GAT_opOQ_noBG4_TSS$l2fold[match(Merged_all$File.accession, GAT_opOQ_noBG4_TSS$ENCODE_ID)]
Merged_all$noBG4atTSS_pvalue <- GAT_opOQ_noBG4_TSS$pvalue[match(Merged_all$File.accession, GAT_opOQ_noBG4_TSS$ENCODE_ID)]
Merged_all$noBG4atTSS_qvalue <- GAT_opOQ_noBG4_TSS$qvalue[match(Merged_all$File.accession, GAT_opOQ_noBG4_TSS$ENCODE_ID)]
Merged_all$noBG4atTSS_rank <- GAT_opOQ_noBG4_TSS$rank[match(Merged_all$File.accession, GAT_opOQ_noBG4_TSS$ENCODE_ID)]

#BG4 peaks shuffled in DHS (open chromatin)
Merged_all$DHS_observed <- Gat_DHS$observed[match(Merged_all$File.accession, Gat_DHS$ENCODE_ID)]
Merged_all$DHS_expected <- Gat_DHS$expected[match(Merged_all$File.accession, Gat_DHS$ENCODE_ID)]
Merged_all$DHS_CI95low <- Gat_DHS$CI95low[match(Merged_all$File.accession, Gat_DHS$ENCODE_ID)]
Merged_all$DHS_CI95high <- Gat_DHS$CI95high[match(Merged_all$File.accession, Gat_DHS$ENCODE_ID)]
Merged_all$DHS_stddev <- Gat_DHS$stddev[match(Merged_all$File.accession, Gat_DHS$ENCODE_ID)]
Merged_all$DHS_fold <- Gat_DHS$fold[match(Merged_all$File.accession, Gat_DHS$ENCODE_ID)]
Merged_all$DHS_l2fold <- Gat_DHS$l2fold[match(Merged_all$File.accession, Gat_DHS$ENCODE_ID)]
Merged_all$DHS_pvalue <- Gat_DHS$pvalue[match(Merged_all$File.accession, Gat_DHS$ENCODE_ID)]
Merged_all$DHS_qvalue <- Gat_DHS$qvalue[match(Merged_all$File.accession, Gat_DHS$ENCODE_ID)]
Merged_all$DHS_rank <- Gat_DHS$rank[match(Merged_all$File.accession, Gat_DHS$ENCODE_ID)]
Merged_all$DHS_rank <- Gat_DHS$rank[match(Merged_all$File.accession, Gat_DHS$ENCODE_ID)]


Merged_all <- Merged_all[rowSums(is.na(Merged_all)) != ncol(Merged_all), ] # remove empty rows 

# export all data unfiltered
if (savetables)
{
  write.csv(Merged_all, file =paste("GAT_K562_GrichnessAtOQs_Dec2019_Gat-Analysis_all.csv",  sep=""))
}


#====== FILTER  q-values 
# Relevant bed files have already been selected prior to downloading. 
# q-Value gives an idea of the reliabiltiy of the shuffling analysis. Only keep cases where shuffling was successfull in case of 'WithBG' and 'noBG4_atTSS' as these will be relevant data sets.
GAT_FILTERED <- Merged_all[ Merged_all$DHS_qvalue  < qValue_cutoff , ]
GAT_FILTERED <- GAT_FILTERED[ GAT_FILTERED$noBG4atTSS_qvalue < qValue_cutoff , ]

GAT_FILTERED <- GAT_FILTERED[GAT_FILTERED$Peaks > PEAK_NUMBER_Cutoff, ]
#Remove data sets that have been explicitly flagged
GAT_FILTERED <- GAT_FILTERED[!(GAT_FILTERED$File.accession  %in% Explicit_exclusion), ]


# Update relative ranking after removing datasets
GAT_FILTERED <- GAT_FILTERED[order(GAT_FILTERED$noBG4atTSS_fold, decreasing = TRUE), ]
GAT_FILTERED$noBG4atTSS_rank <- seq.int(nrow(GAT_FILTERED))
GAT_FILTERED <- GAT_FILTERED[order(GAT_FILTERED$DHS_fold, decreasing = TRUE), ]
GAT_FILTERED$DHS_rank <- seq.int(nrow(GAT_FILTERED))


if(savetables)
{
  write.csv(GAT_FILTERED, file = "GAT_K562_GrichnessAtOQs_Dec2019_Gat-Analysis_FILTERED.csv", row.names = FALSE)
}



# ===========Trimmed and Condensed version of filterd===================================
# Generate a trimmed  version only containing most relevant parameters

GAT_FILTERED_trim <- GAT_FILTERED[,  c("File.accession", "Experiment.accession", "Experiment.target", "Peaks", "noBG4atTSS_fold", "noBG4atTSS_qvalue", "noBG4atTSS_rank", "DHS_fold", "DHS_qvalue", "DHS_rank")] 

# Ratio of DHS vs noBG4atTSS
GAT_FILTERED_trim$ratio <- (GAT_FILTERED_trim$DHS_fold/GAT_FILTERED_trim$noBG4atTSS_fold)
GAT_FILTERED_trim <- GAT_FILTERED_trim[order(GAT_FILTERED_trim$ratio, decreasing = TRUE), ]
GAT_FILTERED_trim$ratio_rank <- seq.int(nrow(GAT_FILTERED))



if(savetables)
{
  write.csv(GAT_FILTERED, file = "GAT_K562_GrichnessAtOQs_Dec2019_Gat-Analysis_FILTERED_trim.csv", row.names = FALSE)
}



# =====
# Plots
# =======

#for plots remove all the additional tag information: eGFP, FLAG, etc.
GAT_FILTERED_trim$Experiment.target <- gsub("eGFP-", "", GAT_FILTERED_trim$Experiment.target)
GAT_FILTERED_trim$Experiment.target <- gsub("3xFLAG-", "" , GAT_FILTERED_trim$Experiment.target)

# Label G4 associated proteins
Known_G4_proteins <- read.csv(file = "G4IPD_Oct2019.csv", header = F)
temp <- GAT_FILTERED_trim
temp$G4IPD <- ifelse(temp$Experiment.target %in% Known_G4_proteins$V1, "#aad2a5", "grey50")


# Label with transparent background 
gg_grey <- ggplot(temp, aes(x=DHS_fold, y=noBG4atTSS_fold, colour=G4IPD)) + 
  geom_point() + 
  theme_minimal()+
  labs(title="TF enrichment vs G-richness", y="enrichment of OQS at TSS \n without BG4", x="enrichment at endogenous G4s") + 
  theme(axis.title.y = element_text(size = rel(1), vjust =2, angle = 90)) + 
  theme(axis.title.x = element_text(size = rel(1), angle = 0)) + 
  theme(plot.title = element_text(face = "bold", size = rel(1.2))) + 
  coord_fixed() +
  scale_y_continuous(breaks=seq(0,5,2)) +
  scale_x_continuous(breaks=seq(0,11,2)) +
  scale_color_manual(values=c("#639b57", "black")) +
  labs(color = "") +
  theme(legend.position="bottom") +
  geom_label_repel(aes(label=ifelse(
    (GAT_FILTERED_trim$ratio > 2.7 & GAT_FILTERED_trim$DHS_fold > 6.5), as.character(Experiment.target),'')),
    box.padding = 0.4, point.padding = 0.4, size=rel(2.5), segment.color = 'grey80') + 
  geom_label_repel(aes(label=ifelse(
    (GAT_FILTERED_trim$ratio < 0.4 & GAT_FILTERED_trim$noBG4atTSS_fold >1), as.character(Experiment.target) ,'')),
    box.padding = 0.4, point.padding = 0.4, size=rel(2.5), segment.color = 'grey20')+
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA))
plot(gg_grey) 



ggsave("figures/GAT_K562_GrichnessAtOQs_Dec2019_vsBG4.pdf", width = 20/2.54, height = 14/2.54, useDingbats=F)



