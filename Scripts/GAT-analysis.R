
# This R script will read blocks of GAT shuffling results (in whitelisted genome (WL), open chromatin (DHS), OQS and 'open' OQS) combine results with ENCODE meta data


# ======= Libraries and working directory
library(data.table)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(scales)

setwd('ENCODE_Analysis/')


# ==== Script parameters
plotdata <- 1
savetables <- 1


# ===  list of segments that were shuffled in different backgrounds
GAT_segments <- c('BG4_Rerun') 

# ====== Parameters for Filtering the data ====
qValue_cutoff <-  0.05
PEAK_NUMBER_Cutoff <- 500 

Explicit_exclusion <- c("ENCFF782GWS") # ("HNRNPF", "RBM15", "RBM17", "RBM34") Datasets have been flagged in original paper (Cell 2019 (10.1016/j.cell.2019.06.001)), but not yet on ENCODE


# ====== Load meta data
ENCODE_bed_meta <- read.table(file = "K562_Rerun_Jan2019/ENCODE_K562_Jan2019_BED_Meta.tsv", sep = '\t', header = TRUE)
ENCODE_bed_meta$Experiment.target <- gsub("-human", "", ENCODE_bed_meta$Experiment.target)
ENCODE_bed_meta <- ENCODE_bed_meta[rowSums(is.na(ENCODE_bed_meta))  != ncol(ENCODE_bed_meta), ] # Remove empty rows


#read in peak numbers
Peak_Number <- read.table(file="K562_Rerun_Jan2019/GAT_Rerun2019_output/K562_Jan2019_NumberOfPeaks.txt", sep= '\t', header = FALSE)
colnames(Peak_Number)[1] <- "ENCODE_ID"
colnames(Peak_Number)[2] <- "Peaks"
Peak_Number$ENCODE_ID <- sapply(as.vector(Peak_Number$ENCODE_ID), function(x) {gsub("\\..*","",x)}) # strip extension
ENCODE_bed_meta$Peaks <- Peak_Number$Peaks[match(ENCODE_bed_meta$File.accession, Peak_Number$ENCODE_ID)]


# =====================================================================================================================================================
# == Cycle Through blocks of GAT results as defined in GAT_segments
# =====================================================================================================================================================
for (i in 1:length(GAT_segments))
{
  Gat_Seg <- (GAT_segments)[i]
  
  
  # =====================================================================================================================================================
  # Load result tables for four GAT runs using different workplaces (Whitelisted genome, DNAse-hypersensitivity sites, Observed Quadruplexes, OQs at DHS)
  # =====================================================================================================================================================
  
  # ====== shuffle in whitelisted genome  
  Gat_WL <- read.table(file =paste("K562_Rerun_Jan2019/GAT_Rerun2019_output/GAT_K562_", Gat_Seg, "_Jan2019.WL.tsv",  sep="") , sep = '\t', header = TRUE) 
  Gat_WL <- Gat_WL[, 2:11]
  colnames(Gat_WL)[1] <- "ENCODE_ID"
  Gat_WL$ENCODE_ID <- gsub(".cut.bed.gz", "", Gat_WL$ENCODE_ID)
  ###Sort by rank and include and extra column 'rank' to see how high the marker scored
  Gat_WL <- Gat_WL[order(Gat_WL$fold, decreasing = TRUE), ]
  Gat_WL$rank <- seq.int(nrow(Gat_WL))
  
  
  # ===== shuffle in  DHS (open chromatin)
  Gat_DHS <- read.table(file = paste("K562_Rerun_Jan2019/GAT_Rerun2019_output/GAT_K562_", Gat_Seg, "_Jan2019.DHS.tsv",  sep=""), sep = '\t', header = TRUE)
  Gat_DHS <- Gat_DHS[, 2:11]
  colnames(Gat_DHS)[1] <- "ENCODE_ID"
  Gat_DHS$ENCODE_ID <- gsub(".cut.bed.gz", "", Gat_DHS$ENCODE_ID)
  ###Sort by rank and include and extra column 'rank' to see how high the marker scored
  Gat_DHS <- Gat_DHS[order(Gat_DHS$fold, decreasing = TRUE), ]
  Gat_DHS$rank <- seq.int(nrow(Gat_DHS))
  
  
  # ========  shuffle in observed Quadruplexes
  Gat_OQs <- read.table(file = paste("K562_Rerun_Jan2019/GAT_Rerun2019_output/GAT_K562_", Gat_Seg, "_Jan2019.OQS.tsv",  sep=""), sep = '\t', header = TRUE)
  Gat_OQs <- Gat_OQs[, 2:11]
  colnames(Gat_OQs)[1] <- "ENCODE_ID"
  Gat_OQs$ENCODE_ID <- gsub(".cut.bed.gz", "", Gat_OQs$ENCODE_ID)
  ###Sort by rank and include and extra column 'rank' to see how high the marker scored
  Gat_OQs <- Gat_OQs[order(Gat_OQs$fold, decreasing = TRUE), ]
  Gat_OQs$rank <- seq.int(nrow(Gat_OQs))
  
  
  # ==== shuffle in observed Quadruplexes at DHS
  Gat_openOQs <- read.table(file = paste("K562_Rerun_Jan2019/GAT_Rerun2019_output/GAT_K562_", Gat_Seg, "_Jan2019.opOQS.tsv",  sep=""), sep = '\t', header = TRUE)
  Gat_openOQs <- Gat_openOQs[, 2:11]
  colnames(Gat_openOQs)[1] <- "ENCODE_ID"
  Gat_openOQs$ENCODE_ID <- gsub(".cut.bed.gz", "", Gat_openOQs$ENCODE_ID)
  ###Sort by rank and include and extra column 'rank' to see how high the marker scored
  Gat_openOQs <- Gat_openOQs[order(Gat_openOQs$fold, decreasing = TRUE), ]
  Gat_openOQs$rank <- seq.int(nrow(Gat_openOQs))
  
  # ============================================================================================================================================================
  # === Generate a Merged data sheet comprising ENCODE meta data and GAT analysis for the different runs ====
  # ============================================================================================================================================================
  
  Merged_all <- ENCODE_bed_meta
  
  #whitelisted genome
  Merged_all$WL_observed <- Gat_WL$observed[match(Merged_all$File.accession, Gat_WL$ENCODE_ID)]
  Merged_all$WL_expected <- Gat_WL$expected[match(Merged_all$File.accession, Gat_WL$ENCODE_ID)]
  Merged_all$WL_CI95low <- Gat_WL$CI95low[match(Merged_all$File.accession, Gat_WL$ENCODE_ID)]
  Merged_all$WL_CI95high <- Gat_WL$CI95high[match(Merged_all$File.accession, Gat_WL$ENCODE_ID)]
  Merged_all$WL_stddev <- Gat_WL$stddev[match(Merged_all$File.accession, Gat_WL$ENCODE_ID)]
  Merged_all$WL_fold <- Gat_WL$fold[match(Merged_all$File.accession, Gat_WL$ENCODE_ID)]
  Merged_all$WL_l2fold <- Gat_WL$l2fold[match(Merged_all$File.accession, Gat_WL$ENCODE_ID)]
  Merged_all$WL_pvalue <- Gat_WL$pvalue[match(Merged_all$File.accession, Gat_WL$ENCODE_ID)]
  Merged_all$WL_qvalue <- Gat_WL$qvalue[match(Merged_all$File.accession, Gat_WL$ENCODE_ID)]
  Merged_all$WL_rank <- Gat_WL$rank[match(Merged_all$File.accession, Gat_WL$ENCODE_ID)]
  
  #DHS (open chromatin)
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
  
  #observed Quadruplexes
  Merged_all$OQs_observed <- Gat_OQs$observed[match(Merged_all$File.accession, Gat_OQs$ENCODE_ID)]
  Merged_all$OQs_expected <- Gat_OQs$expected[match(Merged_all$File.accession, Gat_OQs$ENCODE_ID)]
  Merged_all$OQs_CI95low <- Gat_OQs$CI95low[match(Merged_all$File.accession, Gat_OQs$ENCODE_ID)]
  Merged_all$OQs_CI95high <- Gat_OQs$CI95high[match(Merged_all$File.accession, Gat_OQs$ENCODE_ID)]
  Merged_all$OQs_stddev <- Gat_OQs$stddev[match(Merged_all$File.accession, Gat_OQs$ENCODE_ID)]
  Merged_all$OQs_fold <- Gat_OQs$fold[match(Merged_all$File.accession, Gat_OQs$ENCODE_ID)]
  Merged_all$OQs_l2fold <- Gat_OQs$l2fold[match(Merged_all$File.accession, Gat_OQs$ENCODE_ID)]
  Merged_all$OQs_pvalue <- Gat_OQs$pvalue[match(Merged_all$File.accession, Gat_OQs$ENCODE_ID)]
  Merged_all$OQs_qvalue <- Gat_OQs$qvalue[match(Merged_all$File.accession, Gat_OQs$ENCODE_ID)]
  Merged_all$OQs_rank <- Gat_OQs$rank[match(Merged_all$File.accession, Gat_OQs$ENCODE_ID)]
  
  #observed Quadruplexes at DHS
  Merged_all$opOQs_observed <- Gat_openOQs$observed[match(Merged_all$File.accession, Gat_openOQs$ENCODE_ID)]
  Merged_all$opOQs_expected <- Gat_openOQs$expected[match(Merged_all$File.accession, Gat_openOQs$ENCODE_ID)]
  Merged_all$opOQs_CI95low <- Gat_openOQs$CI95low[match(Merged_all$File.accession, Gat_openOQs$ENCODE_ID)]
  Merged_all$opOQs_CI95high <- Gat_openOQs$CI95high[match(Merged_all$File.accession, Gat_openOQs$ENCODE_ID)]
  Merged_all$opOQs_stddev <- Gat_openOQs$stddev[match(Merged_all$File.accession, Gat_openOQs$ENCODE_ID)]
  Merged_all$opOQs_fold <- Gat_openOQs$fold[match(Merged_all$File.accession, Gat_openOQs$ENCODE_ID)]
  Merged_all$opOQs_l2fold <- Gat_openOQs$l2fold[match(Merged_all$File.accession, Gat_openOQs$ENCODE_ID)]
  Merged_all$opOQs_pvalue <- Gat_openOQs$pvalue[match(Merged_all$File.accession, Gat_openOQs$ENCODE_ID)]
  Merged_all$opOQs_qvalue <- Gat_openOQs$qvalue[match(Merged_all$File.accession, Gat_openOQs$ENCODE_ID)]
  Merged_all$opOQs_rank <- Gat_openOQs$rank[match(Merged_all$File.accession, Gat_openOQs$ENCODE_ID)]
  
  Merged_all <- Merged_all[rowSums(is.na(Merged_all)) != ncol(Merged_all), ] # remove empty rows 
  
  # export all data unfiltered
  if (savetables)
  {
    write.csv(Merged_all, file =paste("K562_Rerun_Jan2019/GAT_analysis/GAT_K562_", Gat_Seg, "_Jan2019_Gat-Analysis_all.csv",  sep=""))
  }
  
  
  #====== FILTER  q-values 
  # Most filtering has been performed prior to downloading the bed-files. 
  # q-Value gives an idea of the reliabiltiy of the shuffling analysis. Only keep cases where shuffling was successfull in case of open OQS as this will be the reference in the end and is consistant with K562 analysis
  #GAT_FILTERED <- Merged_all[ Merged_all$opOQs_qvalue < qValue_cutoff , ]
  GAT_FILTERED <- Merged_all[Merged_all$Peaks > PEAK_NUMBER_Cutoff, ]
  
  #Remove data sets that have been explicitly flagged
  GAT_FILTERED <- GAT_FILTERED[!(GAT_FILTERED$File.accession  %in% Explicit_exclusion), ]
  
  
  # Update relative ranking after removing datasets
  GAT_FILTERED <- GAT_FILTERED[order(GAT_FILTERED$WL_fold, decreasing = TRUE), ]
  GAT_FILTERED$WL_rank <- seq.int(nrow(GAT_FILTERED))
  GAT_FILTERED <- GAT_FILTERED[order(GAT_FILTERED$DHS_fold, decreasing = TRUE), ]
  GAT_FILTERED$DHS_rank <- seq.int(nrow(GAT_FILTERED))
  GAT_FILTERED <- GAT_FILTERED[order(GAT_FILTERED$OQs_fold, decreasing = TRUE), ]
  GAT_FILTERED$OQs_rank <- seq.int(nrow(GAT_FILTERED))
  GAT_FILTERED <- GAT_FILTERED[order(GAT_FILTERED$opOQs_fold, decreasing = TRUE), ]
  GAT_FILTERED$opOQs_rank <- seq.int(nrow(GAT_FILTERED))
  
  if(savetables)
  {
    write.csv(GAT_FILTERED, file =paste("K562_Rerun_Jan2019/GAT_analysis/GAT_K562_", Gat_Seg, '_', PEAK_NUMBER_Cutoff, "Peaks_Jan2019_Gat-Analysis_filt.csv",  sep=""))
  }
  
  
  # ===========Trimmed and Condensed version of filterd===================================
  # Generate a trimmed  version only containing most relevant parameters
  GAT_FILTERED_trim <- GAT_FILTERED[,  c("File.accession", "Experiment.accession", "Experiment.target", "Peaks", "WL_fold", "WL_qvalue", "WL_rank", "DHS_fold", "DHS_qvalue", "DHS_rank", "OQs_fold", "OQs_qvalue", "OQs_rank", "opOQs_fold", "opOQs_qvalue", "opOQs_rank")] 
  GAT_FILTERED_trim$average_rank <- rowMeans(GAT_FILTERED_trim[, c("WL_rank",  "DHS_rank", "OQs_rank", "opOQs_rank")])
  GAT_FILTERED_trim <- GAT_FILTERED_trim[order(GAT_FILTERED_trim$average_rank, decreasing = FALSE), ]
  
  if(savetables)
  {
    write.csv(GAT_FILTERED_trim, file =paste("K562_Rerun_Jan2019/GAT_analysis/GAT_K562_", Gat_Seg, '_', PEAK_NUMBER_Cutoff, "Peaks_Jan2019_Gat-Analysis_trim.csv",  sep=""))
  }
  
  
  
  # ====== Aggregate by Target
  # Comparing enrichment over different cell lines as different experiments in K562 resulted in different enrichments. Aggregate by target with both by highest and  average fold enrichment per target.
  # Tagged and fusion proteins will be treated as unmodified proteins
  
  # remove information on tag
  ENCODE_bed_meta$Experiment.target <- gsub("-human", "", ENCODE_bed_meta$Experiment.target)
  
  
  
  GAT_FILTERED_trim_by_Target <- GAT_FILTERED_trim[, c("Experiment.target", "WL_fold", "DHS_fold", "OQs_fold", "opOQs_fold")]
  # remove information on tag
  GAT_FILTERED_trim_by_Target$Experiment.target <- gsub("eGFP-", "", GAT_FILTERED_trim_by_Target$Experiment.target)
  GAT_FILTERED_trim_by_Target$Experiment.target <- gsub("3xFLAG-", "" , GAT_FILTERED_trim_by_Target$Experiment.target)
  
  # aggregate and keep maximum enrichment
  GAT_FILTERED_trim_by_Target_max <- aggregate(GAT_FILTERED_trim_by_Target[,2:5], by = list(Exper_Target_List=GAT_FILTERED_trim_by_Target$Experiment.target), FUN = max , simplify = TRUE)
  
  GAT_FILTERED_trim_by_Target_max <- GAT_FILTERED_trim_by_Target_max[order(GAT_FILTERED_trim_by_Target_max$WL_fold, decreasing = TRUE), ]
  GAT_FILTERED_trim_by_Target_max$WL_rank <- seq.int(nrow(GAT_FILTERED_trim_by_Target_max))
  GAT_FILTERED_trim_by_Target_max <- GAT_FILTERED_trim_by_Target_max[order(GAT_FILTERED_trim_by_Target_max$DHS_fold, decreasing = TRUE), ]
  GAT_FILTERED_trim_by_Target_max$DHS_rank <- seq.int(nrow(GAT_FILTERED_trim_by_Target_max))
  GAT_FILTERED_trim_by_Target_max <- GAT_FILTERED_trim_by_Target_max[order(GAT_FILTERED_trim_by_Target_max$OQs_fold, decreasing = TRUE), ]
  GAT_FILTERED_trim_by_Target_max$OQs_rank <- seq.int(nrow(GAT_FILTERED_trim_by_Target_max))
  GAT_FILTERED_trim_by_Target_max <- GAT_FILTERED_trim_by_Target_max[order(GAT_FILTERED_trim_by_Target_max$opOQs_fold, decreasing = TRUE), ]
  GAT_FILTERED_trim_by_Target_max$opOQs_rank <- seq.int(nrow(GAT_FILTERED_trim_by_Target_max))
  
  
  # aggregate using the mean enrichment
  GAT_FILTERED_trim_by_Target_mean <- aggregate(GAT_FILTERED_trim_by_Target[,2:5], by = list(Exper_Target_List=GAT_FILTERED_trim_by_Target$Experiment.target), FUN = mean , simplify = TRUE)
  
  GAT_FILTERED_trim_by_Target_mean <- GAT_FILTERED_trim_by_Target_mean[order(GAT_FILTERED_trim_by_Target_mean$WL_fold, decreasing = TRUE), ]
  GAT_FILTERED_trim_by_Target_mean$WL_rank <- seq.int(nrow(GAT_FILTERED_trim_by_Target_mean))
  GAT_FILTERED_trim_by_Target_mean <- GAT_FILTERED_trim_by_Target_mean[order(GAT_FILTERED_trim_by_Target_mean$DHS_fold, decreasing = TRUE), ]
  GAT_FILTERED_trim_by_Target_mean$DHS_rank <- seq.int(nrow(GAT_FILTERED_trim_by_Target_mean))
  GAT_FILTERED_trim_by_Target_mean <- GAT_FILTERED_trim_by_Target_mean[order(GAT_FILTERED_trim_by_Target_mean$OQs_fold, decreasing = TRUE), ]
  GAT_FILTERED_trim_by_Target_mean$OQs_rank <- seq.int(nrow(GAT_FILTERED_trim_by_Target_mean))
  GAT_FILTERED_trim_by_Target_mean <- GAT_FILTERED_trim_by_Target_mean[order(GAT_FILTERED_trim_by_Target_mean$opOQs_fold, decreasing = TRUE), ]
  GAT_FILTERED_trim_by_Target_mean$opOQs_rank <- seq.int(nrow(GAT_FILTERED_trim_by_Target_mean))
  
  if (savetables)
  {
    write.csv(GAT_FILTERED_trim_by_Target_max, file =paste("K562_Rerun_Jan2019/GAT_analysis/GAT_K562_", Gat_Seg, "_Jan2019_Gat-Analysis_byTarget-max.csv",  sep=""))
    write.csv(GAT_FILTERED_trim_by_Target_mean, file =paste("K562_Rerun_Jan2019/GAT_analysis/GAT_K562_", Gat_Seg, "_Jan2019_Gat-Analysis_byTarget-mean.csv",  sep=""))
  }
  
  # ==========
  # Plot Data  
  # =========
  if (plotdata)
  {
    #for plots remove all the additional tag information: eGFP, FLAG, etc.
    GAT_FILTERED_trim$Experiment.target <- gsub("eGFP-", "", GAT_FILTERED_trim$Experiment.target)
    GAT_FILTERED_trim$Experiment.target <- gsub("3xFLAG-", "" , GAT_FILTERED_trim$Experiment.target)
    
    # assign colors based on known features
    Known_G4_proteins <- read.csv(file = "K562_Rerun_Jan2019/GAT_analysis/G4IPD_Oct2019.csv", header = F)
    temp <- GAT_FILTERED_trim
    temp$G4IPD <- ifelse(temp$Experiment.target %in% Known_G4_proteins$V1, "#aad2a5", "grey50")
    
    # =====
    # =DHS
    # =====
    
    temp_plot <- temp[order(temp$DHS_fold, decreasing = TRUE), ]
    
    # all
    gg <- ggplot(temp_plot, aes(x=as.factor(rev(DHS_rank)), y=DHS_fold)) + 
      geom_bar(stat="identity", fill=rev(temp_plot$G4IPD) ) +
      xlab("") +
      ylab("fold enrichment in DHS") +
      theme_minimal() +
      scale_x_discrete(labels="") +
      geom_hline(yintercept=5, linetype="dashed", color = "grey")+
      geom_hline(yintercept=10, linetype="dashed", color = "grey")+
      coord_flip()
    ggsave(file =paste("K562_Rerun_Jan2019/GAT_analysis/figures/DHS_all_", Gat_Seg, "_Jan2019_Gat-Analysis_byTarget-mean.pdf",  sep=""), width = 6/2.54, height = 22/2.54, limitsize = FALSE)
    
    
    # top 20
    temp_plot <- temp[order(temp$DHS_fold, decreasing = TRUE), ]
    temp_plot <- temp_plot[1:20,]
    
    gg <- ggplot(temp_plot, aes(x=as.factor(rev(DHS_rank)), y=DHS_fold)) + 
      geom_bar(stat="identity",color="black", fill=rev(temp_plot$G4IPD), width=0.8 ) +
      xlab("") +
      ylab("fold enrichment in DHS") +
      theme_minimal() +
      scale_x_discrete(labels=rev(temp_plot$Experiment.target)) +
      coord_flip()
    ggsave(file =paste("K562_Rerun_Jan2019/GAT_analysis/figures/DHS_top20_", Gat_Seg, "_Jan2019_Gat-Analysis_byTarget-mean.pdf",  sep=""), width = 6/2.54, height = 8/2.54, limitsize = FALSE)
    
    
    # bottom 20
    temp_plot <- temp[order(temp$DHS_fold, decreasing = TRUE), ]
    temp_plot <- tail(temp_plot, 20)
    
    gg <- ggplot(temp_plot, aes(x=as.factor(rev(DHS_rank)), y=DHS_fold)) + 
      geom_bar(stat="identity",color="black", fill=rev(temp_plot$G4IPD), width=0.8 ) +
      xlab("") +
      ylab("fold enrichment in DHS") +
      theme_minimal() +
      scale_x_discrete(labels=rev(temp_plot$Experiment.target)) +
      scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
      #annotation_logticks(sides= "l") + # unfortunately, this seems incompatible with coord_flip() -> have to edit this manually in illustrator :/
      coord_flip()
    ggsave(file =paste("K562_Rerun_Jan2019/GAT_analysis/figures/DHS_bottom20_", Gat_Seg, "_Jan2019_Gat-Analysis_byTarget-mean.pdf",  sep=""), width = 6/2.54, height = 8/2.54, limitsize = FALSE)
    
    
    
    
    
    #####################
    ### other orientation
    
    # all
    temp_plot <- temp[order(temp$DHS_fold, decreasing = TRUE), ]
    
    gg <- ggplot(temp_plot, aes(x=as.factor((DHS_rank)), y=DHS_fold)) + 
      geom_bar(stat="identity", fill=(temp_plot$G4IPD) ) +
      xlab("") +
      ylab("fold enrichment in DHS") +
      theme_minimal() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      scale_x_discrete(labels="") +
      geom_hline(yintercept=5, linetype="dashed", color = "grey30", size = 0.3 )+
      geom_hline(yintercept=10, linetype="dashed", color = "grey30", size = 0.3)
    plot(gg)
    
    ggsave(file =paste("K562_Rerun_Jan2019/GAT_analysis/figures/DHS_all_", Gat_Seg, "_Jan2019_Gat-Analysis_byTarget-mean_rotated.pdf",  sep=""), width = 20/2.54, height = 6/2.54, limitsize = FALSE)
    
    # top 20
    temp_plot <- temp[order(temp$DHS_fold, decreasing = TRUE), ]
    temp_plot <- temp_plot[1:20,]
    
    gg <- ggplot(temp_plot, aes(x=as.factor((DHS_rank)), y=DHS_fold)) + 
      geom_bar(stat="identity",color="black", fill=(temp_plot$G4IPD), width=0.8 ) +
      xlab("") +
      ylab("fold enrichment in DHS") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90))+
      scale_x_discrete(labels=(temp_plot$Experiment.target)) 
    plot(gg)
    
    ggsave(file =paste("K562_Rerun_Jan2019/GAT_analysis/figures/DHS_top20_", Gat_Seg, "_Jan2019_Gat-Analysis_byTarget-mean_rotated.pdf",  sep=""), width = 9/2.54, height = 6/2.54, limitsize = FALSE)
    
    
    
    
    # bottom 20
    temp_plot <- temp[order(temp$DHS_fold, decreasing = TRUE), ]
    temp_plot <- tail(temp_plot, 20)
    
    gg <- ggplot(temp_plot, aes(x=as.factor((DHS_rank)), y=DHS_fold)) + 
      annotation_logticks(sides= "l")+
      geom_bar(stat="identity",color="black", fill=(temp_plot$G4IPD), width=0.8 ) +
      xlab("") +
      ylab("fold enrichment in DHS") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90))+
      scale_x_discrete(labels=(temp_plot$Experiment.target)) +
      scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) 
    # unfortunately, this seems incompatible with coord_flip() -> have to edit this manually in illustrator :/
    plot(gg)
    ggsave(file =paste("K562_Rerun_Jan2019/GAT_analysis/figures/DHS_bottom20_", Gat_Seg, "_Jan2019_Gat-Analysis_byTarget-mean_rotated.pdf",  sep=""), width = 9/2.54, height = 6/2.54, limitsize = FALSE)
    
    
    
    
    # =====
    # =opOQS
    # =====
    
    temp_plot <- temp[order(temp$opOQs_fold, decreasing = TRUE), ]
    
    # all
    gg <- ggplot(temp_plot, aes(x=as.factor(rev(opOQs_rank)), y=opOQs_fold)) + 
      geom_bar(stat="identity", fill=rev(temp_plot$G4IPD) ) +
      xlab("") +
      ylab("fold enrichment in opOQs") +
      theme_minimal() +
      scale_x_discrete(labels="") +
      geom_hline(yintercept=5, linetype="dashed", color = "grey")+
      geom_hline(yintercept=10, linetype="dashed", color = "grey")+
      coord_flip()
    ggsave(file =paste("K562_Rerun_Jan2019/GAT_analysis/figures/opOQs_all_", Gat_Seg, "_Jan2019_Gat-Analysis_byTarget-mean.pdf",  sep=""), width = 6/2.54, height = 22/2.54, limitsize = FALSE)
    
    
    # top 20
    temp_plot <- temp[order(temp$opOQs_fold, decreasing = TRUE), ]
    temp_plot <- temp_plot[1:20,]
    
    gg <- ggplot(temp_plot, aes(x=as.factor(rev(opOQs_rank)), y=opOQs_fold)) + 
      geom_bar(stat="identity",color="black", fill=rev(temp_plot$G4IPD), width=0.8 ) +
      xlab("") +
      ylab("fold enrichment in opOQs") +
      theme_minimal() +
      scale_x_discrete(labels=rev(temp_plot$Experiment.target)) +
      coord_flip()
    ggsave(file =paste("K562_Rerun_Jan2019/GAT_analysis/figures/opOQs_top20_", Gat_Seg, "_Jan2019_Gat-Analysis_byTarget-mean.pdf",  sep=""), width = 6/2.54, height = 8/2.54, limitsize = FALSE)
    
    
    # bottom 20
    temp_plot <- temp[order(temp$opOQs_fold, decreasing = TRUE), ]
    temp_plot <- tail(temp_plot, 20)
    
    gg <- ggplot(temp_plot, aes(x=as.factor(rev(opOQs_rank)), y=opOQs_fold)) + 
      geom_bar(stat="identity",color="black", fill=rev(temp_plot$G4IPD), width=0.8 ) +
      xlab("") +
      ylab("fold enrichment in opOQs") +
      theme_minimal() +
      scale_x_discrete(labels=rev(temp_plot$Experiment.target)) +
      scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
      #annotation_logticks(sides= "l") + # unfortunately, this seems incompatible with coord_flip() -> have to edit this manually in illustrator :/
      coord_flip()
    ggsave(file =paste("K562_Rerun_Jan2019/GAT_analysis/figures/opOQs_bottom20_", Gat_Seg, "_Jan2019_Gat-Analysis_byTarget-mean.pdf",  sep=""), width = 6/2.54, height = 8/2.54, limitsize = FALSE)
    
    
    # =====
    # = WL
    # =====
    
    temp_plot <- temp[order(temp$WL_fold, decreasing = TRUE), ]
    
    # all
    gg <- ggplot(temp_plot, aes(x=as.factor(rev(WL_rank)), y=WL_fold)) + 
      geom_bar(stat="identity", fill=rev(temp_plot$G4IPD) ) +
      xlab("") +
      ylab("fold enrichment in WL") +
      theme_minimal() +
      scale_x_discrete(labels="") +
      geom_hline(yintercept=100, linetype="dashed", color = "grey")+
      geom_hline(yintercept=200, linetype="dashed", color = "grey")+
      geom_hline(yintercept=300, linetype="dashed", color = "grey")+
      coord_flip()
    ggsave(file =paste("K562_Rerun_Jan2019/GAT_analysis/figures/WL_all_", Gat_Seg, "_Jan2019_Gat-Analysis_byTarget-mean.pdf",  sep=""), width = 6/2.54, height = 22/2.54, limitsize = FALSE)
    
    
    # top 20
    temp_plot <- temp[order(temp$WL_fold, decreasing = TRUE), ]
    temp_plot <- temp_plot[1:20,]
    
    gg <- ggplot(temp_plot, aes(x=as.factor(rev(WL_rank)), y=WL_fold)) + 
      geom_bar(stat="identity",color="black", fill=rev(temp_plot$G4IPD), width=0.8 ) +
      xlab("") +
      ylab("fold enrichment in WL") +
      theme_minimal() +
      scale_x_discrete(labels=rev(temp_plot$Experiment.target)) +
      coord_flip()
    ggsave(file =paste("K562_Rerun_Jan2019/GAT_analysis/figures/WL_top20_", Gat_Seg, "_Jan2019_Gat-Analysis_byTarget-mean.pdf",  sep=""), width = 6/2.54, height = 8/2.54, limitsize = FALSE)
    
    
    # bottom 20
    temp_plot <- temp[order(temp$WL_fold, decreasing = TRUE), ]
    temp_plot <- tail(temp_plot, 20)
    
    gg <- ggplot(temp_plot, aes(x=as.factor(rev(WL_rank)), y=WL_fold)) + 
      geom_bar(stat="identity",color="black", fill=rev(temp_plot$G4IPD), width=0.8 ) +
      xlab("") +
      ylab("fold enrichment in WL") +
      theme_minimal() +
      scale_x_discrete(labels=rev(temp_plot$Experiment.target)) +
      scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
      #annotation_logticks(sides= "l") + # unfortunately, this seems incompatible with coord_flip() -> have to edit this manually in illustrator :/
      coord_flip()
    ggsave(file =paste("K562_Rerun_Jan2019/GAT_analysis/figures/WL_bottom20_", Gat_Seg, "_Jan2019_Gat-Analysis_byTarget-mean.pdf",  sep=""), width = 6/2.54, height = 8/2.54, limitsize = FALSE)
    
  }

}
