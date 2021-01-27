## TF ChIP-seq signal at endogenous G4s

Plot the  TF ChIP-seq signal around endogenous G4s (BG4 ChIP-seq peaks) and control sites (potential G4s in open chromatin promoters that do not form endogenous G4s).
In addition, investigate the distribution of TF binding sites around endogenous G4s.

#### Generate BigWig files
Obtain one TF ChIP-seq bam-file from ENCODE for each of the [selected experiments in K562](https://github.com/sblab-bioinformatics/G4-vs-TFs/blob/6b1da1acae7f35bc39d2b77c853ed67bb7bce86f/Genomic_association_testing.md#selection-of-relevant-meta-data-sets)

```R
library(tidyverse)

# ====== Load meta data
ENCODE_K562_meta <- read_tsv(file = "ENCODE_K562_Jan2019_all_Meta.tsv", col_names = T)
BED_considered <- read_csv(file ="GAT_analysis/GAT_K562_BG4_100Peaks_Jan2019_Gat-Analysis_trim.csv")

# keep bam files from hg19 that were used for considered Bed files
# this still contains multiple replicates and potentially archived files; prefer low number replicates and "released" over "archived"; keep only one bam file each

ENC_filt <-  ENCODE_K562_meta %>% filter(`File format`=='bam', `Output type`=='alignments', Assembly=='hg19', `Experiment accession` %in%  BED_considered$Experiment.accession) %>%
  arrange(`Experiment accession`, desc(`File Status`), `Biological replicate(s)`) %>% 
  group_by(`Experiment accession`) %>% slice(1)  %>% ungroup() %>% 
  select(`File download URL`) %>%
  write_tsv("bam_download/20201220_LINKS_K562_bamfiles_1replicate.txt", col_names = F)

```
Download files

```bash 
sbatch --mem 4G -o DL.out -e DL.err -J bam_DL --wrap "xargs -n 1 curl -O -L < 20201220_LINKS_K562_bamfiles_1replicate.txt"

``` 
Generate bigWig-files normalized to sequencing depth generated using the *bamCoverage* function in *deeptools*

```bash 
# Index bam
for FILE in *.bam
do
sbatch -o index.%j.log -J BamIndex --mem=5G \
--wrap "samtools index $FILE"
done

# Coverage files at 10 nt bin size
for FILE in *.bam
do
bname=`basename $FILE .bam`
sbatch -o log/%j.$f.log -J BCov --mem=8G \
--wrap "bamCoverage -b $FILE \
-o BigWig_10nt/$bname.bs10.bl.RPKM.bw \
--binSize 10 \
--blackListFileName reference_data/hg19/hg19.blacklist_merge1000nt.bed \
--normalizeUsing RPKM"
done

```

#### Signal plots 

Generate signal profile plots 400bp around BG4 ChIP-seq peak centers and respective control regions (potential G4s in open chromatin promoters that do not form endogenous G4s) using functions *computeMatrix* and *plotProfile*.

```bash 
cd BigWig_10nt

for f in *.bs10.bl.RPKM.bw
do
sbatch -o logs/${f%%.bs10.bl.RPKM.bw}.log -J DT_Profile --mem=5G \
--wrap "computeMatrix reference-point \
--referencePoint center \
-b 400 -a 400 \
-S \
$f \
-R \
reference_data/Quadruplex/K562_async_rep1-3.mult.5of8-OQS_Stranded.bed \
reference_data/openOQs_noBG4_Stranded_1kbupTSS.bed \
--skipZeros \
-o ../ChIP_Profile/ENCODE2019_bam_1rep_2020/${f%%.bs10.bl.RPKM.bw}.mat.gz && 
plotProfile -m ../ChIP_Profile/ENCODE2019_bam_1rep_2020/${f%%.bs10.bl.RPKM.bw}.mat.gz \
-out ../ChIP_Profile/ENCODE2019_bam_1rep_2020/${f%%.bs10.bl.RPKM.bw}.pdf \
--refPointLabel PeakCenter \
--regionsLabel G4 OQ \
--outFileNameData ../ChIP_Profile/ENCODE2019_bam_1rep_2020/${f%%.bs10.bl.RPKM.bw}_profile.tab \
--colors green grey \
--yMin 0 \
--plotHeight 6 \
--plotWidth 6" 
done

```

#### Distribution of TF signal maxima around BG4 ChIP-seq peaks

Binned average signal intensity for each TF experiment are saved int *.tab. Extract the information on signal around BG4 ChIP-seq sites (3rd line in each tab file):

```bash 
cd ../ChIP_Profile/ENCODE2019_bam_1rep_2020/

for f in *.tab
do
bname=`basename $f _profile.tab`
sed -n 3p $f > ../../Signal_around_BG4/$bname.tab
done
```

Plot data in using a custom R script

```R
library(tidyverse)
library(readr)
library(ggsignif)
library(ggrepel)

#################################
# <<<<<<  G4s  >>>>>>>>>
#############################
# Fetch data
#generate master file of average TF signal around BG4 peaks
Avr_Signal_at_G4 <- list.files(path="Signal_around_BG4/", full.names = TRUE) %>% lapply(read_tsv) %>% bind_rows 

# reformat and keep relevant columns
Avr_Signal_at_G4 <- Avr_Signal_at_G4 %>% separate(col = "bins", into = c("ENCODE_bam_ID", "B"), sep = ".b")  %>% select(-2, -3)

# find  position and value of maximum count (there are 1-80 '10nt'-bins representing the -400bp to +400bp window)
Max_signal <- Avr_Signal_at_G4 %>% 
  pivot_longer(`1`:`80`, names_to = "position", values_to = "cnt") %>% 
  group_by(ENCODE_bam_ID) %>% 
  slice(which.max(cnt)) %>%
  mutate(position = as.integer(position)*10 -400)    # col1-80 represent -400 to 400 bp; therefore, pos = x*10 -400


# Integrate information on genomic enrichment from GAT analysis 
ENCODE_K562_meta_bam <- read_tsv(file = "ENCODE_K562_Jan2019_BAM_Meta.tsv", col_names = T)
GAT_K562 <- read_csv(file = "GAT_K562_BG4_100Peaks_Jan2019_Gat-Analysis_trim.csv", col_names = T)

Max_signal$ENCODE_exp_ID <- ENCODE_K562_meta_bam$`Experiment accession`[match(Max_signal$ENCODE_bam_ID, ENCODE_K562_meta_bam$`File accession`)]
Max_signal$target <- GAT_K562$Experiment.target[match(Max_signal$ENCODE_exp_ID, GAT_K562$Experiment.accession)]
Max_signal$GAT_DHS <- GAT_K562$DHS_fold[match(Max_signal$ENCODE_exp_ID, GAT_K562$Experiment.accession)]


##################
#### Plot data ###
##################

Signal_Plot <- Max_signal

# Visual inspection suggests that there are a few noisy, low intensity data sets 
# these can easily filtered using a cutoff for intensity
Signal_Plot <- Signal_Plot %>% filter(cnt > 5)


# Remove histone marks and PolR2
Signal_Plot <- Signal_Plot %>% filter(!grepl('H2A|H3|H4|POLR', target))

# Remove  protein-tag information  "eGFP", "3xFLAG" etc
Signal_Plot  <- Signal_Plot %>% mutate(target = str_remove_all(target, "eGFP-|3xFLAG-"))

# Highlight TFs within a 20 bp window around center
Signal_Plot <-  Signal_Plot %>% mutate( In_20bp_window = ifelse(position >= -20 & position <=20, "blue", "black"))

# check how many TFs bind within a 20, 40 or 60 bp window
Signal_Plot %>% filter(position == 0 ) %>% nrow() #7
Signal_Plot %>% filter(position >= -10 & position <=10) %>% nrow() #40
Signal_Plot %>% filter(position >= -20 & position <=20) %>% nrow() #100
Signal_Plot %>% filter(position >= -30 & position <=30) %>% nrow() #148
Signal_Plot %>% filter(position >= -40 & position <=40) %>% nrow() #177

# remove TFs that don't have a maximum inside the explored 400bp interval
Signal_Plot  <- Signal_Plot %>% filter(position > -400 & position < 400) 

# Plot G4 association
gg_G4s_OQS <- ggplot(Signal_Plot, aes(x=position, y=GAT_DHS)) + 
  geom_point(color=Signal_Plot$In_20bp_window) + 
  scale_color_gradient2(midpoint = 20, low = "blue", mid = "black", high = "red", space = "Lab" ) +
  theme_minimal() +
  labs(title="TF position vs GAT DHS", y="enrichment at endogenous G4s", x="TF signal maximum relative to G4") +
  geom_label_repel(aes(label=ifelse(
    (GAT_DHS > 5 & position >= -20 & position <= 20), as.character(target),'')),
    box.padding = 0.4, point.padding = 0.4, size=rel(2.5), segment.color = 'grey80', color=Signal_Plot$In_20bp_window) +
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA))

plot(gg_G4s_OQS) 

ggsave("figures/TFPosition_G4_400bp.pdf", width = 17/2.54, height = 14/2.54, useDingbats=F)


``






