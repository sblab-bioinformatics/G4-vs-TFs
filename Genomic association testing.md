
- filter releavant ENCODE files (Script);
- download ENCODE files -> upload ENCODE file

- prepare workspaces (DHS, OQS, opOQS)


# Genomic association testing
 


## Selection and download of relevant ENCODE ChIP-seq bed files
Genomic binding sites for chromatin-associated factors and histone marks (aligned to hg19) were downloaded from [ENCODE](https://www.encodeproject.org/) (Meta-Data was generated on 10.01.2019 ). All available ChIP-seq data sets were considered.
> filter: 'cell line': "K562"  


## Selection of relevant meta dataset
To maximise the robustness of our analyses, where possible we selected  ‘released’ over ‘archived’ data and generally chose ‘optimal idr’ over ‘conservative idr’ and ‘replicated’ peaks. Different ChIP-seq experiments targeting the same factor were treated independently. 

```R
R

# ====== Load meta data
ENCODE_K562_meta <- read.table(file = "ENCODE_K562_Jan2019_all_Meta.tsv", sep = '\t', header = TRUE, fill=T)


# focus on hg19, ChIP-seq and remove revoked data
ENC_filt <- ENCODE_K562_meta
ENC_filt <- ENC_filt[ENC_filt$File.Status!='revoked',]
ENC_filt <- ENC_filt[ENC_filt$Assay=='ChIP-seq',]
ENC_filt <- ENC_filt[ENC_filt$Assembly=='hg19',]


# ==== select bed files
ENC_filt <- ENC_filt[grep('bed', ENC_filt$File.format),]

# prefer 'released' over archived
Released <-  ENC_filt[ENC_filt$File.Status == "released", ] 
Remaining <- Remaining[!(Remaining$Experiment.accession %in% Released$Experiment.accession),] # keep unique experiments that only have an "archived" but not "released" dataset
ENC_filt <- rbind(Released, Remaining) 


# Identify all experiments for which there is an 'optimal idr' available  
IDR_opt_yes <- ENC_filt[ENC_filt$Output.type == "optimal idr thresholded peaks", ]

# these are the ENCODE gold standard peaks. However, some experiments (in particular, histone marks) do not provide these peaks 

Remaining <- ENC_filt[!(ENC_filt$Experiment.accession %in% IDR_opt_yes$Experiment.accession),] # Keep Experimental.acession that are not present in IDR_opt_yes

# Keep files that result from at least 2 biological replictate and remove experiments that contain only one sample:
Remaining <- Remaining[grep(', ', Remaining$Biological.replicate.s.),]

# Keep only those processed by ENCODE consortium
Remaining <- Remaining[Remaining$Lab =="ENCODE Processing Pipeline" ,]
Remaining <- Remaining[rowSums(is.na(Remaining))  != ncol(Remaining), ] # remove empty rows

# Keep the 'conservative-idr'
IDR_conserv <- Remaining[Remaining$Output.type == "conservative idr thresholded peaks", ] # extract conservative IDR
Remaining <- Remaining[!(Remaining$Experiment.accession %in% IDR_conserv$Experiment.accession),]   #keep remaining experiments

# If there are multiple peaks per experiment keep 'replicated peaks' over 'peaks'
Remaining <- Remaining[order(Remaining$Experiment.accession, Remaining$Output.type, decreasing = TRUE),  ] # order such that, files are grouped by Experiment.Accession number and 'replicated' are in the 1st position
Remaining <- Remaining[!duplicated(Remaining$Experiment.accession),]  



# Merge the selected Experiments
ENC_filt <- rbind(IDR_opt_yes, IDR_conserv, Remaining)

#---- Generate list for Download of relevant bed files
write.table(ENC_filt[,'File.download.URL'], "K562_Bed_Jan2019/20190117_LINKS_K562_bedfiles.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#--- Save filtered meta data for GAT analysi 
write.table(ENC_filt,'ENCODE_K562_Jan2019_BED_Meta.tsv', sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

```


### Download  selected bed files

```
cd K562_Bed_Jan2019
sbatch -o Download.%j.log -J ENCODE_download --mem=10000 --wrap "xargs -n 1 curl -O -L < 20190117_LINKS_K562_bedfiles.txt"
```

# GAT shuffling analysis

## Prepare K562 bedfiles

Cut column 1-3 to be compatible with GAT analysis:

```
cd K562_Bed_Jan2019

for FILE in *.bed.gz
do
bname=`basename $FILE .bed.gz`
echo $bname
zcat $FILE | cut -f 1-3 | bedtools sort -i |  gzip > ../K562_Cut_Bed_Jan2019/$bname.cut.bed.gz
done

```

Generate list of bed files for GAT job scripts:
```
cd ../K562_Cut_Bed_Jan2019

for FILE in *.bed.gz
do
echo '--annotations=/scratchb/sblab/spiege01/ENCODE_K562/GAT_Rerun_Jan2019/K562_Cut_Bed_Jan2019/'"$FILE"' \'
done >> K562_Rerun_Jan2019_annotations_list.txt
```

### Produce different workspacede for randomization









Run shuffling on cluster
```
cd Job_files
for Job in *.sh; do sbatch $Job; done
```

# Analysis 
### peak numbers

```
cd /scratchb/sblab/spiege01/ENCODE_K562/GAT_Rerun_Jan2019/K562_Cut_Bed_Jan2019
srun --mem 8g --pty /usr/bin/bash

for FILE in *.bed.gz
do
Bedname=`basename $FILE .bed.gz`
Peaks=`zcat $FILE | wc -l`
echo -e $Bedname'\t'$Peaks >> K562_Jan2019_NumberOfPeaks.txt
echo $FILE
done

exit
```

# GAT analysis in proximity of G4s
Aim: Expand BG4 peaks by 500bp flanking regions and see if that considerably changes the ranking/enrichment. Maybe there is a protein that only binds in the proximity of G4s but not directly.

### Generate expanded BG4 consensuspeaks

```
cd /scratchb/sblab/spiege01/ENCODE_K562/reference_data/Quadruplex/

bedtools slop -i 20180108_K562_async_rep1-3.mult.5of8.bed -b 500 -g /scratchb/sblab/spiege01/ENCODE_K562/reference_data/hg19/hg19.genome | bedtools sort -i - > 20180108_K562_async_rep1-3.mult.5of8_slop_b_500.bed

```

Generated jobfiles accordingly (see scripts)
```
for JOB in *slop*.sh; do sbatch  $JOB; done
```

Combined [R analysis](https://github.com/sblab-bioinformatics/projects/blob/master/20171123_Jochen_ENCODE/Shuffling-Analysis/K562_Update_Jan2019/Scripts/K562_Rerun_Jan2019_GAT_analysis.R) .



# GAT analysis - disentangle G-richness
## Run shuffling on cluster
Assess influence of inherent G-richness of BG4 ChIP sites. A similar analysis had [previously](https://github.com/sblab-bioinformatics/projects/blob/master/20171123_Jochen_ENCODE/Shuffling-Analysis/20180217_GAT_analysis_Disentangle_G-richness_in_BG4_regions.md) been performed for the pilot data sets.
Repeat with updated list ENCODE bed files. See [Scripts](Scripts/) (GAT_opOQs-noBG4_atTSS_Update_2019Jan_DHS.sh, GAT_opOQs-noBG4_Update_2019Jan_DHS.sh, GAT_opOQs-WithBG4_Update_2019Jan_DHS.sh, GAT_opOQs-noBG4_1kbupstreamTSS_Update_2019Dec_DHS.sh)

Run Shuffling on cluster
```
cd /scratchb/sblab/spiege01/ENCODE_K562/GAT_Rerun_Jan2019/Job_files/GAT_Grichness
for JOB in *opOQs*.sh; do sbatch  $JOB; done
```
See [R analysis]. Technically, ranking based on enrichemnt is not representative, as a couple of very lowly enriched factors are even more depleted for NoBG4.






