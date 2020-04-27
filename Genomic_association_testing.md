# Genomic association testing of G4 ChIP-seq high confidence peaks
 
## Selection and download of relevant ENCODE ChIP-seq bed files

Genomic binding sites for chromatin-associated factors and histone marks (aligned to hg19) were downloaded from [ENCODE](https://www.encodeproject.org/) (Meta-Data was generated on 10.01.2019). All available ChIP-seq data sets from cell line K562 were considered.


## Selection of relevant meta-data sets

To maximise the robustness of our analyses, where possible we selected  ‘released’ over ‘archived’ data and generally chose ‘optimal idr’ over ‘conservative idr’ and ‘replicated’ peaks. Different ChIP-seq experiments targeting the same factor were treated independently. 

```R
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

# Keep files that result from at least 2 biological replicates and remove experiments that contain only one sample:
Remaining <- Remaining[grep(', ', Remaining$Biological.replicate.s.),]

# Keep only those processed by ENCODE consortium
Remaining <- Remaining[Remaining$Lab =="ENCODE Processing Pipeline" ,]
Remaining <- Remaining[rowSums(is.na(Remaining))  != ncol(Remaining), ] # remove empty rows

# Keep the 'conservative-idr'
IDR_conserv <- Remaining[Remaining$Output.type == "conservative idr thresholded peaks", ] # extract conservative IDR
Remaining <- Remaining[!(Remaining$Experiment.accession %in% IDR_conserv$Experiment.accession),] # keep remaining experiments

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


## Download  selected bed files

```R
cd K562_Bed_Jan2019/
xargs -n 1 curl -O -L < 20190117_LINKS_K562_bedfiles.txt
```



# GAT shuffling analysis

## Prepare K562 bedfiles

Cut columns 1-3 to be compatible with GAT analysis:

```bash
cd K562_Bed_Jan2019/

for FILE in *.bed.gz
do
bname=`basename $FILE .bed.gz`
echo $bname
zcat $FILE | cut -f 1-3 | bedtools sort -i |  gzip > ../K562_Cut_Bed_Jan2019/$bname.cut.bed.gz
done
```


## Generate list of bed files for GAT job scripts:

```bash
cd ../K562_Cut_Bed_Jan2019

for FILE in *.bed.gz
do
echo '--annotations=/scratchb/sblab/spiege01/ENCODE_K562/GAT_Rerun_Jan2019/K562_Cut_Bed_Jan2019/'"$FILE"' \'
done >> K562_Rerun_Jan2019_annotations_list.txt
```


## peak numbers in each file

```bash
for FILE in *.bed.gz
do
Bedname=`basename $FILE .bed.gz`
Peaks=`zcat $FILE | wc -l`
echo -e $Bedname'\t'$Peaks >> K562_Jan2019_NumberOfPeaks.txt
done
```


## Workspaces

Several different workspaces were considered for randomization using GAT.

- white-listed genome (hg19): [hg19.wgEncodeDukeMapabilityRegionsExcludable.whitelist.bed](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/)

- sites with G4 forming potential (OQS from G4-seq): [G4-Seq_cat_K_PDS_+-strands.bed](G4-ChIP-seq.md#stranded-oqs-map)

- open Chromatin (ENCODE K562 DNase-seq): [DNAse-seq.concatenated_narrow_rep1_and_rep2.bed](https://www.encodeproject.org/experiments/ENCSR000EPC/)

- G4-Seq OQs in open chromatin: [OQs_in_K562_open_chromatin.bed](G4-ChIP-seq.md#g4-control-data-set-oqs-in-open-chromatin-around-tss)

In addition, chrM was not compatible with GAT tool and needed to be removed from workspaces.

```bash
grep -v chrM hg19.wgEncodeDukeMapabilityRegionsExcludable.whitelist.bed > Whitelist_chrM.bed
grep -v chrM G4-Seq_cat_K_PDS_+-strands.bed > OQs_chrM.bed
grep -v chrM DNAse-seq.concatenated_narrow_rep1_and_rep2.bed > DHS_chrM.bed
grep -v chrM OQs_in_K562_open_chromatin.bed > openOQs_chrM.bed
```


## Randomization and statistical analysis

Scripts for indiviual shuffling analysis.

[GAT_K562_ReRun2019_WL.sh](Scripts/GAT_K562_ReRun2019_WL.sh)
[GAT_K562_ReRun2019_OQS.sh](Scripts/GAT_K562_ReRun2019_OQS.sh)
[GAT_K562_ReRun2019_DHS.sh](Scripts/GAT_K562_ReRun2019_DHS.sh)
[GAT_K562_ReRun2019_opOQS.sh](Scripts/GAT_K562_ReRun2019_opOQS.sh)


## Visualization

Results for different randomizations were combined and visualized [using a custom R script](Scripts/GAT-analysis.R).





# HepG2
Essentially, ENCODE data for HepG2 was analysed as described for K562. Briefly:

## Selection and download of relevant ENCODE ChIP-seq bed files

Genomic binding sites for chromatin-associated factors and histone marks (aligned to hg19) were downloaded from [ENCODE](https://www.encodeproject.org/) (Meta-Data was generated on 09.03.2020). All available ChIP-seq data sets from cell line HepG2 were considered.

```R
# ====== Load meta data
ENCODE_HepG2_meta <- read.table(file = "ENCODE_HepG2_March2020_Meta.tsv", sep = '\t', header = TRUE)

ENC_filt <- ENCODE_HepG2_meta[grep('bed broadPeak|bed narrowPeak', ENCODE_HepG2_meta$File.format),]
ENC_filt <- ENC_filt[ENC_filt$File.Status!='revoked',]
ENC_filt <- ENC_filt[ENC_filt$Assay=='ChIP-seq',]
ENC_filt <- ENC_filt[ENC_filt$Assembly=='hg19',]
ENC_filt <- ENC_filt[ENC_filt$Biosample.treatments=='',]

# Several data sets have been 'archived' and 'released' need to remove

# prefer 'released' over archived
Released <-  ENC_filt[ENC_filt$File.Status == "released", ] 
Remaining <- ENC_filt[!(ENC_filt$Experiment.accession %in% Released$Experiment.accession),] # keep unique experiments that only have an "archived" but not "released" dataset

ENC_filt <- rbind(Released, Remaining) 

# Identify all experiments for which there is an 'optimal idr' available  
IDR_opt_yes <- ENC_filt[ENC_filt$Output.type == "optimal IDR thresholded peaks", ]

# these are the ENCODE gold standard peaks. However, some experiments (in particular, histone marks) do not provide these peaks 

Remaining <- ENC_filt[!(ENC_filt$Experiment.accession %in% IDR_opt_yes$Experiment.accession),] # Keep Experimental.acession that are not present in IDR_opt_yes

# Keep files that result from at least 2 biological replicates and remove experiments that contain only one sample:
Remaining <- Remaining[grep(', ', Remaining$Biological.replicate.s.),]

# Keep only those processed by ENCODE consortium
Remaining <- Remaining[Remaining$Lab =="ENCODE Processing Pipeline" ,]
Remaining <- Remaining[rowSums(is.na(Remaining))  != ncol(Remaining), ] # remove empty rows

# Keep the 'conservative-idr'
IDR_conserv <- Remaining[Remaining$Output.type == "conservative IDR thresholded peaks", ] # extract conservative IDR
Remaining <- Remaining[!(Remaining$Experiment.accession %in% IDR_conserv$Experiment.accession),] # keep remaining experiments

# If there are multiple peaks per experiment keep 'replicated peaks' over 'peaks'
Remaining <- Remaining[order(Remaining$Experiment.accession, Remaining$Output.type, decreasing = TRUE),  ] # order such that, files are grouped by Experiment.Accession number and 'replicated' are in the 1st position
Remaining <- Remaining[!duplicated(Remaining$Experiment.accession),]  


# Merge the selected Experiments
ENC_filt <- rbind(IDR_opt_yes, IDR_conserv, Remaining)


#---- Generate list for Download of relevant bed files
write.table(ENC_filt[,'File.download.URL'], "20200309_LINKS_HepG2_bedfiles.NEW.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(ENC_filt,'ENCODE_HepG2_Mar2020_BED_Meta.tsv', sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
```

## Download  selected bed files
```R
cd HepG2_Bed_Jan2019/
xargs -n 1 curl -O -L < ENCODE_HepG2_Mar2020_BED_Meta.tsv
```
# GAT shuffling analysis in HepG2

## Prepare HepG2 bedfiles

Cut columns 1-3 to be compatible with GAT analysis:

```bash
cd HepG2/bedfiles/

for FILE in *.bed.gz
do
bname=`basename $FILE .bed.gz`
echo $bname
zcat $FILE | cut -f 1-3 | bedtools sort -i |  gzip > ../Cut_bedfiles/$bname.cut.bed.gz
done
```


## Generate list of bed files for GAT job scripts:

```bash
cd ../Cut_bedfiles

for FILE in *.bed.gz
do
echo '--annotations=/scratchb/sblab/spiege01/ENCODE_more_celllines/HepG2/Cut_bedfiles/'"$FILE"' \'
done >> HepG2_2020_annotations_list.txt
```


## peak numbers in each file

```bash
for FILE in *.bed.gz
do
Bedname=`basename $FILE .bed.gz`
Peaks=`zcat $FILE | wc -l`
echo -e $Bedname'\t'$Peaks >> HepG2_2020_NumberOfPeaks.txt
done
```

## Workspaces
As in K562, four different workspaces were considered for randomization using GAT:

- white-listed genome (hg19): [hg19.wgEncodeDukeMapabilityRegionsExcludable.whitelist.bed](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/)

- sites with G4 forming potential (OQS from G4-seq): [G4-Seq_cat_K_PDS_+-strands.bed](G4-ChIP-seq.md#stranded-oqs-map)

- open Chromatin (ENCODE HepG2 DNase-seq): [ENCFF571RHF_DHS.bed.gz](https://www.encodeproject.org/experiments/ENCSR000EJV/)

- G4-Seq OQs in open chromatin: [OQs_in_HepG2_open_chromatin.bed]()

chrM was not compatible with GAT tool and needed to be removed from workspaces.

```bash
zcat ENCFF571RHF_DHS.bed.gz | cut -f 1-3 | grep -v chrM - > DHS_chrM.bed
grep -v chrM OQs_in_HepG2_open_chromatin.bed > openOQs_chrM.bed
```

## Randomization and statistical analysis in HepG2

Scripts for indiviual shuffling analysis.

[GAT_HepG2_WL.sh](Scripts/GAT_HepG2_WL.sh)
[GAT_HepG2_DHS.sh](Scripts/GAT_HepG2_DHS.sh)
[GAT_HepG2_OQS.sh](Scripts/GAT_HepG2_OQS.sh)
[GAT_HepG2_opOQS.sh](Scripts/GAT_HepG2_opOQS.sh)


# G4 secondary structure vs primary sequence

G4-ChIP sites are inherently G-rich. Perform association analysis for [control sites](G4-ChIP-seq.md#g4-control-data-set-oqs-in-open-chromatin-around-tss) that do not form G4 structures, but otherwise share very similar genomic features. (open chromatin, G4-forming potential, around TSS). The script used is [GAT_opOQs-noBG4_1kbupstreamTSS__DHS.sh](Scripts/GAT_opOQs-noBG4_1kbupstreamTSS__DHS.sh). Compare enrichment in endogenous G4s to control sites using [custom R script](Scripts/G4structure_vs_sequence.R)














# G4 binding vs B-DNA consensus
Compare TF chromatin occupancy at G4 ChIP sites and TF consensus motifs.

- retrieve TF motifs from JASPAR
- Scan hg19 for TF motif matches 
- Generate bed files containing motif matches
- GAT analysis 


## retrieve TF motifs from JASPAR

```bash
mkdir -p jaspar && cd jaspar
wget http://jaspar.genereg.net/download/CORE/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt
mv JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt JASPAR2020_CORE_vertebrates_non-redundant.meme
grep "MOTIF" JASPAR2020_CORE_vertebrates_non-redundant.meme | wc -l # 746
grep -E "SP1|SP2" JASPAR2020_CORE_vertebrates_non-redundant.meme -A 20

```


## Scan hg19 for TF motif matches

### Generate background model

```bash
cd reference_data/genomes/gencode/Gencode_human/release_19/fasta/
~/sw/meme/meme_5.0.5/libexec/meme-5.0.5/fasta-get-markov GRCh37.p13.genome.fa > GRCh37.p13.genome.bg
cat GRCh37.p13.genome.bg
# 0-order Markov frequencies from file GRCh37.p13.genome.fa
# seqs: 297    min: 4262    max: 249250621    avg: 10891699.3    sum: 3234834689    alph: DNA
# order 0
# A 2.949e-01
# C 2.051e-01
# G 2.051e-01
# T 2.949e-01
```


### fimo

```bash
cd jaspar

ref_bg=reference_data/genomes/gencode/Gencode_human/release_19/fasta/GRCh37.p13.genome.bg
ref_fa=reference_data/genomes/gencode/Gencode_human/release_19/fasta/GRCh37.p13.genome.fa

# each motif individually
grep "MOTIF" JASPAR2020_CORE_vertebrates_non-redundant.meme | while read line
do
  id=`echo $line | cut -d " " -f 2`
  gene_name=`echo $line | cut -d " " -f 3`
  #echo $id, $gene_name
  sbatch -J $id -o $id.log --mem 8G --wrap "~/sw/meme/meme_5.0.5/bin/fimo --bfile $ref_bg --max-stored-scores 10000000 --motif $id --oc $id JASPAR2020_CORE_vertebrates_non-redundant.meme $ref_fa"
done

# all motifs simultaneously
sbatch -J all -o all.log --mem 32G --wrap "~/sw/meme/meme_5.0.5/bin/fimo --bfile $ref_bg --max-stored-scores 10000000 --oc all JASPAR2020_CORE_vertebrates_non-redundant.meme $ref_fa"
# it just takes too long - this was killed
grep "Using motif +" all.log | wc -l # 54, only managed this small amount of motifs
```


### Generate bed files containing motif matches

```bash
srun --mem 32G --pty /usr/bin/bash

cd jaspar

for tsv in `find . -name "fimo.tsv"`
do
  bed=${tsv%.tsv}.bed
#  echo $tsv, $bed
  nohup tail -n +2 $tsv | \
  awk '{ if ($3 ~ "chr") {print $3 "\t" $4 "\t" $5 "\t" $8 "\t" $9 "\t" $6 "\t" $2} }' | \
  bedtools sort -i > $bed &
done

for bed in `find . -name "fimo.bed"`; do wc -l $bed; done

exit
```



## GAT analysis G4 vs JASPAR

### Obtain BG4 peaks, OQs, DNAse-seq regions, ENCODE files and gene promoters

```bash

# BG4 peaks
mkdir bg4 && cd bg4

# BG4
cp /scratchb/sblab/spiege01/ENCODE_K562/reference_data/Quadruplex/20180108_K562_async_rep1-3.mult.5of8.bed .


# OQs
mkdir ../oqs && cd ../oqs
cp /scratchb/sblab/spiege01/ENCODE_K562/reference_data/Quadruplex/G4-Seq_cat_K_PDS_+-strands.bed .


# DNAse-seq regions
mkdir ../dnaseseq && cd ../dnaseseq
cp /scratchb/sblab/spiege01/ENCODE_K562/reference_data/K562/DNAse-seq.concatenated_narrow_rep1_and_rep2.bed .


# ENCODE bed files and metadata
mkdir ../encode && cd ../encode
cp /scratchb/sblab/spiege01/ENCODE_K562/GAT_Rerun_Jan2019/K562_Cut_Bed_Jan2019/*.cut.bed.gz .


rsync -arvuP martin03@10.20.192.94:/mnt/nfs/nas/group_folders/spiege01/projects/20171123_ENCODE_Analysis/K562_Rerun_Jan2019/GAT_analysis/GAT_K562_BG4_Rerun_500Peaks_Jan2019_Gat-Analysis_trim.csv .

rsync -arvuP martin03@10.20.192.94:/mnt/nfs/nas/group_folders/spiege01/projects/20171123_ENCODE_Analysis/K562_Rerun_Jan2019/ENCODE_K562_Jan2019_BED_Meta.tsv .
```

GenomicFeatures hg19 promoter definition:

```r
library(data.table)
library(GenomicFeatures)

# change width
options(width = 250)

# prepare coordinates table
txdb <- makeTxDbFromGFF("/scratcha/sblab/martin03/reference_data/genomes/gencode/Gencode_human/release_19/gtf/gencode.v19.annotation.gtf", format="gtf")

# gene promoters
gene_promoters <- data.table(data.frame(promoters(genes(txdb, columns="gene_id"), upstream=1000, downstream=0)))[, c("seqnames", "start", "end", "gene_id", "strand")][order(seqnames, start)]
gene_promoters[, gene_id := sapply(gene_promoters$gene_id, function(x) strsplit(x, "\\.")[[1]][1])]
gene_promoters[, start := ifelse(start < 0, 0, start)]
gene_promoters[, dummy := "."]
gene_promoters <- gene_promoters[, .(seqnames, start, end, gene_id, dummy, strand)]

write.table(gene_promoters, file = "/scratchb/sblab/martin03/repository/20191217_jochen/20200106/promoters/genepromoters.bed", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
```

Promoters in open chromatin:

```bash
cd promoters
bedtools intersect -a genepromoters.bed -b ../dnaseseq/DNAse-seq.concatenated_narrow_rep1_and_rep2.bed -wa -u > genepromoters_open.bed
wc -l genepromoters_open.bed # 16596
```


### Match encode ids, jaspar ids and gene names:

```python

# JASPAR
jaspar_file = open("jaspar/JASPAR2020_CORE_vertebrates_non-redundant.meme", "r")
jaspar_lines = jaspar_file.readlines()
jaspar_file.close()

jaspar_d = {}

for line in jaspar_lines:
  if line.startswith("MOTIF"):
    fields = line.split()
    jaspar_id = fields[1]
    gene_name = fields[2]
    if gene_name in jaspar_d:
      jaspar_d[gene_name].append(jaspar_id)
    else:
      jaspar_d[gene_name]= [jaspar_id]

len(jaspar_d)


# ENCODE
encode_file = open("encode/ENCODE_K562_Jan2019_BED_Meta.tsv", "r")
encode_lines = encode_file.readlines()
encode_file.close()

encode_gene_jaspar = []

for line in encode_lines[1:]:
  fields = line.split('\t')
  encode_id = fields[0]
  gene_name = fields[12].replace("-human", "")
  gene_name = gene_name.replace("eGFP-", "") if "eGFP-" in gene_name else gene_name
  gene_name = gene_name.replace("phosphoS5", "") if "phosphoS5" in gene_name else gene_name
  gene_name = gene_name.replace("phosphoS2", "") if "phosphoS2" in gene_name else gene_name
  gene_name = gene_name.replace("3xFLAG-", "") if "3xFLAG-" in gene_name else gene_name
  for gene_id in jaspar_d:
    if gene_id == gene_name:
      encode_gene_jaspar.append("\t".join([encode_id, gene_name, jaspar_d[gene_id][0]]))


len(encode_gene_jaspar) # 193

output_file = open("encode/encode_gene_jaspar.tsv", "w")
output_file.write("\n".join(encode_gene_jaspar) + "\n")
output_file.close()
```


### gat-run.py

```bash
srun --mem 32G --pty /usr/bin/bash

mkdir -p gat/{promoters,promoters_open}

# promoters
cat encode/encode_gene_jaspar.tsv | while read line
do
  encode_id=`echo $line | cut -d " " -f 1`
  gene_name=`echo $line | cut -d " " -f 2`
  jaspar_id=`echo $line | cut -d " " -f 3`
  #echo $encode_id, $gene_name, $jaspar_id
  nohup gat-run.py \
  -s <(zcat encode/${encode_id}.cut.bed.gz) \
  -a <(cat bg4/20180108_K562_async_rep1-3.mult.2of8.bed | sed 's/$/&\tbg4_2of8/') \
  -a <(cat bg4/20180108_K562_async_rep1-3.mult.5of8.bed | sed 's/$/&\tbg4_5of8/') \
  -a <(cat bg4/20180108_K562_async_rep1-3.mult.8of8.bed | sed 's/$/&\tbg4_8of8/') \
  -a <(cut -f1-3 jaspar/${jaspar_id}/fimo.bed | sed 's/$/&\tjaspar/') \
  -w <(cut -f1-3 promoters/genepromoters.bed) \
  --ignore-segment-tracks \
  -n 10000 \
  -L gat/promoters/${encode_id}_${gene_name}.log \
  --num-threads=1 > gat/promoters/${encode_id}_${gene_name}.txt &
done

tail gat/promoters/*.log
tail gat/promoters/*.txt | column -t


# promoters open
cat encode/encode_gene_jaspar.tsv | while read line
do
  encode_id=`echo $line | cut -d " " -f 1`
  gene_name=`echo $line | cut -d " " -f 2`
  jaspar_id=`echo $line | cut -d " " -f 3`
  #echo $encode_id, $gene_name, $jaspar_id
  nohup gat-run.py \
  -s <(zcat encode/${encode_id}.cut.bed.gz) \
  -a <(cat bg4/20180108_K562_async_rep1-3.mult.2of8.bed | sed 's/$/&\tbg4_2of8/') \
  -a <(cat bg4/20180108_K562_async_rep1-3.mult.5of8.bed | sed 's/$/&\tbg4_5of8/') \
  -a <(cat bg4/20180108_K562_async_rep1-3.mult.8of8.bed | sed 's/$/&\tbg4_8of8/') \
  -a <(cut -f1-3 jaspar/${jaspar_id}/fimo.bed | sed 's/$/&\tjaspar/') \
  -w <(cut -f1-3 promoters/genepromoters_open.bed) \
  --ignore-segment-tracks \
  -n 10000 \
  -L gat/promoters_open/${encode_id}_${gene_name}.log \
  --num-threads=1 > gat/promoters_open/${encode_id}_${gene_name}.txt &
done

tail gat/promoters_open/*.log
tail gat/promoters_open/*.txt | column -t
```


### Scatterplots

```r
#cd /scratchb/sblab/martin03/repository/20191217_jochen/20200106/gat
#R

library(data.table)
library(ggplot2)
library(ggpubr)
library(ggrepel)


# Enlarge the view width when printing tables
options(width = 400)

#############
# promoters #
#############

# Load data
data <- fread("tableCat.py -i promoters/*.txt -r .txt --firsthdr")

# Collapse data
data_dcast <- dcast(data, file_id ~ annotation, value.var = "fold")

data_dcast[, file_id := sapply(data_dcast$file_id, function(x) unlist(strsplit(x, "_"))[[2]])]

# Label G4 associated proteins
data_knownG4 <- data_dcast
Known_G4_proteins <- fread(file = "G4IPD_Oct2019.csv", header = F)
data_knownG4$G4IPD <- ifelse(data_knownG4$file_id %in% Known_G4_proteins$V1, "#279615", "black")

# Scatterplots
## bg4_5of8
gg <- ggplot(data_knownG4, aes(x = bg4_5of8, y = jaspar, label = file_id)) +
  geom_point(colour=data_knownG4$G4IPD) +
  theme_minimal() +
  #scale_color_manual(values=color_vector) +
  xlab("fold enrichment at endogeneous G4s") +
  ylab("fold enrichment at JASPAR motifs") +
  theme(axis.text = element_text(size=rel(1), color = "black"), axis.title = element_text(size=rel(1), color = "black"), legend.title = element_blank()) +
  scale_y_continuous(breaks=seq(0, 25, 5), limits=c(0,25)) +
  scale_x_continuous(breaks=seq(0, 25, 5)) +
  coord_fixed(ratio = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_label_repel(data = data_dcast[(bg4_5of8 > 17.5 & (bg4_5of8/jaspar > 2.4)) | ((bg4_5of8/jaspar < 0.5) & jaspar > 5)], 
                   aes(label = file_id), force = 2, box.padding = 0.4, point.padding = 0.4, size=rel(2.5), segment.color = 'grey80')

ggsave("figures/20200122_gat_bg4_5of8_jaspar_promoters.pdf",  height = 13/2.54, useDingbats=F)


#### Keep data from shuffling at promoters for export
export_table <- dcast(data, file_id ~ annotation, value.var = "fold")
export_table[, file.accession := sapply(export_table$file_id, function(x) unlist(strsplit(x, "_"))[[1]])]
export_table[, experiment.target := sapply(export_table$file_id, function(x) unlist(strsplit(x, "_"))[[2]])]
export_table <- export_table[,c("file.accession", "experiment.target", "bg4_5of8", "jaspar", "file_id")]


##################
# promoters_open #
##################

# Load data
data <- fread("tableCat.py -i promoters_open/*.txt -r .txt --firsthdr")

# Collapse data
data_dcast <- dcast(data, file_id ~ annotation, value.var = "fold")
data_dcast[, file_id := sapply(data_dcast$file_id, function(x) unlist(strsplit(x, "_"))[[2]])]

# Label G4 associated proteins
data_knownG4 <- data_dcast
Known_G4_proteins <- fread(file = "G4IPD_Oct2019.csv", header = F)
data_knownG4$G4IPD <- ifelse(data_knownG4$file_id %in% Known_G4_proteins$V1, "#279615", "black")

# Scatterplots
## bg4_5of8
gg <- ggplot(data_knownG4, aes(x = bg4_5of8, y = jaspar, label = file_id)) +
  geom_point(colour=data_knownG4$G4IPD) +
  theme_minimal() +
  xlab("fold enrichment at endogeneous G4s") +
  ylab("fold enrichment at JASPAR motifs") +
  theme(axis.text = element_text(size=rel(1), color = "black"), axis.title = element_text(size=rel(1), color = "black"), legend.title = element_blank()) +
  scale_y_continuous(breaks=seq(0, 25, 5), limits=c(0,25)) +
  scale_x_continuous(breaks=seq(0, 25, 5)) +
  coord_fixed(ratio = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_label_repel(data = data_dcast[(bg4_5of8 > 6.5 & (bg4_5of8/jaspar > 2)) | ((bg4_5of8/jaspar < 0.2) & jaspar > 5)], 
                   aes(label = file_id), force = 3, box.padding = 0.4, point.padding = 0.4, size=rel(2.5), segment.color = 'grey80')

ggsave("figures/20200122_gat_bg4_5of8_jaspar_OPEN_promoters.pdf", height = 13/2.54, useDingbats=F)


#######################
# export result table #
#######################

# data from shuffling at open promoters

data_dcast <- dcast(data, file_id ~ annotation, value.var = "fold")
export_table$bg4_5of8_open <- data_dcast$bg4_5of8[match(export_table$file_id, data_dcast$file_id)]
export_table$bg4_jaspar_open <- data_dcast$jaspar[match(export_table$file_id, data_dcast$file_id)]
export_table <- export_table[, -"file_id"]
  
write.csv(export_table, file ="GAT_ENCODE_G4vsJAPSAR.csv")
```





