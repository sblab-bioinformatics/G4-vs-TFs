# G4 vs TFs


## relevant datasets

Download RNA-seq gene quantifications from ENCODE:

```bash
cd rnaseq

# ENCSR000AEM
## ENCFF934YBO
wget https://www.encodeproject.org/files/ENCFF934YBO/@@download/ENCFF934YBO.tsv
## ENCFF515MUX
wget https://www.encodeproject.org/files/ENCFF515MUX/@@download/ENCFF515MUX.tsv

# ENCSR545DKY
## ENCFF186TXT
wget https://www.encodeproject.org/files/ENCFF186TXT/@@download/ENCFF186TXT.tsv
## ENCFF728TIT
wget https://www.encodeproject.org/files/ENCFF728TIT/@@download/ENCFF728TIT.tsv

wc -l ENCF*
   # 58541 ENCFF186TXT.tsv
   # 58541 ENCFF515MUX.tsv
   # 58541 ENCFF728TIT.tsv
   # 58541 ENCFF934YBO.tsv
```

GenomicFeatures hg19 promoter definition:

```r
# srun --mem 16G --pty /usr/bin/bash
# cd rnaseq/

library(data.table)
library(GenomicFeatures)

# change width
options(width = 250)

# prepare coordinates table
txdb <- makeTxDbFromGFF("reference_data/genomes/gencode/Gencode_human/release_19/gtf/gencode.v19.annotation.gtf", format="gtf")

# gene promoters
gene_promoters <- data.table(data.frame(promoters(genes(txdb, columns="gene_id"), upstream=1000, downstream=0)))[, c("seqnames", "start", "end", "gene_id", "strand")][order(seqnames, start)]
gene_promoters[, gene_id := sapply(gene_promoters$gene_id, function(x) strsplit(x, "\\.")[[1]][1])]
gene_promoters[, start := ifelse(start < 0, 0, start)]
gene_promoters[, dummy := "."]
gene_promoters <- gene_promoters[, .(seqnames, start, end, gene_id, dummy, strand)]

# ENCSR000AEM
## ENCFF934YBO
rnaseq_1_1 <- fread("rnaseq/ENCFF934YBO.tsv")
rnaseq_1_1 <- rnaseq_1_1[gene_id %like% 'ENSG'][, c("gene_id", "TPM", "FPKM")]
rnaseq_1_1[, gene_id := sapply(rnaseq_1_1$gene_id, function(x) strsplit(x, "\\.")[[1]][1])]
## ENCFF515MUX
rnaseq_1_2 <- fread("rnaseq/ENCFF515MUX.tsv")
rnaseq_1_2 <- rnaseq_1_2[gene_id %like% 'ENSG'][, c("gene_id", "TPM", "FPKM")]
rnaseq_1_2[, gene_id := sapply(rnaseq_1_2$gene_id, function(x) strsplit(x, "\\.")[[1]][1])]

# ENCSR545DKY
## ENCFF186TXT
rnaseq_2_1 <- fread("rnaseq/ENCFF186TXT.tsv")
rnaseq_2_1 <- rnaseq_2_1[gene_id %like% 'ENSG'][, c("gene_id", "TPM", "FPKM")]
rnaseq_2_1[, gene_id := sapply(rnaseq_2_1$gene_id, function(x) strsplit(x, "\\.")[[1]][1])]
## ENCFF728TIT
rnaseq_2_2 <- fread("rnaseq/ENCFF728TIT.tsv")
rnaseq_2_2 <- rnaseq_2_2[gene_id %like% 'ENSG'][, c("gene_id", "TPM", "FPKM")]
rnaseq_2_2[, gene_id := sapply(rnaseq_2_2$gene_id, function(x) strsplit(x, "\\.")[[1]][1])]

# Merge
setkey(gene_promoters, "gene_id")
setkey(rnaseq_1_1, "gene_id")
setkey(rnaseq_1_2, "gene_id")
setkey(rnaseq_2_1, "gene_id")
setkey(rnaseq_2_2, "gene_id")

gene_promoters <- gene_promoters[rnaseq_1_1,][rnaseq_1_2,][rnaseq_2_1,][rnaseq_2_2,][order(seqnames, start)]

write.table(gene_promoters, file = "rnaseq/genepromoters_rnaseq.bed", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
```



## Defining BG4 in open chromatin

```bash
cd bg4

bedtools intersect -a 20180108_K562_async_rep1-3.mult.2of8.bed -b ../dnaseseq/DNAse-seq.concatenated_narrow_rep1_and_rep2.bed -wa -u > bg4_2of8_open.bed
wc -l 20180108_K562_async_rep1-3.mult.2of8.bed # 19249
wc -l bg4_2of8_open.bed # 16113

bedtools intersect -a 20180108_K562_async_rep1-3.mult.5of8.bed -b ../dnaseseq/DNAse-seq.concatenated_narrow_rep1_and_rep2.bed -wa -u > bg4_5of8_open.bed
wc -l 20180108_K562_async_rep1-3.mult.5of8.bed # 9205
wc -l bg4_5of8_open.bed # 8768

```




## Intersecting BG4, TFs and gene promoters 5of8

- overlap: BG4 5of8 + TF
- tf_only: no BG4 2of8 + TF
- bg4_only: BG4 5of8 + no TF

```bash
cd /scratchb/sblab/martin03/repository/20191217_jochen/20191217

mkdir -p rnaseq/bg4_5of8
mkdir -p figures/bg4_5of8
mkdir -p figures/bg4_5of8_rep1

for bedgz in encode/*.cut.bed.gz
do
  echo $bedgz
  bedtools intersect \
  -a rnaseq/genepromoters_rnaseq.bed \
  -b <(cat <(bedtools intersect -a bg4/20180108_K562_async_rep1-3.mult.5of8.bed -b <(zcat $bedgz) -wa -u | sed 's/$/&\toverlap/') \
  <(bedtools intersect -a bg4/20180108_K562_async_rep1-3.mult.5of8.bed -b <(zcat $bedgz) -v | sed 's/$/&\tbg4_only/') \
  <(bedtools intersect -a <(zcat $bedgz) -b bg4/20180108_K562_async_rep1-3.mult.2of8.bed -v | sed 's/$/&\ttf_only/')) -loj | \
  awk '$16 != -1' | \
  bedtools sort -i - | \
  pigz > rnaseq/bg4_5of8/`basename $bedgz`
done

ls rnaseq/bg4_5of8/*.bed.gz | wc -l # 528
```



## Boxplots BG4 5of8 - compare all replicates

```r
library(data.table)
library(ggplot2)
library(ggpubr)

# Enlarge the view width when printing tables
options(width = 350)

# Load metadata
metadata <- fread("encode/ENCODE_K562_Jan2019_BED_Meta.tsv")

# Iterate through data
for (f in list.files("rnaseq/bg4_5of8/")){
  if (grepl("cut.bed.gz", f)){
    print(f)
    # Load data
    encode_id <- strsplit(f, "\\.")[[1]][1]
    gene_name <- strsplit(metadata[File.accession == encode_id]$Experiment.target, "-")[[1]][1]
    data <- fread(paste("zcat rnaseq/bg4_5of8/", f, sep = ""))
    setnames(data, c("chr_promoter", "start_promoter", "end_promoter", "ensembl_id", "dummy", "strand", "TPM_rep11", "FPKM_rep11", "TPM_rep12", "FPKM_rep12", "TPM_rep21", "FPKM_rep21", "TPM_rep22", "FPKM_rep22", "chr_mark", "start_mark", "end_mark", "type_mark"))
    # Melt
    data_melt <- melt(data, id.vars = c("type_mark"), measure.vars = c("TPM_rep11", "TPM_rep12", "TPM_rep21", "TPM_rep22"), variable.name = "dataset", value.name = "tpm")
    # Boxplot
    labs <- c("rep 1.1", "rep 1.2", "rep 2.1", "rep 2.2")
    names(labs) <- c("TPM_rep11", "TPM_rep12", "TPM_rep21", "TPM_rep22")
    data_melt <- data_melt[, dataset := factor(dataset, levels = names(labs))]
    gg <- ggplot(data = data_melt, aes(x = factor(data_melt$type_mark, levels = c("overlap", "tf_only", "bg4_only")), y = log10(tpm))) +
    geom_boxplot(outlier.shape=NA) +
    ylab(expression("log"[10]*"TPM")) +
    xlab(expression("")) +
    ggtitle(gene_name) +
    theme_bw() +
    scale_x_discrete(labels=c("overlap" = "TF\nBG4 5of8", "tf_only" = "TF\nno BG4 2of8", "bg4_only" = "no TF\nBG4 5of8")) +
    stat_compare_means(comparisons = list( c("overlap", "tf_only"), c("overlap", "bg4_only")), method = "wilcox.test") +
    coord_cartesian(ylim = c(-2, 6)) +
    facet_wrap(~ dataset, labeller = labeller(dataset = labs)) +
    theme(axis.title = element_text(size=16), axis.text = element_text(size=12, color = "black"), strip.text = element_text(size=16, color = "black"), plot.title = element_text(face="bold", size=16, hjust = 0.5))
    ggsave(paste("figures/bg4_5of8/", encode_id, "_", gene_name, ".pdf", sep = ""))
  }
}
```



## Boxplots BG4 5of8 rep1

```r
library(data.table)
library(ggplot2)
library(ggpubr)

# Load metadata
metadata <- fread("encode/ENCODE_K562_Jan2019_BED_Meta.tsv")

# Iterate through data
for (f in list.files("rnaseq/bg4_5of8/")){
  if (grepl("cut.bed.gz", f)){
    print(f)
    # Load data
    encode_id <- strsplit(f, "\\.")[[1]][1]
    gene_name <- strsplit(metadata[File.accession == encode_id]$Experiment.target, "-")[[1]][1]
    data <- fread(paste("zcat rnaseq/bg4_5of8/", f, sep = ""))
    setnames(data, c("chr_promoter", "start_promoter", "end_promoter", "ensembl_id", "dummy", "strand", "TPM_rep11", "FPKM_rep11", "TPM_rep12", "FPKM_rep12", "TPM_rep21", "FPKM_rep21", "TPM_rep22", "FPKM_rep22", "chr_mark", "start_mark", "end_mark", "type_mark"))

    # remove BG4_only category
    data <- data[data$type_mark != "bg4_only",]
        
    # Boxplot
    gg <- ggplot(data = data, aes(x = factor(data$type_mark, levels = c("overlap", "tf_only")), y = log10(TPM_rep11))) +
          geom_boxplot(outlier.shape=NA) +
          ylab(expression("log"[10]*"TPM")) +
          xlab(expression("")) +
          ggtitle(gene_name) +
          theme_classic() +
          scale_x_discrete(labels=c("overlap" = sprintf("TF\nG4\n(%s)", data.table(table(data$type_mark))[V1 == "overlap"]$N), "tf_only" = sprintf("TF\nno G4\n(%s)", data.table(table(data$type_mark))[V1 == "tf_only"]$N))) +
          stat_compare_means(comparisons = list( c("overlap", "tf_only")), method = "wilcox.test", label.y = 4) +
          coord_cartesian(ylim = c(-2, 4.2)) +
          theme(axis.title = element_text(size=16), axis.text = element_text(size=12, color = "black"), strip.text = element_text(size=16, color = "black"), plot.title = element_text(face="bold", size=16, hjust = 0.5))
    ggsave(paste("figures/bg4_5of8_rep1/", encode_id, "_", gene_name, ".pdf", sep = ""), width = 5, height = 7, units = 'cm')
  }
}
```






## Intersecting BG4 at all open gene promoters (compare BG4 5of8 and no 2of8)

```bash

bedtools intersect \
-a <(bedtools intersect -a rnaseq/genepromoters_rnaseq.bed -b bg4/20180108_K562_async_rep1-3.mult.5of8.bed -wa -u | sed 's/$/&\tbg4_5of8/') \
-b dnaseseq/DNAse-seq.concatenated_narrow_rep1_and_rep2.bed \
-wa -u | wc -l # 5977 promoters overlapping with bg4 5of8 in open chromatin

bedtools intersect \
-a <(bedtools intersect -a rnaseq/genepromoters_rnaseq.bed -b bg4/20180108_K562_async_rep1-3.mult.2of8.bed -v | sed 's/$/&\tno_bg4_2of8/') \
-b dnaseseq/DNAse-seq.concatenated_narrow_rep1_and_rep2.bed \
-wa -u | wc -l # 7904 promoters that do not overlap with bg4 2of8 in open chromatin

cat \
<(bedtools intersect \
-a <(bedtools intersect -a rnaseq/genepromoters_rnaseq.bed -b bg4/20180108_K562_async_rep1-3.mult.5of8.bed -wa -u | sed 's/$/&\tbg4_5of8/') \
-b dnaseseq/DNAse-seq.concatenated_narrow_rep1_and_rep2.bed \
-wa -u) \
<(bedtools intersect \
-a <(bedtools intersect -a rnaseq/genepromoters_rnaseq.bed -b bg4/20180108_K562_async_rep1-3.mult.2of8.bed -v | sed 's/$/&\tno_bg4_2of8/') \
-b dnaseseq/DNAse-seq.concatenated_narrow_rep1_and_rep2.bed \
-wa -u) | \
bedtools sort -i > rnaseq/genepromoters_rnaseq_5of8_open.bed
```



## Boxplot 5of8 rep1 all open

```r
library(data.table)
library(ggplot2)
library(ggpubr)


# Load data
data <- fread("rnaseq/genepromoters_rnaseq_5of8_open.bed")
setnames(data, c("chr_promoter", "start_promoter", "end_promoter", "ensembl_id", "dummy", "strand", "TPM_rep11", "FPKM_rep11", "TPM_rep12", "FPKM_rep12", "TPM_rep21", "FPKM_rep21", "TPM_rep22", "FPKM_rep22", "status"))

sum(data$TPM_rep11 == 0) # 4505
sum(data$TPM_rep11 == 0 & data$status == "bg4_5of8") # 682
sum(data$TPM_rep11 == 0 & data$status == "no_bg4_2of8") # 3823

sum(data$TPM_rep11 != 0 & data$status == "bg4_5of8") # 5295
sum(data$TPM_rep11 != 0 & data$status == "no_bg4_2of8") # 4081

# Boxplot
gg <- ggplot(data = data, aes(x = status, y = log10(TPM_rep11))) +
  geom_boxplot(outlier.shape=NA) +
  ylab(expression("log"[10]*"TPM")) +
  xlab(expression("")) +
  theme_classic() +
  scale_x_discrete(labels=c("bg4_5of8" = sprintf("G4\n(%s)", sum(data$TPM_rep11 != 0 & data$status == "bg4_5of8")), "no_bg4_2of8" = sprintf("no G4\n(%s)", sum(data$TPM_rep11 != 0 & data$status == "no_bg4_2of8")))) +
  stat_compare_means(comparisons = list( c("bg4_5of8", "no_bg4_2of8")), method = "wilcox.test", label.y = 4) +
  coord_cartesian(ylim = c(-2, 4.2)) +
  theme(axis.title = element_text(size=16), axis.text = element_text(size=16, color = "black"), strip.text = element_text(size=16, color = "black"), plot.title = element_text(face="bold", size=16, hjust = 0.5))
ggsave("figures/genepromoters_rnaseq_5of8_open.pdf",  width = 5, height = 7, units = 'cm')
```








# Hyper-recruitment plots 5of8 rep1 promoters with and without bg4 in open chromatin

1. Defined two groups: promoters with BG4 5of8 and no BG4 2of8 overlap in open chromatin
2. Annotated the groups of promoters with the encode bed files (excluding histone marks)
3. Plotted the count of encode files overlapping the groups of promoters vs TPM

```bash
cd /scratchb/sblab/martin03/repository/20191217_jochen/20191217

bedtools intersect \
-a <(bedtools intersect -a rnaseq/genepromoters_rnaseq.bed -b bg4/20180108_K562_async_rep1-3.mult.5of8.bed -wa -u | sed 's/$/&\tbg4_5of8/') \
-b dnaseseq/DNAse-seq.concatenated_narrow_rep1_and_rep2.bed \
-wa -u | wc -l # 5977 promoters overlapping with bg4 5of8 in open chromatin

bedtools intersect \
-a <(bedtools intersect -a rnaseq/genepromoters_rnaseq.bed -b bg4/20180108_K562_async_rep1-3.mult.2of8.bed -v | sed 's/$/&\tno_bg4_2of8/') \
-b dnaseseq/DNAse-seq.concatenated_narrow_rep1_and_rep2.bed \
-wa -u | wc -l # 7904 promoters that do not overlap with bg4 2of8 in open chromatin

bedtools annotate -counts \
-i <(cat \
<(bedtools intersect \
-a <(bedtools intersect -a rnaseq/genepromoters_rnaseq.bed -b bg4/20180108_K562_async_rep1-3.mult.5of8.bed -wa -u | sed 's/$/&\tbg4_5of8/') \
-b dnaseseq/DNAse-seq.concatenated_narrow_rep1_and_rep2.bed \
-wa -u) \
<(bedtools intersect \
-a <(bedtools intersect -a rnaseq/genepromoters_rnaseq.bed -b bg4/20180108_K562_async_rep1-3.mult.2of8.bed -v | sed 's/$/&\tno_bg4_2of8/') \
-b dnaseseq/DNAse-seq.concatenated_narrow_rep1_and_rep2.bed \
-wa -u) | \
bedtools sort -i) \
-files <(zcat `ls encode/*.cut.bed.gz | grep -v -E "ENCFF737AMS|ENCFF730VTO|ENCFF126QYP|ENCFF883POK|ENCFF085JHT|ENCFF378OQB|ENCFF676RWX|ENCFF894VEM|ENCFF769WZF|ENCFF624XRN|ENCFF139CKE|ENCFF173ULG|ENCFF127XXD|ENCFF099LMD|ENCFF183UQD|ENCFF498CMP|ENCFF322IFF|ENCFF044JNJ"`) > rnaseq/genepromoters_rnaseq_5of8_recruitment_open.bed
```

Plotting

```r
library(data.table)
library(ggplot2)
library(ggpubr)

# Load data
data <- fread("rnaseq/genepromoters_rnaseq_5of8_recruitment_open.bed")
setnames(data, c("chr_promoter", "start_promoter", "end_promoter", "ensembl_id", "dummy", "strand", "TPM_rep11", "FPKM_rep11", "TPM_rep12", "FPKM_rep12", "TPM_rep21", "FPKM_rep21", "TPM_rep22", "FPKM_rep22", "status", "count"))

# Categorise encode datasets count
data[, count_cut := cut(count, breaks=c(-1, 50, 100, 150, 200, 250, Inf), labels = c("<50", "50-100", "100-150", "150-200", "200-250", ">250"))]
sum(data[status == "bg4_5of8"]$TPM_rep11 == 0) # 682
sum(data[status == "no_bg4_2of8"]$TPM_rep11 == 0) # 3823

# Boxplot bg4_5of8 and no_bg4_2of8
gg <- ggplot(data = data, aes(x = count_cut, y = log10(TPM_rep11), fill = status)) +
  geom_boxplot(outlier.shape=NA, position=position_dodge(0.9), show.legend = FALSE) +
  ylab(expression("log"[10]*"TPM")) +
  xlab(expression("")) +
  theme_classic() +
  scale_fill_manual(values=c("forestgreen", "dimgray"), labels=c("G4", "non-G4"), name="Open promoters") +
  scale_x_discrete(labels=c("<50" = sprintf("<50\n(%s, %s)", sum(data[status == "bg4_5of8"]$TPM_rep11 != 0 & data[status == "bg4_5of8"]$count_cut == "<50"), sum(data[status == "no_bg4_2of8"]$TPM_rep11 != 0 & data[status == "no_bg4_2of8"]$count_cut == "<50")), "50-100" = sprintf("50-100\n(%s, %s)", sum(data[status == "bg4_5of8"]$TPM_rep11 != 0 & data[status == "bg4_5of8"]$count_cut == "50-100"), sum(data[status == "no_bg4_2of8"]$TPM_rep11 != 0 & data[status == "no_bg4_2of8"]$count_cut == "50-100")), "100-150" = sprintf("100-150\n(%s, %s)", sum(data[status == "bg4_5of8"]$TPM_rep11 != 0 & data[status == "bg4_5of8"]$count_cut == "100-150"), sum(data[status == "no_bg4_2of8"]$TPM_rep11 != 0 & data[status == "no_bg4_2of8"]$count_cut == "100-150")), "150-200" = sprintf("150-200\n(%s, %s)", sum(data[status == "bg4_5of8"]$TPM_rep11 != 0 & data[status == "bg4_5of8"]$count_cut == "150-200"), sum(data[status == "no_bg4_2of8"]$TPM_rep11 != 0 & data[status == "no_bg4_2of8"]$count_cut == "150-200")), "200-250" = sprintf("200-250\n(%s, %s)", sum(data[status == "bg4_5of8"]$TPM_rep11 != 0 & data[status == "bg4_5of8"]$count_cut == "200-250"), sum(data[status == "no_bg4_2of8"]$TPM_rep11 != 0 & data[status == "no_bg4_2of8"]$count_cut == "200-250")), ">250" = sprintf(">250\n(%s, %s)", sum(data[status == "bg4_5of8"]$TPM_rep11 != 0 & data[status == "bg4_5of8"]$count_cut == ">250"), sum(data[status == "no_bg4_2of8"]$TPM_rep11 != 0 & data[status == "no_bg4_2of8"]$count_cut == ">250")))) +
  scale_y_continuous(limits=c(-2,4), breaks=seq(-2,4,1))+
  #coord_cartesian(ylim = c(-2, 4)) +
  theme(axis.title = element_text(size=16), axis.text.x = element_text(size=12, color = "black"), axis.text.y = element_text(size=16, color = "black"), legend.title = element_text(size=12, color = "black", face = "bold"), legend.text = element_text(size=12, color = "black"), legend.justification=c(1,0), legend.position=c(0.975,0.025), legend.background = element_rect(size=0.5, linetype="solid", colour ="black"))
ggsave("figures/genepromoters_rnaseq_5of8_recruitment_open_bg4_5of8_no_bg4_2of8.pdf",  width = 16, height = 12, units = 'cm')



Caculate p-values for individual categories:

## <50
wilcox.test(log10(data[status == "bg4_5of8" & TPM_rep11 != 0 & count_cut == "<50"]$TPM_rep11), log10(data[status == "no_bg4_2of8" & TPM_rep11 != 0 & count_cut == "<50"]$TPM_rep11))$p.value # 0.004586793
## 50-100
wilcox.test(log10(data[status == "bg4_5of8" & TPM_rep11 != 0 & count_cut == "50-100"]$TPM_rep11), log10(data[status == "no_bg4_2of8" & TPM_rep11 != 0 & count_cut == "50-100"]$TPM_rep11))$p.value # 0.003005993
## 100-150
wilcox.test(log10(data[status == "bg4_5of8" & TPM_rep11 != 0 & count_cut == "100-150"]$TPM_rep11), log10(data[status == "no_bg4_2of8" & TPM_rep11 != 0 & count_cut == "100-150"]$TPM_rep11))$p.value # 0.003210076
## 150-200
wilcox.test(log10(data[status == "bg4_5of8" & TPM_rep11 != 0 & count_cut == "150-200"]$TPM_rep11), log10(data[status == "no_bg4_2of8" & TPM_rep11 != 0 & count_cut == "150-200"]$TPM_rep11))$p.value # 0.08692111
## 200-250
wilcox.test(log10(data[status == "bg4_5of8" & TPM_rep11 != 0 & count_cut == "200-250"]$TPM_rep11), log10(data[status == "no_bg4_2of8" & TPM_rep11 != 0 & count_cut == "200-250"]$TPM_rep11))$p.value # 0.1807525
## >250
wilcox.test(log10(data[status == "bg4_5of8" & TPM_rep11 != 0 & count_cut == ">250"]$TPM_rep11), log10(data[status == "no_bg4_2of8" & TPM_rep11 != 0 & count_cut == ">250"]$TPM_rep11))$p.value # 0.7266649
```



