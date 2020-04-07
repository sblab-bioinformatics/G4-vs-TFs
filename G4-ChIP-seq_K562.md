# G4-ChIP-seq analysis

For a more detailed description of our in-house G4-ChIP-seq pipeline please refer to our [previous work](https://github.com/sblab-bioinformatics/dna-secondary-struct-chrom-lands/blob/master/Methods.md) or [HepG2 analysis](https://github.com/sblab-bioinformatics/G4-vs-TFs/blob/master/G4-ChIP-seq_HepG2.md)
Briefly, G4-ChIP-seq experiments in K562 cells generated raw fastq files containing sequencing reads that were trimmed with *cutadapt* to remove adapter sequences and low-quality bases (phred quality score < 10). Trimmed reads were aligned to the human reference genome (version hg19) using *BWA* and duplicates were marked and removed using *Picard* and *samtools*, respectively. Peaks were called using *MACS2* (p < 1e-05) and merged considering different replicates using *bedtools multiinter*. Only peaks overlapping in 5 out 8 replicates were considered as high-confidence. 


### Visualize replicate overlap using Intervene

The folder `macs2/` contains the individual peak files generated as described above.

```bash
cd macs2/

intervene pairwise -i \
G4_bio1_a.narrowPeak \
G4_bio1_b.narrowPeak \
G4_bio2_a.narrowPeak \
G4_bio2_b.narrowPeak \
G4_bio2_c.narrowPeak \
G4_bio3_a.narrowPeak \
G4_bio3_b.narrowPeak \
G4_bio3_c.narrowPeak \
--names \
Rep1_a,Rep1_b,Rep2_a,Rep2_b,Rep2_c,Rep3_a,Rep3_b,Rep3_c \
--compute frac --htype tribar \
-o BG4_ChIP_replicates_PW_triang"
```


### G4-ChIP-seq signal distribution around TSS

Folder `bigWig/` contains G4 ChIP-seq bigWig files from G4 ChIP-seq in K562 and HepG2 normalized to sequencing depth generated using the *bamCoverage* function in *deeptools*. Signal profile plots around TSS were generated using functions *computeMatrix* and *plotProfile*.

```bash
sbatch -o G4.log -J DT_Profile --mem=10000 \
--wrap "computeMatrix reference-point \
--referencePoint center \
-b 1000 -a 1000 \
-S \
bigWig/K562_G4_bio1_a.bs50.bl.RPKM.bw \
bigWig/K562_G4_bio2_a.bs50.bl.RPKM.bw \
bigWig/K562_G4_bio3_a.bs50.bl.RPKM.bw \
bigWig/K562_G4_bio1_IP.bs50.bl.RPKM.bw \
bigWig/HepG2_G4ChIP_REP1_a.SLX-19356.i701_i502.bs50.bl.RPKM.bw \
bigWig/HepG2_G4ChIP_REP2_a.SLX-19356.i703_i502.bs50.bl.RPKM.bw \
bigWig/HepG2_G4ChIP_REP3_a.SLX-19356.i705_i502.bs50.bl.RPKM.bw \
bigWig/HepG2_G4ChIP_REP1_IP.SLX-19356.i701_i517.bs50.bl.RPKM.bw \
-R \
/scratchb/sblab/spiege01/ENCODE_K562/reference_data/hg19/ENCODE_transcripts.bed \
--skipZeros \
-o BG4_around_TSS_HepG2-K562.mat.gz && 
plotProfile -m BG4_around_TSS_HepG2-K562.mat.gz \
-out BG4_around_TSS_HepG2-K562.pdf \
--refPointLabel PeakCenter \
--regionsLabel ENCODE_TSS \
--samplesLabel K562_1a K562_2a KK562_3a K562_IP HepG2_1a HepG2_2a HepG2_3a HepG2_IP \
--dpi 300 \
--plotHeight 10 \
--plotWidth 20 \
--perGroup \
--numPlotsPerRow 1" 
```



# Generation of pseudo-stranded G4-ChIP-seq 

G4 ChIP-seq is not stranded. We used the stranded information from the G4-seq datasets (OQS) to infer strand information for G4-ChIP-seq high confidence peaks.


### Stranded OQS map

Download G4-seq bedfiles from GEO [GSE63874](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63874):

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63874/suppl/GSE63874%5FNa%5FK%5Fminus%5Fhits%5Fintersect%2Ebed%2Egz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63874/suppl/GSE63874%5FNa%5FPDS%5Fminus%5Fhits%5Fintersect%2Ebed%2Egz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63874/suppl/GSE63874%5FNa%5FPDS%5Fplus%5Fhits%5Fintersect%2Ebed%2Egz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63874/suppl/GSE63874%5FNa%5FK%5Fplus%5Fhits%5Fintersect%2Ebed%2Egz
```

Concatenate PDS and K OQs on + and - strands:

```bash
zcat GSE63874_Na_K_minus_hits_intersect.bed.gz GSE63874_Na_PDS_minus_hits_intersect.bed.gz | bedtools sort | bedtools merge | bedtools sort | gzip > G4-Seq_K_PDS_minus.bed.gz
zcat GSE63874_Na_K_plus_hits_intersect.bed.gz GSE63874_Na_PDS_plus_hits_intersect.bed.gz | bedtools sort | bedtools merge | bedtools sort | gzip > G4-Seq_K_PDS_plus.bed.gz

zcat G4-Seq_K_PDS_minus.bed.gz | awk '{print $1"\t"$2"\t"$3"\t.\t.\t-\t"}' - > G4-Seq_K_PDS_stranded.bed
zcat G4-Seq_K_PDS_plus.bed.gz | awk '{print $1"\t"$2"\t"$3"\t.\t.\t+\t"}' - >> G4-Seq_K_PDS_stranded.bed

cat G4-Seq_K_PDS_stranded.bed  | sortBed -i | gzip > G4-Seq_K_PDS_stranded.bed.gz
rm G4-Seq_K_PDS_stranded.bed
```

Non-stranded version:

```bash
bedtools merge -i G4-Seq_K_PDS_stranded.bed.gz > G4-Seq_cat_K_PDS_+-strands.bed   
```


### Stranded G4-ChIP-seq

Intersect stranded OQS with G4-ChIP-seq peaks and add strand information.

```bash
bedtools intersect -wa -a 20180108_K562_async_rep1-3.mult.5of8.bed -b <(zcat G4-Seq_K_PDS_plus.bed.gz) | bedtools merge | bedtools sort | awk -F '\t' '{print $1"\t"$2"\t"$3"\t"".""\t"".""\t""+""\t"}' | gzip > K562_async_rep1-3.mult.5of8_plusstrand.bed.gz

bedtools intersect -wa -a 20180108_K562_async_rep1-3.mult.5of8.bed -b <(zcat G4-Seq_K_PDS_minus.bed.gz) | bedtools merge | bedtools sort | awk -F '\t' '{print $1"\t"$2"\t"$3"\t"".""\t"".""\t""-""\t"}' | gzip > K562_async_rep1-3.mult.5of8_minusstrand.bed.gz
```

Remove ambiguous cases, and retain peaks exclusive to a single strand.

```bash
bedtools intersect -v -a <(zcat K562_async_rep1-3.mult.5of8_plusstrand.bed.gz) -b <(zcat K562_async_rep1-3.mult.5of8_minusstrand.bed.gz) | bedtools sort | bedtools merge | gzip > K562_async_rep1-3.mult.5of8_plusstrand_exclusive.bed.gz

zcat K562_async_rep1-3.mult.5of8_plusstrand_exclusive.bed.gz | wc -l # 2161

bedtools intersect -v -a <(zcat K562_async_rep1-3.mult.5of8_minusstrand.bed.gz) -b <(zcat K562_async_rep1-3.mult.5of8_plusstrand.bed.gz) | bedtools sort | bedtools merge | gzip > K562_async_rep1-3.mult.5of8_minusstrand_exclusive.bed.gz

zcat K562_async_rep1-3.mult.5of8_minusstrand_exclusive.bed.gz | wc -l # 2271
```

In summary, out of the 9205 peaks, 2161 peaks are exclusive to the plus strand and 2271 are exclusive to the minus strand, while 3494 BG4 peaks cannot be assigned unambiguously.

Add '+' for the positive strand and '-' for the negative strand, then combine into a single bedfile.

```bash
awk '{print $1"\t"$2"\t"$3"\t.\t.\t-\t"}' K562_async_rep1-3.mult.5of8_minusstrand_exclusive.bed > K562_async_rep1-3.mult.5of8_OQS-Stranded.bed
awk '{print $1"\t"$2"\t"$3"\t.\t.\t+\t"}' K562_async_rep1-3.mult.5of8_plusstrand_exclusive.bed >> K562_async_rep1-3.mult.5of8_OQS-Stranded.bed

cat K562_async_rep1-3.mult.5of8_OQS-Stranded.bed | sortBed -i | gzip > K562_async_rep1-3.mult.5of8_OQS_Stranded.bed.gz
```



# G4 control data set (OQS in open Chromatin around TSS)

Generate a control data set with sites that have the potential to form G4s (G4-seq) in open chromatin (DHS), reflecting the strong enrichemnt of G4s around TSS. Generate two versions with and without strand information. 

```bash
intersectBed -wa -a <(zcat G4-Seq_K_PDS_stranded.bed.gz) -b DNAse-seq.concatenated_narrow_rep1_and_rep2.bed > openOQs_stranded.bed
bedtools merge -i openOQs_stranded.bed > OQs_in_K562_open_chromatin.bed

intersectBed -v -a openOQs_stranded.bed -b 20180108_K562_async_rep1-3.mult.5of8.bed > openOQs_noBG4_stranded.bed

bedtools merge -i openOQs_noBG4_stranded.bed > openOQs_noBG4.bed
```

Use Pavis tool (https://manticore.niehs.nih.gov/pavis2/) with settings "hg19" "known gene". In the derived txt files, select upstream TSS and 5'UTR.

```bash
grep -e UTR -e Upstream openOQS_noBG4_1kbupTSS.txt | cut -f 2,3,4 | sortBed -i > openOQS_noBG4_1kbupTSS.bed
grep -e UTR -e Upstream openOQs_noBG4_stranded_1kbupTSS.txt | cut -f 2,3,4,5,6,7 | sortBed -i > temp.bed
awk '{print $1"\t"$2"\t"$3"\t.\t.\t"$6}' temp.bed > openOQs_noBG4_Stranded_1kbupTSS.bed
rm temp.bed
```



# Features of data sets

## Fragment size distribution in datasets

```R
BG4_consensus <- read.table(file = "BG4-ChIP/20180108_K562_async_rep1-3.mult.5of8.bed", sep = '\t', header = F)
BG4_consensus$size <- BG4_consensus$V3 - BG4_consensus$V2
mean(BG4_consensus$size)
median(BG4_consensus$size)
max(BG4_consensus$size)
min(BG4_consensus$size)

G4_opOQ_TSSupstream1000 <- read.table(file = "BG4-ChIP/openOQS_noBG4_1kbupTSS.bed", sep = '\t', header = F)
G4_opOQ_TSSupstream1000$size <- G4_opOQ_TSSupstream1000$V3 - G4_opOQ_TSSupstream1000$V2
mean(G4_opOQ_TSSupstream1000$size)
median(G4_opOQ_TSSupstream1000$size)
max(G4_opOQ_TSSupstream1000$size)
min(G4_opOQ_TSSupstream1000$size)
```
