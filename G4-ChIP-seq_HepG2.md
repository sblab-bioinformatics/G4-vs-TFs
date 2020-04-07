For a more detailed description of our in-house G4-ChIP-seq pipeline please refer to our [previous work](https://github.com/sblab-bioinformatics/dna-secondary-struct-chrom-lands/blob/master/Methods.md). 

## G4-ChIP-seq analysis
The folder `fastq/` contains the relevant fastq files.

```
mkdir bam
mkdir clean_bam
mkdir macs2_output
mkdir tdf
mkdir bigWig
mkdir logs
mkdir analysis
mkdir analysis/confidence_peaks
mkdir fastq_trimmed
```


### adapter trimming for single-end sequencing

```
cd fastq


for f in *.fq.gz; do sbatch --mem 8G -o $f.out -e $f.err -J cutadapt --wrap "~/.local/bin/cutadapt -m 10 -q 20 -O 3 -a CTGTCTCTTATACACATCT -o ../fastq_trimmed/${f%%.H757JBGXF.s_1.r_1.fq.gz}.trimmed.fq.gz $f"; done
```



### alignment single-end to hg19
Align to hg19 reference genome.

Align:
```
cd ../fastq_trimmed

genome='~/reference_data/hg19/Sequence/Bowtie2Index/bowtie2-build/genome.fa'; for f in *.trimmed.fq.gz; do sbatch --mem 16G -o $f.bwa.out -e $f.bwa.err -J bwa --wrap "~/bwa mem -t 8 -M $genome $f | ~/bin/samtools view -Sb -F2308 -q 10 > ../bam/${f%%.trimmed.fq.gz}.hg19.bam"; done
```


### Sort bam files, mark and remove duplicates
```
for f in *.bam; do sbatch --mem 16G -o $f.out -e $f.err -J sort --wrap "~/bin/samtools sort -@ 8 $f > ${f%%.hg19.bam}.sort.bam"; done
```


Mark and remove duplicates
```
for f in *sort.bam; do sbatch --mem 8G -o $f.out -e $f.err -J dups --wrap "java -Xmx7g -jar B/picard/picard-2.14.0/picard.jar MarkDuplicates I=$f O=../clean_bam/${f%%.sort.bam}.nodup.bam M=${f%%.sort.bam}.md.txt AS=true REMOVE_DUPLICATES=true"; done
```

Some basic dup and size stats
```
grep 'Unknown Library' *txt > dups_info.txt

srun --mem 32G --pty /usr/bin/bash
cd ../clean_bam/
for f in *nodup.bam; do echo $f;echo $f >> all_counts && ~/bin/samtools view -c $f >> all_counts; done
exit
```

**RESULTS**
Duplication rate
```
HepG2_G4ChIP_REP1_a.SLX-19356.i701_i502.md.txt:Unknown Library	38985038	0	0	0	2743345	0	0	0.070369
HepG2_G4ChIP_REP1_b.SLX-19356.i702_i503.md.txt:Unknown Library	38816297	0	0	0	2745336	0	0	0.070726
HepG2_G4ChIP_REP1_c.SLX-19356.i702_i504.md.txt:Unknown Library	29160060	0	0	0	1796675	0	0	0.061614
HepG2_G4ChIP_REP1_IP.SLX-19356.i701_i517.md.txt:Unknown Library	40206231	0	0	0	2191655	0	0	0.05451
HepG2_G4ChIP_REP2_a.SLX-19356.i703_i502.md.txt:Unknown Library	33848936	0	0	0	2223053	0	0	0.065676
HepG2_G4ChIP_REP2_b.SLX-19356.i704_i503.md.txt:Unknown Library	31486074	0	0	0	2063245	0	0	0.065529
HepG2_G4ChIP_REP2_c.SLX-19356.i704_i504.md.txt:Unknown Library	30949057	0	0	0	1863480	0	0	0.060211
HepG2_G4ChIP_REP2_IP.SLX-19356.i703_i517.md.txt:Unknown Library	45271705	0	0	0	2628556	0	0	0.058062
HepG2_G4ChIP_REP3_a.SLX-19356.i705_i502.md.txt:Unknown Library	36811105	0	0	0	2301494	0	0	0.062522
HepG2_G4ChIP_REP3_b.SLX-19356.i706_i503.md.txt:Unknown Library	41854781	0	0	0	2789298	0	0	0.066642
HepG2_G4ChIP_REP3_c.SLX-19356.i706_i504.md.txt:Unknown Library	30203512	0	0	0	1729902	0	0	0.057275
HepG2_G4ChIP_REP3_IP.SLX-19356.i705_i517.md.txt:Unknown Library	48351618	0	0	0	2864954	0	0	0.059252
```
Duplication rate of 5-7%. OK.


Read count per bam file:
```
HepG2_G4ChIP_REP1_a.SLX-19356.i701_i502.nodup.bam 36241693
HepG2_G4ChIP_REP1_b.SLX-19356.i702_i503.nodup.bam 36070961
HepG2_G4ChIP_REP1_c.SLX-19356.i702_i504.nodup.bam 27363385
HepG2_G4ChIP_REP1_IP.SLX-19356.i701_i517.nodup.bam 38014576
HepG2_G4ChIP_REP2_a.SLX-19356.i703_i502.nodup.bam 31625883
HepG2_G4ChIP_REP2_b.SLX-19356.i704_i503.nodup.bam 29422829
HepG2_G4ChIP_REP2_c.SLX-19356.i704_i504.nodup.bam 29085577
HepG2_G4ChIP_REP2_IP.SLX-19356.i703_i517.nodup.bam 42643149
HepG2_G4ChIP_REP3_a.SLX-19356.i705_i502.nodup.bam 34509611
HepG2_G4ChIP_REP3_b.SLX-19356.i706_i503.nodup.bam 39065483
HepG2_G4ChIP_REP3_c.SLX-19356.i706_i504.nodup.bam 28473610
HepG2_G4ChIP_REP3_IP.SLX-19356.i705_i517.nodup.bam 45486664
```
~30 Mio reads in almost all of the cases. 


### create index bam files and generate tdf files for visual inspection
Generate bai files using samtools index
```
for f in *nodup.bam; do sbatch -o %j.$f.log --mem 16G -J index.$f --wrap "~/bin/samtools index -b $f"; done

```

Generate tdfs
```
for f in *nodup.bam; do sbatch -o %j.$f.log --mem 8G -J igv.$f --wrap "~/bin/igvtools/igvtools-2.3.91/igvtools count -w 15 $f ../tdf/${f%%.nodup.bam}.tdf ~/bin/igvtools/igvtools-2.3.91/genomes/hg19.chrom.sizes"; done
```

### create bigWig files
```
for FILE in *nodup.bam
do
bname=`basename $FILE .bam`
sbatch -o %j.$f.log -J BCoverage --mem=8G \
--wrap "bamCoverage -b $FILE \
-o ../bigWig/${FILE%%.nodup.bam}.bs50.bl.RPKM.bw \
--binSize 50 \
--blackListFileName ~/bin//reference_data/hg19/hg19.blacklist_merge1000nt.bed \
--numberOfProcessors max \
--normalizeUsing RPKM"
done
```



### call peaks
(using macs2 2.1.1.20160309)
Call peaks vs input file.

```
for f in *REP1*nodup.bam; do sbatch --mem 8G -o $f.macs2.out -e $f.macs2.err -J macs2 --wrap "macs2 callpeak --name ../macs2_output/${f%%.nodup.bam}.nodup.q005.all -t $f -c HepG2_G4ChIP_REP1_IP.SLX-19356.i701_i517.nodup.bam --format=BAM --gsize 'hs' --bw=300 --qvalue 0.05"; done
for f in *REP2*nodup.bam; do sbatch --mem 8G -o $f.macs2.out -e $f.macs2.err -J macs2 --wrap "macs2 callpeak --name ../macs2_output/${f%%.nodup.bam}.nodup.q005.all -t $f -c HepG2_G4ChIP_REP2_IP.SLX-19356.i703_i517.nodup.bam --format=BAM --gsize 'hs' --bw=300 --qvalue 0.05"; done
for f in *REP3*nodup.bam; do sbatch --mem 8G -o $f.macs2.out -e $f.macs2.err -J macs2 --wrap "macs2 callpeak --name ../macs2_output/${f%%.nodup.bam}.nodup.q005.all -t $f -c HepG2_G4ChIP_REP3_IP.SLX-19356.i705_i517.nodup.bam --format=BAM --gsize 'hs' --bw=300 --qvalue 0.05"; done
```


### Generate high confidence peaks
```
multiIntersectBed -i *.narrowPeak > ../analysis/confidence_peaks/HepG2_async_G4ChIP_rep1-3_multintersect.bed

cd ../analysis/confidence_peaks

for ((X=1; X<=9; X++))
do
FILENAME="HepG2_async_rep1-3.mult.$X""of9.bed"
awk -v VARIABLE="$X" '$4 >= VARIABLE'< HepG2_async_G4ChIP_rep1-3_multintersect.bed | sortBed -i - | mergeBed -i - > $FILENAME
echo $FILENAME
done



wc -l *.bed

  66501 HepG2_async_rep1-3.mult.1of9.bed
   28387 HepG2_async_rep1-3.mult.2of9.bed
   18471 HepG2_async_rep1-3.mult.3of9.bed
   13791 HepG2_async_rep1-3.mult.4of9.bed
   10922 HepG2_async_rep1-3.mult.5of9.bed
    8805 HepG2_async_rep1-3.mult.6of9.bed
    7019 HepG2_async_rep1-3.mult.7of9.bed
    5285 HepG2_async_rep1-3.mult.8of9.bed
    3606 HepG2_async_rep1-3.mult.9of9.bed

```




# Propteries of G4 ChIP peaks

### Size distribution

```R
BG4_consensus <- read.table(file = "HepG2_async_rep1-3.mult.6of9.bed", sep = '\t', header = F)
BG4_consensus$size <- BG4_consensus$V3 - BG4_consensus$V2
mean(BG4_consensus$size)
median(BG4_consensus$size)
max(BG4_consensus$size)
min(BG4_consensus$size)


mean(BG4_consensus$size)
[1] 189.4319
> median(BG4_consensus$size)
[1] 169
> max(BG4_consensus$size)
[1] 2100
> min(BG4_consensus$size)
[1] 1
```



### Overlap with potential G4s 
(observed quadruplex sequences from G4 ChIP)

```
intersectBed -a HepG2_async_rep1-3.mult.6of9.bed -b ~/reference_data/Quadruplex/G4-Seq_cat_K_PDS_+-strands.bed | sortBed -i | mergeBed | wc -l
6894 (78%)


# Venn diagramms
sbatch -J Venn -o pw.log --mem 2g --wrap " intervene venn -i \
HepG2_async_rep1-3.mult.6of9.bed  \
~/reference_data/Quadruplex/G4-Seq_cat_K_PDS_+-strands.bed \
--names \
G4,OQS \
-o BG4_vs_OQS"

```

### Overlap with open chromatin
```
intersectBed -a HepG2_async_rep1-3.mult.6of9.bed -b ../../../../HepG2_reference/ENCFF571RHF_DHS.bed | sortBed -i | mergeBed | wc -l
8430 (96%)


# Venn diagramms
sbatch -J Venn -o pw.log --mem 2g --wrap " intervene venn -i \
HepG2_async_rep1-3.mult.6of9.bed  \
../../../../HepG2_reference/ENCFF571RHF_DHS.bed  \
--names \
G4,DHS \
-o BG4_vs_DHS"
```




### Overlap with K562
Check overlap with K562 cells:

```
intersectBed -a HepG2_async_rep1-3.mult.6of9.bed -b 20180108_K562_async_rep1-3.mult.5of8.bed  | wc -l
4824 (55% oof HepG2 peaks; 52% of 5of8_K562)
```



### Clustering of replicates in K562 and HepG2
Assess similarity of datasets comparing G4 ChIP-seq of HepG2 and K562 in union of G4 ChiP peaks.

```
cd analysis
mkdir deeptools && cd deeptools


cat ../../macs2_output/HepG2_G4ChIP_REP*narrowPeak ~/ENCODE_K562/BG4_features/SLX* | sortBed -i - | mergeBed -d 500 -i  - > G4ChIP_K562_HepG2_union.merged_500.bed


cd ../../bigWig
sbatch -o ../analysis/deeptools/multiBigWigSum.%j.log -J MultBWSum --mem=20000 \
--wrap "multiBigwigSummary BED-file \
-b \
HepG2_G4ChIP_REP1_a.SLX-19356.i701_i502.bs50.bl.RPKM.bw \
HepG2_G4ChIP_REP1_b.SLX-19356.i702_i503.bs50.bl.RPKM.bw \
HepG2_G4ChIP_REP1_c.SLX-19356.i702_i504.bs50.bl.RPKM.bw \
HepG2_G4ChIP_REP2_a.SLX-19356.i703_i502.bs50.bl.RPKM.bw \
HepG2_G4ChIP_REP2_b.SLX-19356.i704_i503.bs50.bl.RPKM.bw \
HepG2_G4ChIP_REP2_c.SLX-19356.i704_i504.bs50.bl.RPKM.bw \
HepG2_G4ChIP_REP3_a.SLX-19356.i705_i502.bs50.bl.RPKM.bw \
HepG2_G4ChIP_REP3_b.SLX-19356.i706_i503.bs50.bl.RPKM.bw \
HepG2_G4ChIP_REP3_c.SLX-19356.i706_i504.bs50.bl.RPKM.bw \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio1_a.bs50.bl.RPKM.bw \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio1_c.bs50.bl.RPKM.bw \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio2_a.bs50.bl.RPKM.bw \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio2_b.bs50.bl.RPKM.bw \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio2_c.bs50.bl.RPKM.bw \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio3_a.bs50.bl.RPKM.bw \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio3_b.bs50.bl.RPKM.bw \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio3_c.bs50.bl.RPKM.bw \
-out ../analysis/deeptools/G4ChIP_K562_HepG2_MultiBWSummary_per_50bin2.npz \
--binSize 50 \
--BED ../analysis/deeptools/G4ChIP_K562_HepG2_union.merged_500.bed \
--numberOfProcessors max \
--outRawCounts ../analysis/deeptools/G4ChIP_K562_HepG2_MultiBWSummary_per_50bin2.tab &&
plotCorrelation \
-in ../analysis/deeptools/G4ChIP_K562_HepG2_MultiBWSummary_per_50bin2.npz \
--whatToPlot heatmap \
--corMethod spearman \
--plotTitle G4ChIP_K562_HepG2_spearman \
--outFileCorMatrix ../analysis/deeptools/G4ChIP_K562_HepG2_MultiBWSummary_per_50bin2.tab \
--plotNumbers \
-o ../analysis/deeptools/G4ChIP_K562_HepG2_MultiBWSummary_per_50bin_spearman.pdf"
```



## Intervene analysis

```
mkdir ../analysis/Intervene
cd  ../../macs2_output/


sbatch -J Intervene -o pw.log --mem 2g --wrap " intervene pairwise -i \
HepG2_G4ChIP_REP1_a.SLX-19356.i701_i502.nodup.q005.all_peaks.narrowPeak \
HepG2_G4ChIP_REP1_b.SLX-19356.i702_i503.nodup.q005.all_peaks.narrowPeak \
HepG2_G4ChIP_REP1_c.SLX-19356.i702_i504.nodup.q005.all_peaks.narrowPeak \
HepG2_G4ChIP_REP2_a.SLX-19356.i703_i502.nodup.q005.all_peaks.narrowPeak \
HepG2_G4ChIP_REP2_b.SLX-19356.i704_i503.nodup.q005.all_peaks.narrowPeak \
HepG2_G4ChIP_REP2_c.SLX-19356.i704_i504.nodup.q005.all_peaks.narrowPeak \
HepG2_G4ChIP_REP3_a.SLX-19356.i705_i502.nodup.q005.all_peaks.narrowPeak \
HepG2_G4ChIP_REP3_b.SLX-19356.i706_i503.nodup.q005.all_peaks.narrowPeak \
HepG2_G4ChIP_REP3_c.SLX-19356.i706_i504.nodup.q005.all_peaks.narrowPeak \
--names \
--figsize 7 4 \
Rep1_a,Rep1_b,Rep1_c,Rep2_a,Rep2_b,Rep2_c,Rep3_a,Rep3_b,Rep3_c \
--compute frac --htype tribar \
-o HepG2_BG4_ChIP_replciates_PW_triang"
```


## PAVIS
Shuffle confidence peask in OQS to calculate enrichment from PAVIS analysis.
```
for f in {1..5}
do
echo run $f
bedtools shuffle \
-excl /scratchb/sblab/spiege01/ENCODE_K562/reference_data/hg19/hg19.blacklist_merge1000nt.bed \
-incl /scratchb/sblab/spiege01/ENCODE_K562/reference_data/Quadruplex/G4-Seq_cat_K_PDS_+-strands.bed \
-maxTries 1000 \
-noOverlapping \
-i HepG2_async_rep1-3.mult.6of9.bed \
-g /scratchb/sblab/spiege01/ENCODE_K562/reference_data/GM/hg19/hg19.genome2 > ../Bedshuffle/HepG2_async_rep1-3.mult.6of9_shuffled_$f.bed
done
```



