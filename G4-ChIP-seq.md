# G4-ChIP-seq analysis
For a more detailed description of our in-house ChIP-seq pipeline please refer to our [previous work](https://github.com/sblab-bioinformatics/dna-secondary-struct-chrom-lands/blob/master/Methods.md).
Briefly, raw fastq reads from G4-ChIP-seq in K562 cells were trimmed with cutadapt to remove adapter sequences and low-quality reads (mapping quality < 10). Reads were aligned to the human genome (version hg19) with BWA and duplicates were removed using Picard. Peaks were called by MACS2 (p < 1e-05). Peaks were merged from different replicates with bedtools multiIntersect. Only peaks overlapping in 5 out 8 replicates were considered high-confidence. 



# Generation of pseudo-stranded G4 ChIP-seq 
## Aim:
G4 ChIP-seq is not stranded. Use stranded information for G4-seq (OQS)to inferr standinformation for G4 ChIP-seq high confidence peaks.

### Stranded OQS map

Download G4-seq bedfiles from GEO (GSE63874):

```
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63874/suppl/GSE63874%5FNa%5FK%5Fminus%5Fhits%5Fintersect%2Ebed%2Egz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63874/suppl/GSE63874%5FNa%5FPDS%5Fminus%5Fhits%5Fintersect%2Ebed%2Egz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63874/suppl/GSE63874%5FNa%5FPDS%5Fplus%5Fhits%5Fintersect%2Ebed%2Egz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63874/suppl/GSE63874%5FNa%5FK%5Fplus%5Fhits%5Fintersect%2Ebed%2Egz
```

Concatinate PDS and K OQs on + and -strand:
```
zcat GSE63874_Na_K_minus_hits_intersect.bed.gz GSE63874_Na_PDS_minus_hits_intersect.bed.gz | bedtools sort | bedtools merge | bedtools sort | gzip > G4-Seq_K_PDS_minus.bed.gz
zcat GSE63874_Na_K_plus_hits_intersect.bed  GSE63874_Na_PDS_plus_hits_intersect.bed | bedtools sort | bedtools merge | bedtools sort | gzip > G4-Seq_K_PDS_plus.bed.gz

zcat G4-Seq_K_PDS_minus.bed.gz | awk '{print $1"\t"$2"\t"$3"\t.\t.\t-\t"}' - > G4-Seq_K_PDS_stranded.bed
zcat G4-Seq_K_PDS_plus.bed.gz | awk '{print $1"\t"$2"\t"$3"\t.\t.\t+\t"}' - >> G4-Seq_K_PDS_stranded.bed
cat G4-Seq_K_PDS_stranded.bed  | sortBed -i | gzip > G4-Seq_K_PDS_stranded.bed.gz
rm G4-Seq_K_PDS_stranded.bed
```
Non-stranded version:
```
bedtools merge -i G4-Seq_K_PDS_stranded.bed.gz > G4-Seq_cat_K_PDS_+-strands.bed   

```


### Stranded G4 ChIP
Intersect stranded OQS with G4 ChIP-seq peaks and add strand information 
```
bedtools intersect -wa -a 20180108_K562_async_rep1-3.mult.5of8.bed -b G4-Seq_K_PDS_plus.bed.gz | bedtools merge | bedtools sort | awk -F '\t' '{print $1"\t"$2"\t"$3"\t"".""\t"".""\t""+""\t"}' | gzip > K562_async_rep1-3.mult.5of8_plusstrand.bed.gz

bedtools intersect -wa -a 20180108_K562_async_rep1-3.mult.5of8.bed -b G4-Seq_K_PDS_minus.bed.gz | bedtools merge | bedtools sort | awk -F '\t' '{print $1"\t"$2"\t"$3"\t"".""\t"".""\t""-""\t"}' | gzip > K562_async_rep1-3.mult.5of8_minusstrand.bed.gz
```

Remove ambigous cases, and retain peaks exclusive to a single strand:
```
bedtools intersect -v -a K562_async_rep1-3.mult.5of8_plusstrand.bed.gz -b K562_async_rep1-3.mult.5of8_minusstrand.bed.gz| bedtools sort | bedtools merge | gzip > K562_async_rep1-3.mult.5of8_plusstrand_exclusive.bed.gz

zcat K562_async_rep1-3.mult.5of8_plusstrand_exclusive.bed.gz | wc -l
2161


bedtools intersect -v -a K562_async_rep1-3.mult.5of8_minusstrand.bed.gz -b K562_async_rep1-3.mult.5of8_plusstrand.bed.gz | bedtools sort | bedtools merge | gzip > K562_async_rep1-3.mult.5of8_minusstrand_exclusive.bed.gz

zcat K562_async_rep1-3.mult.5of8_minusstrand_exclusive.bed.gz | wc -l
2271
```
> In Summary, out of the 9205 peaks, 2161 peaks are exclusive to the plus strand and 2271 are exclusive to the minus strand, while 3494 BG4 peaks cannot be assinged unambigously.


Add '+' for positive strand and '-' for negative strand, then combine into a single bedfile.
```
awk '{print $1"\t"$2"\t"$3"\t.\t.\t-\t"}' K562_async_rep1-3.mult.5of8_minusstrand_exclusive.bed > K562_async_rep1-3.mult.5of8_OQS-Stranded.bed
awk '{print $1"\t"$2"\t"$3"\t.\t.\t+\t"}' K562_async_rep1-3.mult.5of8_plusstrand_exclusive.bed >> K562_async_rep1-3.mult.5of8_OQS-Stranded.bed
cat K562_async_rep1-3.mult.5of8_OQS-Stranded.bed  | sortBed -i | gzip > K562_async_rep1-3.mult.5of8_OQS_Stranded.bed.gz

```



# G4 control data set (OQS in open Chromatin around TSS)
Generate a control data set with sites that have the potential to form G4s (G-seq) in open chromatin (DHS), reflecting the strong enrichemnt of G4s around TSS. Generate two versions and without strand information. 

```

intersectBed -wa -a G4-Seq_K_PDS_stranded.bed.gz -b DNAse-seq.concatenated_narrow_rep1_and_rep2.bed > openOQs_stranded.bed

intersectBed -v -a openOQs_stranded.bed -b 20180108_K562_async_rep1-3.mult.5of8.bed > openOQs_noBG4_stranded.bed

bedtools merge -i openOQs_noBG4_stranded.bed > openOQs_noBG4.bed

```


Use Pavis tool (https://manticore.niehs.nih.gov/pavis2/) with settings "hg19" "known gene".
In the derived txt files, select upstream TSS and 5'UTR
```

grep -e UTR -e Upstream openOQS_noBG4_1kbupTSS.txt | cut -f 2,3,4 | sortBed -i > openOQS_noBG4_1kbupTSS.bed
grep -e UTR -e Upstream openOQs_noBG4_stranded_1kbupTSS.txt | cut -f 2,3,4,5,6,7 | sortBed -i > temp.bed
awk '{print $1"\t"$2"\t"$3"\t.\t.\t"$6}' temp.bed > openOQs_noBG4_Stranded_1kbupTSS.bed
rm temp.bed
```



