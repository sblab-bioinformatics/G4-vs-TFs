
| Sample | Index | SLX-Number |	
| --- | --- | --- | 
| K562_rep2_F_input_Input_DMSO_i709-i501 | D709dna-D501dna | SLX-15865 |
| K562_rep1_A1_SP2_50uMPDS_i710-i502 | D710dna-D502dna | SLX-15865 |
| K562_rep1_A2_SP2_50uMPDS_i711-i503 | D711dna-D503dna | SLX-15865 |
| K562_rep1_A3_SP2_50uMPDS_i712-i504 | D712dna-D504dna | SLX-15865 |
| K562_rep1_B1_SP2_DMSO_i701-i505 | D701dna-D505dna | SLX-15865 |
| K562_rep1_B2_SP2_DMSO_i702-i506 | D702dna-D506dna | SLX-15865 |
| K562_rep1_B3_SP2_DMSO_i703-i507 | D703dna-D507dna | SLX-15865 |
| K562_rep1_B_input_Input_DMSO_i710-i506 | D710dna-D506dna | SLX-15865 |



### Generate relevant folders for analysis

```bash
mkdir fastq
mkdir fastq2
mkdir fastq/fastqc
mkdir fastQC
mkdir bam
mkdir merged_bam
mkdir macs2_output
mkdir tdf
mkdir bigWig
mkdir logs
mkdir analysis
```



### Download data from basespace and batch rename

Download from basespace. As I have resequenced the same pool fq files of both runs will be named the same. Therefore I have to transfer fq files into 2 different folders fastq and fastq2 and merge by name at a later stage. Ignore unindexed reads for the moment.

```bash

for f in K562_NativeChIP-136472337/FASTQ_Generation_2019-07-18_00_16_44Z-189977835/*/*.gz 
do
fq_basename=`basename $f`
mv $f fastq/$fq_basename
done


for f in K562_NativeChIP-136472337/FASTQ_Generation_2019-07-19_02_25_01Z-190049955/*/*.gz 
do
fq_basename=`basename $f`
mv $f fastq2/$fq_basename
done


```

### Quality check
```bash
cd fastq

for fq in *.fastq.gz
do
  bname=${fq%.fastq.gz}
  sbatch -J $bname -o fastqc/$bname.log --mem 4G --wrap "/home/bioinformatics/software/fastqc/fastqc-v0.11.5/fastqc --noextract --nogroup -q -o fastqc $fq"
done
```

Libraries look  quite unbalanced,  overall FastQC of the reads that made it past the filter is fine. Save html files remove rest to save space.
```bash
cd fastqc
for f in *.html
do
mv $f ../../fastQC/$f
done

cd ..
rm -r fastqc
```



### adapter trimming for single-end sequencing 
(ignore lost of failed libraries and unindexed reads)

```
cd fastq
mkdir fastq_trimmed

for f in SLX-15865*.gz; do sbatch --mem 8G -o $f.out -e $f.err -J cutadapt --wrap "~/.local/bin/cutadapt -m 10 -q 20 -O 3 -a AGATCGGAAGAGC -g GCTCTTCCGATCT -o fastq_trimmed/${f%%.fastq.gz}.trimmed.fq.gz $f"; done


cd ../fastq2
mkdir fastq_trimmed

for f in SLX-15865*.gz; do sbatch --mem 8G -o $f.out -e $f.err -J cutadapt --wrap "~/.local/bin/cutadapt -m 10 -q 20 -O 3 -a AGATCGGAAGAGC -g GCTCTTCCGATCT -o fastq_trimmed/${f%%.fastq.gz}.trimmed.fq.gz $f"; done
```



### alignment single-end to hg19
Align to hg19 using Giovannis reference genome. Giovanni's tools have be removed by IT, so I switch to the bioinformatics core moving forward (same version number 0.7.15, so this should not be a problem).

Align:
```
cd ../fastq/fastq_trimmed

genome='/scratchb/sblab/spiege01/ENCODE_K562/reference_data/GM/hg19/Sequence/Bowtie2Index/bowtie2-build/genome.fa'; for f in *.trimmed.fq.gz; do sbatch --mem 16G -o $f.bwa.out -e $f.bwa.err -J bwa --wrap "/home/bioinformatics/software/bwa/bwa-0.7.15/bwa mem -t 8 -M $genome $f | /Users/spiege01/bin/samtools view -Sb -F2308 -q 10 > /scratchb/sblab/spiege01/ENCODE_K562/Native_ChIP/SLX-18311_nChIP_SP2_FUS_vsPDS_rep1-2/bam/${f%%.trimmed.fq.gz}.hg19.bam"; done

cd ../../fastq2/fastq_trimmed
genome='/scratchb/sblab/spiege01/ENCODE_K562/reference_data/GM/hg19/Sequence/Bowtie2Index/bowtie2-build/genome.fa'; for f in *.trimmed.fq.gz; do sbatch --mem 16G -o $f.bwa.out -e $f.bwa.err -J bwa --wrap "/home/bioinformatics/software/bwa/bwa-0.7.15/bwa mem -t 8 -M $genome $f | /Users/spiege01/bin/samtools view -Sb -F2308 -q 10 > /scratchb/sblab/spiege01/ENCODE_K562/Native_ChIP/SLX-18311_nChIP_SP2_FUS_vsPDS_rep1-2/bam/${f%%.trimmed.fq.gz}b.hg19.bam"; done


```


Merge files of same conditions
```
cd ../../bam

for f in *_L001_R1_001.hg19.bam
do
  sbatch --mem 8G -o $f.out -e $f.err -J merge --wrap "samtools merge -@8 ${f%%_L001_R1_001.*}.hg19.merged.bam ${f%%_L001_R1_001*}_L00*_R1_001*.hg19.bam"
done
```

Rename bam files and move to merged_bam folder. (Use basespace Samplesheet to generate command lines in excel -> make new table and then use Atom to replace tabs with space, etc.)

```
mv SLX-15865_S10.hg19.merged.bam K562_rep1_A3_SP2_50uMPDS_i712-i504_SLX-18311-2.hg19.merged.bam
mv SLX-15865_S11.hg19.merged.bam K562_rep2_F4_SP2_DMSO_i708-i508_SLX-18311-2.hg19.merged.bam
mv SLX-15865_S12.hg19.merged.bam K562_rep2_A2_SP2_50uMPDS_i702-i502_SLX-18311-2.hg19.merged.bam
mv SLX-15865_S13.hg19.merged.bam K562_rep2_A4_SP2_50uMPDS_i704-i504_SLX-18311-2.hg19.merged.bam
mv SLX-15865_S14.hg19.merged.bam K562_rep1_A2_SP2_50uMPDS_i711-i503_SLX-18311-2.hg19.merged.bam
mv SLX-15865_S15.hg19.merged.bam K562_rep1_B2_SP2_DMSO_i702-i506_SLX-18311-2.hg19.merged.bam
mv SLX-15865_S16.hg19.merged.bam K562_rep2_A1_SP2_50uMPDS_i701-i501_SLX-18311-2.hg19.merged.bam
mv SLX-15865_S17.hg19.merged.bam K562_rep1_B_input_Input_DMSO_i710-i506_SLX-18311-2.hg19.merged.bam
mv SLX-15865_S18.hg19.merged.bam K562_rep2_A3_SP2_50uMPDS_i703-i503_SLX-18311-2.hg19.merged.bam
mv SLX-15865_S19.hg19.merged.bam K562_rep1_D2_FUS_DMSO_i708-i504_SLX-18311-2.hg19.merged.bam
mv SLX-15865_S1.hg19.merged.bam K562_rep1_D3_FUS_DMSO_i709-i505_SLX-18311-2.hg19.merged.bam
mv SLX-15865_S20.hg19.merged.bam K562_rep2_F2_SP2_DMSO_i706-i506_SLX-18311-2.hg19.merged.bam
mv SLX-15865_S21.hg19.merged.bam K562_rep2_F_input_Input_DMSO_i709-i501_SLX-18311-2.hg19.merged.bam
mv SLX-15865_S22.hg19.merged.bam K562_rep1_B1_SP2_DMSO_i701-i505_SLX-18311-2.hg19.merged.bam
mv SLX-15865_S2.hg19.merged.bam K562_rep1_A1_SP2_50uMPDS_i710-i502_SLX-18311-2.hg19.merged.bam
mv SLX-15865_S3.hg19.merged.bam K562_rep2_F1_SP2_DMSO_i705-i505_SLX-18311-2.hg19.merged.bam
mv SLX-15865_S4.hg19.merged.bam K562_rep2_F3_SP2_DMSO_i707-i507_SLX-18311-2.hg19.merged.bam
mv SLX-15865_S5.hg19.merged.bam K562_rep1_D1_FUS_DMSO_i707-i503_SLX-18311-2.hg19.merged.bam
mv SLX-15865_S6.hg19.merged.bam K562_rep1_C1_FUS_50uMPDS_i704-i508_SLX-18311-2.hg19.merged.bam
mv SLX-15865_S7.hg19.merged.bam K562_rep1_B3_SP2_DMSO_i703-i507_SLX-18311-2.hg19.merged.bam
mv SLX-15865_S8.hg19.merged.bam K562_rep1_C3_FUS_50uMPDS_i706-i502_SLX-18311-2.hg19.merged.bam
mv SLX-15865_S9.hg19.merged.bam K562_rep1_C2_FUS_50uMPDS_i705-i501_SLX-18311.hg19.merged.bam

for f in *hg19.merged.bam; do mv $f ../merged_bam/$f; done
```


### Sort bam files, mark and remove duplicates
```
cd ../merged_bam 

for f in *.bam; do sbatch --mem 16G -o $f.out -e $f.err -J sort --wrap "/Users/spiege01/bin/samtools sort -@ 8 $f > ${f%%.merged.bam}.sort.bam"; done
```


Mark and remove duplicates
```
mkdir nodup
for f in *sort.bam; do sbatch --mem 8G -o $f.out -e $f.err -J dups --wrap "java -Xmx7g -jar /home/bioinformatics/software/picard/picard-2.14.0/picard.jar MarkDuplicates I=$f O=nodup/${f%%.sort.bam}.nodup.bam M=${f%%.sort.bam}.md.txt AS=true REMOVE_DUPLICATES=true"; done
```



Some basic dup and size stats
```
grep 'Unknown Library' *txt > dups_info.txt

cd nodup/
for f in *nodup.bam; do echo $f;echo $f >> all_counts && /Users/spiege01/bin/samtools view -c $f >> all_counts; done
```

**RESULTS**
Duplication rate
```
K562_rep1_A1_SP2_50uMPDS_i710-i502_SLX-18311-2.hg19.md.txt:Unknown Library	13715644	0	0	0	914500	0	0	0.066676	
K562_rep1_A2_SP2_50uMPDS_i711-i503_SLX-18311-2.hg19.md.txt:Unknown Library	18273296	0	0	0	616515	0	0	0.033739	
K562_rep1_A3_SP2_50uMPDS_i712-i504_SLX-18311-2.hg19.md.txt:Unknown Library	18249566	0	0	0	1075344	0	0	0.058924	
K562_rep1_B1_SP2_DMSO_i701-i505_SLX-18311-2.hg19.md.txt:Unknown Library	11308114	0	0	0	481355	0	0	0.042567	
K562_rep1_B2_SP2_DMSO_i702-i506_SLX-18311-2.hg19.md.txt:Unknown Library	18595385	0	0	0	1929675	0	0	0.103772	
K562_rep1_B3_SP2_DMSO_i703-i507_SLX-18311-2.hg19.md.txt:Unknown Library	12838571	0	0	0	713414	0	0	0.055568	
K562_rep1_B_input_Input_DMSO_i710-i506_SLX-18311-2.hg19.md.txt:Unknown Library	23549966	0	0	0	768359	0	0	0.032627	
K562_rep1_C1_FUS_50uMPDS_i704-i508_SLX-18311-2.hg19.md.txt:Unknown Library	19045020	0	0	0	1181916	0	0	0.062059	
K562_rep1_C2_FUS_50uMPDS_i705-i501_SLX-18311.hg19.md.txt:Unknown Library	20697598	0	0	0	3127237	0	0	0.151092	
K562_rep1_C3_FUS_50uMPDS_i706-i502_SLX-18311-2.hg19.md.txt:Unknown Library	21045287	0	0	0	2191891	0	0	0.104151	
K562_rep1_D1_FUS_DMSO_i707-i503_SLX-18311-2.hg19.md.txt:Unknown Library	13479461	0	0	0	541943	0	0	0.040205	
K562_rep1_D2_FUS_DMSO_i708-i504_SLX-18311-2.hg19.md.txt:Unknown Library	9457822	0	0	0	391243	0	0	0.041367	
K562_rep1_D3_FUS_DMSO_i709-i505_SLX-18311-2.hg19.md.txt:Unknown Library	12194706	0	0	0	1058728	0	0	0.086819	
K562_rep2_A1_SP2_50uMPDS_i701-i501_SLX-18311-2.hg19.md.txt:Unknown Library	29728889	0	0	0	1058596	0	0	0.035608	
K562_rep2_A2_SP2_50uMPDS_i702-i502_SLX-18311-2.hg19.md.txt:Unknown Library	33877182	0	0	0	1741975	0	0	0.05142	
K562_rep2_A3_SP2_50uMPDS_i703-i503_SLX-18311-2.hg19.md.txt:Unknown Library	24335660	0	0	0	842027	0	0	0.034601	
K562_rep2_A4_SP2_50uMPDS_i704-i504_SLX-18311-2.hg19.md.txt:Unknown Library	36249376	0	0	0	1578358	0	0	0.043542	
K562_rep2_F1_SP2_DMSO_i705-i505_SLX-18311-2.hg19.md.txt:Unknown Library	24558154	0	0	0	2175258	0	0	0.088576	
K562_rep2_F2_SP2_DMSO_i706-i506_SLX-18311-2.hg19.md.txt:Unknown Library	11685373	0	0	0	2052273	0	0	0.175628	
K562_rep2_F3_SP2_DMSO_i707-i507_SLX-18311-2.hg19.md.txt:Unknown Library	27821354	0	0	0	1445039	0	0	0.05194	
K562_rep2_F4_SP2_DMSO_i708-i508_SLX-18311-2.hg19.md.txt:Unknown Library	18620264	0	0	0	943806	0	0	0.050687	
K562_rep2_F_input_Input_DMSO_i709-i501_SLX-18311-2.hg19.md.txt:Unknown Library	27129941	0	0	0	1070461	0	0	0.039457
```
Low number of duplicates only F2_SP_DMSO has 17%.


Read count per bam file:
```
K562_rep1_A1_SP2_50uMPDS_i710-i502_SLX-18311-2.hg19.nodup.bam 12801144
K562_rep1_A2_SP2_50uMPDS_i711-i503_SLX-18311-2.hg19.nodup.bam 17656781
K562_rep1_A3_SP2_50uMPDS_i712-i504_SLX-18311-2.hg19.nodup.bam 17174222
K562_rep1_B1_SP2_DMSO_i701-i505_SLX-18311-2.hg19.nodup.bam 10826759
K562_rep1_B2_SP2_DMSO_i702-i506_SLX-18311-2.hg19.nodup.bam 16665710
K562_rep1_B3_SP2_DMSO_i703-i507_SLX-18311-2.hg19.nodup.bam 12125157
K562_rep1_B_input_Input_DMSO_i710-i506_SLX-18311-2.hg19.nodup.bam 22781607
K562_rep1_C1_FUS_50uMPDS_i704-i508_SLX-18311-2.hg19.nodup.bam 17863104
K562_rep1_C2_FUS_50uMPDS_i705-i501_SLX-18311.hg19.nodup.bam 17570361
K562_rep1_C3_FUS_50uMPDS_i706-i502_SLX-18311-2.hg19.nodup.bam 18853396
K562_rep1_D1_FUS_DMSO_i707-i503_SLX-18311-2.hg19.nodup.bam 12937518
K562_rep1_D2_FUS_DMSO_i708-i504_SLX-18311-2.hg19.nodup.bam 9066579
K562_rep1_D3_FUS_DMSO_i709-i505_SLX-18311-2.hg19.nodup.bam 11135978
K562_rep2_A1_SP2_50uMPDS_i701-i501_SLX-18311-2.hg19.nodup.bam 28670293
K562_rep2_A2_SP2_50uMPDS_i702-i502_SLX-18311-2.hg19.nodup.bam 32135207
K562_rep2_A3_SP2_50uMPDS_i703-i503_SLX-18311-2.hg19.nodup.bam 23493633
K562_rep2_A4_SP2_50uMPDS_i704-i504_SLX-18311-2.hg19.nodup.bam 34671018
K562_rep2_F1_SP2_DMSO_i705-i505_SLX-18311-2.hg19.nodup.bam 22382896
K562_rep2_F2_SP2_DMSO_i706-i506_SLX-18311-2.hg19.nodup.bam 9633100
K562_rep2_F3_SP2_DMSO_i707-i507_SLX-18311-2.hg19.nodup.bam 26376315
K562_rep2_F4_SP2_DMSO_i708-i508_SLX-18311-2.hg19.nodup.bam 17676458
K562_rep2_F_input_Input_DMSO_i709-i501_SLX-18311-2.hg19.nodup.bam 26059480



```
~20 Mio reads in most of the cases. Not a lot, but this should still be enough for robust peak calling and DBA. Can still eliminate libraries with low read count.

### create index bam files and generate tdf files for visual inspection
Generate bai files using samtools index
```
for f in *nodup.bam; do sbatch -o %j.$f.log --mem 16G -J index.$f --wrap "/home/bioinformatics/software/samtools/samtools-1.6/bin/samtools index -b $f"; done

```

Generate tdfs
```

for f in *nodup.bam; do sbatch -o %j.$f.log --mem 8G -J igv.$f --wrap "/Users/martin03/sw/igvtools/igvtools-2.3.91/igvtools count -w 15 $f ../../tdf/${f%%.nodup.bam}.tdf /Users/martin03/sw/igvtools/igvtools-2.3.91/genomes/hg19.chrom.sizes"; done
```

### create bigWig files
```
for FILE in *nodup.bam 
do
bname=`basename $FILE .bam`
sbatch -o %j.$f.log -J BCoverage --mem=8G \
--wrap "bamCoverage -b $FILE \
-o ../../bigWig/${FILE%%.nodup.bam}.bs50.bl.RPKM.bw \
--binSize 50 \
--blackListFileName /scratchb/sblab/spiege01/ENCODE_K562/reference_data/hg19/hg19.blacklist_merge1000nt.bed \
--numberOfProcessors max \
--normalizeUsing RPKM"
done
```
> SP2: Visual inspection suggests that native ChIPs show strong overlap with X-ChIP. Seems like displacement was extremely efficient in rep2 (or PDS replicates failed), and also strong in rep1. FUS does not show any peaks. Probably need to increase antibody input as for other replicates.



### call peaks 
(using macs2 2.1.1.20160309)
Call peaks vs input file. 

```
for f in *rep1*nodup.bam; do sbatch --mem 8G -o $f.macs2.out -e $f.macs2.err -J macs2 --wrap "macs2 callpeak --name ../../macs2_output/${f%%.nodup.bam}.nodup.q005.all -t $f -c K562_rep1_B_input_Input_DMSO_i710-i506_SLX-18311-2.hg19.nodup.bam --format=BAM --gsize 'hs' --bw=300 --qvalue 0.05"; done

for f in *rep2*nodup.bam; do sbatch --mem 8G -o $f.macs2.out -e $f.macs2.err -J macs2 --wrap "macs2 callpeak --name ../../macs2_output/${f%%.nodup.bam}.nodup.q005.all -t $f -c K562_rep2_F_input_Input_DMSO_i709-i501_SLX-18311-2.hg19.nodup.bam --format=BAM --gsize 'hs' --bw=300 --qvalue 0.05"; done

```



Results:
```
cd ../../macs2_output/
wc -l *.narrowPeak

1114 K562_rep1_A1_SP2_50uMPDS_i710-i502_SLX-18311-2.hg19.nodup.q005.all_peaks.narrowPeak
 311 K562_rep1_A2_SP2_50uMPDS_i711-i503_SLX-18311-2.hg19.nodup.q005.all_peaks.narrowPeak
1453 K562_rep1_A3_SP2_50uMPDS_i712-i504_SLX-18311-2.hg19.nodup.q005.all_peaks.narrowPeak
2045 K562_rep1_B1_SP2_DMSO_i701-i505_SLX-18311-2.hg19.nodup.q005.all_peaks.narrowPeak
3360 K562_rep1_B2_SP2_DMSO_i702-i506_SLX-18311-2.hg19.nodup.q005.all_peaks.narrowPeak
1987 K562_rep1_B3_SP2_DMSO_i703-i507_SLX-18311-2.hg19.nodup.q005.all_peaks.narrowPeak
  50 K562_rep1_C1_FUS_50uMPDS_i704-i508_SLX-18311-2.hg19.nodup.q005.all_peaks.narrowPeak
 107 K562_rep1_C2_FUS_50uMPDS_i705-i501_SLX-18311.hg19.nodup.q005.all_peaks.narrowPeak
  49 K562_rep1_C3_FUS_50uMPDS_i706-i502_SLX-18311-2.hg19.nodup.q005.all_peaks.narrowPeak
  55 K562_rep1_D1_FUS_DMSO_i707-i503_SLX-18311-2.hg19.nodup.q005.all_peaks.narrowPeak
  32 K562_rep1_D2_FUS_DMSO_i708-i504_SLX-18311-2.hg19.nodup.q005.all_peaks.narrowPeak
 114 K562_rep1_D3_FUS_DMSO_i709-i505_SLX-18311-2.hg19.nodup.q005.all_peaks.narrowPeak
   0 K562_rep2_A1_SP2_50uMPDS_i701-i501_SLX-18311-2.hg19.nodup.q005.all_peaks.narrowPeak
   1 K562_rep2_A2_SP2_50uMPDS_i702-i502_SLX-18311-2.hg19.nodup.q005.all_peaks.narrowPeak
   0 K562_rep2_A3_SP2_50uMPDS_i703-i503_SLX-18311-2.hg19.nodup.q005.all_peaks.narrowPeak
   0 K562_rep2_A4_SP2_50uMPDS_i704-i504_SLX-18311-2.hg19.nodup.q005.all_peaks.narrowPeak
3032 K562_rep2_F1_SP2_DMSO_i705-i505_SLX-18311-2.hg19.nodup.q005.all_peaks.narrowPeak
4890 K562_rep2_F2_SP2_DMSO_i706-i506_SLX-18311-2.hg19.nodup.q005.all_peaks.narrowPeak
2070 K562_rep2_F3_SP2_DMSO_i707-i507_SLX-18311-2.hg19.nodup.q005.all_peaks.narrowPeak
3590 K562_rep2_F4_SP2_DMSO_i708-i508_SLX-18311-2.hg19.nodup.q005.all_peaks.narrowPeak

```







> Peak numbers look a lot better after reseqencing. FUS is not useful. Used more anti-FUS in 3rd replicate and this seems to help a lot. Focus on SP2 from here on.



# QC native ChIP samples
## Compare native and X-ChIP
### Overlap of native ChIP with X-ChIP peaks for DMSO
For more statistical power should use different peakcaller (e.g. PePr). 

```
#peak numbers
# rep 1
multiIntersectBed -i *rep1*SP2_DMSO*.narrowPeak | awk '{if ($4 > 1) print $1 "\t" $2 "\t" $3}' | sortBed | mergeBed > ../analysis/K562_REP1_SP2_DMSO_SLX-18311-2_Mult_2out3.bed 
wc -l K562_REP1_SP2_DMSO_SLX-18311-2_Mult_2out3.bed 2234

# rep 1&2 
multiIntersectBed -i *SP2_DMSO*.narrowPeak | awk '{if ($4 > 3) print $1 "\t" $2 "\t" $3}' | sortBed | mergeBed > ../analysis/K562_SP2_DMSO_SLX-18311-2_Mult_4out7.bed 

wc -l K562_SP2_DMSO_SLX-18311-2_Mult_4out7.bed 2686


cd ../analysis
for f in *.bed; do intersectBed -a $f -b ../macs2_output/K562_SP2_DMSO_SLX-18311-2_Mult_4out7.bed | wc -l; echo out of; wc -l $f; done


1546 (58%) out of 3124 ENCFF002CMO.cut.bed
1884 (70%) out of 7554 K562_ReChIP_SP2-BG4_biorep1-3.SLX-16769.hg19.nodup.q005.all_peaks_2out3.bed
1122 (43%) out of 1444 K562_SP2_DMSO_3h_mult_2out2.bed
```
If I take my own SP2 consensus peak set, 70% of native peaks overlap.  This should be good enough.

Test overlapp with BG4:
```
intersectBed -a K562_SP2_DMSO_SLX-18311-2_Mult_4out7.bed -b /scratchb/sblab/spiege01/ENCODE_K562/reference_data/Quadruplex/20180108_K562_async_rep1-3.mult.5of8.bed  | wc -l
1956 (73%)
```
**73% of native SP2 peaks overlap with BG4 which is similar to what BG4 technical replicates overlap. This suggests that nearly all of the SP2 sites are BG4 related and the rest might be lost in noise!**




### Perform similar analysis for PDS-treated sample.

```
#peak numbers
# rep 1
multiIntersectBed -i *rep1*SP2_50uMPDS*.narrowPeak | awk '{if ($4 > 1) print $1 "\t" $2 "\t" $3}' | sortBed | mergeBed > ../analysis/K562_REP1_SP2_50uMPDS_SLX-18311-2_Mult_2out3.bed 
wc -l K562_REP1_SP2_50uMPDS_SLX-18311-2_Mult_2out3.bed  986

# rep 1&2 
multiIntersectBed -i *SP2_50uMPDS*.narrowPeak | awk '{if ($4 > 3) print $1 "\t" $2 "\t" $3}' | sortBed | mergeBed > ../analysis/K562_SP2_50uMPDS_SLX-18311-2_Mult_4out7.bed 

wc -l K562_SP2_50uMPDS_SLX-18311-2_Mult_4out7.bed  1


cd ../analysis
for f in *.bed; do intersectBed -a $f -b K562_REP1_SP2_50uMPDS_SLX-18311-2_Mult_2out3.bed | wc -l; echo out of; wc -l $f; done

750 out of 3124 ENCFF002CMO.cut.bed
980 out of 2234 K562_REP1_SP2_DMSO_SLX-18311-2_Mult_2out3.bed

847 out of 7554 K562_ReChIP_SP2-BG4_biorep1-3.SLX-16769.hg19.nodup.q005.all_peaks_2out3.bed
614 out of 1444 K562_SP2_DMSO_3h_mult_2out2.bed


```
Test overlapp with BG4:
```
intersectBed -a K562_REP1_SP2_50uMPDS_SLX-18311-2_Mult_2out3.bed -b /scratchb/sblab/spiege01/ENCODE_K562/reference_data/Quadruplex/20180108_K562_async_rep1-3.mult.5of8.bed  | wc -l
718 (73%)
```
**Upon PDS treatment number of peaks drop from 2234 to 987 (44%) of which 980 (99%) are a subset of the DMSO peaks. 750 (75%) overlapp with ENCODE data and 718 (73%) overlap with BG4.**

### Venn diagrams:
Generate venn diagrams comparing overlap of replicates
```bash
intervene venn -i *rep1*SP2_DMSO*.narrowPeak \
--names \
rep1_DMSO_1,rep1_DMSO_2,rep1_DMSO_3  \
-o ../analysis/Intervene/Overlap_rep1_SP2_DMSO

intervene venn -i *rep1*SP2_50uM*.narrowPeak \
--names \
rep1_50uMPDS_1,rep1_50uMPDS_2,rep1_50uMPDS_3  \
-o ../analysis/Intervene/Overlap_rep1_SP2_50uMPDS

```
Redraw with R/Illustrator:
```R
library(VennDiagram)

# Venn diagram for PDS (numbers from intervene)
venn.plot <- draw.triple.venn(
    area1 = 1114,
  area2 = 311,
  area3 = 1453,
  n12 = 302,
  n13 = 980,
  n23 =302,
  n123 = 301,
  #scaled = T,
  fill = c("red3", "forestgreen", "deepskyblue3"),
  cex = 1.5,
  fontfamily = "sans",
  category = c("rep1", "rep2", "rep3"),
  cat.dist = c(0.08, 0.08, 0.04),
  cat.cex = 1.5,
  cat.col = c("red3", "forestgreen", "deepskyblue3"),
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  scaled = FALSE,
  print.mode = c("raw"),
  sigdigs = 2,
  margin = 0.075)


pdf("nChIP_SP2_rep1_50uMPDS.pdf", width = 20/2.54, height = 20/2.54)
g <- grid.draw(venn.plot)
dev.off()
```


## Deeptools analysis

### Clustering of replicates
Assess similarity of datasets in X-ChIP DMSO consensus. Prepare a number of alternative peak sets for DBA: Rep1&2 (in PDS&DMSO or DMSO); Rep1 (in PDS & DMSO)
```
cd analysis
mkdir deeptools

cat ../macs2_output/*SP2_DMSO*narrowPeak | sortBed -i - | mergeBed -i - > K562_SP2_DMSO_SLX-18311-2_union.merged.bed
mergeBed -d 500 -i  K562_SP2_DMSO_SLX-18311-2_union.merged.bed > K562_SP2_DMSO_SLX-18311-2_union.merged_500.bed

cat ../macs2_output/*SP2*narrowPeak | sortBed -i - | mergeBed -i - > K562_SP2_SLX-18311-2_union.merged.bed
mergeBed -d 500 -i  K562_SP2_SLX-18311-2_union.merged.bed > K562_SP2_SLX-18311-2_union.merged_500.bed

#remove blacklisted peaks:
intersectBed -v -a K562_SP2_SLX-18311-2_union.merged_500.bed -b /scratchb/sblab/spiege01/ENCODE_K562/reference_data/hg19/hg19.blacklist_merge1000nt.bed > K562_SP2_SLX-18311-2_union.merged_500_BL.bed

# SP2 REP1 only 
cat ../macs2_output/*rep1*SP2*narrowPeak | sortBed -i - | mergeBed -i - > K562_SP2_REP1_SLX-18311-2_union.merged.bed
mergeBed -d 500 -i  K562_SP2_REP1_SLX-18311-2_union.merged.bed > K562_SP2_REP1_SLX-18311-2_union.merged_500.bed
intersectBed -v -a K562_SP2_REP1_SLX-18311-2_union.merged_500.bed -b /scratchb/sblab/spiege01/ENCODE_K562/reference_data/hg19/hg19.blacklist_merge1000nt.bed > K562_SP2_REP1_SLX-18311-2_union.merged_500_BL.bed






cd ../bigWig

sbatch -o ../analysis/deeptools/multiBigWigSum.%j.log -J MultBWSum --mem=20000 \
--wrap "multiBigwigSummary BED-file \
-b \
K562_rep1_A1_SP2_50uMPDS_i710-i502_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep1_A2_SP2_50uMPDS_i711-i503_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep1_A3_SP2_50uMPDS_i712-i504_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep1_B1_SP2_DMSO_i701-i505_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep1_B2_SP2_DMSO_i702-i506_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep1_B3_SP2_DMSO_i703-i507_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep2_A1_SP2_50uMPDS_i701-i501_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep2_A2_SP2_50uMPDS_i702-i502_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep2_A3_SP2_50uMPDS_i703-i503_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep2_A4_SP2_50uMPDS_i704-i504_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep2_F1_SP2_DMSO_i705-i505_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep2_F2_SP2_DMSO_i706-i506_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep2_F3_SP2_DMSO_i707-i507_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep2_F4_SP2_DMSO_i708-i508_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
-out ../analysis/deeptools/K562_SP2_SLX-18311-2_MultiBWSummary_per_50bin2.npz \
--binSize 50 \
--BED ../analysis/ K562_SP2_SLX-18311-2_union.merged_500_BL.bed \
--numberOfProcessors max \
--outRawCounts ../analysis/deeptools/K562_SP2_SLX-18311-2_MultiBWSummary_per_50bin2.tab"

cd ../analysis/deeptools


sbatch -o PlotCorrelation.%j.log -J PCA --mem=1000 \
--wrap "plotCorrelation \
-in K562_SP2_SLX-18311-2_MultiBWSummary_per_50bin2.npz \
--whatToPlot heatmap \
--corMethod spearman \
--plotTitle SP2_nativeChIP_rep1-2_PDS_bed_spearman \
--outFileCorMatrix K562_SP2_SLX-18311-2_MultiBWSummary_per_50bin2.tab \
--plotNumbers \
-o K562_SP2_SLX-18311-2_MultiBWSummary_per_50bin2_spearman1.pdf"
```
In rep1 two PDS replicates fairly similar to  DMSO, but still cluster together, while 1 other rep1_A2 is fairly different. In rep2 all the PDS are very different from DMSO. 
**Not sure if low peak number represents 'failed' ChIP or actual succesful displacement...** However, spearman correlation suggests, that rep2_PDS is not realiable. Need another sequencing run.




### Plot Coverage
```
cd ../../merged_bam/nodup

sbatch -o ../../analysis/deeptools/plotCoverage.%j.log -J PlotCov --mem=20000 \
--wrap "plotCoverage -b \
K562_rep1_A1_SP2_50uMPDS_i710-i502_SLX-18311-2.hg19.nodup.bam \
K562_rep1_A2_SP2_50uMPDS_i711-i503_SLX-18311-2.hg19.nodup.bam \
K562_rep1_A3_SP2_50uMPDS_i712-i504_SLX-18311-2.hg19.nodup.bam \
K562_rep1_B1_SP2_DMSO_i701-i505_SLX-18311-2.hg19.nodup.bam \
K562_rep1_B2_SP2_DMSO_i702-i506_SLX-18311-2.hg19.nodup.bam \
K562_rep1_B3_SP2_DMSO_i703-i507_SLX-18311-2.hg19.nodup.bam \
K562_rep2_A1_SP2_50uMPDS_i701-i501_SLX-18311-2.hg19.nodup.bam \
K562_rep2_A2_SP2_50uMPDS_i702-i502_SLX-18311-2.hg19.nodup.bam \
K562_rep2_A3_SP2_50uMPDS_i703-i503_SLX-18311-2.hg19.nodup.bam \
K562_rep2_A4_SP2_50uMPDS_i704-i504_SLX-18311-2.hg19.nodup.bam \
K562_rep2_F1_SP2_DMSO_i705-i505_SLX-18311-2.hg19.nodup.bam \
K562_rep2_F2_SP2_DMSO_i706-i506_SLX-18311-2.hg19.nodup.bam \
K562_rep2_F3_SP2_DMSO_i707-i507_SLX-18311-2.hg19.nodup.bam \
K562_rep2_F4_SP2_DMSO_i708-i508_SLX-18311-2.hg19.nodup.bam \
-o ../../analysis/deeptools/K562_SP2_SLX-18311-2_coverage.png \
--skipZeros \
-n 1000000 \
-p max"  
```
Reasonable coverage.


### Assess ChIP strength - plotFingerprint

```
sbatch -o ../../analysis/deeptools/plotCoverage.%j.log -J PlotCov --mem=20000 \
--wrap "plotFingerprint -b \
nohup plotFingerprint -b \
K562_rep1_A1_SP2_50uMPDS_i710-i502_SLX-18311-2.hg19.nodup.bam \
K562_rep1_A2_SP2_50uMPDS_i711-i503_SLX-18311-2.hg19.nodup.bam \
K562_rep1_A3_SP2_50uMPDS_i712-i504_SLX-18311-2.hg19.nodup.bam \
K562_rep1_B1_SP2_DMSO_i701-i505_SLX-18311-2.hg19.nodup.bam \
K562_rep1_B2_SP2_DMSO_i702-i506_SLX-18311-2.hg19.nodup.bam \
K562_rep1_B3_SP2_DMSO_i703-i507_SLX-18311-2.hg19.nodup.bam \
K562_rep2_A1_SP2_50uMPDS_i701-i501_SLX-18311-2.hg19.nodup.bam \
K562_rep2_A2_SP2_50uMPDS_i702-i502_SLX-18311-2.hg19.nodup.bam \
K562_rep2_A3_SP2_50uMPDS_i703-i503_SLX-18311-2.hg19.nodup.bam \
K562_rep2_A4_SP2_50uMPDS_i704-i504_SLX-18311-2.hg19.nodup.bam \
K562_rep2_F1_SP2_DMSO_i705-i505_SLX-18311-2.hg19.nodup.bam \
K562_rep2_F2_SP2_DMSO_i706-i506_SLX-18311-2.hg19.nodup.bam \
K562_rep2_F3_SP2_DMSO_i707-i507_SLX-18311-2.hg19.nodup.bam \
K562_rep2_F4_SP2_DMSO_i708-i508_SLX-18311-2.hg19.nodup.bam \
K562_rep1_B_input_Input_DMSO_i710-i506_SLX-18311-2.hg19.nodup.bam \
K562_rep2_F_input_Input_DMSO_i709-i501_SLX-18311-2.hg19.nodup.bam \
-plot ../../analysis/deeptools/K562_SP2_SLX-18311-2_fingerprint.png \
-p max" 

```


# Compare native ChIP SP2 DMSO vs 50 µM PDS
Visual inspection and deeptools QC suggests that the experiment may have worked excellently. Already quite clear that changes in  rep1 are not nearly as pronounced. Will need a 3rd replicate to assess, what is 'real'.



## DBA analysis
### prepare peaks for DBA
Use the union of native SP2 ChIP peaks. Try both, rep1&rep2 together and individually.


## Peak coverage calculation (for DBA) (on cluster)
Run DBA in both union of peaks for all conditions (IN_union) or in previous consensus peaks (in_SLX16706_7)

``` bash
cd analysis/
mkdir DBA
mkdir DBA/coverage_peaks


cd ../../merged_bam/nodup/

for f in *SP2*nodup.bam ; do sbatch --mem 50G -o $f.out -e $f.err -J coverage --wrap "bedtools coverage -a ../../analysis/K562_SP2_SLX-18311-2_union.merged_500.bed  -b $f -counts > ../../analysis/DBA/coverage_peaks/In_union_nativePeaks.${f%%.nodup.bam}.bedgraph"; done


cd ../../analysis/DBA/coverage_peaks/

for f in In_union*bedgraph
do
  cut -f 4 $f > tmp.$f
done
paste ../../K562_SP2_SLX-18311-2_union.merged_500.bed tmp.In_union* > In_union_NativePeaks_SP2_SLX-18311-2.bedgraph.tab 

# added headers manually using excel
```
Save the binding affinity matrix to sblab-server for analysis in R.

#### Rerun Counst removing blacklisted files
**Update**: Seemed like DBA was influenced by some very high occupancy reagions. Rerun with blacklisted bed file. In addition, test rep1 only as rep2_PDS doesn't show any peaks -> not sure if the experiment worked extremely well or not at all.



``` bash
cd analysis/
mkdir DBA
mkdir DBA/coverage_peaks


cd ../../merged_bam/nodup/

for f in *SP2*nodup.bam ; do sbatch --mem 50G -o $f.out -e $f.err -J coverage --wrap "bedtools coverage -a ../../analysis/K562_SP2_SLX-18311-2_union.merged_500_BL.bed -b $f -counts > ../../analysis/DBA/coverage_peaks/In_union_nativePeaks.BL.${f%%.nodup.bam}.bedgraph"; done


cd ../../analysis/DBA/coverage_peaks/

for f in In_union*BL*_i*bedgraph
do
  cut -f 4 $f > tmp.$f
done
paste ../../K562_SP2_SLX-18311-2_union.merged_500_BL.bed tmp.In_union* > In_union_BL_NativePeaks_SP2_SLX-18311-2.bedgraph.tab 




printf Chr'\t'Start'\t'End'\t' > testheader.tab

for f in In_union*BL*_i*bedgraph
do
f=${f%%_SLX*} # strip everything after including _SLX
f=${f##*K562_} # strip everything until including K562_
printf $f'\t'
done >> testheader.tab
printf '\n' >> testheader.tab
cat testheader.tab In_union_BL_NativePeaks_SP2_SLX-18311-2.bedgraph.tab   > In_union_BL_NativePeaks_SP2_SLX-18311-2.bedgraph.head.tab 


rm tmp.*
rm In_union_BL_NativePeaks_SP2_SLX-18311-2.bedgraph.tab 
rm testheader.tab
```


#### Rep1 only
``` bash

cd ../../merged_bam/nodup/

for f in *rep1*SP2*_i*nodup.bam ; do sbatch --mem 20G -o $f.out -e $f.err -J coverage --wrap "bedtools coverage -a ../../analysis/K562_SP2_REP1_SLX-18311-2_union.merged_500_BL.bed -b $f -counts > ../../analysis/DBA/coverage_peaks/In_union_REP1_nativePeaks.BL.${f%%.nodup.bam}.bedgraph"; done


cd ../../analysis/DBA/coverage_peaks/

rm tmp.*
for f in In_union_REP1*bedgraph
do
  cut -f 4 $f > tmp.$f
done
paste ../../K562_SP2_REP1_SLX-18311-2_union.merged_500_BL.bed  tmp.In_union* > In_union_REP1_BL_NativePeaks_SP2_SLX-18311-2.bedgraph.tab 




printf Chr'\t'Start'\t'End'\t' > testheader.tab

for f in In_union_REP1*_i*bedgraph
do
f=${f%%_SLX*} # strip everything after including _SLX
f=${f##*K562_} # strip everything until including K562_
printf $f'\t'
done >> testheader.tab
printf '\n' >> testheader.tab
cat testheader.tab In_union_REP1_BL_NativePeaks_SP2_SLX-18311-2.bedgraph.tab > In_union_REP1_BL_NativePeaks_SP2_SLX-18311-2.bedgraph.head.tab 


rm tmp.*
rm In_union_REP1_BL_NativePeaks_SP2_SLX-18311-2.bedgraph.tab 
rm testheader.tab
```








## Count reads (for libsize in diff. binding analysis, which is needed to normalize depth)
See above, read counts per bam file. (12801144, 17656781, 17174222, 10826759, 16665710, 12125157, 28670293, 32135207, 23493633, 34671018, 22382896, 9633100, 26376315,17676458)

## R differential binding
Copy folder to sblabserver and then run Rstudio on server. 

``` R
library(data.table)
library(edgeR)
library(reshape2)
setwd('/Users/js2260/mnt/sblabserver/spiege01/projects/20171123_ENCODE_Analysis/NativeChIP/SLX-18311-2_nChIP_FUS_SP2_vsPDS_rep1-2/analysis/DBA')


# ======================
# === CONTROL Parameters
# ======================

savetables <- 1
saveMAplot <- 1
#define all groups of experiments that you want to have compared to one another:
all_exp_groups <- c('DMSO', '50uMPDS')
relevant_replicates <- c('rep1', 'rep2')




# ====================
# == Input files
# ====================

binding_affinity_matrix <- read.table('In_union_BL_NativePeaks_SP2_SLX-18311-2.bedgraph.head.tab',stringsAsFactors=F,header=T,sep='\t', fill=T)
binding_affinity_matrix <- read.table('In_union_BL_NativePeaks_SP2_SLX-18311-2.bedgraph.head.tab',stringsAsFactors=F,header=T,sep='\t', fill=T)[, 1:17]
colnames(binding_affinity_matrix)
# [1] "Chr"                                "Start"                              "End"                                "K562_rep1_A1_SP2_50uMPDS_i710.i502"
# [5] "K562_rep1_A2_SP2_50uMPDS_i711.i503" "K562_rep1_A3_SP2_50uMPDS_i712.i504" "K562_rep1_B1_SP2_DMSO_i701.i505"    "K562_rep1_B2_SP2_DMSO_i702.i506"   
# [9] "K562_rep1_B3_SP2_DMSO_i703.i507"    "K562_rep2_A1_SP2_50uMPDS_i701.i501" "K562_rep2_A2_SP2_50uMPDS_i702.i502" "K562_rep2_A3_SP2_50uMPDS_i703.i503"
# [13] "K562_rep2_A4_SP2_50uMPDS_i704.i504" "K562_rep2_F1_SP2_DMSO_i705.i505"    "K562_rep2_F2_SP2_DMSO_i706.i506"    "K562_rep2_F3_SP2_DMSO_i707.i507"   
# [17] "K562_rep2_F4_SP2_DMSO_i708.i508"    

libSize_all <- c(12801144, 17656781, 17174222, 10826759, 16665710, 12125157, 28670293, 32135207, 23493633, 34671018, 22382896, 9633100, 26376315,17676458)
# Lib size (reads per bam file, see above)


# ================================================================
# = List of all possible combinations that will be analyzed in DBA
# ===============================================================

all_combinations <- combn(all_exp_groups, 2) #matrix with all possible combinations


# ================================================================
# = differential binding analyis with edgeR
# ===============================================================


# loop that performs DBA for all selected combinations
if(saveMAplot)
  {
pdf(file='SLX18311-2_DBA_SP2_native-ChIP_PDS-displacement_rep1&2.pdf', 8, 8)
for (i in 1:ncol(all_combinations))
{
 
exp_group1 <- all_combinations[1,i]
exp_group2 <- all_combinations[2,i]

y <- binding_affinity_matrix  
row.names(y) <- paste(binding_affinity_matrix$Chr,':',binding_affinity_matrix$Start,'-', binding_affinity_matrix$End,sep='')
y <- y[,4:length(y)] #transfer loci into row names

relevant_columns <- c(grep(exp_group1, colnames(y)), grep(exp_group2, colnames(y)))
y <- y[relevant_columns] #keep relevant columns
libSize <- libSize_all[relevant_columns] # deterimne relevant libery sizes

cnt <- y 
cnt$locus <- row.names(y) 
number_replicates_exp_group1 <- length(grep(exp_group1, colnames(y)))
number_replicates_exp_group2 <- length(grep(exp_group2, colnames(y)))
group <- factor(c(rep(exp_group1, number_replicates_exp_group1), rep(exp_group2, number_replicates_exp_group2)))

y<- DGEList(counts=y, group=group)
stopifnot(rownames(y$samples) == names(libSize))
y$samples$lib.size<- libSize
y<- calcNormFactors(y, method= 'none')
y<- estimateDisp(y)
y<- estimateCommonDisp(y)
y<- estimateTagwiseDisp(y)

et<- exactTest(y, pair= as.vector(unique(y$samples$group)))   #read in as vector; levels would order alphabetically and might flip samples

detable<- data.frame(topTags(et, n= Inf)$table)
detable$locus<- rownames(detable)
detable<- data.table(detable)
detable<- merge(detable, cnt, by= 'locus')

#  MA plot
pal<- colorRampPalette(c("white", "lightblue", "yellow", "red"), space = "Lab")
par(las= 1, mgp= c(1.75, 0.5, 0), bty= 'l', mar= c(3, 3, 3, 0.5))
smoothScatter(x= detable$logCPM, y= detable$logFC, xlab= 'logCPM', ylab= 'logFC',
              main= paste("Differential binding", exp_group1, 'vs', exp_group2), colramp= pal, col= 'blue')
lines(loess.smooth(x= detable$logCPM, y= detable$logFC, span= 0.1), lwd= 2, col= 'grey60')
abline(h= 0, col= 'grey30')
points(x= detable$logCPM, y= detable$logFC, col= ifelse(detable$FDR < 0.05, '#FF000080', 'transparent'), cex= 0.5, pch= '.')
mtext(side= 3, line= -1.2, text= sprintf('FDR < 0.05: %s', nrow(detable[FDR < 0.05 & logFC > 0])), adj= 1)
mtext(side= 1, line= -1.2, text= sprintf('FDR < 0.05: %s', nrow(detable[FDR < 0.05 & logFC < 0])), adj= 1)
grid(col= 'grey50')

if(savetables)
{
write.table(detable, paste("Differential binding", exp_group1, 'vs', exp_group2, '.tab'), row.names= FALSE, col.names= TRUE, sep= '\t', quote= FALSE)
}


grid()
}
dev.off()
}
 


#### Rep2 only

  pdf(file='SLX18311-2_DBA_SP2_native-ChIP_PDS-displacement_rep2.pdf', 8, 8)
  for (i in 1:ncol(all_combinations))
  {
    
    exp_group1 <- all_combinations[1,i]
    exp_group2 <- all_combinations[2,i]
    y <- binding_affinity_matrix  
    row.names(y) <- paste(binding_affinity_matrix$Chr,':',binding_affinity_matrix$Start,'-', binding_affinity_matrix$End,sep='')
    y <- y[,4:length(y)] #transfer loci into row names
    
    relevant_columns <- c(grep(glob2rx('*rep2*DMSO*'), colnames(y)), grep(glob2rx('*rep2*PDS*'), colnames(y)))
    y <- y[relevant_columns] #keep relevant columns
    libSize <- libSize_all[relevant_columns] # deterimne relevant libery sizes
    
    cnt <- y 
    cnt$locus <- row.names(y) 
    number_replicates_exp_group1 <- length(grep(exp_group1, colnames(y)))
    number_replicates_exp_group2 <- length(grep(exp_group2, colnames(y)))
    group <- factor(c(rep(exp_group1, number_replicates_exp_group1), rep(exp_group2, number_replicates_exp_group2)))
    
    y<- DGEList(counts=y, group=group)
    stopifnot(rownames(y$samples) == names(libSize))
    y$samples$lib.size<- libSize
    y<- calcNormFactors(y, method= 'none')
    y<- estimateDisp(y)
    y<- estimateCommonDisp(y)
    y<- estimateTagwiseDisp(y)
    
    et<- exactTest(y, pair= as.vector(unique(y$samples$group)))   #read in as vector; levels would order alphabetically and might flip samples
    
    detable<- data.frame(topTags(et, n= Inf)$table)
    detable$locus<- rownames(detable)
    detable<- data.table(detable)
    detable<- merge(detable, cnt, by= 'locus')
    
    #  MA plot
    pal<- colorRampPalette(c("white", "lightblue", "yellow", "red"), space = "Lab")
    par(las= 1, mgp= c(1.75, 0.5, 0), bty= 'l', mar= c(3, 3, 3, 0.5))
    smoothScatter(x= detable$logCPM, y= detable$logFC, xlab= 'logCPM', ylab= 'logFC',
                  main= paste("Differential binding", exp_group1, 'vs', exp_group2), colramp= pal, col= 'blue')
    lines(loess.smooth(x= detable$logCPM, y= detable$logFC, span= 0.1), lwd= 2, col= 'grey60')
    abline(h= 0, col= 'grey30')
    points(x= detable$logCPM, y= detable$logFC, col= ifelse(detable$FDR < 0.05, '#FF000080', 'transparent'), cex= 0.5, pch= '.')
    mtext(side= 3, line= -1.2, text= sprintf('FDR < 0.05: %s', nrow(detable[FDR < 0.05 & logFC > 0])), adj= 1)
    mtext(side= 1, line= -1.2, text= sprintf('FDR < 0.05: %s', nrow(detable[FDR < 0.05 & logFC < 0])), adj= 1)
    grid(col= 'grey50')
    
    if(savetables)
    {
      write.table(detable, paste("Differential binding", exp_group1, 'vs', exp_group2, '.tab'), row.names= FALSE, col.names= TRUE, sep= '\t', quote= FALSE)
    }
    
    
    grid()
  }
  dev.off()

  #### Rep1 only
  # loop that performs DBA for all selected combinations
  pdf(file='SLX18311-2_DBA_SP2_native-ChIP_PDS-displacement_rep1.pdf', 8, 8)
  for (i in 1:ncol(all_combinations))
  {
    
    exp_group1 <- all_combinations[1,i]
    exp_group2 <- all_combinations[2,i]
    y <- binding_affinity_matrix  
    row.names(y) <- paste(binding_affinity_matrix$Chr,':',binding_affinity_matrix$Start,'-', binding_affinity_matrix$End,sep='')
    y <- y[,4:length(y)] #transfer loci into row names
    
    relevant_columns <- c(grep(glob2rx('*rep1*DMSO*'), colnames(y)), grep(glob2rx('*rep1*PDS*'), colnames(y)))
    y <- y[relevant_columns] #keep relevant columns
    libSize <- libSize_all[relevant_columns] # deterimne relevant libery sizes
    
    cnt <- y 
    cnt$locus <- row.names(y) 
    number_replicates_exp_group1 <- length(grep(exp_group1, colnames(y)))
    number_replicates_exp_group2 <- length(grep(exp_group2, colnames(y)))
    group <- factor(c(rep(exp_group1, number_replicates_exp_group1), rep(exp_group2, number_replicates_exp_group2)))
    
    y<- DGEList(counts=y, group=group)
    stopifnot(rownames(y$samples) == names(libSize))
    y$samples$lib.size<- libSize
    y<- calcNormFactors(y, method= 'none')
    y<- estimateDisp(y)
    y<- estimateCommonDisp(y)
    y<- estimateTagwiseDisp(y)
    
    et<- exactTest(y, pair= as.vector(unique(y$samples$group)))   #read in as vector; levels would order alphabetically and might flip samples
    
    detable<- data.frame(topTags(et, n= Inf)$table)
    detable$locus<- rownames(detable)
    detable<- data.table(detable)
    detable<- merge(detable, cnt, by= 'locus')
    
    #  MA plot
    pal<- colorRampPalette(c("white", "lightblue", "yellow", "red"), space = "Lab")
    par(las= 1, mgp= c(1.75, 0.5, 0), bty= 'l', mar= c(3, 3, 3, 0.5))
    smoothScatter(x= detable$logCPM, y= detable$logFC, xlab= 'logCPM', ylab= 'logFC',
                  main= paste("Differential binding", exp_group1, 'vs', exp_group2), colramp= pal, col= 'blue')
    lines(loess.smooth(x= detable$logCPM, y= detable$logFC, span= 0.1), lwd= 2, col= 'grey60')
    abline(h= 0, col= 'grey30')
    points(x= detable$logCPM, y= detable$logFC, col= ifelse(detable$FDR < 0.05, '#FF000080', 'transparent'), cex= 0.5, pch= '.')
    mtext(side= 3, line= -1.2, text= sprintf('FDR < 0.05: %s', nrow(detable[FDR < 0.05 & logFC > 0])), adj= 1)
    mtext(side= 1, line= -1.2, text= sprintf('FDR < 0.05: %s', nrow(detable[FDR < 0.05 & logFC < 0])), adj= 1)
    grid(col= 'grey50')
    
    if(savetables)
    {
      write.table(detable, paste("Differential binding", exp_group1, 'vs', exp_group2, '.tab'), row.names= FALSE, col.names= TRUE, sep= '\t', quote= FALSE)
    }
    
    
    grid()
  }
  dev.off()

  

```

The DBA results for [rep1](Files/SLX-18311/SLX18311-2_DBA_SP2_native-ChIP_PDS-displacement_rep1.pdf) [rep2](Files/SLX-18311/SLX18311-2_DBA_SP2_native-ChIP_PDS-displacement_rep2.pdf) both show significant changes upon PDS treatment. I will need a 3rd replicate to assess, whether rep2 worked particularly well or PDS treatment killed everything in this experiment.


### Check Signal intensity in SP2 peak reagions:


```
cd analysis

SP2_DMSO_PDS+-BG4_rep1.mat.gz
K562_SP2_DMSO_SLX-18311-2_Mult_4out7_NO_BG4.bed
K562_SP2_DMSO_SLX-18311-2_Mult_4out7_WITH_BG4.bed


intersectBed -a K562_SP2_DMSO_SLX-18311-2_Mult_4out7.bed  -b /scratchb/sblab/spiege01/ENCODE_K562/reference_data/Quadruplex/20180108_K562_async_rep1-3.mult.5of8.bed > signal/K562_SP2_DMSO_SLX-18311-2_Mult_4out7_WITH_BG4.bed

intersectBed -v -a K562_SP2_DMSO_SLX-18311-2_Mult_4out7.bed  -b /scratchb/sblab/spiege01/ENCODE_K562/reference_data/Quadruplex/20180108_K562_async_rep1-3.mult.5of8.bed > signal/K562_SP2_DMSO_SLX-18311-2_Mult_4out7_NO_BG4.bed


cd ../../bigWig

sbatch -o Plotprofile.%j.log -J DT_Profile --mem=20000 \
--wrap "computeMatrix reference-point \
--referencePoint center \
-b 2000 -a 2000 \
-S \
K562_rep1_A1_SP2_50uMPDS_i710-i502_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep1_A2_SP2_50uMPDS_i711-i503_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep1_A3_SP2_50uMPDS_i712-i504_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep1_B1_SP2_DMSO_i701-i505_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep1_B2_SP2_DMSO_i702-i506_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep1_B3_SP2_DMSO_i703-i507_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep2_A1_SP2_50uMPDS_i701-i501_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep2_A2_SP2_50uMPDS_i702-i502_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep2_A3_SP2_50uMPDS_i703-i503_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep2_A4_SP2_50uMPDS_i704-i504_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep2_F1_SP2_DMSO_i705-i505_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep2_F2_SP2_DMSO_i706-i506_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep2_F3_SP2_DMSO_i707-i507_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep2_F4_SP2_DMSO_i708-i508_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
-R \
../../analysis/signal/K562_SP2_DMSO_SLX-18311-2_Mult_4out7_NO_BG4.bed \
../../analysis/signal/K562_SP2_DMSO_SLX-18311-2_Mult_4out7_WITH_BG4.bed \
--skipZeros \
-o SP2_DMSO_PDS+-BG4_rep1-2.mat.gz"


sbatch -o Plotprofile.%j.log -J DT_Profile --mem=1000 \
--wrap "plotProfile -m SP2_DMSO_PDS+-BG4_rep1-2.mat.gz \
-out SP2_DMSO_PDS+-BG4_rep1-2.pdf \
--outFileNameData SP2_DMSO_PDS+-BG4_rep1-2.tab \
--dpi 500 \
--plotHeight 15 \
--plotWidth 20 \
--perGroup \
--refPointLabel Center \
--yMin 0 \
--samplesLabel PDS_rep1_1 PDS_rep1_2 PDS_rep1_3 DMSO_rep1_1 DMSO_rep1_2 DMSO_rep1_3 PDS_rep2_1 PDS_rep2_2 PDS_rep3_3 PDS_rep3_4 DMSO_rep2_1 DMSO_rep2_2 DMSO_rep2_3 DMSO_rep2_4 \
--plotTitle Signalstrength_union_+-BG4 \
--numPlotsPerRow 2"




### Also look at rep1 only

sbatch -o Plotprofile.%j.log -J DT_Profile --mem=20000 \
--wrap "computeMatrix reference-point \
--referencePoint center \
-b 2000 -a 2000 \
-S \
K562_rep1_A1_SP2_50uMPDS_i710-i502_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep1_A2_SP2_50uMPDS_i711-i503_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep1_A3_SP2_50uMPDS_i712-i504_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep1_B1_SP2_DMSO_i701-i505_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep1_B2_SP2_DMSO_i702-i506_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
K562_rep1_B3_SP2_DMSO_i703-i507_SLX-18311-2.hg19.bs50.bl.RPKM.bw \
-R \
../../analysis/signal/K562_SP2_DMSO_SLX-18311-2_Mult_4out7_NO_BG4.bed \
../../analysis/signal/K562_SP2_DMSO_SLX-18311-2_Mult_4out7_WITH_BG4.bed \
--skipZeros \
-o SP2_DMSO_PDS+-BG4_rep1.mat.gz"


sbatch -o Plotprofile.%j.log -J DT_Profile --mem=1000 \
--wrap "plotProfile -m SP2_DMSO_PDS+-BG4_rep1.mat.gz \
-out SP2_DMSO_PDS+-BG4_rep1.pdf \
--outFileNameData SP2_DMSO_PDS+-BG4_rep1.tab \
--dpi 500 \
--plotHeight 15 \
--plotWidth 20 \
--perGroup \
--refPointLabel Center \
--yMin 0 \
--samplesLabel PDS_rep1_1 PDS_rep1_2 PDS_rep1_3 DMSO_rep1_1 DMSO_rep1_2 DMSO_rep1_3 PDS_rep2_1 PDS_rep2_2 PDS_rep3_3 PDS_rep3_4 DMSO_rep2_1 DMSO_rep2_2 DMSO_rep2_3 DMSO_rep2_4 \
--plotTitle Signalstrength_union_+-BG4_rep1 \
--numPlotsPerRow 2"

```

Looking at [all replicates](Files/SLX-18311/SP2_DMSO_PDS%2B-BG4.pdf) and  [rep1](Files/SLX-18311/SP2_DMSO_PDS%2B-BG4_rep1.pdf) there seems barely any difference between +/-BG4 sites. On the contrary, there is a massive difference between 4/6 PDS samples. Still, another biol repl. will have to confirm that this is real and not just a drop in ChIP quality (e.g. testing FOXA1)

The low difference between +/-BG4 may result from the fact that -BG4 regions still are G4s but were not called peaks in the original G4 ChIP analysis.

Try to compare G4 ChIP-signal in these regions:
```
cd ../analysis/signal

sbatch -o Plotprofile.%j.log -J DT_Profile --mem=20000 \
--wrap "computeMatrix reference-point \
--referencePoint center \
-b 2000 -a 2000 \
-S \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio1_a.bs50.bl.RPKM.bw \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio1_c.bs50.bl.RPKM.bw \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio2_a.bs50.bl.RPKM.bw \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio2_b.bs50.bl.RPKM.bw \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio2_c.bs50.bl.RPKM.bw \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio3_a.bs50.bl.RPKM.bw \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio3_b.bs50.bl.RPKM.bw \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio3_c.bs50.bl.RPKM.bw \
-R \
K562_SP2_DMSO_SLX-18311-2_Mult_4out7_WITH_BG4.bed \
K562_SP2_DMSO_SLX-18311-2_Mult_4out7_NO_BG4.bed \
--skipZeros \
-o BG4-ChIP_in_SP2_DMSO4out7+-BG4.mat.gz"


sbatch -o Plotprofile.%j.log -J DT_Profile --mem=1000 \
--wrap "plotProfile -m BG4-ChIP_in_SP2_DMSO4out7+-BG4.mat.gz \
-out BG4-ChIP_in_SP2_DMSO4out7+-BG4.pdf \
--outFileNameData BG4-ChIP_in_SP2_DMSO4out7+-BG4.tab \
--dpi 500 \
--plotHeight 15 \
--plotWidth 20 \
--perGroup \
--refPointLabel Center \
--yMin 0 \
--samplesLabel BG4_rep1a BG4_rep1c BG4_rep2a BG4_rep2b BG4_rep2c BG4_rep3a BG4_rep3b BG4_rep3c \
--plotTitle Signalstrength_BG4-ChIP_in_SP2-nChIP+-G4 \
--numPlotsPerRow 2"

```
The [G4 ChIP intensity in SP2 nChIP peaks +-BG4](Files/SLX-18311/BG4-ChIP_in_SP2_DMSO4out7%2B-BG4.pdf) reveals a ~2-fold stronger BG4 occupancy in the '+BG4' regions, however, also substantial signal in the '-BG4' peaks suggesting that these might also be G4 sites. (Lowever BG4 occupancy could be explained by either less stable G4s or more masking by other proteins).

Overall, the results looks very promising. PDS seems to promote  a genome-wide drop of signal. Now I need to show that a regular TF (FOXA1) is not affected by 60 uM PDS (and that iPDS does not cause the same effect).



### Rerun DBA
Clean up previous analysis:
- perform DBA for rep1&rep2 in respective union of peaks using a design matrix and general linear model
- perform DBA for rep1 ONLY ind REP1 peaks! (might be a conservative approach of presenting the data, trying not to overestimate effects from rep2)
- generate top and bottom candidates

```
library(data.table)
library(edgeR)
library(reshape2)
setwd('/Users/js2260/mnt/sblabserver/spiege01/projects/20171123_ENCODE_Analysis/NativeChIP/SLX-18311-2_nChIP_FUS_SP2_vsPDS_rep1-2/analysis/DBA')


# ======================
# === CONTROL Parameters
# ======================

savetables <- 1
saveMAplot <- 1
#define all groups of experiments that you want to have compared to one another:

# ====================
# == Input files
# ====================

binding_affinity_matrix <- read.table('In_union_BL_NativePeaks_SP2_SLX-18311-2.bedgraph.head.tab',stringsAsFactors=F,header=T,sep='\t', fill=T)
binding_affinity_matrix <- read.table('In_union_BL_NativePeaks_SP2_SLX-18311-2.bedgraph.head.tab',stringsAsFactors=F,header=T,sep='\t', fill=T)[, 1:17]
colnames(binding_affinity_matrix)
# [1] "Chr"                                "Start"                              "End"                                "K562_rep1_A1_SP2_50uMPDS_i710.i502"
# [5] "K562_rep1_A2_SP2_50uMPDS_i711.i503" "K562_rep1_A3_SP2_50uMPDS_i712.i504" "K562_rep1_B1_SP2_DMSO_i701.i505"    "K562_rep1_B2_SP2_DMSO_i702.i506"   
# [9] "K562_rep1_B3_SP2_DMSO_i703.i507"    "K562_rep2_A1_SP2_50uMPDS_i701.i501" "K562_rep2_A2_SP2_50uMPDS_i702.i502" "K562_rep2_A3_SP2_50uMPDS_i703.i503"
# [13] "K562_rep2_A4_SP2_50uMPDS_i704.i504" "K562_rep2_F1_SP2_DMSO_i705.i505"    "K562_rep2_F2_SP2_DMSO_i706.i506"    "K562_rep2_F3_SP2_DMSO_i707.i507"   
# [17] "K562_rep2_F4_SP2_DMSO_i708.i508"    

libSize_all <- c(12801144, 17656781, 17174222, 10826759, 16665710, 12125157, 28670293, 32135207, 23493633, 34671018, 22382896, 9633100, 26376315,17676458)
# Lib size (reads per bam file, see above)




# ================================================================
# = differential binding analyis with edgeR
# ===============================================================


###############
# Perform GLM analysis on rep1 and rep2 in union of peaks (no BL)
#####

  y <- binding_affinity_matrix  
  row.names(y) <- paste(binding_affinity_matrix$Chr,':',binding_affinity_matrix$Start,'-', binding_affinity_matrix$End,sep='')
  y <- as.data.frame(y[,4:17]) #transfer loci into row names
  libSize <- libSize_all # deterimne relevant libery sizes
  
  cnt <- binding_affinity_matrix 
  cnt$locus <- paste(binding_affinity_matrix$Chr,':',binding_affinity_matrix$Start,'-', binding_affinity_matrix$End,sep='') 


  # build design matrix
  TreatType <- factor(c(rep('50uMPDS',3), rep('DMSO',3), rep('50uMPDS',4), rep('DMSO',4)))
  rep_day <- factor(c(rep('rep1',6),rep('rep2',8)))
  design_matrix <- model.matrix(~rep_day+TreatType)
  #design_matrix <- model.matrix(~0+TreatType+rep_day)
 
  
  # Calculate normalization factors
  y <- DGEList(counts = y,group= TreatType)

  y$samples$lib.size<- libSize
  y <- calcNormFactors(y, method = "none")
  
  # Estimate dispertion
  y <- estimateDisp(y, design_matrix)
  
  
  
  #GLM
  fit <- glmQLFit(y, design_matrix)
  qlf_PDS_DMSO <- glmQLFTest(fit, contrast = c(0,0,-1))

  detable<- data.frame(topTags(qlf_PDS_DMSO, n= Inf)$table)
  detable$locus<- rownames(detable)
  detable<- data.table(detable)
  detable<- merge(detable, cnt, by= 'locus')
  
  #  MA plot
  pdf('Differential binding_GLM_PDSvsDMSO_rep12.pdf', 4,4)
  pal<- colorRampPalette(c("white", "lightblue", "yellow", "red"), space = "Lab")
  par(las= 1, mgp= c(1.75, 0.5, 0), bty= 'l', mar= c(3, 3, 3, 0.5))
  smoothScatter(x= detable$logCPM, y= detable$logFC, xlab= 'logCPM', ylab= 'logFC',
                main= "Differential binding_DMSOvsPDS", colramp= pal, col= 'blue')
  lines(loess.smooth(x= detable$logCPM, y= detable$logFC, span= 0.1), lwd= 2, col= 'grey60')
  abline(h= 0, col= 'grey30')
  points(x= detable$logCPM, y= detable$logFC, col= ifelse(detable$FDR < 0.05, '#FF000080', 'transparent'), cex= 0.5, pch= '.')
  mtext(side= 3, line= -1.2, text= sprintf('FDR < 0.05: %s', nrow(detable[FDR < 0.05 & logFC > 0])), adj= 1)
  mtext(side= 1, line= -1.2, text= sprintf('FDR < 0.05: %s', nrow(detable[FDR < 0.05 & logFC < 0])), adj= 1)
  grid(col= 'grey50')  
  dev.off()
  if(savetables)
  {
    write.table(detable, "Differential binding_GLM_PDSvsDMSO_rep12.tab", row.names= FALSE, col.names= TRUE, sep= '\t', quote= FALSE)
  }
  
  # generate list of top changing sites and non chaning (cpm minimum >)
  detable <- detable[order(detable$FDR, decreasing = F), ]
  detable <- detable[detable$FDR < 0.05, ]
  top100 <- detable[1:100, c("Chr", "Start", "End")]
  top100_cpm <- detable[detable$logCPM > 0 ,]
  top100_cpm <- top100_cpm[1:100, c("Chr", "Start", "End")]
  test$Start <- top100$locus

  
  
  
  
  
  
  
  
  
  
  
  
  
######################
## REP1 only
######################

  
  binding_affinity_matrix <- read.table('In_union_REP1_BL_NativePeaks_SP2_SLX-18311-2.bedgraph.head.tab',stringsAsFactors=F,header=T,sep='\t', fill=T)
  binding_affinity_matrix <- binding_affinity_matrix[, 1:9]
  colnames(binding_affinity_matrix)
  #[1] "Chr"                           "Start"                         "End"                          
  #[4] "rep1_A1_SP2_50uMPDS_i710.i502" "rep1_A2_SP2_50uMPDS_i711.i503" "rep1_A3_SP2_50uMPDS_i712.i504"
  #[7] "rep1_B1_SP2_DMSO_i701.i505"    "rep1_B2_SP2_DMSO_i702.i506"    "rep1_B3_SP2_DMSO_i703.i507"   
  
  libSize_all <- c(12801144, 17656781, 17174222, 10826759, 16665710, 12125157)
  # Lib size (reads per bam file, see above)
  
  
  
  # ================================================================
  # = differential binding analyis with edgeR
  # ===============================================================
  
  
  ###############
  # Perform GLM analysis on rep1  union of rep PDS/DMSO peaks (no BL)
  #####

  
  y <- binding_affinity_matrix  
  row.names(y) <- paste(binding_affinity_matrix$Chr,':',binding_affinity_matrix$Start,'-', binding_affinity_matrix$End,sep='')
  y <- as.data.frame(y[,4:length(y)]) #transfer loci into row names
  libSize <- libSize_all # deterimne relevant libery sizes
  
  cnt <- binding_affinity_matrix 
  cnt$locus <- paste(binding_affinity_matrix$Chr,':',binding_affinity_matrix$Start,'-', binding_affinity_matrix$End,sep='') 
  
  
  # build design matrix
  TreatType <- factor(c(rep('50uMPDS',3), rep('DMSO',3)))
  design_matrix <- model.matrix(~0+TreatType)
  

  # Calculate normalization factors
  y <- DGEList(counts = y,group= TreatType)
  
  y$samples$lib.size<- libSize
  y <- calcNormFactors(y, method = "none")
  
  # Estimate dispertion
  y <- estimateDisp(y, design_matrix)
  
  #GLM
  my.contrasts <- makeContrasts(PDSvsDMSO = TreatType50uMPDS - TreatTypeDMSO, levels=design_matrix )
  fit <- glmQLFit(y, design_matrix)
  lrt_PDSvsDMSO <- glmLRT(fit, contrast=my.contrasts[, "PDSvsDMSO"])
  
  detable<- data.frame(topTags(lrt_PDSvsDMSO, n= Inf)$table)
  detable$locus<- rownames(detable)
  detable<- data.table(detable)
  detable<- merge(detable, cnt, by= 'locus')
  
  #  MA plot
  pdf('Differential binding_GLM_DMSOvsPDS_REP1.pdf', 4,4)
  pal<- colorRampPalette(c("white", "lightblue", "yellow", "red"), space = "Lab")
  par(las= 1, mgp= c(1.75, 0.5, 0), bty= 'l', mar= c(3, 3, 3, 0.5))
  smoothScatter(x= detable$logCPM, y= detable$logFC, xlab= 'logCPM', ylab= 'logFC',
                main= "Differential binding_DMSOvsPDS", colramp= pal, col= 'blue')
  lines(loess.smooth(x= detable$logCPM, y= detable$logFC, span= 0.1), lwd= 2, col= 'grey60')
  abline(h= 0, col= 'grey30')
  points(x= detable$logCPM, y= detable$logFC, col= ifelse(detable$FDR < 0.05, '#FF000080', 'transparent'), cex= 0.5, pch= '.')
  mtext(side= 3, line= -1.2, text= sprintf('FDR < 0.05: %s', nrow(detable[FDR < 0.05 & logFC > 0])), adj= 1)
  mtext(side= 1, line= -1.2, text= sprintf('FDR < 0.05: %s', nrow(detable[FDR < 0.05 & logFC < 0])), adj= 1)
  grid(col= 'grey50')  
  dev.off()
  if(savetables)
  {
    write.table(detable, "Differential binding_GLM_DMSOvsPDS_REP1.tab", row.names= FALSE, col.names= TRUE, sep= '\t', quote= FALSE)
  }
  
  
  
  
  # generate list of top changing sites and non changing (changing: lowest logFC ; not chaning: loweste absolute value of logFC ) 
  subset <- detable
  top100 <- subset[subset$FDR < 0.05, ]
  top100 <- top100[order(top100$logFC, decreasing = F), ]
  top100 <- top100[1:100, c("Chr", "Start", "End")]
  
  top400 <- subset[subset$FDR < 0.05, ]
  top400 <- top400[order(top400$logFC, decreasing = F), ]
  top400 <- top400[1:400, c("Chr", "Start", "End")]
  
  
  bottom100 <- subset[order(abs(subset$logFC), decreasing = F), ]
  bottom100 <- bottom100[1:100,c("Chr", "Start", "End")]

  bottom400 <- subset[order(abs(subset$logFC), decreasing = F), ]
  bottom400 <- bottom400[1:400,c("Chr", "Start", "End")]
  
  NoChange <- subset[abs(subset$logFC) < 0.5,  c("Chr", "Start", "End")] 
  
  if(savetables)
  {
    write.table(top100, "DBA_GLM_DMSOvsPDS_REP1_top100.bed", row.names= F, col.names= F, sep= '\t', quote= F)
    write.table(top400, "DBA_GLM_DMSOvsPDS_REP1_top400.bed", row.names= F, col.names= F, sep= '\t', quote= F)
    write.table(bottom100, "DBA_GLM_DMSOvsPDS_REP1_bottom100.bed", row.names= F, col.names= F, sep= '\t', quote= F)
    write.table(bottom400, "DBA_GLM_DMSOvsPDS_REP1_bottom400.bed", row.names= F, col.names= F, sep= '\t', quote= F)
    write.table(NoChange, "DBA_GLM_DMSOvsPDS_REP1_NoChange.bed", row.names= F, col.names= F, sep= '\t', quote= F)    
  }
  

  # generate files at regions with higher signal intensity
  subset <- detable[detable$logCPM > 0]
  top100 <- subset[subset$FDR < 0.05, ]
  top100 <- top100[order(top100$logFC, decreasing = F), ]
  top100 <- top100[1:100, c("Chr", "Start", "End")]
  
  top400 <- subset[subset$FDR < 0.05, ]
  top400 <- top400[order(top400$logFC, decreasing = F), ]
  top400 <- top400[1:400, c("Chr", "Start", "End")]
  
  
  bottom100 <- subset[order(abs(subset$logFC), decreasing = F), ]
  bottom100 <- bottom100[1:100,c("Chr", "Start", "End")]
  
  bottom400 <- subset[order(abs(subset$logFC), decreasing = F), ]
  bottom400 <- bottom400[1:400,c("Chr", "Start", "End")]
  
  NoChange <- subset[abs(subset$logFC) < 0.5,  c("Chr", "Start", "End")] 
  
  
  if(savetables)
  {
    write.table(top100, "DBA_GLM_DMSOvsPDS_REP1_top100_cpm.bed", row.names= F, col.names= F, sep= '\t', quote= F)
    write.table(top400, "DBA_GLM_DMSOvsPDS_REP1_top400_cpm.bed", row.names= F, col.names= F, sep= '\t', quote= F)
    write.table(bottom100, "DBA_GLM_DMSOvsPDS_REP1_bottom100_cpm.bed", row.names= F, col.names= F, sep= '\t', quote= F)
    write.table(bottom400, "DBA_GLM_DMSOvsPDS_REP1_bottom400_cpm.bed", row.names= F, col.names= F, sep= '\t', quote= F)
    write.table(NoChange, "DBA_GLM_DMSOvsPDS_REP1_NoChange_cpm.bed", row.names= F, col.names= F, sep= '\t', quote= F)
  }
 
  
  

 
  
  

```


#### Analysis of changing vs constant regions:

Transfer top/botom changing SP2 nChIP sites from sblabserver to cluster. Check for overlap with BG4 and plot BG4 signal strength in respective sites.

Overlap with BG4
```
cd /scratchb/sblab/spiege01/ENCODE_K562/Native_ChIP/SLX-18311_nChIP_SP2_FUS_vsPDS_rep1-2/analysis/DBA/changing_sites

for f in *bed
do
printf $f' '
intersectBed -a $f -b /scratchb/sblab/spiege01/ENCODE_K562/reference_data/Quadruplex/20180108_K562_async_rep1-3.mult.5of8.bed  | wc -l
printf "out of "
wc -l $f 
done


DBA_GLM_DMSOvsPDS_REP1_bottom100.bed 56 out of 100 DBA_GLM_DMSOvsPDS_REP1_bottom100.bed
DBA_GLM_DMSOvsPDS_REP1_bottom100_cpm.bed 67 out of 100 DBA_GLM_DMSOvsPDS_REP1_bottom100_cpm.bed
DBA_GLM_DMSOvsPDS_REP1_bottom400.bed 235 out of 400 DBA_GLM_DMSOvsPDS_REP1_bottom400.bed
DBA_GLM_DMSOvsPDS_REP1_bottom400_cpm.bed 303 out of 400 DBA_GLM_DMSOvsPDS_REP1_bottom400_cpm.bed
DBA_GLM_DMSOvsPDS_REP1_NoChange.bed 66 out of 121 DBA_GLM_DMSOvsPDS_REP1_NoChange.bed
DBA_GLM_DMSOvsPDS_REP1_NoChange_cpm.bed 16 out of 22 DBA_GLM_DMSOvsPDS_REP1_NoChange_cpm.bed
DBA_GLM_DMSOvsPDS_REP1_top100.bed 58 out of 100 DBA_GLM_DMSOvsPDS_REP1_top100.bed
DBA_GLM_DMSOvsPDS_REP1_top100_cpm.bed 89 out of 100 DBA_GLM_DMSOvsPDS_REP1_top100_cpm.bed
DBA_GLM_DMSOvsPDS_REP1_top400.bed 251 out of 400 DBA_GLM_DMSOvsPDS_REP1_top400.bed
DBA_GLM_DMSOvsPDS_REP1_top400_cpm.bed 321 out of 400 DBA_GLM_DMSOvsPDS_REP1_top400_cpm.bed

```

Plot Signal strenght of BG4 at these sites.

```
cd /scratchb/sblab/spiege01/ENCODE_K562/Native_ChIP/SLX-18311_nChIP_SP2_FUS_vsPDS_rep1-2/analysis/DBA/changing_sites

sbatch -o Plotprofile.%j.log -J DT_Profile --mem=20000 \
--wrap "computeMatrix reference-point \
--referencePoint center \
-b 300 -a 300 \
-S \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio1_a.bs50.bl.RPKM.bw \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio1_c.bs50.bl.RPKM.bw \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio2_a.bs50.bl.RPKM.bw \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio2_b.bs50.bl.RPKM.bw \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio2_c.bs50.bl.RPKM.bw \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio3_a.bs50.bl.RPKM.bw \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio3_b.bs50.bl.RPKM.bw \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio3_c.bs50.bl.RPKM.bw \
-R \
DBA_GLM_DMSOvsPDS_REP1_bottom100.bed \
DBA_GLM_DMSOvsPDS_REP1_bottom400.bed \
DBA_GLM_DMSOvsPDS_REP1_top100.bed \
DBA_GLM_DMSOvsPDS_REP1_top400.bed \
--skipZeros \
-o BG4-ChIP_in_SP2_changing_REP1.mat.gz &&
plotProfile -m BG4-ChIP_in_SP2_changing_REP1.mat.gz \
-out BG4-ChIP_in_SP2_changing_REP1.pdf \
--outFileNameData BG4-ChIP_in_SP2_changing_REP1.tab \
--dpi 500 \
--plotHeight 15 \
--plotWidth 20 \
--perGroup \
--refPointLabel Center \
--yMin 0 \
--samplesLabel BG4_rep1a BG4_rep1c BG4_rep2a BG4_rep2b BG4_rep2c BG4_rep3a BG4_rep3b BG4_rep3c \
--plotTitle Signalstrength_BG4-ChIP_in_nChIP_SP2_chaning \
--numPlotsPerRow 2"


sbatch -o Plotprofile.%j.log -J DT_Profile --mem=20000 \
--wrap "computeMatrix reference-point \
--referencePoint center \
-b 300 -a 300 \
-S \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio1_a.bs50.bl.RPKM.bw \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio1_c.bs50.bl.RPKM.bw \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio2_a.bs50.bl.RPKM.bw \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio2_b.bs50.bl.RPKM.bw \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio2_c.bs50.bl.RPKM.bw \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio3_a.bs50.bl.RPKM.bw \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio3_b.bs50.bl.RPKM.bw \
/scratchb/sblab/spiege01/ENCODE_K562/Signal_based_Analysis/bigWig/G4_bio3_c.bs50.bl.RPKM.bw \
-R \
DBA_GLM_DMSOvsPDS_REP1_top100_cpm.bed \
DBA_GLM_DMSOvsPDS_REP1_bottom100_cpm.bed \
DBA_GLM_DMSOvsPDS_REP1_top400_cpm.bed \
DBA_GLM_DMSOvsPDS_REP1_bottom400_cpm.bed \
--skipZeros \
-o BG4-ChIP_in_SP2_changing_cpm_REP1.mat.gz &&
plotProfile -m BG4-ChIP_in_SP2_changing_cpm_REP1.mat.gz \
-out BG4-ChIP_in_SP2_changing_cpm_REP1.pdf \
--outFileNameData BG4-ChIP_in_SP2_changing_cpm_REP1.tab \
--dpi 500 \
--plotHeight 15 \
--plotWidth 20 \
--perGroup \
--refPointLabel Center \
--yMin 0 \
--samplesLabel BG4_rep1a BG4_rep1c BG4_rep2a BG4_rep2b BG4_rep2c BG4_rep3a BG4_rep3b BG4_rep3c \
--plotTitle Signalstrength_BG4-ChIP_in_nChIP_SP2_chaning \
--numPlotsPerRow 2"

```
