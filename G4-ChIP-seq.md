## G4-ChIP-seq analysis
Raw fastq reads from G4-ChIP-seq in K562 cells were trimmed with cutadapt to remove adapter sequences and low-quality reads (mapping quality < 10). Reads were aligned to the human genome (version hg19) with BWA and duplicates were removed using Picard. Peaks were called by MACS2 (p < 10-5). Peaks were merged from different replicates with bedtools multiIntersect. Only peaks overlapping in 3 out 5 replicates were considered high-confidence. 

For a more detailed description of our in-house ChIP-seq pipeline please refer to our [previous work](https://github.com/sblab-bioinformatics/dna-secondary-struct-chrom-lands/blob/master/Methods.md).




