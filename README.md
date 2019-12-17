
# Transcription factors are routinely recruited to G-quadruplex DNA secondary structures in chromatin 

This repository contains the data access and computational analyses for the methods developed in our *submitted manuscript*.

### Data

- All in-house sequencing data have been deposited in the NCBI GEO database under accession number [XXX](). 

- Additional replicates for G4 ChIP-seq in K562 cells have previously been deposited in NCBI GEO [GSE107690](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107690).

- G4 forming sequences in genomic (OQS), single-stranded DNA observed via G4-seq are available in NCBI GEO [GSE63874](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63874).

- DNase-seq data for open chromatin maps and ChIP-seq maps of transcription factors and histone marks in K562 cells were downloaded from ENCODE [ENCSR000EOT](https://www.encodeproject.org/experiments/ENCSR000EOT/). Files were mapped to hg19 in Jan 2019. The meta data of the considered ChIP-seq files can be found [here]().

- RNA-seq data for K562 cells were downloaded from ENCODE [ENCSR000AEM](https://www.encodeproject.org/experiments/ENCSR000AEM/) and [ENCSR545DKY](https://www.encodeproject.org/experiments/ENCSR545DKY/)


### Code

- [**BG4 ChIP-seq analysis and G4 control data set**](G4-ChIP-seq.md)

- [**Genomic association testing**](Genomic_association_testing.md)

- [**SP2 native ChIP-seq and differential binding**](SP2_native_ChIP-seq_and_differential_binding.md)

- [**Expression analysis**](Transcript_analysis.md)



### Tools 

|Tool           | version                         | URL                                                           |
| ------------- |:--------------------------------| --------------------------------------------------------------|
| cutadapt      | 1.16                            |http://cutadapt.readthedocs.io/en/stable/installation.html     |
| BWA           | 0.7.15                          |http://bio-bwa.sourceforge.net/                                |
| Picard        | 2.14.0                          |http://broadinstitute.github.io/picard                         |
| MACS          | 2.1.1                           |http://liulab.dfci.harvard.edu/MACS/                           |
| Bedtools      | 2.26.0                          |http://bedtools.readthedocs.io/en/latest/content/overview.html |
| samtools      | 1.6                             |http://www.htslib.org/download/                                |
| Deeptools     | 3.1.2                           |https://deeptools.readthedocs.io/                              |
| PAVIS         | Accessed Jan 2019 to Dec 2019   |https://manticore.niehs.nih.gov/pavis2/                        |
| MEME-ChIP     | Accessed Jan 2019 to Dec 2019   |http://meme-suite.org/tools/meme-chip                          |
| GAT           | 1.3.5                           |https://gat.readthedocs.io                                     |
| Intervene     | 0.6.4                           |https://github.com/asntech/intervene                           |
