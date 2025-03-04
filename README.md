# Moroidin-Transcriptomics
Softwares: sratoolkit v2.10.9, trimgalore v0.6.7, seqkit v2.3.0, orfipy v0.0.4, spades v3.15.5, trinity v2.15.5, megahit v1.2.9, blast-plus v2.16.0, stringtie v2.2.1, samtools v1.21, star v2.7.11a

The following scripts were submitted to a computational cluster via SLURM.


# seqkit-grep search of QLLVW motifs in unassembled data
Unassembled RNA-seq data was searched for the presence of the stephanotic acid core peptide motif QLLVW with seqkit from 6frame-translated raw read data as follows:
1. SRA-download – see transcriptome assembly
Raw RNA-seq datasets were downloaded as described under SRA download with sratoolkit (v2.10.9).
2. Trimming
Raw RNA-seq datasets were trimmed by TrimGalore with default settings with the following script in batch mode (see instructions in Transcriptome assembly):
