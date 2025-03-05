# Moroidin-Transcriptomics
Softwares: 

sratoolkit v2.10.9

trimgalore v0.6.7

seqkit v2.3.0

orfipy v0.0.4

spades v3.15.5

trinity v2.15.5

megahit v1.2.9

blast-plus v2.16.0

stringtie v2.2.1

samtools v1.21

star v2.7.11a

The following scripts were submitted to a computational cluster via SLURM.

# Seqkit-grep search of QLLVW motifs in unassembled data
Unassembled RNA-seq data was searched for the presence of the stephanotic acid core peptide motif QLLVW with seqkit from 6frame-translated raw read data as follows:
1. SRA-download – see transcriptome assembly
Raw RNA-seq datasets were downloaded as described under SRA download with sratoolkit (v2.10.9).
2. Trimming
Raw RNA-seq datasets were trimmed by TrimGalore with default settings with the following script in batch mode (see instructions in Transcriptome assembly):
```
#!/bin/bash
#SBATCH --job-name=trimgalore
#SBATCH --account=your_account
#SBATCH --partition=standard 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=7g
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END
#SBATCH --output=./trimgalore-%j
module load Bioinformatics
module load trimgalore/0.6.7-ztb2tpz
trim_galore --cores 4 --paired ./SRA#_1.fastq ./ SRA#_2.fastq
```

3. Seqkit-remove duplicates
Seqkit (v 2.3.0) remove duplicate command was used to remove duplicate reads in raw RNA-seq datasets.

    a.Generate directories for trimmed fwd reads and trimmed rev reads of trimmed NCBI SRA fastq-files
```
mkdir input_data_trimmed_1
mkdir input_data_trimmed_2
```

b. Move trimmed fwd reads to input_data_trimmed_1/ directory
```
mv *_1.fq /path/to/input_data_trimmed_1/
```
