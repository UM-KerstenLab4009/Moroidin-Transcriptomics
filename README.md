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

a. Generate directories for trimmed fwd reads and trimmed rev reads of trimmed NCBI SRA fastq-files
```
mkdir input_data_trimmed_1
mkdir input_data_trimmed_2
```

b. Move trimmed fwd reads to input_data_trimmed_1/ directory
```
mv *_1.fq /path/to/input_data_trimmed_1/
```

c. Move trimmed rev reads to input_data_trimmed_2/ directory
```
mv *_2.fq /path/to/input_data_trimmed_2/
```

d. Run batch seqkit-rmdup script:
```
#!/bin/bash
#SBATCH --job-name=remove-duplicates
#SBATCH --account=your_account
#SBATCH --partition=standard 
#SBATCH --array=1-(insert-number-of-datasets-in-input_data_trimmed_1_directory)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --time=24:00:00
#SBATCH --mem=48g
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END
#SBATCH --output=./rmdup-%j 
module load Bioinformatics
module load seqkit
file1=$(ls ./input_data_trimmed_1/ | sed -n ${SLURM_ARRAY_TASK_ID}p)
file2=$(ls ./input_data_trimmed_2/ | sed -n ${SLURM_ARRAY_TASK_ID}p)
seqkit rmdup -s ./input_data_trimmed_1/${file1} -o ./${file1}
seqkit rmdup -s ./input_data_trimmed_2/${file2} -o ./${file2}
```

4. Seqkit-fastq-to-fasta-conversion
Seqkit (v 2.3.0) fastq-to-fasta command was used to remove duplicate reads in raw RNA-seq datasets.

a. Generate directories for fwd reads and rev reads of NCBI SRA fastq-files:
```
mkdir input_data_rmdup_1
mkdir input_data_rmdup_2
```

b. Move remove-duplicate fwd reads to input_data_rmdup_1/ directory
```
mv *_1.fq /path/to/input_data_rmdup_1/
```

c. Move remove-duplicate rev reads to input_data_rmdup_2/ directory
```
mv *_2.fq /path/to/input_data_rmdup_2/
```

d. Run batch seqkit-rmdup script:
```
#!/bin/bash
#SBATCH --job-name=fastq2fasta
#SBATCH --account=your_account
#SBATCH --partition=standard 
#SBATCH --array=1-(insert-number-of-datasets-in-input_data_rmdup_1_directory)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=7g
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END
#SBATCH --output=./fq2fa-%j
module load Bioinformatics
module load seqkit
file1=$(ls ./input_data_rmdup_1/ | sed -n ${SLURM_ARRAY_TASK_ID}p)
file2=$(ls ./input_data_rmdup_2/ | sed -n ${SLURM_ARRAY_TASK_ID}p)
seqkit fq2fa ./input_data_rmdup_1/${file1} -o ./${file1}.fasta
seqkit fq2fa ./input_data_rmdup_2/${file2} -o ./${file2}.fasta
```

e. Combination of unassembled fasta files
Combine unassembled fasta files into one fasta file in their directory with the following command:
```
cat *.fasta > unassembled_data_search.fasta
```

5. Orfipy frame translation
Install orfipy:
```
pip install orfipy
```
The unassembled RNA-seq read fasta file was translated in 6 frames with orfipy (v0.0.4). The minimum bp length was set to 90 bp to include short read lengths in the raw RNA-seq data search.

```
#!/bin/bash
#SBATCH --job-name=orfipy
#SBATCH --account=your_account
#SBATCH --partition=standard 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --time=24:00:00
#SBATCH --mem=48g
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END
#SBATCH --output=./orfipy-%j
orfipy --pep unassembled_data_search.pep --min 90 --between-stops unassembled_data_search.fasta
```

6. Seqkit-grep-QLLVW-search
Seqkit (v 2.3.0) grep command was used to search the unassembled, 6frame-translated RNA-seq data for the core peptide motif QLLVW with the command:
```
#!/bin/bash
#SBATCH --job-name=seqkit-grep-QLLVW
#SBATCH --account=your_account
#SBATCH --partition=standard 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --time=24:00:00
#SBATCH --mem=48g
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END
#SBATCH --output=./seqkit-QLLVW-%j
module load Bioinformatics
module load seqkit
awk '/^>/ {print $1; next} {print}' unassembled_data_search.pep > cleaned_unassembled_data_search.pep
seqkit grep -s -r -p "QLLVW" cleaned_unassembled_data_search.pep > seqkit_grep_hits.pep
```

7. Extract SRA codes and count QLLVW reads per SRA code
SRA codes of unassembled raw RNA-seq datasets with reads encoding QLLVW motifs were extracted with the awk command in the following script:
```
#!/bin/bash
#SBATCH --job-name=SRA-codes
#SBATCH --account=your_account
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=7g
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END
#SBATCH --output=./sra-extract-%j
awk '/^>/ { split($1, a, "."); print substr(a[1], 2) }' seqkit_grep_hits.pep | sort | uniq -c > seqkit_grep_hits_sra_codes.txt
```

# Transcriptome assembly – single dataset
1. SRA-download
The following scripts were submitted to a computational cluster via SLURM.
Configure your cluster for sratools commands by running the command:
```
./vdb-config -i
```

Target SRA datasets were primarily paired-ended data. To download a single raw RNA-seq dataset from the NCBI sequence read archive (SRA), run the script:
```
#!/bin/bash
#SBATCH --job-name=sra-download
#SBATCH --account=your_account
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=10g
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END
#SBATCH --output=./sratools-%j
source /etc/profile.d/http_proxy.sh
module load Bioinformatics
module load sratoolkit/2.10.9-udmejx7
fasterq-dump SRR8782583 --split-files
```
2. Trimming

To trim one raw RNA-seq dataset with TrimGalore default settings, run the script:
```
#!/bin/bash
#SBATCH --job-name=trimgalore
#SBATCH –account=your_account
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

3. Transcriptome assembly
   
Three assembler softwares were used for de novo transcriptome assembly from TrimGalore-trimmed RNA-seq datasets

a. SPAdes – paired-ended datasets
```
#!/bin/bash
#SBATCH --job-name=SPAdes
#SBATCH --account=your_account
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=48g
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END
#SBATCH --output=./SPAdes-%j
module load Bioinformatics
module load spades/3.15.5-jhe6qq2
spades.py --rna -1 ./SRA#_1_val_1.fq -2 ./SRA#_2_val_2.fq -o spades_SRA#
cd spades_SRA#
mv transcripts.fasta spades_SRA#.fasta
```

b. SPAdes – single-ended datasets
```
#!/bin/bash
#SBATCH --job-name=SPAdes
#SBATCH --account=your_account
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=48g
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END
#SBATCH --output=./SPAdes-%j
module load Bioinformatics
module load spades/3.15.5-jhe6qq2
spades.py --rna -s ./SRA#_trimmed.fq -o spades_SRA#
cd spades_SRA#
mv transcripts.fasta spades_SRA#.fasta
```

c. Trinity
```
#!/bin/bash
#SBATCH --job-name=Trinity
#SBATCH --account=your_account
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=48g
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END
#SBATCH --output=./Trinity-%j
module load Bioinformatics
module load trinity/2.15.1
LOG_FILE="trinity_log.txt"
TRINITY_PARAMS="--SS_lib_type RF --max_memory 180G --CPU 8 --trimmomatic --seqType fq --left SRA#_1_val_1.fq --right SRA#_2_val_2.fq --min_glue 2 --min_kmer_cov 1 --full_cleanup --no_normalize_reads" 
Trinity $TRINITY_PARAMS
```

d. MEGAHIT
```
#!/bin/bash
#SBATCH --job-name=megahit
#SBATCH --account=your_account
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=48g
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END
#SBATCH --output=./megahit-%j
module load MEGAHIT/1.2.9
MEGAHIT -1 ./SRA#_1_val_1.fq -2 ./SRA#_2_val_2.fq -o megahit_ SRA#
```

# Transcriptome assembly – multiple datasets

1. 
