# Moroidin-Transcriptomics

**Softwares** 
- **sratoolkit v2.10.9**:  https://github.com/ncbi/sra-tools
- **trimgalore v0.6.7**: https://github.com/FelixKrueger/TrimGalore
- **seqkit v2.3.0**: https://github.com/shenwei356/seqkit/releases
- **orfipy v0.0.4**: https://github.com/urmi-21/orfipy
- **spades v3.15.5**: https://github.com/ablab/spades
- **trinity v2.15.5**: https://github.com/trinityrnaseq/trinityrnaseq
- **megahit v1.2.9**: https://github.com/voutcn/megahit
- **blast-plus v2.16.0**: https://github.com/ncbi/blast_plus_docs
- **stringtie v2.2.1**: https://github.com/gpertea/stringtie
- **samtools v1.21**: https://github.com/samtools/samtools
- **star v2.7.11a**: https://github.com/alexdobin/STAR
- **gffread v0.12.7**: https://github.com/gpertea/gffread
  
The following scripts were submitted to a high performance computing cluster via SLURM.

# Transcriptome assembly – single dataset
**1. SRA-download**
- Configure your cluster for sratools commands by running the command:
```
./vdb-config -i
```
  - Target SRA datasets were primarily paired-ended data. To download a single raw RNA-seq dataset from the NCBI sequence read archive (SRA), run the script:
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
**2. Trimming**
  - To trim one raw RNA-seq dataset with TrimGalore default settings, run the script:
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
**3. Transcriptome assembly** - *Three assembler softwares were used for de novo transcriptome assembly from TrimGalore-trimmed RNA-seq datasets.*
- **a. SPAdes – paired-end datasets**
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
- **b. SPAdes – single-end datasets**
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
- **c. Trinity**
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
- **d. MEGAHIT**
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
**1. Batch SRA-download**
- For large-scale transcriptome mining, SRA datasets were downloaded in batches of 100 datasets with the following script:
```
#!/bin/bash
#SBATCH --job-name=sra-download          
#SBATCH --account=your_account            
#SBATCH --partition=standard            
#SBATCH --nodes=1                             
#SBATCH --ntasks=1                               
#SBATCH --cpus-per-task=8                
#SBATCH --time=24:00:00                         
#SBATCH --mem=48g                               
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END                 
#SBATCH --output=./sratools-%j
source /etc/profile.d/http_proxy.sh
module load Bioinformatics
module load sratoolkit/2.10.9-udmejx7
sed -i 's/\r$//' SRA.txt & 
cat SRA.txt |parallel -j 10 xargs -n 25 -P 0 fasterq-dump --split-files --outdir "/path-to-directory" & jobs -l
wait
printf "\n...done\n\n"
```
**2. Batch trimming**
- Generate directories for fwd reads and rev reads of NCBI SRA fastq-files:
```
mkdir
input_data_1
```
```
mkdir
input_data_2
```
- Move fwd reads to input_data_1/ directory:
```
mv *_1.fastq /path/to/input_data_1/
```
- Move rev reads to input_data_2/ directory:
```
mv *_2.fastq /path/to/input_data_2/
```
- Run batch TrimGalore-trimming script:
```
#!/bin/bash
#SBATCH --job-name=trimgalore
#SBATCH --account=your_account
#SBATCH --partition=standard 
#SBATCH --array=1-(insert-number-of-datasets-in-input_data_1_directory)
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
file1=$(ls ./input_data_1/ | sed -n ${SLURM_ARRAY_TASK_ID}p)
file2=$(ls ./input_data_2/ | sed -n ${SLURM_ARRAY_TASK_ID}p)
trim_galore --cores 4 --paired ./input_data_1/${file1} ./input_data_2/${file2}
```
**3. Batch transcriptome assembly (SPAdes)**
- Generate directories for trimmed fwd reads and trimmed rev reads of NCBI SRA fastq-files:
```
mkdir input_data_trimmed_1
```
```
mkdir input_data_trimmed_2
```
- Move trimmed fwd reads to input_data_trimmed_1/ directory:
```
mv *_1.fq /path/to/input_data_trimmed_1/
```
- Move trimmed rev reads to input_data_trimmed_2/ directory:
```
mv *_2.fq /path/to/input_data_trimmed_2/
```
- Run SPAdes batch assembly:
```
#!/bin/bash
#SBATCH --job-name=SPAdes
#SBATCH --account=your_account 
#SBATCH --partition=standard
#SBATCH --array=1-(insert-number-of-datasets-in-input_data_trimmed_1_directory)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --time=24:00:00
#SBATCH --mem=48g
#SBATCH --mail-user=your@mail.com
#SBATCH --mail-type=END
#SBATCH --output=./SPAdes-%j
module load Bioinformatics
module load spades/3.15.5-jhe6qq2
file1=$(ls ./input_data_trimmed_1/ | sed -n ${SLURM_ARRAY_TASK_ID}p)
file2=$(ls ./input_data_trimmed_2/ | sed -n ${SLURM_ARRAY_TASK_ID}p)
spades.py --rna -1 ./input_data_trimmed_1/${file1} -2 ./input_data_trimmed_2/${file2} -o spades_$file1
cd spades_$file1
mv transcripts.fasta /path-to-directory/spades_$file1\.fasta
```
# Sequenceserver-based BLAST search and burpitide prediction
**1. BLAST database formatting** - *Addition of transcriptome assemblies to Sequenceserver requires reduction of fasta headers to less than 51 letters. Below are several scripts for assembler-specific databases formatting of assemblies for Sequenceserver addition.*
- **a. SPAdes**
```
#!/bin/bash
#SBATCH --job-name=rename
#SBATCH --account=your_account
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=10g
#SBATCH --mail-user=your@mail.com
#SBATCH --mail-type=END
#SBATCH --output=./rename-%j
for i in *fasta; do n="${i%.fasta}"; sed -i.bak "s/>[^_]\+/>$n/" $i; done
cat *fasta
sed -i 's/_length//g' *fasta
sed -i 's/_cov//g' *fasta
sed -i 's/_g//g' *fasta
sed -i 's/_i//g' *fasta
```
- **b. MEGAHIT**
```
#!/bin/bash
#SBATCH --job-name=rename
#SBATCH --account=your_account
#SBATCH --partition=standard
##SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=10g
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END 
#SBATCH --output=./rename-%j 
for i in *fa; do n="${i%.fa}"; sed -i.bak "s/>[^_]\+/>$n/" $i; done
cat *fa
sed -i 's/ flag=/_/g' *fa
sed -i 's/ multi=/_/g' *fa
sed -i 's/ len=/_/g' *fa
```
- **c. Trinity**
```
#!/bin/bash
#SBATCH --job-name=rename
#SBATCH --account=your_account
#SBATCH --partition=standard 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=10g
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END
#SBATCH --output=./rename-%j
for i in *fa; do n="${i%.fa}"; sed -i.bak "s/>[^_]\+/>$n/" $i; done
cat *fa
sed -i 's/len.*]/ /g' *.fa
```
**2. Sequenceserver search and search-hit.fasta download**
- For large-scale search of BURP domain transcripts, 100 de novo transcriptomes were combined into databases and searched via tblastn (v2.16.0+) for homologs of Sali3-2 as a BURP-domain-containing protein query (GenBank ID AAB66369) on Sequenceserver (v3.1.0) with the following BLAST parameters: evalue 1e-05, matrix BLOSUM62, gap-open 11, gap-extend 1, filter L, max, max_target_seqs 5000. tblastn hits from all databases were combined into a fasta-file.
```
>AAB66369.1 Sali3-2 [Glycine max]
MEFRCSVISFTILFSLALAGESHVHASLPEEDYWEAVWPNTPIPTALRDVLKPLPAGVEIDQLPKQIDDT
QYPKTFFYKEDLHPGKTMKVQFTKRPYAQPYGVYTWLTDIKDTSKEGYSFEEICIKKEAFEGEEKFCAKS
LGTVIGFAISKLGKNIQVLSSSFVNKQEQYTVEGVQNLGDKAVMCHGLNFRTAVFYCHKVRETTAFVVPL
VAGDGTKTQALAVCHSDTSGMNHHILHELMGVDPGTNPVCHFLGSKAILWVPNISMDTAYQTNVVV
```
**3. Orfipy**
- Search hits from Sequenceserver were 6-frame translated with orfipy with 450 bp minimum open reading frames due to the size of the target BURP protein domain of >200 amino acids (>600 bp).
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
orfipy --pep sequenceserver-hits.pep --min 450 --between-stops sequenceserver-hits.fasta
```
**4. RepeatFinder analysis** - *RepeatFinder was run as a standalone version on a computational cluster.*
- Install RepeatFinder:
```
git clone https://github.com/FlorisdeWaal/repeatfinder_standalone
```
- Copy sequenceserver-hits.pep file into the standalone directory:
```
cd repeatfinder_standalone/
cp /your-directory/sequenceserver-hits.pep .
```
- Run RepeatFinder via python2:
```
python2 run_rf.py -i sequenceserver-hits.pep -o sequenceserver-hits.html
```
- The resulting html-file includes predicted stephanotic acid-type burpitide cyclases based on the ‘QLxxW’ motif assignment. To search for other core peptide patterns, open the file known_motifs.txt in the repeatfinder_standalone/ directory:
```
nano known_motifs.txt
```
- and add a target core peptide motif, for example for moroidins:
```
>moroidin pattern
QL..W..H
```
- Preloaded core peptides in known_motifs.txt are:
```
>Lyciumin pattern
QP.....W
>Legumenin pattern
QP...Y.W
>bicyclic cyclopeptide alkaloid pattern
.L.Y..Y
>Cca-like peptide pattern
QI..W
>stephanotic acid pattern
QL..W
```

# Commandline-PHI-BLAST search and burpitide prediction
**1. Commandline BLAST search**  - *An alternative to Sequenceserver-based BLAST search and RepeatFinder-based burpitide prediction is commandline PHI-BLAST. It was applied to de novo assembled plant transcriptomes (SPAdes assembly) combined as a single fasta-file.*
- **a.	Generate transcriptome input file from multiple transcriptome assembly fasta files in a directory:**
```
cat *.fasta > all.fasta
```
- **b.	Generate query.faa file:**
```
touch query.faa
nano query.faa
```
- **Copy target protein sequence (moroidin cyclase KjaBURP), for searching stephanotic acid-type burpitides:**
```
>QIG55799.1 BURP domain protein [Kerria japonica]
MACRLSLIFAFLCLTLVACHAALSPQEVYWNSVFPQTPMPKTLSALVQPAAKNFIRYKKVDDGQTQDIDV
AADNQLLVWRGHVAIDDDAAADNQLLVWRGHVAIDDDDAAADNQLLVWRGHVAIHDDAAADNQLLVWRAH 
VANDDVDARNLLRKDPSRVLVFLEKDVHPGKTMKYSLIRSLPKSATFLPRNTAESIPFSSSKLLEILIQF 
SVQPKSVEANAMTEAILKCEVPAMRGEAKYCATSLESMIDFVTSRLGRNIRAISTEVEEGATHVQNYTIY 
HGVKKLTDKKVITCHRLRYPYVVFYCHELENTSIYMVPLKGADGTNAKAITTCHEDTSEWDPKSFVLQLL 
KVKPGTDPVCHFLSESDVVWVSNHGTYKPA
```
- **c. Generate phi_pattern.txt file:**
```
touch phi_pattern.txt
nano phi_pattern.txt
```
- **Add the following text for searching a core peptide ‘QLxxW’ (where x is any proteinogenic amino acid) and save:**
```
PA QQL-x(2)-W
```
- **d. PHI-BLAST search**
  - Install orfipy before PHI-BLAST search as described above. The same orfipy translation parameters were applied for PHI-BLAST search as for Sequenceserver-based              tblastn search of burpitide cyclase sequences (i.e. minimum of 450 bp open reading frame length, translation between stop codons).
```
#!/bin/bash
#SBATCH --job-name=phiblast-QLxxW
#SBATCH --account=your_account
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --time=24:00:00
#SBATCH --mem=48g
#SBATCH --mail-user=your@mail.com
#SBATCH --mail-type=END
#SBATCH --output=./phiblast-%j
module load Bioinformatics
module load blast-plus/2.16.0
module load seqkit
orfipy --pep all.pep --min 450 --between-stops all.fasta
cp orfipy_all.fasta_out/all.pep .
awk '/^>/ {print $1; next} {print}' all.pep > cleaned_all.pep
makeblastdb -in cleaned_all.pep -dbtype prot -out cleaned_all_db
psiblast -db cleaned_all_db -query query.faa -out cleaned_all_phiblast.txt -phi_pattern phi_pattern.txt
grep ">" cleaned_all_phiblast.txt | awk '{print $1}' | tr -d '>' > hits.txt
seqkit subseq cleaned_all.pep hits.txt > phiblast_all_sequences.pep
```
  - For searching other core peptide motifs, please change the phi_pattern.txt file according to the ___ manual and use a query burpitide cyclase protein sequence which includes the target core peptide motif.


# Genome-guided transcriptome assembly
**1.	SRA-download**   
- Download SRA RNA-seq datasets as specified above for single datasets.
  
**2. Trimming**
- Trim RNA-seq datasets as specified above for single datasets.
  
**3.	STAR alignment**
- Download the RefSeq genome assembly fasta file (e.g. genome.fasta) and genome annotation gtf/gff3 file (e.g. genome.gtf) from NCBI Genome database or a genome-specific      repository to a personal computer and transfer the files to a computational cluster in the directory for star alignment. In the following STAR alignment script, specify    the -–sjdbOverhang parameter as the read length of the target SRA RNA-seq dataset minus 1 (e.g. if the read length is 151 bp, specify 150 as -–sjdbOverhang parameter).
```
#!/bin/bash
#SBATCH --job-name=star
#SBATCH --account=your_account
#SBATCH --partition=standard 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --time=24:00:00
#SBATCH --mem=48g
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END
#SBATCH --output=./star-align-%j
module load Bioinformatics
module load star/2.7.11a-hdp2onj
STAR --runThreadN 7 --runMode genomeGenerate --genomeDir /path-to-directory/ --genomeFastaFiles ./genomic.fasta --sjdbGTFfile ./genome.gtf --sjdbOverhang 150
```
- For large genomes, more memory might be required for STAR alignment. For example, STAR alignment of the wheat genome (16 Gbp) required 422 GB memory and STAR alignment of the Nicotiana tabacum genome (4.5 Gbp) required 139 GB memory.
  
**4. STAR mapping**
- The following script was used for STAR mapping of paired-ended trimmed RNA-seq datasets:
```
#!/bin/bash
#SBATCH --job-name=star-map
#SBATCH --account=your account       
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --time=24:00:00
#SBATCH --mem=48g
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END
#SBATCH --output=./star-map-%j
module load Bioinformatics
module load star/2.7.11a-hdp2onj
STAR --runThreadN 7 --runMode alignReads --genomeDir /path-to-directory/ --readFilesIn ./SRA#_1_val_1.fq ./SRA#_2_val_2.fq
```
- The resulting output files include an Aligned.out.sam file. For large genomes, more memory might be required for STAR mapping. For example, STAR mapping of the wheat genome (16 Gbp) required 163 GB memory and STAR alignment of the Nicotiana tabacum genome (4.5 Gbp) required 63 GB memory.

**5.	BAM file generation**
- The following script formats the Aligned.out.sam file to an Aligned.out.bam which can be used as an input file for StringTie assembly.
```
#!/bin/bash
#SBATCH --job-name=samtools
#SBATCH --account=your_account
#SBATCH --partition=standard 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=7g
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END
#SBATCH --output=./samtools-%j 
module load Bioinformatics
module load samtools/1.21
samtools sort -o Aligned.out.bam Aligned.out.sam
```

**6. Transcriptome assembly**
- **a. StringTie**
  - StringTie (v2.2.1) was applied for genome-guided transcriptome assembly with the Aligned.out.bam file as follows:
```
#!/bin/bash
#SBATCH --job-name=stringtie
#SBATCH --account=your_account
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --time=24:00:00
#SBATCH --mem=48g
#SBATCH --mail-user=your@mail.com
#SBATCH --mail-type=END
#SBATCH --output=./stringtie-%j
module load Bioinformatics
module load stringtie
stringtie -o stringtie-SRA#.gtf Aligned.out.bam
```
         The output gtf-file was converted to a fasta-file for subsequent BLAST searches (described under ‘Sequenceserver-based BLAST search and burpitide prediction’ or             ‘Commandline-PHI-BLAST search and burpitide prediction’) with gffread. Install gffread as specified under https://github.com/gpertea/gffread. Run gffread                    conversion of gtf.file to fasta.file with command:
```
gffread -w stringtie-SRA#.fa -g /path-to-directory/stringtie-SRA#.gtf
```
      b. Trinity

         Trinity (v2.15.1) was applied for genome-guided transcriptome assembly with the Aligned.out.bam file as follows:  
```
#!/bin/bash
#SBATCH --job-name=trinity
#SBATCH --account=your_account
#SBATCH --partition=standard 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=36:00:00
#SBATCH --mem=48g
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END
#SBATCH --output=./trinity-%j
module load Bioinformatics
module load trinity/2.15.1
LOG_FILE="trinity_log.txt"
TRINITY_PARAMS="--SS_lib_type RF --max_memory 48G --CPU 8 --genome_guided_bam Aligned.out.bam --genome_guided_max_intron 50000 --output trinity-SRA#" 
Trinity $TRINITY_PARAMS
```

# PHI-BLAST search of KjaBURP and QLLVW motif in unassembled RNA-seq data
Unassembled RNA-seq data was searched by PHI-BLAST (v2.16.0) with fused moroidin cyclase KjaBURP as a query (GenBank QIG55799.1) in the default setting for the presence of the stephanotic acid core peptide motif QLLVW as follows:

1. **SRA-download**
      - Raw RNA-seq datasets were downloaded as described under SRA download (transcriptome assembly) with sratoolkit (v2.10.9). 100 RNA-seq datasets were
   downloaded per directory with a total number of 11 directories.
3. Trimming

   Raw RNA-seq datasets were trimmed by TrimGalore with default settings with the following script in batch mode (see instructions in Transcriptome assembly).
4. Combine paired reads
  
   To reduce file size, paired-end reads were combined into one file.

   Generate directories for fwd reads
```
mkdir input_data_1
mkdir input_data_2
```
   Move fwd reads to input_data_1/ directory:
```
mv *_1.fq /path/to/input_data_1/
```
   Move rev reads to input_data_2/ directory:
```
mv *_2.fq /path/to/input_data_2/
```

   Combine paired-end fastq files (example script for 100 paired-end files): 
```
#!/bin/bash
#SBATCH --job-name=cat-fastq
#SBATCH --account=your_account
#SBATCH --partition=standard
#SBATCH --array=1-100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=7g
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END
#SBATCH --output=./cat-%j
file1=$(ls ./input_data_1/ | sed -n ${SLURM_ARRAY_TASK_ID}p)
file2=$(ls ./input_data_2/ | sed -n ${SLURM_ARRAY_TASK_ID}p)
cat ./input_data_1/${file1} ./input_data_2/${file2} > ${file1}
```

4. Remove duplicate reads with seqkit2

The seqkit2 (v2.3.0) remove duplicate command was used to remove duplicate reads in combined, trimmed raw RNA-seq datasets.

4a. Generate directory for combined trimmed reads of NCBI SRA fastq-files:
```
mkdir input_data_rmdup
```

4b. Move combined trimmed reads to input_data directory:
```
mv *.fastq > input_data_rmdup
```

      c. Run the following script to remove duplicates:
```
#!/bin/bash
#SBATCH --job-name=remove-duplicates
#SBATCH --account=your_account
#SBATCH --partition=standard
#SBATCH --array=1-100 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=7g
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END
#SBATCH --output=./rmdup-%j
module load Bioinformatics
module load seqkit
file1=$(ls ./input_data/ | sed -n ${SLURM_ARRAY_TASK_ID}p)
seqkit rmdup -s ./input_data/${file1} -o ./${file1}
```

5. Seqkit-fastq-to-fasta-conversion

The seqkit2 (v2.3.0) fastq-to-fasta command was used to convert fastq reads to fasta format.
      a. Generate directories for fwd reads and rev reads of NCBI SRA fastq-files:
```
mkdir input_data_fq2fa
```

      b. Move fastq reads to input_data/ directory:
```
mv *.fq input_data_fq2fa/
```

      c. Run batch seqkit-fq2fa script:
```
#!/bin/bash
#SBATCH --job-name=fastq2fasta
#SBATCH --account=your_account
#SBATCH --partition=standard 
#SBATCH --array=1-100 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=7g
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END
#SBATCH --output=./fq2fa-%j
module load Bioinformatics
module load seqkit
file1=$(ls ./input_data_fq2fa/ | sed -n ${SLURM_ARRAY_TASK_ID}p)
seqkit fq2fa ./input_data_fq2fa/${file1} -o ./${file1}.fasta
```

6. Combine fasta files

Fasta files were combined in one directory (usually 100 fasta files were combined) to a single file.
```
#!/bin/bash
#SBATCH --job-name=cat-fasta
#SBATCH --account=rkersten1
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=7g
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END
#SBATCH --output=./cat-%j
cat *.fasta > unassembled_data.fasta
```

7. Fasta header truncation

To reduce the size of the unassembled data fasta file in one directory, fasta headers in the file were reduced with the following script:
```
#!/bin/bash
#SBATCH --job-name=cleaning-fasta
#SBATCH --account=your_account                                    
#SBATCH --partition=standard 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --time=24:00:00
#SBATCH --mem=48g
#SBATCH --mail-user=rkersten@umich.edu
#SBATCH --mail-type=END
#SBATCH --output=./cleaning-%j
module load Bioinformatics
awk '/^>/ {print $1; next} {print}' unassembled_data.fasta > cleaned_unassembled_data.fasta
```

8. 6-frame translation with orfipy

The unassembled RNA-seq read fasta file was translated in 6 frames with orfipy (v0.0.4). The minimum bp length was set to 90 bp to include short read lengths in the raw RNA-seq data search. The resulting minimum protein sequence length was 30 amino acids which is the recommended minimum input length for PHI-BLAST.

Install orfipy:
```
pip install orfipy
```
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
orfipy --pep unassembled_data.pep --min 90 --between-stops cleaned_unassembled_data.fasta
```
Move the pep file to the original directory:
```
mv ./orfipy_cleaned_unassembled_data.fasta_out/unassembled_data.pep .
```

9. Generation of BLAST database
    
      a. Split fasta file

         To generate a searchable BLAST database from the large unassembled data, the final fasta file was split into files of 1000000000 lines (~60GB).
```
#!/bin/bash
#SBATCH --job-name=split
#SBATCH --account=your_account
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=7g
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END
#SBATCH --output=./split-%j
split -l 1000000000 unassembled_data.pep small_unassembled_
```
         Resulting files were formatted to a .pep file.

      b. Makeblastdb

         An input_data directory was created:
```
mkdir input_data
```
         and split pep. files were moved to the input_data directory.
```
mv small* input_data/ 
```
         BLAST databases of the corresponding split pep. files were generated with the following script:
```
#!/bin/bash
#SBATCH --job-name=makeblastdb
#SBATCH --account=your_account
#SBATCH --partition=standard 
#SBATCH --array=1-(number of .pep files in input_data/)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --time=24:00:00
#SBATCH --mem=48g
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END
#SBATCH --output=./phiblast-%j
module load Bioinformatics
module load blast-plus/2.16.0
file1=$(ls ./input_data/ | sed -n ${SLURM_ARRAY_TASK_ID}p)
makeblastdb -in ./input_data800/${file1} -dbtype prot -out ${file1}_db
```

10. PHI-BLAST search

      a. Generate a query file:
```
touch query.faa
```
         open the query.faa file:
```
nano query.faa
```
         and copy KjaBURP protein sequence into it:
```
>QIG55799.1 BURP domain protein [Kerria japonica]
MACRLSLIFAFLCLTLVACHAALSPQEVYWNSVFPQTPMPKTLSALVQPAAKNFIRYKKVDDGQTQDIDV
AADNQLLVWRGHVAIDDDAAADNQLLVWRGHVAIDDDDAAADNQLLVWRGHVAIHDDAAADNQLLVWRAH 
VANDDVDARNLLRKDPSRVLVFLEKDVHPGKTMKYSLIRSLPKSATFLPRNTAESIPFSSSKLLEILIQF 
SVQPKSVEANAMTEAILKCEVPAMRGEAKYCATSLESMIDFVTSRLGRNIRAISTEVEEGATHVQNYTIY 
HGVKKLTDKKVITCHRLRYPYVVFYCHELENTSIYMVPLKGADGTNAKAITTCHEDTSEWDPKSFVLQLL 
KVKPGTDPVCHFLSESDVVWVSNHGTYKPA
```

      b. Generate a phi_pattern.txt file:
```
touch phi_pattern.txt
```
         open the phi_pattern.txt file:
```
nano phi_pattern.txt
```
         and copy the stephanotic acid QLLVW phi_pattern into it. The additional Q is required for successful PHI-BLAST search of the QLLVW motif.
```
PA QQLLVW
```
      c. Run PHI-BLAST

         The following script was run for PHI-BLAST of the unassembled BLAST databases:
```
#!/bin/bash
#SBATCH --job-name=phiblast-QLLVW-all-1
#SBATCH --account=your_account
#SBATCH --partition=standard
#SBATCH --array=1 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --time=24:00:00
#SBATCH --mem=48g
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=END
#SBATCH --output=./phiblast-%j
module load Bioinformatics
module load blast-plus/2.16.0
psiblast -db ./100/small_100_aa.pep_db -query query.faa -out small_100_aa.pep_db_phiblast.txt -phi_pattern phi_pattern.txt
```
phiblast.txt were analyzed by RepeatFinder for SRA datasets with KjaBURP search hits with QLLVW core peptides.

# Seqkit-grep search of QLLVW motifs in unassembled data
Unassembled RNA-seq data was searched for the presence of the stephanotic acid core peptide motif QLLVW with seqkit from 6frame-translated raw read data as follows:
      a.  Seqkit-grep-QLLVW-search

         Seqkit (v2.3.0) grep command was used to search unassembled, 6-frame-translated RNA-seq data for the core peptide motif QLLVW with the following command.       
         unassembled, 6frame-translated RNA-seq data was generated as for PHI-BLAST search of unassembled RNA-seq data (steps 1-8).
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
      b. Extract SRA codes and count QLLVW reads per SRA code

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
