#!/bin/bash
#$ -l h_rt=96:00:00
#$ -l mem=12G
#$ -l rmem=12G
#$ -P alloteropsis
#$ -q alloteropsis.q
#$ -pe openmp 8
#$ -v OMP_NUM_THREADS=8
#$ -j y

#####################################################################################
#       Script Name:    GenomeAnnotation.sh
#       Description:    Annotate gene models
#       Author:         LPereiraG
#       Last updated:   24/11/2022
#####################################################################################

source /usr/local/extras/Genomics/.bashrc

#### Directories and files
wd=/mnt/fastdata/bo1lpg/annotation/Stipagrostis
SRR=${wd}/SRR.list

#### Parameters
sample='Stipagrostis'

#### First, get gene models using RNA-seq available for the target species

#### Step 1: Download RNA-seq data
cd ${wd}
mkdir data
cd data
cat ${SRR} | while read line; do fastq-dump --split-files ${line}; done

#### Step 2: Filter and trim reads
# Copy TruSeq.fa file containing the Illumina adapters in this folder
mkdir trimmed
cat ${SRR} | while read line; do trimmomatic PE -phred33 ${line}_1.fastq ${line}_2.fastq \
  trimmed/${line}_1_paired.fastq trimmed/${line}_1_unpaired.fastq \
  trimmed/${line}_2_paired.fastq trimmed/${line}_2_unpaired.fastq \
  ILLUMINACLIP:TruSeq.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:50; done

#### Step 3: Align reads with splicer-aware aligner
cd ${wd}
mkdir ref-genome
cd ref-genome
cp /mnt/fastdata/bo1lpg/Stipagrostis/A08_results/${sample}_nuclear_4n-softMasked.fa ${sample}_nuclear_4n-softMasked.fa
hisat2-build -f ${sample}_nuclear_4n-softMasked.fa ${sample}
cd ..
cat ${SRR} | while read line; do hisat2 -x ref-genomes/Stipagrostis -1 data/trimmed/${line}_1_paired.fastq \
  -2 data/trimmed/${line}_2_paired.fastq -p 8 -q --met-file ${sample}_${line}.stats \
  | samtools sort -O BAM > ${sample}_${line}.sorted.bam; samtools flagstat ${sample}_${line}.sorted.bam; \
  samtools view -b -q 40 -f 2 -F 12 ${sample}_${line}.sorted.bam > ${sample}_${line}.sorted.clean.bam; done

#### Step 4: Run braker2 with RNA-seq data
module load apps/python/conda
source activate /shared/dunning_lab/Shared/conda_env/braker2
export PATH="/shared/dunning_lab/Shared/programs/GeneMark/gmes_linux_64_4/ProtHint/bin:$PATH"
export PATH="/shared/dunning_lab/Shared/programs/GeneMark/gmes_linux_64_4:$PATH"
export PATH="/shared/dunning_lab/Shared/programs/GeneMark/gmes_linux_64_4/ProtHint/dependencies:$PATH"
# Need to write down the list of bam files separated by commas
braker.pl --genome ref-genomes/${sample}_nuclear_4n-softMasked.fa --bam Stipagrostis_SRR15729218.sorted.clean.bam,Stipagrostis_SRR6221622.sorted.clean.bam,Stipagrostis_SRR6221624.sorted.clean.bam,Stipagrostis_SRR6221636.sorted.clean.bam,Stipagrostis_SRR6221637.sorted.clean.bam \
    --softmasking --cores 8 --workingdir braker_RNA

#### Second, use a protein database to annotate gene models in the target species

#### Step 5: Download and prepare the proteins
cd ${wd}
mkdir prot-database
cd prot-database
wget https://v100.orthodb.org/download/odb10_plants_fasta.tar.gz
tar xvf odb10_plants_fasta.tar.gz
cat plants/Rawdata/* > ../../proteins.fasta

#### Step 6: Run ProtHint to obtain hints file, needed by braker2
python /shared/dunning_lab/Shared/programs/GeneMark/gmes_linux_64_4/ProtHint/bin/prothint.py \
  ${wd}/ref-genome/${sample}_nuclear_4n-softMasked.fa proteins.fasta \
  --workdir ${wd}/prot-database --threads 8

#### Step 7: run braker2 with protein hints data from all plants
cd ${wd}
braker.pl --genome ${wd}/ref-genome/${sample}_nuclear_4n-softMasked.fa --hints ${wd}/prot-database/prothint_augustus.gff \
  --softmasking --cores 8 --workingdir braker_prot

#### Third, combine results from both analyses to get a final list of transcripts

#### Step 8: Combine transcripts with TSEBRA
tsebra=/shared/dunning_lab/Shared/programs/TSEBRA/bin/tsebra.py
cd ${wd}
${tsebra} -g braker_RNA/augustus.hints.gtf,braker_prot/augustus.hints.gtf -c ${config} \
  -e braker_RNA/hintsfile.gff,braker_prot/hintsfile.gff -o braker_combined.gtf

#### Step 9: Generate common files for genome annotations
gffread -o ${sample}_genes.gff3 braker_combined.gtf
gffread -g ${wd}/ref-genome/${sample}_nuclear_4n-softMasked.fa -x ${sample}_CDS.fa ${sample}_genes.gff3
gffread -g ${wd}/ref-genome/${sample}_nuclear_4n-softMasked.fa -y ${sample}_aa.fa ${sample}_genes.gff3
