#!/bin/bash
#$ -l h_rt=96:00:00
#$ -l mem=12G
#$ -l rmem=12G
#$ -P alloteropsis
#$ -q alloteropsis.q
#$ -pe openmp 6
#$ -v OMP_NUM_THREADS=6
#$ -j y

#####################################################################################
#       Script Name:    OrganelleAssembly.sh
#       Description:    Fish organelle HiFi reads and assemble them
#       Author:         LPereiraG
#       Last updated:   24/11/2022
#####################################################################################

source /usr/local/extras/Genomics/.bashrc

#### Directories and files
wd=/mnt/fastdata/bo1lpg/Stipagrostis
mitochondria=${wd}/data
chloroplast=${wd}/data/Stipagrostis-hirtigluma-chloroplast.fa
reads=/shared/dunning_lab/User/bo1lpg-backup/Stipagrostis/data/Stipagrostis_reads.fasta

#### Parameters
sample='Stipagrostis'

#### Step 1: Fish chloroplast reads
cd ${wd}
mkdir A01_results
cd A01_results
module load apps/python/conda
source activate /shared/dunning_lab/Shared/conda_env/blasr
blasr ${reads} ${chloroplast} -m 1 --header --minAlnLength 5000 --minPctSimilarity 99 --nproc 6 >> Blast_raw
cat Blast_raw | cut -f 1 -d " " | sort | uniq > reads_to_fish_chloroplast.txt
sed -i 's/\/0_[0-9]\{5\}$//g' reads_to_fish_chloroplast.txt
sed -i 's/\/0_[0-9]\{4\}$//g' reads_to_fish_chloroplast.txt
source /usr/local/extras/Genomics/.bashrc
seqtk subseq ${reads} reads_to_fish_chloroplast.txt >> chloroplast_reads.fa

#### Step 2: Assemble plastome
/shared/dunning_lab/Shared/programs/hifiasm-0.16.1/./hifiasm -o ${sample}_chloroplast.asm -t 6 chloroplast_reads.fa
awk '/^S/{print ">"$2;print $3}' ${sample}_chloroplast.asm.bp.p_ctg.gfa > ${sample}_chloroplast.asm.bp.p_ctg.fa

#### Step 3: Circularize plastome
source activate /shared/dunning_lab/Shared/conda_env/circlator
circlator clean ${sample}_chloroplast.asm.bp.p_ctg.fa ${sample}_chloroplast_clean

#### OPTIONAL: ONLY IF THE RESULTS ARE NOT ACCURATE
#### Subsample 1/3 of the reads and repeat assembly to see if it improves
mkdir subsampling
cd subsampling
num=$(grep ">" ${wd}/A01_results/chloroplast_reads.fa | wc -l)
sub=$(expr $num / 3)
seqtk sample -s100 ../chloroplast_reads.fa ${sub} > subsample_chloroplast.fa
/shared/dunning_lab/Shared/programs/hifiasm-0.16.1/./hifiasm -o ${sample}_chloroplast_subsample.asm -t 6 subsample_chloroplast.fa
awk '/^S/{print ">"$2;print $3}' ${sample}_chloroplast_subsample.asm.bp.p_ctg.gfa > ${sample}_chloroplast_subsample.asm.bp.p_ctg.fa
source activate /shared/dunning_lab/Shared/conda_env/circlator
circlator clean ${sample}_chloroplast_subsample.asm.bp.p_ctg.fa ${sample}_chloroplast_subsample_clean

#### Step 4: Get a similar mitome
mkdir ${wd}/A02_results
cd ${wd}/A02_results
apptainer exec --bind ${wd} /usr/local/packages/singularity/images/mitohifi/mitohifi.sif findMitoReference.py \
  --species "Stipagrostis hirtigluma" --outfolder ${wd}/A02_results -t mitochondrion

#### Step 5: Run the main script to assemble mitome
cp ${reads} ${mitochondria}
apptainer exec --bind ${wd} /usr/local/packages/singularity/images/mitohifi/mitohifi.sif mitohifi.py \
 -r ${mitochondria}/${sample}_reads.fasta -f ${wd}/A02_results/*.fasta -g ${wd}/A02_results/*.gb \
 -t 6 -a plant

#### Step 6: Clean hifiasm results to keep only unique sequences
cd reads_mapping_and_assembly
source activate /shared/dunning_lab/Shared/conda_env/circlator
circlator clean hifiasm.contigs.fasta hifiasm.contigs.clean

#### OPTIONAL: ONLY IF THE RESULTS ARE NOT ACCURATE
#### Subsample 1/3 of the reads and repeat assembly to see if it improves
mkdir subsampling
cd subsampling
num=$(grep ">" ${wd}/A02_results/reads_mapping_and_assembly/gbk.HiFiMapped.bam.fasta | wc -l) 
sub=$(expr $num / 3)
seqtk sample -s100 ../gbk.HiFiMapped.bam.fasta ${sub} > subsample_mitochondria.fa
/shared/dunning_lab/Shared/programs/hifiasm-0.16.1/./hifiasm --primary -f 0 -o ${sample}_mitochondria_subsample.asm -t 6 subsample_mitochondria.fa
awk '/^S/{print ">"$2;print $3}' ${sample}_mitochondria_subsample.asm.p_ctg.gfa > ${sample}_mitochondria_subsample.asm.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' ${sample}_mitochondria_subsample.asm.a_ctg.gfa > ${sample}_mitochondria_subsample.asm.a_ctg.fa
cat ${sample}_mitochondria_subsample.asm.p_ctg.fa ${sample}_mitochondria_subsample.asm.a_ctg.fa > ${sample}_mitochondria.contigs.fa
source activate /shared/dunning_lab/Shared/conda_env/circlator
circlator clean ${sample}_mitochondria.contigs.fa ${sample}_mitochondria.clean.fa
