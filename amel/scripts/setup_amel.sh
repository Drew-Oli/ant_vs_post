#!/bin/sh

#######################
#Amel mellifera set up#
#######################

#log#
exec &> setup_amel.log

#get raw_reads files and concatenate lane replicates#
mkdir -p /Volumes/archive/deardenlab/drew_oliphant/projects/ant_v_post/amel/data/raw_reads
read_dir=/Volumes/archive/deardenlab/drew_oliphant/projects/ant_v_post/amel/data/raw_reads
zcat /Volumes/archive/deardenlab/HTS_raw_sequencing_reads/NZGL02506qmp_flies\&ant_post/NZGL02506_fastq/Anterior_posterior/CBBUJANXX-2506-11-21-01_S110_L00*_R1_001.fastq.gz > $read_dir/amel_ant1_r1.fq
zcat /Volumes/archive/deardenlab/HTS_raw_sequencing_reads/NZGL02506qmp_flies\&ant_post/NZGL02506_fastq/Anterior_posterior/CBBUJANXX-2506-11-21-01_S110_L00*_R2_001.fastq.gz > $read_dir/amel_ant1_r2.fq
zcat /Volumes/archive/deardenlab/HTS_raw_sequencing_reads/NZGL02506qmp_flies\&ant_post/NZGL02506_fastq/Anterior_posterior/CBBUJANXX-2506-12-21-01_S111_L00*_R1_001.fastq.gz > $read_dir/amel_ant2_r1.fq
zcat /Volumes/archive/deardenlab/HTS_raw_sequencing_reads/NZGL02506qmp_flies\&ant_post/NZGL02506_fastq/Anterior_posterior/CBBUJANXX-2506-12-21-01_S111_L00*_R2_001.fastq.gz > $read_dir/amel_ant2_r2.fq
zcat /Volumes/archive/deardenlab/HTS_raw_sequencing_reads/NZGL02506qmp_flies\&ant_post/NZGL02506_fastq/Anterior_posterior/CBBUJANXX-2506-13-21-01_S112_L00*_R1_001.fastq.gz > $read_dir/amel_ant3_r1.fq
zcat /Volumes/archive/deardenlab/HTS_raw_sequencing_reads/NZGL02506qmp_flies\&ant_post/NZGL02506_fastq/Anterior_posterior/CBBUJANXX-2506-13-21-01_S112_L00*_R2_001.fastq.gz > $read_dir/amel_ant3_r2.fq
zcat /Volumes/archive/deardenlab/HTS_raw_sequencing_reads/NZGL02506qmp_flies\&ant_post/NZGL02506_fastq/Anterior_posterior/CBBUJANXX-2506-14-21-01_S113_L00*_R1_001.fastq.gz > $read_dir/amel_post2_r1.fq
zcat /Volumes/archive/deardenlab/HTS_raw_sequencing_reads/NZGL02506qmp_flies\&ant_post/NZGL02506_fastq/Anterior_posterior/CBBUJANXX-2506-14-21-01_S113_L00*_R2_001.fastq.gz > $read_dir/amel_post2_r2.fq
zcat /Volumes/archive/deardenlab/HTS_raw_sequencing_reads/NZGL02506qmp_flies\&ant_post/NZGL02506_fastq/Anterior_posterior/CBBUJANXX-2506-15-21-01_S114_L00*_R1_001.fastq.gz > $read_dir/amel_post3_r1.fq
zcat /Volumes/archive/deardenlab/HTS_raw_sequencing_reads/NZGL02506qmp_flies\&ant_post/NZGL02506_fastq/Anterior_posterior/CBBUJANXX-2506-15-21-01_S114_L00*_R2_001.fastq.gz > $read_dir/amel_post3_r2.fq
zcat /Volumes/archive/deardenlab/HTS_raw_sequencing_reads/NZGL02506qmp_flies\&ant_post/NZGL02506_fastq/Anterior_posterior/CBBUJANXX-2506-16-21-01_S115_L00*_R1_001.fastq.gz > $read_dir/amel_post7_r1.fq
zcat /Volumes/archive/deardenlab/HTS_raw_sequencing_reads/NZGL02506qmp_flies\&ant_post/NZGL02506_fastq/Anterior_posterior/CBBUJANXX-2506-16-21-01_S115_L00*_R2_001.fastq.gz > $read_dir/amel_post7_r2.fq

#make sample.txt#
mkdir -p /Volumes/archive/deardenlab/drew_oliphant/projects/ant_v_post/amel/output/02_salmon
output_dir=/Volumes/archive/deardenlab/drew_oliphant/projects/ant_v_post/amel/output
echo "
sample  group
amel_ant1       anterior
amel_ant2       anterior
amel_ant3       anterior
amel_post2      posterior
amel_post3      posterior
amel_post7      posterior" > $output_dir/samples.txt

#download reference data#
mkdir -p /Volumes/archive/deardenlab/drew_oliphant/projects/ant_v_post/amel/data/ref
ref_dir=/Volumes/archive/deardenlab/drew_oliphant/projects/ant_v_post/amel/data/ref
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz -P $ref_dir
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.gff.gz -P $ref_dir
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_rna.fna.gz -P $ref_dir
gunzip $ref_dir/*.gz
singularity run -B /Volumes/ docker://quay.io/biocontainers/eggnog-mapper:2.0.1--py_1 download_eggnog_data.py -y -f --data_dir $ref_dir

#make virtual environment in which to run snakemake#
mkdir -p /Volumes/archive/deardenlab/drew_oliphant/projects/ant_v_post/amel/venv
venv_dir=/Volumes/archive/deardenlab/drew_oliphant/projects/ant_v_post/amel/venv
python3 -m venv $venv_dir
source $venv_dir/bin/activate
pip3 install --upgrade pip
pip3 install snakemake
deactivate
