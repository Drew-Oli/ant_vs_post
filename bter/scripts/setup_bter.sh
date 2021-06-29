#!/bin/sh

###################
#Bombus terrestris#
###################

#log#
exec &> bter.log

# get read data#
mkdir -p /Volumes/archive/deardenlab/drew_oliphant/projects/ant_v_post/bter/data/raw_reads
read_dir=/Volumes/archive/deardenlab/drew_oliphant/projects/ant_v_post/bter/data/raw_reads
zcat /Volumes/archive/deardenlab/HTS_raw_sequencing_reads/NZGL02506qmp_flies\&ant_post/NZGL02506_fastq/Anterior_posterior/CBBUJANXX-2506-23-21-01_S122_L00*_R1_001.fastq.gz > $read_dir/bter_ant1_r1.fq
zcat /Volumes/archive/deardenlab/HTS_raw_sequencing_reads/NZGL02506qmp_flies\&ant_post/NZGL02506_fastq/Anterior_posterior/CBBUJANXX-2506-23-21-01_S122_L00*_R2_001.fastq.gz > $read_dir/bter_ant1_r2.fq
zcat /Volumes/archive/deardenlab/HTS_raw_sequencing_reads/NZGL02506qmp_flies\&ant_post/NZGL02506_fastq/Anterior_posterior/CBBUJANXX-2506-24-21-01_S123_L00*_R1_001.fastq.gz > $read_dir/bter_ant2_r1.fq
zcat /Volumes/archive/deardenlab/HTS_raw_sequencing_reads/NZGL02506qmp_flies\&ant_post/NZGL02506_fastq/Anterior_posterior/CBBUJANXX-2506-24-21-01_S123_L00*_R2_001.fastq.gz > $read_dir/bter_ant2_r2.fq
zcat /Volumes/archive/deardenlab/HTS_raw_sequencing_reads/NZGL02506qmp_flies\&ant_post/NZGL02506_fastq/Anterior_posterior/CBBUJANXX-2506-25-21-01_S124_L00*_R1_001.fastq.gz > $read_dir/bter_ant3_r1.fq
zcat /Volumes/archive/deardenlab/HTS_raw_sequencing_reads/NZGL02506qmp_flies\&ant_post/NZGL02506_fastq/Anterior_posterior/CBBUJANXX-2506-25-21-01_S124_L00*_R2_001.fastq.gz > $read_dir/bter_ant3_r2.fq
zcat /Volumes/archive/deardenlab/HTS_raw_sequencing_reads/NZGL02506qmp_flies\&ant_post/NZGL02506_fastq/Anterior_posterior/CBBUJANXX-2506-26-21-01_S125_L00*_R1_001.fastq.gz > $read_dir/bter_post1_r1.fq
zcat /Volumes/archive/deardenlab/HTS_raw_sequencing_reads/NZGL02506qmp_flies\&ant_post/NZGL02506_fastq/Anterior_posterior/CBBUJANXX-2506-26-21-01_S125_L00*_R2_001.fastq.gz > $read_dir/bter_post1_r2.fq
zcat /Volumes/archive/deardenlab/HTS_raw_sequencing_reads/NZGL02506qmp_flies\&ant_post/NZGL02506_fastq/Anterior_posterior/CBBUJANXX-2506-27-21-01_S126_L00*_R1_001.fastq.gz > $read_dir/bter_post2_r1.fq
zcat /Volumes/archive/deardenlab/HTS_raw_sequencing_reads/NZGL02506qmp_flies\&ant_post/NZGL02506_fastq/Anterior_posterior/CBBUJANXX-2506-27-21-01_S126_L00*_R2_001.fastq.gz > $read_dir/bter_post2_r2.fq
zcat /Volumes/archive/deardenlab/HTS_raw_sequencing_reads/NZGL02506qmp_flies\&ant_post/NZGL02506_fastq/Anterior_posterior/CBBUJANXX-2506-28-21-01_S127_L00*_R1_001.fastq.gz > $read_dir/bter_post3_r1.fq
zcat /Volumes/archive/deardenlab/HTS_raw_sequencing_reads/NZGL02506qmp_flies\&ant_post/NZGL02506_fastq/Anterior_posterior/CBBUJANXX-2506-28-21-01_S127_L00*_R2_001.fastq.gz > $read_dir/bter_post3_r2.fq

#make sample.txt#
mkdir -p /Volumes/archive/deardenlab/drew_oliphant/projects/ant_v_post/bter/output/02_salmon
output_dir=/Volumes/archive/deardenlab/drew_oliphant/projects/ant_v_post/bter/output
echo "
sample	group
bter_ant1	anterior
bter_ant2	anterior
bter_ant3	anterior
bter_post1	posterior
bter_post2	posterior
bter_post3	posterior" > $output_dir/samples.txt


#download reference data#
mkdir -p /Volumes/archive/deardenlab/drew_oliphant/projects/ant_v_post/bter/data/ref
ref_dir=/Volumes/archive/deardenlab/drew_oliphant/projects/ant_v_post/bter/data/ref
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/214/255/GCF_000214255.1_Bter_1.0/GCF_000214255.1_Bter_1.0_genomic.fna.gz -P $ref_dir
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/214/255/GCF_000214255.1_Bter_1.0/GCF_000214255.1_Bter_1.0_genomic.gff.gz -P $ref_dir
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/214/255/GCF_000214255.1_Bter_1.0/GCF_000214255.1_Bter_1.0_rna.fna.gz -P $ref_dir
gunzip $ref_dir/*.gz
singularity run -B /Volumes/ docker://quay.io/biocontainers/eggnog-mapper:2.0.1--py_1 download_eggnog_data.py -y -f --data_dir $ref_dir

#make virtual environment in which to run snakemake#
mkdir -p /Volumes/archive/deardenlab/drew_oliphant/projects/ant_v_post/bter/venv
venv_dir=/Volumes/archive/deardenlab/drew_oliphant/projects/ant_v_post/bter/venv
python3 -m venv $venv_dir
source $venv_dir/bin/activate
pip3 install --upgrade pip
pip3 install snakemake
deactivate
