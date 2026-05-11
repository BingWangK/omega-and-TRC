#!/bin/bash

CONDA_PATH=$(conda info | grep -i 'base environment' | awk '{print $4}')
source $CONDA_PATH/etc/profile.d/conda.sh


mkdir 1_DEMULTIPLEXING
conda activate bcl2fastq2
mv ./*NB* ./0_RUNFOLDER
cd ./0_RUNFOLDER
bcl2fastq \
-R ./ \
-o ../ \
--minimum-trimmed-read-length 19 \
--mask-short-adapter-reads 19 \
--no-lane-splitting
cd ../
mv Reports Stats 1_DEMULTIPLEXING
mv 1_DEMULTIPLEXING 0_RUNFOLDER
mv Undetermined*.fastq.gz 0_RUNFOLDER

filenames_R1=($(find ./ -maxdepth 1 -type f -name "*R1*.fastq.gz" | while read F; do basename $F | grep -oP '.*(?=_S\d+)'; done | sort | uniq))
for i in "${filenames_R1[@]}"; do mv "$i"*R1*.fastq.gz "$i".fastq.gz; done