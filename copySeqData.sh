#!/bin/bash


cd $SLURM_SUBMIT_DIR

# Copy fastq sequences to current directory
mkdir rawData
cd rawData

find /rfs/VasCeGenS_RNAseq/Shared/HUVEC_RNASeq_100_CambCollab/ -name "*.fastq.gz" -type f -exec rsync -avi --update --ignore-existing --progress {} . \; &



find . -type f -exec chmod 0644 {} \;