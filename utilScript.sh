#!/bin/bash

# Create necessary directories
mkdir -p rawReads
mkdir -p fastqcResults
mkdir -p mappedReads
mkdir -p salmonQuant
mkdir -p readCount

# Copy analyses scripts to the current directory
# Check if there's any folder with "_Scripts" and remove it
if [ -n "$(find . -maxdepth 1 -type d -name '*_Scripts' -print -quit)" ]; then
    rm -r *_Scripts
fi


cp -r /scratch/vasccell/cs806/huvecRNASEQ/huvecRNASEQ_Scripts/ .

# Rename script directory
wdir1=$PWD
wdir2=$(basename $PWD)
SF="${wdir2}_Scripts"
mv huvecRNASEQ_Scripts $SF


# Update file paths in scripts
todel="/scratch/vasccell/cs806/huvecRNASEQ"
cd $SF
sed -i "s|$todel|$wdir1|g" *
sed -i "s|huvecRNASEQ_Scripts|$SF|g" *


# Check if raw reads are available
echo "Are the raw reads available in the rawReads directory?"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) break;;
        No ) exit;;
    esac
done

# Edit and submit jobs
gedit salmonQuant.job
sbatch salmonQuant.job

bash qc_starmap.sh


# TO DO ####
# - Archive .cram to /rfs with project name
