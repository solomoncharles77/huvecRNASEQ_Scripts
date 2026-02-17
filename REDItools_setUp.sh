
# # installation
# git clone https://github.com/BioinfoUNIBA/REDItools
# cd REDItools
# module purge
# ## conda create -n py27 python=2.7 anaconda
# conda activate py27
# python setup.py install
# cd ..


module purge
source ~/miniforge3/bin/activate py27


# Format REDIportal known sites for REDItools
zcat editDB/TABLE1_hg38_v3.txt.gz | awk '{OFS ="\t"; sub(/^chr/, "", $2)} {print $2,$3,$6}' > temp.tab
tail -n +2 temp.tab > editDB/knownEditingSites_TABLE1_hg38_v3.tab
rm temp.tab

REDItoolKnown.py -i testBAM/A1070.sorted.bam \
  -f ../exprPhenoData/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
  -l editDB/knownEditingSites_TABLE1_hg38_v3.tab.gz \
  -o rediKnownTest

# Run REDItoolsKnown
while IFS= read -r sample; do
        echo ${sample}
        REDItoolKnown.py -i bamFiles/${sample}.sorted.bam \
        -f /scratch/vasccell/cs806/exprPhenoData/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
        -l editDB/knownEditingSites_TABLE1_hg38_v3.tab.gz \
        -o rediKnownTest/${sample} \
        -t 8
        
done < /scratch/vasccell/cs806/rnaEditing/rnaEditing_scripts/baiList.txt
