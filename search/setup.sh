# Set up
conda create --name blast -y

conda activate blast
mamba install -c bioconda blast=2.10.1 -y
conda deactivate

# Set up blast database
conda activate blast
dir=/mnt/e/lncrna/
mkdir -p $dir/blast_db
makeblastdb \
-in $dir/data/lncrna_exon_seqs.fasta \
-input_type fasta \
-out $dir/blast_db/lncrna_exons \
-title lncrna_exons \
-parse_seqids \
-blastdb_version 5 \
-logfile $dir/blast_db/lncrna_exons.log \
-dbtype nucl
conda deactivate