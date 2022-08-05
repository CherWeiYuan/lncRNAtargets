"""
Script to create fasta file from lncRNA gtf
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from gtfparse import read_gtf
import sys
from tqdm import tqdm

# Filter non-exons of lncRNA from gtf
df = read_gtf("genome/gencode.v41.long_noncoding_RNAs.gtf")
df = df[df["feature"]=="exon"]

# Load genome sequence
def parse_genome(genome_fasta):
    """
    Input
        genome_fasta: directory and name of fasta file
    Output
        genome_dict: dictionary of key (chromosome) and value (chromosomal 
                     sequence)
    Note 
        chrom.description is used instead of chrom.id to get full fasta header
        Otherwise, spaces in fasta header will lead to its truncation
    """
    ref_genome = SeqIO.parse(genome_fasta, "fasta")
    genome_dict = {}
    for contig in ref_genome:
        genome_dict[str(contig.name)] = str(contig.seq)
    if genome_dict:
        pass
    else:
        sys.exit("Empty genome dictionary. Check if input fasta is empty.")
    return genome_dict

genome_dict = parse_genome("genome/GRCh38.primary_assembly.genome.fa")

# Get fasta of lncRNA
# Remember to store name of gene in header
# Remember to reverse complement minus strand so all fasta seq in 5' -> 3'
record_list = []
genes_recorded = []
pbar = tqdm(total = len(df))
for index, row in df.iterrows():
    # Get lncRNA metadata
    gene_name = df.loc[index, "gene_name"]
    exon_id = df.loc[index, "exon_id"]
    
    if gene_name not in genes_recorded:
        # Get lncRNA seq
        chrom = df.loc[index, "seqname"]
        start = df.loc[index, "start"]
        end = df.loc[index, "end"]
        seq = genome_dict[chrom][start: end + 1]
        strand = df.loc[index, "strand"]
        exon_num = df.loc[index, "exon_number"]
        # We need lncRNA in 3' -> 5' direction
        if strand == "+":
            record = SeqRecord(Seq(seq), 
                               id = f"{gene_name}_{exon_id}_{exon_num}")
        elif strand == "-":
            record = SeqRecord(Seq(seq).reverse_complement(), 
                               id = f"{gene_name}_{exon_id}_{exon_num}")
        else:
            sys.exit("Unknown strand")
        
    genes_recorded.append(gene_name)
    record_list.append(record)
    pbar.update()

# Write out fasta
with open("data/lncrna_exon_seqs.fasta", "w") as output_handle:
    SeqIO.write(record_list, output_handle, "fasta")    

