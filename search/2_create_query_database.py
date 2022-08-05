"""
Script to create query gene database from the monogenic disease compendium
"""

from Bio import Seq
from Bio import SeqIO
from gtfparse import read_gtf
from numpy import nan
import pandas as pd
import sys

# Get majority genes
df = pd.read_csv("data/disease_compendium.csv", encoding = "ISO-8859-1")
major_genes = df.majority_gene.to_list()

# Prepare gtf for majority genes
gtf = read_gtf("genome/gencode.v41.annotation.gtf")
gtf = gtf[gtf.feature == "gene"]

keep_index = []
for index, row in gtf.iterrows():
    gene_name = gtf.loc[index, "gene_name"]
    if gene_name in major_genes:
        keep_index.append(index)
gtf = gtf.loc[keep_index]

# Update gtf with seq
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

gtf["seq"] = nan
for index, row in gtf.iterrows():
    chrom = gtf.loc[index, "seqname"]
    start = gtf.loc[index, "start"]
    end = gtf.loc[index, "end"]
    strand = gtf.loc[index, "strand"]
    if strand == "+":
        seq = genome_dict[chrom][start: end + 1]
    elif strand == "-":
        seq = str(Seq(genome_dict[chrom][start: end + 1]).reverse_complement())
    else:
        sys.exit("Unknown strand")
    gtf.loc[index, "seq"] = seq

# Merge gtf with disease compendium
gtf = gtf.rename(columns = {"gene_name": "majority_gene"})
df = pd.merge(gtf, df, on = "majority_gene", how = "left")

df.to_csv("data/monogenic_disease_genes.csv", index = False)


