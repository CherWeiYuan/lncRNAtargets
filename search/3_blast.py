"""
Script to blast monogenic disease genes against lncRNA exons
"""

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from datetime import datetime
from sys import exit

def date_time():
    """
    Returns dd/mm/YY H:M:S
    """
    now = datetime.now().strftime("DD%dMM%mYY%Y,%H:%M:%S")
    return now

def run_blast(query_seq, blast_db, threads, strand = "minus", evalue = 0.05, 
              word_size = 11):
    """
    Run blast client and output XML file in temp folder

    Blast default is minus strand since all mRNA input is plus strand 5' -> 3',
    and all lncRNA exon is plus strand 5' -> 3', thereby if minus strand of 
    mRNA aligns to lncRNA exon plus strand, it means lncRNA exon binds to
    mRNA plus strand
    """
    # Determine which blast task to use by length of query seq
    if len(query_seq) <= 30:
        blast_mode = "blastn-short"
    else:
        blast_mode = "blastn"

    # BLAST set-up
    # Note max_target_seqs is 500 by default
    time = date_time()
    cline = NcbiblastnCommandline(
        db = blast_db,
        word_size = word_size,      # Small word size since probe is short
        task = blast_mode,          # BLASTN/ BLASTN short used for sensitivity
        strand = strand,            # Query strand(s) to search against database
        out = f"blast_temp/{time}blastn.xml",  # Output file name
        outfmt = 5,                 # Output as XML format (outfmt=5)
        max_hsps = 10,              # 
        evalue = evalue,              
        num_threads = threads)      # Number of threads to use

    # Run BLAST which output XML
    stdout, stderr = cline(stdin = query_seq)
    return f"blast_temp/{time}blastn.xml"

def parse_blast(probe_name, probe_df, xml,
                low_threshold_id, low_threshold_len,
                mid_threshold_id, mid_threshold_len,
                high_threshold_id, high_threshold_len):
    """
    Reads blast XML output in temp folder and store alignment information in 
    dataframe
    
    Parameters
    ----------
    probe_name : TYPE str. Probe name
    probe_df : TYPE Pandas dataframe. Contains at least two columns 'Seq' and 
               'Name' of probes in 5' -> 3' direction. 
               IMPORTANT: Probe name must be unique if not programme will stop
    xml : TYPE str. XML output from run_blast() function.
    {low/ mid/ high}_threshold_id : TYPE int. Minimum percent identity to  
                                    consider alignment as a valid hit
    {low/ mid/ high}_threshold_len : TYPE int. Minimum alignment length to  
                                     consider alignment as a valid hit
    Returns
    -------
    hsps_list : TYPE list. List of [chromosome, subject start, subject end, 
                                    alignment length, percent identity]
    
    Notes
    -----
    1. All positions of blast outputs are 1-based but all
       positions in alignments are converted to 0-based
    2. All stored positions from blast in alignments refer to 
       the same strand of query/ subject sequence (plus strand; 
       thereby bed file must also give coordinates using plus strand)
    3. BLAST hits (in terms of HSPs) exceeding high threshold id and len will
       have their blast output statistics stored for future reference:
       probe name, alignment title, 
       HSP subject start, HSP subject end,
       HSP alignment length, percent identity
    """
    # Check if probe names are unique (no duplicates)
    unique_probe_names = list(probe_df['Name'])
    if len(unique_probe_names) != len(probe_df.index):
        exit('ERROR: Probe name in dataframe is not unique')

    # Store alignment information if off-target count is below threshold_len
    result_handle = open(xml)
    blast_records = NCBIXML.read(result_handle)
    hsps_list = []
    num_alignment_high = 0
    num_alignment_mid = 0
    num_alignment_low = 0
    num_alignment_stratified = 0
    for alignment in blast_records.alignments:
        for hsp in alignment.hsps:
            align_length = hsp.align_length
            percent_id = (hsp.identities/ align_length)*100
            if percent_id >= high_threshold_id and \
               align_length >= high_threshold_len:
                # Note all blast output coordinates are 1-based but converted
                # to 0-based
                # Note structure of hsp_list: 
                # [probe name, chromosome, start, end, 
                # alignment length, percent identity]
                num_alignment_high += 1
                if hsp.sbjct_start-1 <= hsp.sbjct_end-1:
                    hsps_list.append([probe_name, alignment.title, 
                                      hsp.sbjct_start-1, hsp.sbjct_end-1,
                                      align_length, percent_id]) 
                elif hsp.sbjct_start-1 > hsp.sbjct_end-1:
                    hsps_list.append([probe_name, alignment.title, 
                                      hsp.sbjct_end-1, hsp.sbjct_start-1,
                                      align_length, percent_id]) 
                else:
                    exit('BLAST ERROR: Unexpected subject strand')
            # The parameters below generate results close to IDT's algorithm
            # for reference
            if percent_id >= mid_threshold_id and \
               align_length >= mid_threshold_len:
                num_alignment_mid += 1
            # Lenient thresholds help detect partial alignment to 
            # repetitive elements
            if percent_id >= low_threshold_id and \
               align_length >= low_threshold_len:
                num_alignment_low += 1   
            # Extremely lenient and comprehensive partial alignments
            if percent_id >= 50 and align_length >= 100:
                num_alignment_stratified += 1
            elif percent_id >= 60 and align_length >= 80:
                num_alignment_stratified += 1
            elif percent_id >= 80 and align_length >= 60:
                num_alignment_stratified += 1
            elif percent_id >= 50 and align_length >= 95:
                num_alignment_stratified += 1

    # Find index of unique probe name
    probe_index = probe_df.index[probe_df["Name"] == probe_name].tolist()[0]
    probe_df.loc[probe_index, "hitCountHighThreshold"] = int(num_alignment_high)
    probe_df.loc[probe_index, "hitCountMidThreshold"] = int(num_alignment_mid)
    probe_df.loc[probe_index, "hitCountLowThreshold"] = int(num_alignment_low)
    probe_df.loc[probe_index, "hitCountStratified"] = int(num_alignment_stratified)

    # Remove duplicate lists
    hsps_list = set(map(tuple, hsps_list))
    hsps_list = list(map(list, hsps_list))
    return (probe_df, hsps_list, num_alignment_high)
