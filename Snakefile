#!/usr/bin/env python

"""
Author: K. Chen
Affiliation: NCI/NIH
Aim: Python-based isomiR quantification and analysis pipeline.
Run: snakemake
"""

import glob
import csv
import collections as co
import numpy as np
import pandas as pd
import difflib
import mmap
import time
import logging
from os.path import join
from os.path import splitext
from itertools import takewhile
from decimal import *
from Bio import SeqIO

###############################################################################

TIMESTAMP = time.strftime('%d-%b-%Y@%I:%M:%S%p')
SAMPLES = [os.path.basename(f) for f in glob.glob('data/*')]

###############################################################################

def has_substitution_3p(len_trim, seq_end, consensus_end):
    if len_trim > 3:
        right = min(len(seq_end), len(consensus_end))
        left = max(0, len(seq_end) - len_trim + 1)
        if left < right and all(seq_end[i] == consensus_end[i] or consensus_end[i] == 'N' for i in range(left, right)):
            return True
    return False

def has_substitution_5p(seq_end, consensus_end):
    for i in range(0, min(len(seq_end), len(consensus_end))):
        if (seq_end[len(seq_end) - i - 1] !=
            consensus_end[len(consensus_end) - i - 1]) and consensus_end[len(consensus_end) - i - 1] != 'N':
            return True
    return False

def calc_trimming(seq_end, consensus_end):
    for i in range(0, min(len(seq_end), len(consensus_end))):
        if seq_end[i] != consensus_end[i] and consensus_end[i] != 'N':
            return (len(consensus_end) - i)
    if len(seq_end) < len(consensus_end):
        return len(consensus_end) - len(seq_end)
    return 0

def calc_trimming_5p(seq_end, consensus_end):
    for i in range(0, min(len(seq_end), len(consensus_end))):
        if seq_end[len(seq_end) - i - 1] != consensus_end[len(
                consensus_end) - i - 1] and consensus_end[len(
                consensus_end) - i - 1] != 'N':
            return (len(consensus_end) - i)
    if len(seq_end) < len(consensus_end):
        return len(consensus_end) - len(seq_end)
    return 0

def calc_tailing(seq_end, consensus_end, trim_len):
    return (len(seq_end) - len(consensus_end) + trim_len)

def get_tailing_seq(seq, tail_len):
    if tail_len is not 0:
        return seq[-tail_len:]
    else:
        return '-'

def get_tailing_seq_5p(seq, tail_len):
    if tail_len is not 0:
        return seq[:tail_len]
    else:
        return '-'

def find_in_file(file, string):
    try:
        with open(file, 'r') as f, \
            mmap.mmap(f.fileno(), 0,
                      access=mmap.ACCESS_READ) as s:
            if s.find(bytes(string, 'utf-8')) != -1:
                return True
    except FileNotFoundError:
        return False
    else:
        return False

def motif_consensus_to_dict(file):
    ordered_dict = co.OrderedDict()
    handle = open(file)
    seq_records = list(SeqIO.parse(handle, 'fasta'))
    for seq_record in seq_records:
        mirna = seq_record.description.split()[0]
        motif = seq_record.description.split()[1]
        consensus = str(seq_record.seq)
        if motif in ordered_dict:
            logging.info('Motif ' + motif + ' encountered more than once in ' + file + ': only first occurrence is kept')
            continue
            raise Exception("\n************************************\n" +
                "DUPLICATE MOTIFS FOUND IN '" + file + "'\n" +
                "FIRST INSTANCE: " + motif + "\n" +
                "PLEASE CHECK YOUR MOTIF CONSENSUS FILE\n" +
                "AND FIX OR DELETE BEFORE RERUNNING\n" +
                "************************************\n")
        ordered_dict[motif] = [mirna, consensus]
    return ordered_dict

def other_motifs_pulled_seq(dict_mirna_consensus, line, motif):
    annotation = ""
    list_motifs = list(dict_mirna_consensus.keys())
    matching_motifs = [
        x for x in list_motifs if contains_motif(line, chop_motif(x)) if x != motif]
    for matched_motif in matching_motifs:
        annotation += dict_mirna_consensus[matched_motif][0]
        annotation += " "
    return annotation

def chop_motif(motif):
    n_pos = str.find(motif, 'N')
    if n_pos==-1:
        return([motif])
    n_num = 0
    i = n_pos
    while i<len(motif) and motif[i]=='N':
        i+=1
        n_num+=1
    return([motif[:n_pos], n_num] + chop_motif(motif[n_pos+n_num:]))

def starts_with(line, motif):
    if len(motif) == 1:
        return(str.startswith(line, motif[0]))
    if not str.startswith(line, motif[0]) or len(motif[0])+motif[1] > len(line):
        return(False)
    return(starts_with(line[len(motif[0])+motif[1]:], motif[2:]))

def find_motif(line, motif):
    if len(motif) == 1:
        return(str.find(line, motif[0]))
    motif_pos = str.find(line, motif[0])
    if motif_pos == -1 or motif_pos + len(motif[0]) + motif[1] > len(line):
        return(-1)
    if starts_with(line[motif_pos + len(motif[0]) + motif[1]:], motif[2:]):
        return(motif_pos)
    subline_pos = find_motif(line[motif_pos+1:], motif)
    if subline_pos == -1:
        return(-1)
    return(motif_pos + 1 + subline_pos)
    
def contains_motif(line, motif):
    if find_motif(line, motif) == -1:
        return(False)
    return(True)

def chunkify_file(infilepath, delim=">"):
    with open(infilepath) as infile:
        answer = []
        tinfile = iter(infile)
        while 1:
            try:
                chunk = [next(tinfile)]
                chunk.extend(takewhile(lambda line: not line.startswith(">"), tinfile))
                answer.append(chunk)
            except StopIteration:
                break
    return answer

configfile:
    'config.yaml'

rule all:
    input:
        expand('results/{A}.isomir.tsv', A=SAMPLES),
        expand('results/{A}.isomir.sequence_info.tsv', A=SAMPLES),
        expand('results/{A}.isomir.expression.tsv', A=SAMPLES)

rule collapse_fastq:
    input:
        'data/{A}'
    output:
        'collapsed/{A}.collapsed'
    shell:
        'awk "NR%4==2" {input} | sort -S1900G | uniq -c > {output}'

rule analyze_isomir:
    input:
        motif_consensus = config['motif_consensus_file'],
        collapsed_fasta = 'collapsed/{A}.collapsed',
        input_files = 'data/{A}'
    output:
        'results/{A}.isomir.tsv',
        'results/{A}.isomir.sequence_info.tsv',
        'results/{A}.isomir.nucleotide_dist.tsv'
    log:
        os.path.join("logs/", TIMESTAMP)
    run:
        # SECTION | SETUP LOGGER AND LOG CONFIG ###############################
        logging.basicConfig(
            filename=log[0],
            level=logging.DEBUG,
            format='%(levelname)s: %(message)s')
        if os.stat(log[0]).st_size == 0:  # per-run log
            for key, val in config.items():
                logging.info(str(key) + ':\t' + str(val))
            config_logged = True
        logging.debug("Start sample: " + input.collapsed_fasta)

        # SECTION | SETUP MIRNA CONSENSUS DICT #####################################
        dict_mirna_consensus = motif_consensus_to_dict(input.motif_consensus)

        # with open(str(input.motif_txt), 'rt') as txt:
        #     total_reads_in_sample = txt.readline().split(": ")[1]
        #     total_mapped_reads = 0
        #     for ln in txt:
        #         if ln.startswith("total-reads"):
        #             total_mapped_reads += int(ln.split("\t")[1])

        # SECTION | MIRNA LOOP ################################################
        first_ds = True
        first_dsi = True
        first_dnd = True
        for motif, value in dict_mirna_consensus.items():
            motif_choped = chop_motif(motif)
            mirna = value[0]
            consensus = value[1]
            table_out = []
            freq_nt = co.defaultdict(lambda: co.defaultdict(int))
            logging.debug("Start motif: " + motif + ' mirna: ' + mirna)

        # SECTION | GENERATE SINGLE STATISTICS ################################
            with open(str(input.collapsed_fasta), "rt") as sample:
                # calculate total reads
                pulled_lines = []
                for line in sample:
                    reads = line.rpartition(' ')[0]
                    if contains_motif(line, motif_choped):
                        pulled_lines.append(line)

                for line in pulled_lines:
                    num_reads = int(line.rpartition(' ')[0])
                    seq = line.rpartition(' ')[2].rstrip()

                    # ascertain sequences pulled in by several miRNA motifs
                    has_other = other_motifs_pulled_seq(
                        dict_mirna_consensus, line, motif)

                    # sequence manipulations
                    consensus_index_3p = find_motif(
                        consensus, motif_choped) + len(motif)
                    seq_index_3p = find_motif(seq, motif_choped) + len(motif)
                    consensus_end_3p = consensus[consensus_index_3p:]
                    seq_end_3p = seq[seq_index_3p:]
                    consensus_index_5p = find_motif(consensus, motif_choped)
                    seq_index_5p = find_motif(seq, motif_choped)
                    consensus_end_5p = consensus[:consensus_index_5p]
                    seq_end_5p = seq[:seq_index_5p]

                    # calculation of single-statistics
                    ratio = float(num_reads)
                    len_read = len(seq)
                    len_trim = calc_trimming(seq_end_3p,
                        consensus_end_3p)
                    len_tail = calc_tailing(
                        seq_end_3p, consensus_end_3p, len_trim)
                    seq_tail = get_tailing_seq(seq, len_tail)
                    len_trim_5p = calc_trimming_5p(
                        seq_end_5p, consensus_end_5p)
                    vari_5p = max(len_trim_5p,
                        calc_tailing(
                            seq_end_5p, consensus_end_5p, len_trim_5p))

                    # option for not counting same sequence multiple times
                    if (config['destructive_motif_pull'] and
                       len(has_other) > 0 and  # seq maps to multiple miR
                       find_in_file(output[0], seq)):  # check prev
                           logging.warning(
                           'Skipped ' + seq + ' ' + mirna)
                    elif has_substitution_5p(seq_end_5p, consensus_end_5p):
                        logging.warning(
                        'Skipped (5p substitution) ' + seq + ' ' + mirna)
                    elif has_substitution_3p(len_trim, seq_end_3p,
                        consensus_end_3p):
                        logging.warning(
                        'Skipped (3p sequencing error) ' + seq + ' ' +
                            str(len_trim) + ' ' + mirna)
                    else:
                        # calculation of nt frequencies at each position
                        nt_offset = seq_index_5p - consensus_index_5p
                        for index, nt in enumerate(seq):
                            freq_nt[nt][index - nt_offset] += num_reads
                            for nt2 in ['A', 'C', 'G', 'T']:
                                if nt2!=nt and (freq_nt[nt]==None or freq_nt[nt2][index - nt_offset] == None):
                                    freq_nt[nt2][index - nt_offset] = 0
                        # add to display queue
                        table_out.append([mirna, seq, len_read, num_reads,
                                ratio, len_trim, len_tail, seq_tail,
                                vari_5p, has_other])

    # SECTION | MOVE STATISTICS INTO DATAFRAME ############################
                df = pd.DataFrame(table_out,
                              columns=["MIRNA",
                                       "SEQUENCE",
                                       "LEN_READ",
                                       "READS",
                                       "RATIO",
                                       "LEN_TRIM",
                                       "LEN_TAIL",
                                       "SEQ_TAIL",
                                       "VAR_5P",
                                       "MATCH"])

                if len(table_out) > 0:
                    df.sort_values(by="READS", ascending=0, inplace=1)

                # calculate total reads
                total_reads = float(df['READS'].sum())
                if total_reads == 0:
                    logging.info("Mirna " + mirna + " skipped: no supporting/matched reads")
                    continue

                # calculate ratio
                df['RATIO'] = df['READS'].apply(lambda x: round(100*x/total_reads, 2))
                df2 = pd.DataFrame(freq_nt)
                df2['MIRNA'] = mirna
                if "N" not in df2:
                    df2["N"] = 0.0
                else:
                    df2["N"]=df2["N"].fillna(0.0)
                if len(table_out) > 0:
                    df2['READS'] = df2.sum(axis=1)
                    df2.loc[:, "A":"T"] = df2.loc[
                        :, "A":"T"].div(df2["READS"], axis=0)
                    df2 = np.round(df2, decimals=4)
                    df2.index.name = 'NT_POSITION'
                df2.set_index('MIRNA', append=True, inplace=True)
                df2 = df2.swaplevel(0, 1)
                df2 = df2[["A", "C", "G", "T", "N", "READS"]]

    # SECTION | GENERATE SUMMARY STATISTICS ###############################
                # calculate 5' fidelity score
                vals_vari_5p = df['VAR_5P'].tolist()
                vals_reads = df['READS'].tolist()
                fidelity = 0
                for vari_5p, read in zip(vals_vari_5p, vals_reads):
                    fidelity += (vari_5p * read)
                fidelity = round((fidelity / total_reads), 4)

                # calculate individual nt tailing ratios
                vals_tail_seq = df['SEQ_TAIL'].tolist()
                array_nt = [0, 0, 0, 0]  # A C G T
                for seq_tail, read in zip(vals_tail_seq, vals_reads):
                    count_nt = co.Counter(seq_tail)
                    array_nt[0] += (count_nt['A'] * read)
                    array_nt[1] += (count_nt['C'] * read)
                    array_nt[2] += (count_nt['G'] * read)
                    array_nt[3] += (count_nt['T'] * read)
                total_nt_tail = sum(array_nt)
                ratio_a_tailing = "nan"
                ratio_c_tailing = "nan"
                ratio_g_tailing = "nan"
                ratio_t_tailing = "nan"
                if total_nt_tail > 0:
                    ratio_a_tailing = str(100*round(array_nt[0] / total_nt_tail, 4))
                    ratio_c_tailing = str(100*round(array_nt[1] / total_nt_tail, 4))
                    ratio_g_tailing = str(100*round(array_nt[2] / total_nt_tail, 4))
                    ratio_t_tailing = str(100*round(array_nt[3] / total_nt_tail, 4))

                # calculate ratios of trimmed/tailed sequences
                ratio_seq_trim_only = 0
                ratio_seq_tail_only = 0
                ratio_seq_trim_and_tail = 0
                vals_len_trim = df['LEN_TRIM'].tolist()
                vals_len_tail = df['LEN_TAIL'].tolist()
                for len_trim, len_tail, read in zip(vals_len_trim, vals_len_tail, vals_reads):
                    read_ratio = 100*round(read / total_reads, 4)
                    if len_trim > 0 and len_tail == 0:
                        ratio_seq_trim_only += read_ratio
                    if len_tail > 0 and len_trim == 0:
                        ratio_seq_tail_only += read_ratio
                    if len_trim > 0 and len_tail > 0:
                        ratio_seq_trim_and_tail += read_ratio

    # SECTION | FILTER OUT DISPLAYED READS AND ADD % SIGN
                # only show reads above threshold
                df = df[(df.RATIO > config['min_ratio']) |
                    (df.READS > config['min_read'])]
                df['RATIO'] = df['RATIO'].apply(lambda x: str(x) + '%')

                # calculate number of isomirs above threshold
                total_isomirs = df.shape[0]
                if (df['LEN_TRIM'].iloc[0] == 0 and df['LEN_TAIL'].iloc[0] == 0):
                    total_isomirs -= 1

    # SECTION | DISPLAY HEADER AND SUMMARY STATISTICS #################
                if config['display_summary']:
                    with open(output[0], 'a') as out:
                        total_reads_in_sample = sum(1 for line in open(input.input_files))/4
                        summary_out = [[mirna,
                                        motif,
                                        consensus,
                                        str(int(total_reads)),
                                        str(total_isomirs),
                                        str(fidelity),
                                        ratio_a_tailing,
                                        ratio_c_tailing,
                                        ratio_g_tailing,
                                        ratio_t_tailing,
                                        str(ratio_seq_trim_only),
                                        str(ratio_seq_trim_only + ratio_seq_trim_and_tail),
                                        str(ratio_seq_tail_only),
                                        str(ratio_seq_tail_only + ratio_seq_trim_and_tail),
                                        str(ratio_seq_trim_and_tail),
                                        str(int(total_reads_in_sample))]]
                        df3 = pd.DataFrame(summary_out,
                                           columns = ["MIRNA",
                                                      "MOTIF",
                                                      "CONSENSUS",
                                                      "TOTAL READS",
                                                      "TOTAL ISOMIRS",
                                                      "FIDELITY 5P",
                                                      "A TAILING",
                                                      "C TAILING",
                                                      "G TAILING",
                                                      "T TAILING",
                                                      "SEQUENCE TRIMMING ONLY",
                                                      "SEQUENCE TRIMMING",
                                                      "SEQUENCE TAILING ONLY",
                                                      "SEQUENCE TAILING",
                                                      "SEQUENCE TRIMMING AND TAILING",
                                                      "TOTAL READS IN SAMPLE"])
                        if first_ds:
                            df3.to_csv(out, sep='\t', index = False)
                            first_ds = False
                        else:
                            df3.to_csv(out, sep='\t', index=False, header = False)
                else:
                    open(output[0], 'a').close()
    # SECTION | DISPLAY SEQUENCES AND SINGLE STATISTICS ###############
                if config['display_sequence_info']:
                    with open(output[1], 'a') as out:
                        if first_dsi:
                            df.to_csv(out, sep='\t', index = False)
                            first_dsi = False
                        else:
                            df.to_csv(out, sep='\t', index=False, header = False)
                else:
                    open(output[1], 'a').close()
                if config['display_nucleotide_dist']:
                    with open(output[2], 'a') as out:
                        if first_dnd:
                            df2.to_csv(out, sep='\t')
                            first_dnd = False
                        else:
                            df2.to_csv(out, sep='\t', header = False)
                else:
                    open(output[2], 'a').close()

rule cpm_normalize_motifs:
    input:
        'results/{A}.isomir.sequence_info.tsv',
        'results/{A}.isomir.tsv'
    output:
        'results/{A}.isomir.expression.tsv'
    run:
        first_exp = True
        if config['display_sequence_info'] and config['display_summary']:
            df = pd.read_csv(input[0], delimiter="\t", header = 0)
            total_reads = pd.read_csv(input[1], delimiter="\t", header = 0, nrows = 2)['TOTAL READS IN SAMPLE'].iloc[0]
            df['CPM'] = df['READS'] / df['READS'].sum() * float(10^6) # TPM is also len-norm
            df['RPKM'] = df['READS'] / df['LEN_READ'] / total_reads * float(10^9)
            with open(output[0], 'a') as out:
                df.to_csv(out, sep='\t', index=False, header=True)
        else:
            open(output[0], 'a').close()
