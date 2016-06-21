#!/usr/bin/env python

""" A pipeline to analyze
raw miRNA-seq files for
isomiR composition.
Currently takes in .fastq_ready
files in /data folder and motif.txt"""

import glob
import csv
import collections as co
import numpy as np
import pandas as pd
import difflib
from os.path import join
from os.path import splitext
from decimal import *
from Bio import SeqIO

###############################################################################

MOTIF_CONSENSUS = 'motif-consensus.fa'

SAMPLES = [os.path.basename(f) for f in glob.glob('data/*.fastq_ready')]

###############################################################################


def calc_trimming_5p(seq_end, consensus_end, seq):
    for i in range(0, min(len(seq_end), len(consensus_end))):
        if seq_end[len(seq_end) - i - 1] != consensus_end[len(
                consensus_end) - i - 1]:
            return (len(consensus_end) - i)
    if len(seq_end) < len(consensus_end):
        return len(consensus_end) - len(seq_end)
    return 0


def calc_trimming(seq_end, consensus_end):
    for i in range(0, min(len(seq_end), len(consensus_end))):
        if seq_end[i] != consensus_end[i]:
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


def matches(large_string, query_string, threshold):
    words = large_string.split()
    for word in words:
        s = difflib.SequenceMatcher(None, word, query_string)
        match = ''.join(word[i:i + n]
                        for i, j, n in s.get_matching_blocks() if n)
        if len(match) / float(len(query_string)) >= threshold:
            yield match

##########################################################################

configfile:
    "config.yaml"

rule all:
    input:
        expand("results/{A}.results.txt", A=SAMPLES)

rule collapse_fastq:
    input:
        'data/{A}'
    output:
        'data/collapsed/{A}.collapsed'
    shell:
        "awk 'NR%4==2' {input} | sort -S1900G | uniq -c > {output}"

rule analyze_isomir:
    input:
        motif_consensus = MOTIF_CONSENSUS,
        collapsed_fasta = 'data/collapsed/{A}.collapsed'
    output:
        "results/{A}.results.txt"
    run:
        # SECTION | SETUP #####################################################
        # open input miRNA info file and input info into dict
        mirna_dict = {}
        handle = open('motif-consensus.fa')
        seq_records = list(SeqIO.parse(handle, 'fasta'))
        for seq_record in seq_records:
            mirna = seq_record.description.split()[0]
            motif = seq_record.description.split()[1]
            consensus = str(seq_record.seq)
            mirna_dict[motif] = [mirna, consensus]

        # SECTION | MIRNA LOOP ################################################
        for motif, value in mirna_dict.items():
            mirna = value[0]
            consensus = value[1]
            table_out = []
            freq_nt = co.defaultdict(lambda: co.defaultdict(int))

        # SECTION | GENERATE SINGLE STATISTICS ################################
            with open(str(input.collapsed_fasta), "rt") as sample:
                # calculate total reads
                total_reads = 0
                total_sequences = 0
                for line in sample:
                    reads = line.rpartition(' ')[0]
#                    if config['fuzzy_motif']:  # if fuzzy motif matching
#                        if len(list(matches(motif, line, 0.8))) > 0:
#                            motif_con = list(matches(motif, line, 0.8))[0]
#                            total_reads += int(reads)
#                            total_sequences += 1
                    if motif in line:  # regular motif pull
                        total_reads += int(reads)
                        total_sequences += 1

                if total_reads > 0:
                    sample.seek(0)
                    gen = (z for z in sample if motif in z)
                    for line in gen:
                        num_reads = int(line.rpartition(' ')[0])
                        seq = line.rpartition(' ')[2].rstrip()

                        # sequence manipulations
                        consensus_index_3p = str.find(
                            consensus, motif) + len(motif)
                        seq_index_3p = str.find(seq, motif) + len(motif)
                        consensus_end_3p = consensus[consensus_index_3p:]
                        seq_end_3p = seq[seq_index_3p:]
                        consensus_index_5p = str.find(
                            consensus, motif)
                        seq_index_5p = str.find(seq, motif)
                        consensus_end_5p = consensus[
                            :consensus_index_5p]
                        seq_end_5p = seq[:seq_index_5p]

                        # calculation of single-statistics
                        ratio = round(
                            float(num_reads) / float(total_reads), 4)
                        len_read = len(seq)
                        len_trim = calc_trimming(seq_end_3p, consensus_end_3p)
                        len_tail = calc_tailing(
                            seq_end_3p, consensus_end_3p, len_trim)
                        seq_tail = get_tailing_seq(seq, len_tail)
                        len_trim_5p = calc_trimming_5p(
                            seq_end_5p, consensus_end_5p, seq)
                        vari_5p = max(
                            len_trim_5p,
                            calc_tailing(
                                seq_end_5p, consensus_end_5p, len_trim_5p))
                        annotation = ""

                        # calculation of nt frequencies at each position
                        nt_offset = seq_index_5p - consensus_index_5p
                        for index, nt in enumerate(seq):
                            freq_nt[nt][index - nt_offset] += num_reads

                        # ascertain sequences pulled in by several miRNA motifs
                        list_motifs = list(mirna_dict.keys())
                        matching_motifs = [
                            x for x in list_motifs if x in line if x != motif]
                        for matched_motif in matching_motifs:
                            annotation += mirna_dict[matched_motif][0]
                            annotation += " "

                        if ((ratio > config['min_ratio']) or
                                (num_reads > config['min_read'])):
                            table_out.append([seq, len_read, num_reads, ratio,
                                              len_trim, len_tail, seq_tail,
                                              vari_5p, annotation])

        # SECTION | MOVE STATISTICS INTO DATAFRAME ############################
                    df = pd.DataFrame(table_out,
                                      columns=["SEQUENCE",
                                               "READ_LENGTH",
                                               "MIRNA_READS",
                                               "ABUNDANCE_RATIO",
                                               "TRIM_LENGTH",
                                               "TAIL_LENGTH",
                                               "TAIL_SEQUENCE",
                                               "VARIATION_5P",
                                               "ANNOTATION"])
                    df.sort_values(by="MIRNA_READS", ascending=0, inplace=1)

                    df2 = pd.DataFrame(freq_nt).fillna(value=0)
                    df2['READS'] = df2.sum(axis=1)
                    df2.loc[:, "A":"T"] = df2.loc[
                        :, "A":"T"].div(df2["READS"], axis=0)
                    df2 = np.round(df2, decimals=4)
                    df2.index.name = 'NT_POSITION'

        # SECTION | GENERATE SUMMARY STATISTICS ###############################
                    # calculate 5' fidelity score
                    vals_vari_5p = df['VARIATION_5P'].tolist()
                    vals_reads = df['MIRNA_READS'].tolist()
                    fidelity = 0
                    for vari_5p, read in zip(vals_vari_5p, vals_reads):
                        fidelity += (vari_5p * read)
                    fidelity = round((fidelity / total_reads), 4)

                    # calculate individual nt tailing ratios
                    vals_tail_seq = df['TAIL_SEQUENCE'].tolist()
                    array_nt = [0, 0, 0, 0]  # A T C G
                    for seq_tail, read in zip(vals_tail_seq, vals_reads):
                        count_nt = co.Counter(seq_tail)
                        array_nt[0] += (count_nt['A'] * read)
                        array_nt[1] += (count_nt['T'] * read)
                        array_nt[2] += (count_nt['C'] * read)
                        array_nt[3] += (count_nt['G'] * read)
                    total_nt_tailing = sum(array_nt)
                    ratio_nt_tailing = (
                        'A:' + str(round(array_nt[0] / total_nt_tailing, 4)) +
                        '/T:' + str(round(array_nt[1] / total_nt_tailing, 4)) +
                        '/C:' + str(round(array_nt[2] / total_nt_tailing, 4)) +
                        '/G:' + str(round(array_nt[3] / total_nt_tailing, 4)))

                    # calculate ratios of trimmed/tailed sequences
                    vals_len_trim = df['TRIM_LENGTH'].tolist()
                    vals_len_tail = df['TAIL_LENGTH'].tolist()
                    ratio_seq_trim = round((len(
                        [x for x in vals_len_trim if x != 0]) * num_reads) /
                        total_reads, 4)
                    ratio_seq_tail = round((len(
                        [x for x in vals_len_tail if x != 0]) * num_reads) /
                        total_reads, 4)

        # SECTION | DISPLAY HEADER AND SUMMARY STATISTICS #################
                    with open(output[0], 'a') as out:
                        out.write('==========================\n' +
                                  mirna + '\n')
                        out.write('**Motif: ' + motif + '\n')
                        out.write('**Consensus: ' + consensus + '\n')
                        out.write('**Total-Reads: ' + str(total_reads) + '\n')
                        out.write('**Fidelity-5p: ' + str(fidelity) + '\n')
                        out.write('**NT-Tailing: ' + ratio_nt_tailing + '\n')
                        out.write('**Sequence-Trimming: ' +
                                  str(ratio_seq_trim) + '\n')
                        out.write('**Sequence-Tailing: ' +
                                  str(ratio_seq_tail) + '\n')

        # SECTION | DISPLAY SEQUENCES AND SINGLE STATISTICS ###############
                        out.write('\n[Sequence-Information]\n')
                        df.to_csv(out, sep='\t', index=False)
                        out.write('\n[Nucleotide-Distribution]\n')
                        df2.to_csv(out, sep='\t')
