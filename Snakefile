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
import mmap
import time
import logging
from os.path import join
from os.path import splitext
from decimal import *
from Bio import SeqIO

###############################################################################

MOTIF_CONSENSUS = 'motif-consensus.fa'

SAMPLES = [os.path.basename(f) for f in glob.glob('data/*.fastq_ready')]

TIMESTR = time.strftime("%Y%m%d-%H%M%S")

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

    # def matches(large_string, query_string, threshold):
    #    words = large_string.split()
    #    for word in words:
    #        s = difflib.SequenceMatcher(None, word, query_string)
    #        match = ''.join(word[i:i + n]
    #                        for i, j, n in s.get_matching_blocks() if n)
    #      if len(match) / float(len(query_string)) >= threshold:
    #            yield match

    ###################################################################

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
    log:
        "logs/{A}%s", TIMESTR
    run:
        # SECTION | DISPLAY CONFIGURATION OPTIONS #########################
        with open(output[0], 'a') as out:
            out.write('QuagmiR-generated@' + TIMESTR + '\n')
            for key, val in config.items():
                out.write(str(key) + ':\t' + str(val) + '\n')
            out.write('\n')

        # SECTION | LOG CONFIG ################################################
        logging.basicConfig(
            filename=log[0],
            level=logging.DEBUG,
            format='%(asctime)s %(message)s')
        logging.info('----')
        for key, val in config.items():
            logging.info(key)
            logging.info(val)

        # SECTION | SETUP MIRNA INFO DICT #####################################
        dict_mirna = co.OrderedDict()
        handle = open(input.motif_consensus)
        seq_records = list(SeqIO.parse(handle, 'fasta'))
        for seq_record in seq_records:
            mirna = seq_record.description.split()[0]
            motif = seq_record.description.split()[1]
            consensus = str(seq_record.seq)
            dict_mirna[motif] = [mirna, consensus]

        # SECTION | MIRNA LOOP ################################################
        for motif, value in dict_mirna.items():
            mirna = value[0]
            consensus = value[1]
            table_out = []
            freq_nt = co.defaultdict(lambda: co.defaultdict(int))

        # SECTION | GENERATE SINGLE STATISTICS ################################
            with open(str(input.collapsed_fasta), "rt") as sample:
                # calculate total reads
                total_reads = 0
                pulled_lines = []
                for line in sample:
                    reads = line.rpartition(' ')[0]
#                    if config['fuzzy_motif']:  # if fuzzy motif matching
#                        if len(list(matches(motif, line, 0.8))) > 0:
#                            motif_con = list(matches(motif, line, 0.8))[0]
#                            total_reads += int(reads)
#                            total_sequences += 1
                    if motif in line:
                        pulled_lines.append(line)
                        total_reads += int(reads)

                if total_reads > 0:
                    for line in pulled_lines:
                        num_reads = int(line.rpartition(' ')[0])
                        seq = line.rpartition(' ')[2].rstrip()

                        # ascertain sequences pulled in by several miRNA motifs
                        annotation = ""
                        list_motifs = list(dict_mirna.keys())
                        matching_motifs = [
                            x for x in list_motifs if x in line if x != motif]
                        for matched_motif in matching_motifs:
                            annotation += dict_mirna[matched_motif][0]
                            annotation += " "

                        # option for not counting same sequence multiple times
                        if config['destructive_motif_pull']:
                            if len(annotation) > 0:  # maps to multiple miR
                                if find_in_file(output[0], seq):  # check prev

                                    continue

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

                        # calculation of nt frequencies at each position
                        nt_offset = seq_index_5p - consensus_index_5p
                        for index, nt in enumerate(seq):
                            freq_nt[nt][index - nt_offset] += num_reads

                        # if sequence has minimum reads, add to display queue
                        if ((ratio > config['min_ratio']) or
                                (num_reads > config['min_read'])):
                            table_out.append([seq, len_read, num_reads, ratio,
                                              len_trim, len_tail, seq_tail,
                                              vari_5p, annotation])

        # SECTION | MOVE STATISTICS INTO DATAFRAME ############################
                    df = pd.DataFrame(table_out,
                                      columns=["SEQUENCE",
                                               "LEN_READ",
                                               "READS",
                                               "RATIO",
                                               "LEN_TRIM",
                                               "LEN_TAIL",
                                               "SEQ_TAIL",
                                               "VAR_5P",
                                               "MATCH"])
                    df.sort_values(by="READS", ascending=0, inplace=1)

                    df2 = pd.DataFrame(freq_nt).fillna(value=0)
                    df2['READS'] = df2.sum(axis=1)
                    df2.loc[:, "A":"T"] = df2.loc[
                        :, "A":"T"].div(df2["READS"], axis=0)
                    df2 = np.round(df2, decimals=4)
                    df2.index.name = 'NT_POSITION'

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
                    array_nt = [0, 0, 0, 0]  # A T C G
                    for seq_tail, read in zip(vals_tail_seq, vals_reads):
                        count_nt = co.Counter(seq_tail)
                        array_nt[0] += (count_nt['A'] * read)
                        array_nt[1] += (count_nt['C'] * read)
                        array_nt[2] += (count_nt['G'] * read)
                        array_nt[3] += (count_nt['T'] * read)
                    total_nt_tail = sum(array_nt)
                    ratio_nt_tailing = "nan"
                    if total_nt_tail > 0:
                        ratio_nt_tailing = (
                            str(round(array_nt[0] / total_nt_tail, 4)) +
                            '\t' + str(round(array_nt[1] / total_nt_tail, 4)) +
                            '\t' + str(round(array_nt[2] / total_nt_tail, 4)) +
                            '\t' + str(round(array_nt[3] / total_nt_tail, 4)))

                    # calculate ratios of trimmed/tailed sequences
                    ratio_seq_trim = 0
                    ratio_seq_tail = 0
                    vals_len_trim = df['LEN_TRIM'].tolist()
                    vals_len_tail = df['LEN_TAIL'].tolist()
                    for len_trim, read in zip(vals_len_trim, vals_reads):
                        if len_trim > 0:
                            ratio_seq_trim += round(read / total_reads, 4)
                    for len_tail, read in zip(vals_len_tail, vals_reads):
                        if len_tail > 0:
                            ratio_seq_tail += round(read / total_reads, 4)

        # SECTION | DISPLAY HEADER AND SUMMARY STATISTICS #################
                    with open(output[0], 'a') as out:
                        if config['display_summary']:
                            out.write('==========================\n' +
                                      mirna + '\n')
                            out.write('**Motif:\t' + motif + '\n')
                            out.write('**Consensus:\t' + consensus + '\n')
                            out.write('**Total-Reads:\t' +
                                      str(total_reads) + '\n')
                            out.write('**Fidelity-5p:\t' +
                                      str(fidelity) + '\n')
                            out.write('**ACGT-Tailing:\t' +
                                      ratio_nt_tailing + '\n')
                            out.write('**Sequence-Trimming:\t' +
                                      str(ratio_seq_trim) + '\n')
                            out.write('**Sequence-Tailing:\t' +
                                      str(ratio_seq_tail) + '\n')

        # SECTION | DISPLAY SEQUENCES AND SINGLE STATISTICS ###############
                        if config['display_sequence_info']:
                            out.write('\n[Sequence-Information]\n')
                            df.to_csv(out, sep='\t', index=False)
                        if config['display_nucleotide_dist']:
                            out.write('\n[Nucleotide-Distribution]\n')
                            df2.to_csv(out, sep='\t')
