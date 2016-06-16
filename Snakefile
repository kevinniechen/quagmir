#!/usr/bin/env python

""" A pipeline to analyze
raw miRNA-seq files for
isomiR composition.
Currently takes in .fastq_ready
files in /data folder and motif.txt"""

import glob
import csv
import numpy as np
import pandas as pd
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

    ##########################################################################

rule all:
    input:
        expand("results/{A}.results.txt", A=SAMPLES)

rule collapse_fastq:
    input:
        'data/{A}'
    output:
        'data/{A}.collapsed'
    shell:
        "awk 'NR%4==2' {input} | sort -S1900G | uniq -c > {output}"

rule analyze_isomir:
    input:
        motif_consensus = MOTIF_CONSENSUS,
        collapsed_fasta = 'data/{A}.collapsed'
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

        # SECTION | GENERATE SINGLE STATISTICS ################################
            with open(str(input.collapsed_fasta), "rt") as sample:
                # calculate total reads
                total_reads = 0
                for line in sample:
                    reads = line.rpartition(' ')[0]
                    if motif in line:
                        total_reads += int(reads)

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
                        percent = round(
                            100 * float(num_reads) / float(total_reads), 2)
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

                        # ascertain sequences pulled in by several miRNA motifs
                        list_motifs = list(mirna_dict.keys())
                        matching_motifs = [
                            x for x in list_motifs if x in line if x != motif]
                        for matched_motif in matching_motifs:
                            annotation += mirna_dict[matched_motif][0]
                            annotation += " "

                        table_out.append([seq, num_reads, percent, len_read,
                                          len_trim, len_tail, seq_tail,
                                          vari_5p, annotation])
                # create pandas dataframe containing output
                df = pd.DataFrame(table_out,
                                  columns=["SEQUENCE",
                                           "MIRNA_READS",
                                           "PERCENT",
                                           "READ_LENGTH",
                                           "TRIM_LENGTH",
                                           "TAIL_LENGTH",
                                           "TAIL_SEQUENCE",
                                           "VARIATION_5P",
                                           "ANNOTATION"])

        # SECTION | GENERATE SUMMARY STATISTICS ###############################
                # calculate 5' fidelity score
                vari_5p_vals = [table_out[i][7] for i in range(len(table_out))]
                percent_vals = [table_out[i][2] for i in range(len(table_out))]
                fidelity = 0
                for vari_5p, percent in zip(vari_5p_vals, percent_vals):
                    fidelity += round((vari_5p * percent), 2)

        # SECTION | DISPLAY HEADER AND SUMMARY STATISTICS #####################
                with open(output[0], 'a') as out:
                    out.write('===========================================\n' +
                              mirna + '\n')
                    out.write('**Motif: ' + motif + '\n')
                    out.write('**Consensus: ' + consensus + '\n')
                    out.write('**Fidelity(5p): ' + str(fidelity) + '\n')

        # SECTION | DISPLAY SEQUENCES AND SINGLE STATISTICS ###################
                    df.to_csv(out, sep='\t', index=False)
