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
from decimal import *
from Bio import SeqIO

###############################################################################

TIMESTAMP = time.strftime('%d-%b-%Y@%I:%M:%S%p')
SAMPLES = [os.path.basename(f) for f in glob.glob('data/*.fastq_ready')]

###############################################################################


def has_substitution_3p(len_trim, seq_end, consensus_end):
    min_len = min(len(seq_end), len(consensus_end))
    if len_trim > 3:
        if (seq_end[min_len - len_trim + 1:min_len] ==
            consensus_end[min_len - len_trim + 1:min_len]):
            return True
    return False


def has_substitution_5p(seq_end, consensus_end):
    for i in range(0, min(len(seq_end), len(consensus_end))):
        if (seq_end[len(seq_end) - i - 1] !=
            consensus_end[len(consensus_end) - i - 1]):
            return True
    return False

def calc_trimming(seq_end, consensus_end):
    for i in range(0, min(len(seq_end), len(consensus_end))):
        if seq_end[i] != consensus_end[i]:
            return (len(consensus_end) - i)
    if len(seq_end) < len(consensus_end):
        return len(consensus_end) - len(seq_end)
    return 0


def calc_trimming_5p(seq_end, consensus_end):
    for i in range(0, min(len(seq_end), len(consensus_end))):
        if seq_end[len(seq_end) - i - 1] != consensus_end[len(
                consensus_end) - i - 1]:
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
            raise Exception("\n************************************\n" +
                "DUPLICATE MOTIFS FOUND IN '" + file + "'\n" +
                "PLEASE CHECK YOUR MOTIF CONSENSUS FILE\n" +
                "AND FIX OR DELETE BEFORE RERUNNING\n" +
                "************************************\n")
        ordered_dict[motif] = [mirna, consensus]
    return ordered_dict


def other_motifs_pulled_seq(dict_mirna_consensus, line, motif):
    annotation = ""
    list_motifs = list(dict_mirna_consensus.keys())
    matching_motifs = [
        x for x in list_motifs if x in line if x != motif]
    for matched_motif in matching_motifs:
        annotation += dict_mirna_consensus[matched_motif][0]
        annotation += " "
    return annotation


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
    'config.yaml'

rule all:
    input:
        expand('results/tabular/{A}.isomir.tsv', A=SAMPLES),
        expand('results/text/{A}.mirna.txt', A=SAMPLES)

rule collapse_fastq:
    input:
        'data/{A}'
    output:
        'data/collapsed/{A}.collapsed'
    shell:
        'awk "NR%4==2" {input} | sort -S1900G | uniq -c > {output}'

rule analyze_isomir:
    input:
        motif_consensus = config['motif_consensus_file'],
        collapsed_fasta = 'data/collapsed/{A}.collapsed'
    output:
        'results/tabular/{A}.isomir.tsv'
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
        with open(output[0], 'a') as out:
            out.write(TIMESTAMP + '\n')

        # SECTION | SETUP MIRNA CONSENSUS DICT #####################################
        dict_mirna_consensus = motif_consensus_to_dict(input.motif_consensus)

        # SECTION | MIRNA LOOP ################################################
        for motif, value in dict_mirna_consensus.items():
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
                    if motif in line:
                        pulled_lines.append(line)

                for line in pulled_lines:
                    num_reads = int(line.rpartition(' ')[0])
                    seq = line.rpartition(' ')[2].rstrip()

                    # ascertain sequences pulled in by several miRNA motifs
                    has_other = other_motifs_pulled_seq(
                        dict_mirna_consensus, line, motif)

                    # sequence manipulations
                    consensus_index_3p = str.find(
                        consensus, motif) + len(motif)
                    seq_index_3p = str.find(seq, motif) + len(motif)
                    consensus_end_3p = consensus[consensus_index_3p:]
                    seq_end_3p = seq[seq_index_3p:]
                    consensus_index_5p = str.find(consensus, motif)
                    seq_index_5p = str.find(seq, motif)
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
                        # if sequence has minimum reads, add to display queue
                        if (ratio > config['min_ratio'] or
                            num_reads > config['min_read']):
                            table_out.append([seq, len_read, num_reads,
                                ratio, len_trim, len_tail, seq_tail,
                                vari_5p, has_other])

    # SECTION | MOVE STATISTICS INTO DATAFRAME ############################
                if len(table_out):
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
                # calculate total reads
                total_reads = float(df['READS'].sum())
                if total_reads == 0:
                    raise Exception("\n************************************\n" +
                        "NO MATCHED READS FOUND IN '" + input.collapsed_fasta + "'\n" +
                        "PLEASE CHECK YOUR FASTQ_READY FILE\n" +
                        "AND FIX OR DELETE BEFORE RERUNNING\n" +
                        "************************************\n")

                # calculate ratio
                df['RATIO'] = df['READS'].apply(lambda x: str(round(100*x/total_reads, 2))+'%')

                # calculate number of isomirs (visible)
                total_isomirs = df.shape[0]
                if df['LEN_TRIM'].iloc[0] == 0 and df['LEN_TAIL'].iloc[0] == 0:
                    total_isomirs -= 1

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

    # SECTION | DISPLAY HEADER AND SUMMARY STATISTICS #################
                with open(output[0], 'a') as out:
                    if config['display_summary']:
                        out.write('==========================\n' +
                                  '>' + mirna + '\n')
                        out.write('motif:\t' + motif + '\n')
                        out.write('consensus:\t' + consensus + '\n')
                        out.write('total-reads:\t' +
                                  str(int(total_reads)) + '\n')
                        out.write('total-isomirs:\t' +
                                   str(total_isomirs) + '\n')
                        out.write('fidelity-5p:\t' +
                                  str(fidelity) + '\n')
                        out.write('A-tailing:\t' +
                                  ratio_a_tailing + '%\n')
                        out.write('C-tailing:\t' +
                                  ratio_c_tailing + '%\n')
                        out.write('G-tailing:\t' +
                                  ratio_g_tailing + '%\n')
                        out.write('T-tailing:\t' +
                                  ratio_t_tailing + '%\n')
                        out.write('sequence-trimming-only:\t' +
                                  str(ratio_seq_trim_only) + '%\n')
                        out.write('sequence-trimming:\t' +
                                  str(ratio_seq_trim_only + ratio_seq_trim_and_tail) + '%\n')
                        out.write('sequence-tailing-only:\t' +
                                  str(ratio_seq_tail_only) + '%\n')
                        out.write('sequence-tailing:\t' +
                                  str(ratio_seq_tail_only + ratio_seq_trim_and_tail) + '%\n')
                        out.write('sequence-trimming-and-tailing:\t' +
                                  str(ratio_seq_trim_and_tail) + '%\n')
                        out.write('\n')

    # SECTION | DISPLAY SEQUENCES AND SINGLE STATISTICS ###############
                    if config['display_sequence_info']:
                        out.write('[sequence-information]\n')
                        df.to_csv(out, sep='\t', index=False)
                        out.write('\n')
                    if config['display_nucleotide_dist']:
                        out.write('[nucleotide-distribution]\n')
                        df2.to_csv(out, sep='\t')
                        out.write('\n')

rule summarize_fastq:
    input:
        'data/{A}',
        'results/tabular/{A}.isomir.tsv'
    output:
        'results/text/{A}.mirna.txt'
    shell:
        """
        total_reads=$(($(cat {input[0]} | wc -l) / 4))
        echo "Total miR reads in sample: $total_reads" >> {output}
        grep ">hsa\|total-reads" {input[1]} | cat >> {output}
        echo "\n" >> {output}
        """
