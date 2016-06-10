#!/usr/bin/env python

""" A pipeline to analyze
raw miRNA-seq files for
isomiR composition.
Currently takes in .fastq_ready
files in /data folder and motif.txt"""

import glob
import csv
from os.path import join
from os.path import splitext
from decimal import *
from Bio import SeqIO

###############################################################################

MOTIF_CONSENSUS = 'motif-consensus.fa'

SAMPLES = [os.path.basename(f) for f in glob.glob('data/*.fastq_ready')]

###############################################################################


def calc_trimming_5p(seq_end, consensus_end, seq):
    trimming_threshold = 5
    for i in range(0, min(len(seq_end), len(consensus_end))):
        if seq_end[len(seq_end) - i - 1] != consensus_end[len(
                consensus_end) - i - 1]:
            # snp check
            if (len(consensus_end) - i) > trimming_threshold:
                return '0-SNP'
            else:
                return (len(consensus_end) - i)
    if len(seq_end) < len(consensus_end):
        return len(consensus_end) - len(seq_end)
    return 0


def calc_trimming(seq_end, consensus_end):
    trimming_threshold = 5
    for i in range(0, min(len(seq_end), len(consensus_end))):
        if seq_end[i] != consensus_end[i]:
            # snp check
            if (len(consensus_end) - i) > trimming_threshold:
                return '0-SNP'
            else:
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
        "report.txt"

rule collapse_fastq:
    input:
        expand('data/{A}', A=SAMPLES)
    output:
        expand('data/{A}.collapsed', A=SAMPLES)
    shell:
        "awk 'NR%4==2' {input} | sort -S1900G | uniq -c > {output}"

rule analyze_isomir:
    input:
        motif_consensus = MOTIF_CONSENSUS,
        collapsed_fasta = expand('data/{A}.collapsed', A=SAMPLES)
    output:
        "report.txt"
    run:
        input_sequences = SeqIO.parse(open('motif-consensus.fa'), 'fasta')
        for seq_record in input_sequences:
            mirna = seq_record.description.split()[0]
            motif = seq_record.description.split()[1]
            consensus = str(seq_record.seq)

            # going in mirna-seq collapsed file, searching for single motif
            with open(str(input.collapsed_fasta), "rt") as f:
                # calculate total reads of miR
                total_reads = 0
                variation = 0
                for line in f:
                    reads = line.rpartition(' ')[0]
                    seq = line.rpartition(' ')[2].rstrip()
                    if motif in line:
                        total_reads += int(reads)
                f.seek(0)

                with open(output[0], 'a') as out:
                    out.write('===========================================\n' +
                              mirna + '\n')
                    out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                              "SEQUENCE",
                              "MIRNA_READS",
                              "PERCENTAGE",
                              "TRIM_LENGTH",
                              "TAIL_LENGTH",
                              "TAIL_SEQUENCE",
                              "VARIATION_5P"))

                for line in f:
                    reads = line.rpartition(' ')[0]
                    seq = line.rpartition(' ')[2].rstrip()
                    if motif in line:
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

                        # calculation of output values
                        percentage = '{0:.2f}%'.format(
                            100 * Decimal(reads) / Decimal(total_reads))
                        trim_len = calc_trimming(seq_end_3p, consensus_end_3p)
                        tail_len = calc_tailing(
                            seq_end_3p, consensus_end_3p, trim_len)
                        tail_seq = get_tailing_seq(seq, tail_len)
                        trim_len_5p = calc_trimming_5p(
                            seq_end_5p, consensus_end_5p, seq)
                        variation_5p = max(
                            trim_len_5p,
                            calc_tailing(
                                seq_end_5p, consensus_end_5p, trim_len_5p))

                        # display
                        with open(output[0], 'a') as out:
                            out.write(
                                '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                                    seq,
                                    reads,
                                    percentage,
                                    trim_len,
                                    tail_len,
                                    tail_seq,
                                    variation_5p))

                        # calculate variation score (5' end variation)
                        variation += ((Decimal(reads) *
                                       variation_5p) / total_reads)

                with open(output[0], 'a') as out:
                    out.write('***5P-Variation: ' +
                              '{0:.4f}'.format(variation) + '\n' +
                              '***Motif: ' + motif + '\n' +
                              '***Consensus: ' + consensus + '\n')
