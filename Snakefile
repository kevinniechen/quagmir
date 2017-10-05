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
from weighted_levenshtein import lev
from os.path import join
from os.path import splitext
from itertools import takewhile
from decimal import *
from Bio import SeqIO

###############################################################################

TIMESTAMP = time.strftime('%d-%b-%Y@%I:%M:%S%p')
SAMPLES = [os.path.basename(f) for f in glob.glob('data/*')]

###############################################################################

comparison = dict(A=['A'], C=['C'], G=['G'], T=['T'], R=['A', 'G'],
                  Y=['C', 'T'], S=['G', 'C'], W=['A', 'T'], K=['G', 'T'],
                  M=['A', 'C'], B=['C', 'G', 'T'], D=['A', 'G', 'T'],
                  H=['A', 'C', 'T'], V=['A', 'C', 'G'])


#class mystr(str):
#    def __contains__(self, value):
#        if 'R' in value:
#            return super().__contains__('A') or super().__contains__('G')
#        return super().__contains__(value)
#class MyChars:
#    __slots__ = ['N', 'A']
#R = 78
#A = 65
#G = 71

def compare_chars(fastq, consensus):

#    l = ord(fastq)
#    r = ord(consensus)
#    if not r ^ R:
#        if not (l ^ A and l ^ G ):
#            return True
#        return False

    if fastq == 'N' or consensus == 'N':
        return True
    return fastq in comparison[consensus]


#class mystr(str):
#    def __contains__(self, fastq):
#        for i in range(len(fastq)):
#            if not compare_chars(fastq[i], self[i]):
#                return False
#        return True


def compare_strings(fastq, consensus):
#    consensus = mystr(consensus)
#    return fastq in consensus
    if len(fastq) != len(consensus):
        return False
    for i in range(len(fastq)):
        if not compare_chars(fastq[i], consensus[i]):
            return False
    return True


def startswith_string(str1, str2):
    return compare_strings(str1[:len(str2)], str2)


def find_in_string(str1, str2):
    if startswith_string(str1, str2):
        return 0
    if len(str1) > len(str2):
        recursive_find = find_in_string(str1[1:], str2)
        if recursive_find == -1:
            return -1
        return recursive_find + 1
    return -1


def has_substitution_3p(len_trim, seq_end, consensus_end):
    if len_trim > 3:
        right = min(len(seq_end), len(consensus_end))
        left = max(0, len(seq_end) - len_trim + 1)
        if left < right and compare_strings(seq_end[left:right], consensus_end[left:right]):
            return False
    return True


def has_substitution_5p(seq_end, consensus_end):
    for i in range(0, min(len(seq_end), len(consensus_end))):
        if not compare_chars(seq_end[len(seq_end) - i - 1], consensus_end[len(consensus_end) - i - 1]):
            return True
    return False


def calc_trimming(seq_end, consensus_end):
    for i in range(0, min(len(seq_end), len(consensus_end))):
        if not compare_chars(seq_end[i], consensus_end[i]):
            return len(consensus_end) - i
    if len(seq_end) < len(consensus_end):
        return len(consensus_end) - len(seq_end)
    return 0


def calc_trimming_5p(seq_end, consensus_end):
    for i in range(0, min(len(seq_end), len(consensus_end))):
        if not compare_chars(seq_end[len(seq_end) - i - 1], consensus_end[len(consensus_end) - i - 1]):
            return len(consensus_end) - i
    if len(seq_end) < len(consensus_end):
        return len(consensus_end) - len(seq_end)
    return 0


def calc_tailing(seq_end, consensus_end, trim_len):
    return len(seq_end) - len(consensus_end) + trim_len


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
        ordered_dict[motif] = [mirna, consensus]
    return ordered_dict


def other_motifs_pulled_seq(dict_mirna_consensus, line, motif):
    annotation = ""
    list_motifs = list(dict_mirna_consensus.keys())
    matching_motifs = [x for x in list_motifs if contains_motif(line, x) if x != motif]
    for matched_motif in matching_motifs:
        annotation += dict_mirna_consensus[matched_motif][0]
        annotation += " "
    return annotation


def starts_with(line, motif):
    if len(motif) == 1:
        return str.startswith(line, motif[0])
    if not str.startswith(line, motif[0]) or len(motif[0])+motif[1] > len(line):
        return False
    return starts_with(line[len(motif[0])+motif[1]:], motif[2:])


def contains_motif(line, motif):
    if find_in_string(line, motif) == -1:
        return False
    return True


def chunkify_file(infilepath, delim=">"):
    with open(infilepath) as infile:
        answer = []
        tinfile = iter(infile)
        while 1:
            try:
                chunk = [next(tinfile)]
                chunk.extend(takewhile(lambda line: not line.startswith(delim), tinfile))
                answer.append(chunk)
            except StopIteration:
                break
    return answer


def input_name(prefix, sufix):
    answer = []
    for sample in SAMPLES:
        answer.append(prefix + sample + sufix)
    return answer


def sample_name(input, prefix, sufix):
    return input[len(prefix):-len(sufix)]


configfile:
    'config.yaml'


rule all:
    input:
        expand('results/{A}.isomir.tsv', A=SAMPLES),
        expand('results/{A}.isomir.sequence_info.tsv', A=SAMPLES),
        expand('group_results/' + config['group_output_name'] + '.isomir.tsv'),
        expand('group_results/' + config['group_output_name'] + '.isomir.sequence_info.tsv'),
        expand('group_results/' + config['group_output_name'] + '.isomir.nucleotide_dist.tsv')


rule collapse_fastq:
    input:
        'data/{A}'
    output:
        'collapsed/{A}.collapsed'
    shell:
        'awk "NR%4==2" {input} | sort -S1900G | uniq -c > {output}'


rule analyze_isomir:
    input:
        collapsed_fasta='collapsed/{A}.collapsed',
        input_files='data/{A}',
        motif_consensus=config['motif_consensus_file']
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
        # SECTION | SETUP MIRNA CONSENSUS DICT ################################
        dict_mirna_consensus = motif_consensus_to_dict(input.motif_consensus)

        # SECTION | MIRNA LOOP ################################################
        first_ds = True
        first_dsi = True
        first_dnd = True
        open(output[0], 'a').close()
        open(output[1], 'a').close()
        open(output[2], 'a').close()
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
                    if contains_motif(line, motif):
                        pulled_lines.append(line)

                for line in pulled_lines:
                    num_reads = int(line.rpartition(' ')[0])
                    seq = line.rpartition(' ')[2].rstrip()

                    # ascertain sequences pulled in by several miRNA motifs
                    has_other = other_motifs_pulled_seq(
                        dict_mirna_consensus, line, motif)

                    # sequence manipulations
                    consensus_index_3p = find_in_string(
                        consensus, motif) + len(motif)
                    seq_index_3p = find_in_string(seq, motif) + len(motif)
                    consensus_end_3p = consensus[consensus_index_3p:]
                    seq_end_3p = seq[seq_index_3p:]
                    consensus_index_5p = find_in_string(consensus, motif)
                    seq_index_5p = find_in_string(seq, motif)
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
                    if config['destructive_motif_pull'] and len(has_other) > 0 and find_in_file(output[0], seq):
                        logging.warning('Skipped ' + seq + ' ' + mirna)
                    elif config['filter_out_5p'] and has_substitution_5p(seq_end_5p, consensus_end_5p):
                        logging.warning(
                        'Skipped (5p substitution) ' + seq + ' ' + mirna)
                    elif config['filter_out_3p'] and has_substitution_3p(len_trim, seq_end_3p, consensus_end_3p):
                        logging.warning(
                        'Skipped (3p sequencing error) ' + seq + ' ' +
                            str(len_trim) + ' ' + mirna)
                    else:
                        # calculation of nt frequencies at each position
                        nt_offset = seq_index_5p - consensus_index_5p
                        for index, nt in enumerate(seq):
                            freq_nt[nt][index - nt_offset] += num_reads
                            for nt2 in ['A', 'C', 'G', 'T']:
                                if nt2!=nt and (
                                        freq_nt[nt] is None or freq_nt[nt2][index - nt_offset] is None):
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
                for base in comparison:
                    if base not in df2:
                        df2[base] = 0.0
                    else:
                        df2[base] = df2[base].fillna(0.0)
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
                    df2 = df2[["A", "C", "G", "T", "R", "Y", "S", "W", "K", "M",
                               "B", "D", "H", "V", "N", "READS"]]

    # SECTION | GENERATE SUMMARY STATISTICS ###################################
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

    # SECTION | CALCULATE CPM AND RPKM OF READS ###############################
                total_reads_in_sample = sum(1 for line in open(input.input_files)) / 4
                df = pd.read_csv(input[0], delimiter="\t", header=0)
                df['CPM'] = df['READS'] / df['READS'].sum() * float(
                    10 ** 6)  # TPM is also len-norm
                df['RPKM'] = df['READS'] / df[
                    'LEN_READ'] / total_reads_in_sample * float(10 ** 9)

    # SECTION | DISPLAY HEADER AND SUMMARY STATISTICS #########################
                if config['display_summary']:
                    with open(output[0], 'a') as out:
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

    # SECTION | DISPLAY SEQUENCES AND SINGLE STATISTICS ###############
                if config['display_sequence_info']:
                    with open(output[1], 'a') as out:
                        if first_dsi:
                            df.to_csv(out, sep='\t', index = False)
                            first_dsi = False
                        else:
                            df.to_csv(out, sep='\t', index=False, header = False)
                if config['display_nucleotide_dist']:
                    with open(output[2], 'a') as out:
                        if first_dnd:
                            df2.to_csv(out, sep='\t')
                            first_dnd = False
                        else:
                            df2.to_csv(out, sep='\t', header = False)


rule group_outputs:
    input:
        input_name("results/", ".isomir.tsv"),
        input_name("results/", ".isomir.sequence_info.tsv"),
        input_name("results/", ".isomir.nucleotide_dist.tsv")
    output:
        'group_results/' + config['group_output_name'] + '.isomir.tsv',
        'group_results/' + config['group_output_name'] + '.isomir.sequence_info.tsv',
        'group_results/' + config['group_output_name'] + '.isomir.nucleotide_dist.tsv'
    run:
        prefix = "results/"
        sufix = [".isomir.tsv",
                 ".isomir.sequence_info.tsv",
                 ".isomir.nucleotide_dict.tsv"]
        condition = [config['display_group_output'] and config['display_summary'],
                     config['display_group_output'] and config['display_sequence_info'],
                     config['display_group_output'] and config['display_nucleotide_dist']]

        if config['display_distance_metric']:
            motifs = motif_consensus_to_dict(config['motif_consensus_file'])
            mirna_consensus = {v[0]:v[1] for k, v in motifs.items()}

        for i in range(3):
            open(output[i], 'a').close()
            with open(input[i]) as input_file:
                first_char = input_file.read(1)
            if condition[i] and first_char:
                with open(output[i], 'a') as out:
                    df = pd.read_csv(input[i*len(SAMPLES)], delimiter="\t", header = 0)
                    df['SAMPLE'] = sample_name(input[i*len(SAMPLES)], prefix, sufix[i])
                    cols = df.columns.tolist()
                    cols = cols[-1:] + cols[:-1]
                    df = df[cols]
                    for inp in input[i*len(SAMPLES) + 1 : (i+1)*len(SAMPLES)]:
                        df_rest = pd.read_csv(inp, delimiter="\t", header = 0)
                        df_rest['SAMPLE'] = sample_name(inp, prefix, sufix[i])
                        cols = df_rest.columns.tolist()
                        cols = cols[-1:] + cols[:-1]
                        df_rest = df_rest[cols]
                        df = df.append(df_rest)

                    # add distance metrics to seq_info or expression matrix
                    if config['display_distance_metric'] and (i==1 or i==2):
                        consensus = df["MIRNA"].apply(lambda x: mirna_consensus[x])
                        unique_seqs = set(zip(df["SEQUENCE"], consensus))
                        edit_distances = {}

                        # set scores from config
                        delete_costs = np.ones(128, dtype=np.float64)
                        insert_costs = np.ones(128, dtype=np.float64)
                        substitute_costs = np.ones((128, 128), dtype=np.float64)
                        if config['deletion_score']:
                            delete_costs[ord('A')] = config['deletion_score']
                            delete_costs[ord('C')] = config['deletion_score']
                            delete_costs[ord('G')] = config['deletion_score']
                            delete_costs[ord('T')] = config['deletion_score']
                        if config['insertion_score']:
                            insert_costs[ord('A')] = config['insertion_score']
                            insert_costs[ord('C')] = config['insertion_score']
                            insert_costs[ord('G')] = config['insertion_score']
                            insert_costs[ord('T')] = config['insertion_score']
                        if config['substitution_AG']:
                            substitute_costs[ord('A'), ord('G')] = config['substitution_AG']
                        if config['substitution_GA']:
                            substitute_costs[ord('G'), ord('A')] = config['substitution_GA']
                        if config['substitution_CT']:
                            substitute_costs[ord('C'), ord('T')] = config['substitution_CT']
                        if config['substitution_TC']:
                            substitute_costs[ord('T'), ord('C')] = config['substitution_TC']
                        if config['substitution_AT']:
                            substitute_costs[ord('A'), ord('T')] = config['substitution_AT']
                        if config['substitution_TA']:
                            substitute_costs[ord('T'), ord('A')] = config['substitution_TA']
                        if config['substitution_AC']:
                            substitute_costs[ord('A'), ord('C')] = config['substitution_AC']
                        if config['substitution_CA']:
                            substitute_costs[ord('C'), ord('A')] = config['substitution_CA']
                        if config['substitution_GC']:
                            substitute_costs[ord('G'), ord('C')] = config['substitution_GC']
                        if config['substitution_CG']:
                            substitute_costs[ord('C'), ord('G')] = config['substitution_CG']
                        if config['substitution_GT']:
                            substitute_costs[ord('G'), ord('T')] = config['substitution_GT']
                        if config['substitution_TG']:
                            substitute_costs[ord('T'), ord('G')] = config['substitution_TG']
                        for base in comparison.keys():
                            substitute_costs[ord('N'), ord(base)] = 0
                        for base in ['A', 'C', 'G', 'T']:
                            substitute_costs[ord(base), ord('N')] = 0
                        for base in ['A', 'C', 'G', 'T']:
                            ord_base = ord(base)
                            for base1 in comparison.keys():
                                if base1 not in ['A', 'C', 'G', 'T', 'N']:
                                    substitute_costs[ord_base, ord(base1)] = min(
                                        [substitute_costs[ord_base, ord(base2)] for base2 in
                                         comparison[base1]])
                        for pair in unique_seqs:
                            edit_distances[pair] = lev(pair[0], pair[1],
                                delete_costs = delete_costs,
                                substitute_costs = substitute_costs,
                                insert_costs = insert_costs)

                        df['DISTANCE'] = df[['SEQUENCE', 'MIRNA']].apply(lambda row:
                                         edit_distances[(row[0], mirna_consensus[row[1]])],
                                         axis=1)
                    df.to_csv(out, sep='\t', index=False, header=True)