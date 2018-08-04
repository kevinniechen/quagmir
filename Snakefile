#!/usr/bin/env python

"""
Author: K. Chen
Affiliation: NCI/NIH
Aim: Python-based isomiR quantification and analysis pipeline.
Run: snakemake
"""


import glob
import collections as co
import numpy as np
import pandas as pd
import time
import logging
from weighted_levenshtein import lev
from Bio import SeqIO

###############################################################################

TIMESTAMP = time.strftime('%d-%b-%Y@%I:%M:%S%p')
SAMPLES = [os.path.basename(f) for f in glob.glob('data/*')]

###############################################################################

base_ords = (A, C, G, T, R, Y, S, W, K, M, B, D, H, V,
         N) = (65, 67, 71, 84, 82, 89, 83, 87, 75, 77, 66, 68, 72, 86, 78)
base_ords = set(base_ords)
fastq_ords = {A, C, G, T, N}
acgt = ["A", "C", "G", "T"]
acgtn = ["A", "C", "G", "T", "N"]
acgt_ords = {A, C, G, T}
fastq_bases = {"A", "C", "G", "T", "N"}
all_bases = {"A", "C", "G", "T", "N", "R", "Y", "S", "W", "K", "M", "B", "D",
             "H", "V"}
ambiguous_letters = {"N", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"}
comparison_ord = dict([(A, {A}), (C, {C}), (G, {G}), (T, {T}), (R, {A, G}),
                       (Y, {C, T}), (S, {C, G}), (W, {A, T}), (K, {G, T}),
                       (M, {A, C}), (B, {C, G, T}), (D, {A, G, T}),
                       (H, {A, C, T}), (V, {A, C, G})])
comparison = dict([('A', {'A', 'N'}), ('C', {'C', 'N'}), ('G', {'G', 'N'}),
                   ('T', {'T', 'N'}), ('R', {'A', 'G', 'N'}),
                   ('Y', {'C', 'T', 'N'}), ('S', {'C', 'G', 'N'}),
                   ('W', {'A', 'T', 'N'}), ('K', {'G', 'T', 'N'}),
                   ('M', {'A', 'C', 'N'}), ('B', {'C', 'G', 'T', 'N'}),
                   ('D', {'A', 'G', 'T', 'N'}), ('H', {'A', 'C', 'T', 'N'}),
                   ('V', {'A', 'C', 'G', 'N'}),
                   ('N', {'A', 'C', 'G', 'T', 'N', 'R', 'Y', 'S', 'W', 'K',
                          'M', 'B', 'D', 'H', 'V', 'N'})])
codes = {"AAAAA":"BB", "AAAAC":"BD", "AAAAG":"B0", "AAAAT":"BE", "AAACA":"BF",
         "AAACC":"B1", "AAACG":"BH", "AAACT":"BI", "AAAGA":"B2", "AAAGC":"BJ",
         "AAAGG":"BK", "AAAGT":"B3", "AAATA":"BL", "AAATC":"BM", "AAATG":"B4",
         "AAATT":"BN", "AACAA":"BO", "AACAC":"B5", "AACAG":"BP", "AACAT":"BQ",
         "AACCA":"B6", "AACCC":"BR", "AACCG":"BS", "AACCT":"B7", "AACGA":"BU",
         "AACGC":"BV", "AACGG":"B8", "AACGT":"BW", "AACTA":"BX", "AACTC":"B9",
         "AACTG":"BY", "AACTT":"BZ", "AAGAA":"DB", "AAGAC":"DD", "AAGAG":"D0",
         "AAGAT":"DE", "AAGCA":"DF", "AAGCC":"D1", "AAGCG":"DH", "AAGCT":"DI",
         "AAGGA":"D2", "AAGGC":"DJ", "AAGGG":"DK", "AAGGT":"D3", "AAGTA":"DL",
         "AAGTC":"DM", "AAGTG":"D4", "AAGTT":"DN", "AATAA":"DO", "AATAC":"D5",
         "AATAG":"DP", "AATAT":"DQ", "AATCA":"D6", "AATCC":"DR", "AATCG":"DS",
         "AATCT":"D7", "AATGA":"DU", "AATGC":"DV", "AATGG":"D8", "AATGT":"DW",
         "AATTA":"DX", "AATTC":"D9", "AATTG":"DY", "AATTT":"DZ", "ACAAA":"0B",
         "ACAAC":"0D", "ACAAG":"00", "ACAAT":"0E", "ACACA":"0F", "ACACC":"01",
         "ACACG":"0H", "ACACT":"0I", "ACAGA":"02", "ACAGC":"0J", "ACAGG":"0K",
         "ACAGT":"03", "ACATA":"0L", "ACATC":"0M", "ACATG":"04", "ACATT":"0N",
         "ACCAA":"0O", "ACCAC":"05", "ACCAG":"0P", "ACCAT":"0Q", "ACCCA":"06",
         "ACCCC":"0R", "ACCCG":"0S", "ACCCT":"07", "ACCGA":"0U", "ACCGC":"0V",
         "ACCGG":"08", "ACCGT":"0W", "ACCTA":"0X", "ACCTC":"09", "ACCTG":"0Y",
         "ACCTT":"0Z", "ACGAA":"EB", "ACGAC":"ED", "ACGAG":"E0", "ACGAT":"EE",
         "ACGCA":"EF", "ACGCC":"E1", "ACGCG":"EH", "ACGCT":"EI", "ACGGA":"E2",
         "ACGGC":"EJ", "ACGGG":"EK", "ACGGT":"E3", "ACGTA":"EL", "ACGTC":"EM",
         "ACGTG":"E4", "ACGTT":"EN", "ACTAA":"EO", "ACTAC":"E5", "ACTAG":"EP",
         "ACTAT":"EQ", "ACTCA":"E6", "ACTCC":"ER", "ACTCG":"ES", "ACTCT":"E7",
         "ACTGA":"EU", "ACTGC":"EV", "ACTGG":"E8", "ACTGT":"EW", "ACTTA":"EX",
         "ACTTC":"E9", "ACTTG":"EY", "ACTTT":"EZ", "AGAAA":"FB", "AGAAC":"FD",
         "AGAAG":"F0", "AGAAT":"FE", "AGACA":"FF", "AGACC":"F1", "AGACG":"FH",
         "AGACT":"FI", "AGAGA":"F2", "AGAGC":"FJ", "AGAGG":"FK", "AGAGT":"F3",
         "AGATA":"FL", "AGATC":"FM", "AGATG":"F4", "AGATT":"FN", "AGCAA":"FO",
         "AGCAC":"F5", "AGCAG":"FP", "AGCAT":"FQ", "AGCCA":"F6", "AGCCC":"FR",
         "AGCCG":"FS", "AGCCT":"F7", "AGCGA":"FU", "AGCGC":"FV", "AGCGG":"F8",
         "AGCGT":"FW", "AGCTA":"FX", "AGCTC":"F9", "AGCTG":"FY", "AGCTT":"FZ",
         "AGGAA":"1B", "AGGAC":"1D", "AGGAG":"10", "AGGAT":"1E", "AGGCA":"1F",
         "AGGCC":"11", "AGGCG":"1H", "AGGCT":"1I", "AGGGA":"12", "AGGGC":"1J",
         "AGGGG":"1K", "AGGGT":"13", "AGGTA":"1L", "AGGTC":"1M", "AGGTG":"14",
         "AGGTT":"1N", "AGTAA":"1O", "AGTAC":"15", "AGTAG":"1P", "AGTAT":"1Q",
         "AGTCA":"16", "AGTCC":"1R", "AGTCG":"1S", "AGTCT":"17", "AGTGA":"1U",
         "AGTGC":"1V", "AGTGG":"18", "AGTGT":"1W", "AGTTA":"1X", "AGTTC":"19",
         "AGTTG":"1Y", "AGTTT":"1Z", "ATAAA":"HB", "ATAAC":"HD", "ATAAG":"H0",
         "ATAAT":"HE", "ATACA":"HF", "ATACC":"H1", "ATACG":"HH", "ATACT":"HI",
         "ATAGA":"H2", "ATAGC":"HJ", "ATAGG":"HK", "ATAGT":"H3", "ATATA":"HL",
         "ATATC":"HM", "ATATG":"H4", "ATATT":"HN", "ATCAA":"HO", "ATCAC":"H5",
         "ATCAG":"HP", "ATCAT":"HQ", "ATCCA":"H6", "ATCCC":"HR", "ATCCG":"HS",
         "ATCCT":"H7", "ATCGA":"HU", "ATCGC":"HV", "ATCGG":"H8", "ATCGT":"HW",
         "ATCTA":"HX", "ATCTC":"H9", "ATCTG":"HY", "ATCTT":"HZ", "ATGAA":"IB",
         "ATGAC":"ID", "ATGAG":"I0", "ATGAT":"IE", "ATGCA":"IF", "ATGCC":"I1",
         "ATGCG":"IH", "ATGCT":"II", "ATGGA":"I2", "ATGGC":"IJ", "ATGGG":"IK",
         "ATGGT":"I3", "ATGTA":"IL", "ATGTC":"IM", "ATGTG":"I4", "ATGTT":"IN",
         "ATTAA":"IO", "ATTAC":"I5", "ATTAG":"IP", "ATTAT":"IQ", "ATTCA":"I6",
         "ATTCC":"IR", "ATTCG":"IS", "ATTCT":"I7", "ATTGA":"IU", "ATTGC":"IV",
         "ATTGG":"I8", "ATTGT":"IW", "ATTTA":"IX", "ATTTC":"I9", "ATTTG":"IY",
         "ATTTT":"IZ", "CAAAA":"2B", "CAAAC":"2D", "CAAAG":"20", "CAAAT":"2E",
         "CAACA":"2F", "CAACC":"21", "CAACG":"2H", "CAACT":"2I", "CAAGA":"22",
         "CAAGC":"2J", "CAAGG":"2K", "CAAGT":"23", "CAATA":"2L", "CAATC":"2M",
         "CAATG":"24", "CAATT":"2N", "CACAA":"2O", "CACAC":"25", "CACAG":"2P",
         "CACAT":"2Q", "CACCA":"26", "CACCC":"2R", "CACCG":"2S", "CACCT":"27",
         "CACGA":"2U", "CACGC":"2V", "CACGG":"28", "CACGT":"2W", "CACTA":"2X",
         "CACTC":"29", "CACTG":"2Y", "CACTT":"2Z", "CAGAA":"JB", "CAGAC":"JD",
         "CAGAG":"J0", "CAGAT":"JE", "CAGCA":"JF", "CAGCC":"J1", "CAGCG":"JH",
         "CAGCT":"JI", "CAGGA":"J2", "CAGGC":"JJ", "CAGGG":"JK", "CAGGT":"J3",
         "CAGTA":"JL", "CAGTC":"JM", "CAGTG":"J4", "CAGTT":"JN", "CATAA":"JO",
         "CATAC":"J5", "CATAG":"JP", "CATAT":"JQ", "CATCA":"J6", "CATCC":"JR",
         "CATCG":"JS", "CATCT":"J7", "CATGA":"JU", "CATGC":"JV", "CATGG":"J8",
         "CATGT":"JW", "CATTA":"JX", "CATTC":"J9", "CATTG":"JY", "CATTT":"JZ",
         "CCAAA":"KB", "CCAAC":"KD", "CCAAG":"K0", "CCAAT":"KE", "CCACA":"KF",
         "CCACC":"K1", "CCACG":"KH", "CCACT":"KI", "CCAGA":"K2", "CCAGC":"KJ",
         "CCAGG":"KK", "CCAGT":"K3", "CCATA":"KL", "CCATC":"KM", "CCATG":"K4",
         "CCATT":"KN", "CCCAA":"KO", "CCCAC":"K5", "CCCAG":"KP", "CCCAT":"KQ",
         "CCCCA":"K6", "CCCCC":"KR", "CCCCG":"KS", "CCCCT":"K7", "CCCGA":"KU",
         "CCCGC":"KV", "CCCGG":"K8", "CCCGT":"KW", "CCCTA":"KX", "CCCTC":"K9",
         "CCCTG":"KY", "CCCTT":"KZ", "CCGAA":"3B", "CCGAC":"3D", "CCGAG":"30",
         "CCGAT":"3E", "CCGCA":"3F", "CCGCC":"31", "CCGCG":"3H", "CCGCT":"3I",
         "CCGGA":"32", "CCGGC":"3J", "CCGGG":"3K", "CCGGT":"33", "CCGTA":"3L",
         "CCGTC":"3M", "CCGTG":"34", "CCGTT":"3N", "CCTAA":"3O", "CCTAC":"35",
         "CCTAG":"3P", "CCTAT":"3Q", "CCTCA":"36", "CCTCC":"3R", "CCTCG":"3S",
         "CCTCT":"37", "CCTGA":"3U", "CCTGC":"3V", "CCTGG":"38", "CCTGT":"3W",
         "CCTTA":"3X", "CCTTC":"39", "CCTTG":"3Y", "CCTTT":"3Z", "CGAAA":"LB",
         "CGAAC":"LD", "CGAAG":"L0", "CGAAT":"LE", "CGACA":"LF", "CGACC":"L1",
         "CGACG":"LH", "CGACT":"LI", "CGAGA":"L2", "CGAGC":"LJ", "CGAGG":"LK",
         "CGAGT":"L3", "CGATA":"LL", "CGATC":"LM", "CGATG":"L4", "CGATT":"LN",
         "CGCAA":"LO", "CGCAC":"L5", "CGCAG":"LP", "CGCAT":"LQ", "CGCCA":"L6",
         "CGCCC":"LR", "CGCCG":"LS", "CGCCT":"L7", "CGCGA":"LU", "CGCGC":"LV",
         "CGCGG":"L8", "CGCGT":"LW", "CGCTA":"LX", "CGCTC":"L9", "CGCTG":"LY",
         "CGCTT":"LZ", "CGGAA":"MB", "CGGAC":"MD", "CGGAG":"M0", "CGGAT":"ME",
         "CGGCA":"MF", "CGGCC":"M1", "CGGCG":"MH", "CGGCT":"MI", "CGGGA":"M2",
         "CGGGC":"MJ", "CGGGG":"MK", "CGGGT":"M3", "CGGTA":"ML", "CGGTC":"MM",
         "CGGTG":"M4", "CGGTT":"MN", "CGTAA":"MO", "CGTAC":"M5", "CGTAG":"MP",
         "CGTAT":"MQ", "CGTCA":"M6", "CGTCC":"MR", "CGTCG":"MS", "CGTCT":"M7",
         "CGTGA":"MU", "CGTGC":"MV", "CGTGG":"M8", "CGTGT":"MW", "CGTTA":"MX",
         "CGTTC":"M9", "CGTTG":"MY", "CGTTT":"MZ", "CTAAA":"4B", "CTAAC":"4D",
         "CTAAG":"40", "CTAAT":"4E", "CTACA":"4F", "CTACC":"41", "CTACG":"4H",
         "CTACT":"4I", "CTAGA":"42", "CTAGC":"4J", "CTAGG":"4K", "CTAGT":"43",
         "CTATA":"4L", "CTATC":"4M", "CTATG":"44", "CTATT":"4N", "CTCAA":"4O",
         "CTCAC":"45", "CTCAG":"4P", "CTCAT":"4Q", "CTCCA":"46", "CTCCC":"4R",
         "CTCCG":"4S", "CTCCT":"47", "CTCGA":"4U", "CTCGC":"4V", "CTCGG":"48",
         "CTCGT":"4W", "CTCTA":"4X", "CTCTC":"49", "CTCTG":"4Y", "CTCTT":"4Z",
         "CTGAA":"NB", "CTGAC":"ND", "CTGAG":"N0", "CTGAT":"NE", "CTGCA":"NF",
         "CTGCC":"N1", "CTGCG":"NH", "CTGCT":"NI", "CTGGA":"N2", "CTGGC":"NJ",
         "CTGGG":"NK", "CTGGT":"N3", "CTGTA":"NL", "CTGTC":"NM", "CTGTG":"N4",
         "CTGTT":"NN", "CTTAA":"NO", "CTTAC":"N5", "CTTAG":"NP", "CTTAT":"NQ",
         "CTTCA":"N6", "CTTCC":"NR", "CTTCG":"NS", "CTTCT":"N7", "CTTGA":"NU",
         "CTTGC":"NV", "CTTGG":"N8", "CTTGT":"NW", "CTTTA":"NX", "CTTTC":"N9",
         "CTTTG":"NY", "CTTTT":"NZ", "GAAAA":"OB", "GAAAC":"OD", "GAAAG":"O0",
         "GAAAT":"OE", "GAACA":"OF", "GAACC":"O1", "GAACG":"OH", "GAACT":"OI",
         "GAAGA":"O2", "GAAGC":"OJ", "GAAGG":"OK", "GAAGT":"O3", "GAATA":"OL",
         "GAATC":"OM", "GAATG":"O4", "GAATT":"ON", "GACAA":"OO", "GACAC":"O5",
         "GACAG":"OP", "GACAT":"OQ", "GACCA":"O6", "GACCC":"OR", "GACCG":"OS",
         "GACCT":"O7", "GACGA":"OU", "GACGC":"OV", "GACGG":"O8", "GACGT":"OW",
         "GACTA":"OX", "GACTC":"O9", "GACTG":"OY", "GACTT":"OZ", "GAGAA":"5B",
         "GAGAC":"5D", "GAGAG":"50", "GAGAT":"5E", "GAGCA":"5F", "GAGCC":"51",
         "GAGCG":"5H", "GAGCT":"5I", "GAGGA":"52", "GAGGC":"5J", "GAGGG":"5K",
         "GAGGT":"53", "GAGTA":"5L", "GAGTC":"5M", "GAGTG":"54", "GAGTT":"5N",
         "GATAA":"5O", "GATAC":"55", "GATAG":"5P", "GATAT":"5Q", "GATCA":"56",
         "GATCC":"5R", "GATCG":"5S", "GATCT":"57", "GATGA":"5U", "GATGC":"5V",
         "GATGG":"58", "GATGT":"5W", "GATTA":"5X", "GATTC":"59", "GATTG":"5Y",
         "GATTT":"5Z", "GCAAA":"PB", "GCAAC":"PD", "GCAAG":"P0", "GCAAT":"PE",
         "GCACA":"PF", "GCACC":"P1", "GCACG":"PH", "GCACT":"PI", "GCAGA":"P2",
         "GCAGC":"PJ", "GCAGG":"PK", "GCAGT":"P3", "GCATA":"PL", "GCATC":"PM",
         "GCATG":"P4", "GCATT":"PN", "GCCAA":"PO", "GCCAC":"P5", "GCCAG":"PP",
         "GCCAT":"PQ", "GCCCA":"P6", "GCCCC":"PR", "GCCCG":"PS", "GCCCT":"P7",
         "GCCGA":"PU", "GCCGC":"PV", "GCCGG":"P8", "GCCGT":"PW", "GCCTA":"PX",
         "GCCTC":"P9", "GCCTG":"PY", "GCCTT":"PZ", "GCGAA":"QB", "GCGAC":"QD",
         "GCGAG":"Q0", "GCGAT":"QE", "GCGCA":"QF", "GCGCC":"Q1", "GCGCG":"QH",
         "GCGCT":"QI", "GCGGA":"Q2", "GCGGC":"QJ", "GCGGG":"QK", "GCGGT":"Q3",
         "GCGTA":"QL", "GCGTC":"QM", "GCGTG":"Q4", "GCGTT":"QN", "GCTAA":"QO",
         "GCTAC":"Q5", "GCTAG":"QP", "GCTAT":"QQ", "GCTCA":"Q6", "GCTCC":"QR",
         "GCTCG":"QS", "GCTCT":"Q7", "GCTGA":"QU", "GCTGC":"QV", "GCTGG":"Q8",
         "GCTGT":"QW", "GCTTA":"QX", "GCTTC":"Q9", "GCTTG":"QY", "GCTTT":"QZ",
         "GGAAA":"6B", "GGAAC":"6D", "GGAAG":"60", "GGAAT":"6E", "GGACA":"6F",
         "GGACC":"61", "GGACG":"6H", "GGACT":"6I", "GGAGA":"62", "GGAGC":"6J",
         "GGAGG":"6K", "GGAGT":"63", "GGATA":"6L", "GGATC":"6M", "GGATG":"64",
         "GGATT":"6N", "GGCAA":"6O", "GGCAC":"65", "GGCAG":"6P", "GGCAT":"6Q",
         "GGCCA":"66", "GGCCC":"6R", "GGCCG":"6S", "GGCCT":"67", "GGCGA":"6U",
         "GGCGC":"6V", "GGCGG":"68", "GGCGT":"6W", "GGCTA":"6X", "GGCTC":"69",
         "GGCTG":"6Y", "GGCTT":"6Z", "GGGAA":"RB", "GGGAC":"RD", "GGGAG":"R0",
         "GGGAT":"RE", "GGGCA":"RF", "GGGCC":"R1", "GGGCG":"RH", "GGGCT":"RI",
         "GGGGA":"R2", "GGGGC":"RJ", "GGGGG":"RK", "GGGGT":"R3", "GGGTA":"RL",
         "GGGTC":"RM", "GGGTG":"R4", "GGGTT":"RN", "GGTAA":"RO", "GGTAC":"R5",
         "GGTAG":"RP", "GGTAT":"RQ", "GGTCA":"R6", "GGTCC":"RR", "GGTCG":"RS",
         "GGTCT":"R7", "GGTGA":"RU", "GGTGC":"RV", "GGTGG":"R8", "GGTGT":"RW",
         "GGTTA":"RX", "GGTTC":"R9", "GGTTG":"RY", "GGTTT":"RZ", "GTAAA":"SB",
         "GTAAC":"SD", "GTAAG":"S0", "GTAAT":"SE", "GTACA":"SF", "GTACC":"S1",
         "GTACG":"SH", "GTACT":"SI", "GTAGA":"S2", "GTAGC":"SJ", "GTAGG":"SK",
         "GTAGT":"S3", "GTATA":"SL", "GTATC":"SM", "GTATG":"S4", "GTATT":"SN",
         "GTCAA":"SO", "GTCAC":"S5", "GTCAG":"SP", "GTCAT":"SQ", "GTCCA":"S6",
         "GTCCC":"SR", "GTCCG":"SS", "GTCCT":"S7", "GTCGA":"SU", "GTCGC":"SV",
         "GTCGG":"S8", "GTCGT":"SW", "GTCTA":"SX", "GTCTC":"S9", "GTCTG":"SY",
         "GTCTT":"SZ", "GTGAA":"7B", "GTGAC":"7D", "GTGAG":"70", "GTGAT":"7E",
         "GTGCA":"7F", "GTGCC":"71", "GTGCG":"7H", "GTGCT":"7I", "GTGGA":"72",
         "GTGGC":"7J", "GTGGG":"7K", "GTGGT":"73", "GTGTA":"7L", "GTGTC":"7M",
         "GTGTG":"74", "GTGTT":"7N", "GTTAA":"7O", "GTTAC":"75", "GTTAG":"7P",
         "GTTAT":"7Q", "GTTCA":"76", "GTTCC":"7R", "GTTCG":"7S", "GTTCT":"77",
         "GTTGA":"7U", "GTTGC":"7V", "GTTGG":"78", "GTTGT":"7W", "GTTTA":"7X",
         "GTTTC":"79", "GTTTG":"7Y", "GTTTT":"7Z", "TAAAA":"UB", "TAAAC":"UD",
         "TAAAG":"U0", "TAAAT":"UE", "TAACA":"UF", "TAACC":"U1", "TAACG":"UH",
         "TAACT":"UI", "TAAGA":"U2", "TAAGC":"UJ", "TAAGG":"UK", "TAAGT":"U3",
         "TAATA":"UL", "TAATC":"UM", "TAATG":"U4", "TAATT":"UN", "TACAA":"UO",
         "TACAC":"U5", "TACAG":"UP", "TACAT":"UQ", "TACCA":"U6", "TACCC":"UR",
         "TACCG":"US", "TACCT":"U7", "TACGA":"UU", "TACGC":"UV", "TACGG":"U8",
         "TACGT":"UW", "TACTA":"UX", "TACTC":"U9", "TACTG":"UY", "TACTT":"UZ",
         "TAGAA":"VB", "TAGAC":"VD", "TAGAG":"V0", "TAGAT":"VE", "TAGCA":"VF",
         "TAGCC":"V1", "TAGCG":"VH", "TAGCT":"VI", "TAGGA":"V2", "TAGGC":"VJ",
         "TAGGG":"VK", "TAGGT":"V3", "TAGTA":"VL", "TAGTC":"VM", "TAGTG":"V4",
         "TAGTT":"VN", "TATAA":"VO", "TATAC":"V5", "TATAG":"VP", "TATAT":"VQ",
         "TATCA":"V6", "TATCC":"VR", "TATCG":"VS", "TATCT":"V7", "TATGA":"VU",
         "TATGC":"VV", "TATGG":"V8", "TATGT":"VW", "TATTA":"VX", "TATTC":"V9",
         "TATTG":"VY", "TATTT":"VZ", "TCAAA":"8B", "TCAAC":"8D", "TCAAG":"80",
         "TCAAT":"8E", "TCACA":"8F", "TCACC":"81", "TCACG":"8H", "TCACT":"8I",
         "TCAGA":"82", "TCAGC":"8J", "TCAGG":"8K", "TCAGT":"83", "TCATA":"8L",
         "TCATC":"8M", "TCATG":"84", "TCATT":"8N", "TCCAA":"8O", "TCCAC":"85",
         "TCCAG":"8P", "TCCAT":"8Q", "TCCCA":"86", "TCCCC":"8R", "TCCCG":"8S",
         "TCCCT":"87", "TCCGA":"8U", "TCCGC":"8V", "TCCGG":"88", "TCCGT":"8W",
         "TCCTA":"8X", "TCCTC":"89", "TCCTG":"8Y", "TCCTT":"8Z", "TCGAA":"WB",
         "TCGAC":"WD", "TCGAG":"W0", "TCGAT":"WE", "TCGCA":"WF", "TCGCC":"W1",
         "TCGCG":"WH", "TCGCT":"WI", "TCGGA":"W2", "TCGGC":"WJ", "TCGGG":"WK",
         "TCGGT":"W3", "TCGTA":"WL", "TCGTC":"WM", "TCGTG":"W4", "TCGTT":"WN",
         "TCTAA":"WO", "TCTAC":"W5", "TCTAG":"WP", "TCTAT":"WQ", "TCTCA":"W6",
         "TCTCC":"WR", "TCTCG":"WS", "TCTCT":"W7", "TCTGA":"WU", "TCTGC":"WV",
         "TCTGG":"W8", "TCTGT":"WW", "TCTTA":"WX", "TCTTC":"W9", "TCTTG":"WY",
         "TCTTT":"WZ", "TGAAA":"XB", "TGAAC":"XD", "TGAAG":"X0", "TGAAT":"XE",
         "TGACA":"XF", "TGACC":"X1", "TGACG":"XH", "TGACT":"XI", "TGAGA":"X2",
         "TGAGC":"XJ", "TGAGG":"XK", "TGAGT":"X3", "TGATA":"XL", "TGATC":"XM",
         "TGATG":"X4", "TGATT":"XN", "TGCAA":"XO", "TGCAC":"X5", "TGCAG":"XP",
         "TGCAT":"XQ", "TGCCA":"X6", "TGCCC":"XR", "TGCCG":"XS", "TGCCT":"X7",
         "TGCGA":"XU", "TGCGC":"XV", "TGCGG":"X8", "TGCGT":"XW", "TGCTA":"XX",
         "TGCTC":"X9", "TGCTG":"XY", "TGCTT":"XZ", "TGGAA":"9B", "TGGAC":"9D",
         "TGGAG":"90", "TGGAT":"9E", "TGGCA":"9F", "TGGCC":"91", "TGGCG":"9H",
         "TGGCT":"9I", "TGGGA":"92", "TGGGC":"9J", "TGGGG":"9K", "TGGGT":"93",
         "TGGTA":"9L", "TGGTC":"9M", "TGGTG":"94", "TGGTT":"9N", "TGTAA":"9O",
         "TGTAC":"95", "TGTAG":"9P", "TGTAT":"9Q", "TGTCA":"96", "TGTCC":"9R",
         "TGTCG":"9S", "TGTCT":"97", "TGTGA":"9U", "TGTGC":"9V", "TGTGG":"98",
         "TGTGT":"9W", "TGTTA":"9X", "TGTTC":"99", "TGTTG":"9Y", "TGTTT":"9Z",
         "TTAAA":"YB", "TTAAC":"YD", "TTAAG":"Y0", "TTAAT":"YE", "TTACA":"YF",
         "TTACC":"Y1", "TTACG":"YH", "TTACT":"YI", "TTAGA":"Y2", "TTAGC":"YJ",
         "TTAGG":"YK", "TTAGT":"Y3", "TTATA":"YL", "TTATC":"YM", "TTATG":"Y4",
         "TTATT":"YN", "TTCAA":"YO", "TTCAC":"Y5", "TTCAG":"YP", "TTCAT":"YQ",
         "TTCCA":"Y6", "TTCCC":"YR", "TTCCG":"YS", "TTCCT":"Y7", "TTCGA":"YU",
         "TTCGC":"YV", "TTCGG":"Y8", "TTCGT":"YW", "TTCTA":"YX", "TTCTC":"Y9",
         "TTCTG":"YY", "TTCTT":"YZ", "TTGAA":"ZB", "TTGAC":"ZD", "TTGAG":"Z0",
         "TTGAT":"ZE", "TTGCA":"ZF", "TTGCC":"Z1", "TTGCG":"ZH", "TTGCT":"ZI",
         "TTGGA":"Z2", "TTGGC":"ZJ", "TTGGG":"ZK", "TTGGT":"Z3", "TTGTA":"ZL",
         "TTGTC":"ZM", "TTGTG":"Z4", "TTGTT":"ZN", "TTTAA":"ZO", "TTTAC":"Z5",
         "TTTAG":"ZP", "TTTAT":"ZQ", "TTTCA":"Z6", "TTTCC":"ZR", "TTTCG":"ZS",
         "TTTCT":"Z7", "TTTGA":"ZU", "TTTGC":"ZV", "TTTGG":"Z8", "TTTGT":"ZW",
         "TTTTA":"ZX", "TTTTC":"Z9", "TTTTG":"ZY", "TTTTT":"ZZ", "A":"B",
         "C":"D", "G":"0", "T":"E", "AA":"F", "AC":"1", "AG":"H", "AT":"I",
         "CA":"2", "CC":"J", "CG":"K", "CT":"3", "GA":"L", "GC":"M", "GG":"4",
         "GT":"N", "TA":"O", "TC":"5", "TG":"P", "TT":"Q", "AAA":"6",
         "AAC":"R", "AAG":"S", "AAT":"7", "ACA":"U", "ACC":"V", "ACG":"8",
         "ACT":"W", "AGA":"X", "AGC":"9", "AGG":"Y", "AGT":"Z", "ATA":"DB",
         "ATC":"DD", "ATG":"D0", "ATT":"DE", "CAA":"DF", "CAC":"D1",
         "CAG":"DH", "CAT":"DI", "CCA":"D2", "CCC":"DJ", "CCG":"DK",
         "CCT":"D3", "CGA":"DL", "CGC":"DM", "CGG":"D4", "CGT":"DN",
         "CTA":"DO", "CTC":"D5", "CTG":"DP", "CTT":"DQ", "GAA":"D6",
         "GAC":"DR", "GAG":"DS", "GAT":"D7", "GCA":"DU", "GCC":"DV",
         "GCG":"D8", "GCT":"DW", "GGA":"DX", "GGC":"D9", "GGG":"DY",
         "GGT":"DZ", "GTA":"0B", "GTC":"0D", "GTG":"00", "GTT":"0E",
         "TAA":"0F", "TAC":"01", "TAG":"0H", "TAT":"0I", "TCA":"02",
         "TCC":"0J", "TCG":"0K", "TCT":"03", "TGA":"0L", "TGC":"0M",
         "TGG":"04", "TGT":"0N", "TTA":"0O", "TTC":"05", "TTG":"0P",
         "TTT":"0Q", "AAAA":"06", "AAAC":"0R", "AAAG":"0S", "AAAT":"07",
         "AACA":"0U", "AACC":"0V", "AACG":"08", "AACT":"0W", "AAGA":"0X",
         "AAGC":"09", "AAGG":"0Y", "AAGT":"0Z", "AATA":"EB", "AATC":"ED",
         "AATG":"E0", "AATT":"EE", "ACAA":"EF", "ACAC":"E1", "ACAG":"EH",
         "ACAT":"EI", "ACCA":"E2", "ACCC":"EJ", "ACCG":"EK", "ACCT":"E3",
         "ACGA":"EL", "ACGC":"EM", "ACGG":"E4", "ACGT":"EN", "ACTA":"EO",
         "ACTC":"E5", "ACTG":"EP", "ACTT":"EQ", "AGAA":"E6", "AGAC":"ER",
         "AGAG":"ES", "AGAT":"E7", "AGCA":"EU", "AGCC":"EV", "AGCG":"E8",
         "AGCT":"EW", "AGGA":"EX", "AGGC":"E9", "AGGG":"EY", "AGGT":"EZ",
         "AGTA":"FB", "AGTC":"FD", "AGTG":"F0", "AGTT":"FE", "ATAA":"FF",
         "ATAC":"F1", "ATAG":"FH", "ATAT":"FI", "ATCA":"F2", "ATCC":"FJ",
         "ATCG":"FK", "ATCT":"F3", "ATGA":"FL", "ATGC":"FM", "ATGG":"F4",
         "ATGT":"FN", "ATTA":"FO", "ATTC":"F5", "ATTG":"FP", "ATTT":"FQ",
         "CAAA":"F6", "CAAC":"FR", "CAAG":"FS", "CAAT":"F7", "CACA":"FU",
         "CACC":"FV", "CACG":"F8", "CACT":"FW", "CAGA":"FX", "CAGC":"F9",
         "CAGG":"FY", "CAGT":"FZ", "CATA":"1B", "CATC":"1D", "CATG":"10",
         "CATT":"1E", "CCAA":"1F", "CCAC":"11", "CCAG":"1H", "CCAT":"1I",
         "CCCA":"12", "CCCC":"1J", "CCCG":"1K", "CCCT":"13", "CCGA":"1L",
         "CCGC":"1M", "CCGG":"14", "CCGT":"1N", "CCTA":"1O", "CCTC":"15",
         "CCTG":"1P", "CCTT":"1Q", "CGAA":"16", "CGAC":"1R", "CGAG":"1S",
         "CGAT":"17", "CGCA":"1U", "CGCC":"1V", "CGCG":"18", "CGCT":"1W",
         "CGGA":"1X", "CGGC":"19", "CGGG":"1Y", "CGGT":"1Z", "CGTA":"HB",
         "CGTC":"HD", "CGTG":"H0", "CGTT":"HE", "CTAA":"HF", "CTAC":"H1",
         "CTAG":"HH", "CTAT":"HI", "CTCA":"H2", "CTCC":"HJ", "CTCG":"HK",
         "CTCT":"H3", "CTGA":"HL", "CTGC":"HM", "CTGG":"H4", "CTGT":"HN",
         "CTTA":"HO", "CTTC":"H5", "CTTG":"HP", "CTTT":"HQ", "GAAA":"H6",
         "GAAC":"HR", "GAAG":"HS", "GAAT":"H7", "GACA":"HU", "GACC":"HV",
         "GACG":"H8", "GACT":"HW", "GAGA":"HX", "GAGC":"H9", "GAGG":"HY",
         "GAGT":"HZ", "GATA":"IB", "GATC":"ID", "GATG":"I0", "GATT":"IE",
         "GCAA":"IF", "GCAC":"I1", "GCAG":"IH", "GCAT":"II", "GCCA":"I2",
         "GCCC":"IJ", "GCCG":"IK", "GCCT":"I3", "GCGA":"IL", "GCGC":"IM",
         "GCGG":"I4", "GCGT":"IN", "GCTA":"IO", "GCTC":"I5", "GCTG":"IP",
         "GCTT":"IQ", "GGAA":"I6", "GGAC":"IR", "GGAG":"IS", "GGAT":"I7",
         "GGCA":"IU", "GGCC":"IV", "GGCG":"I8", "GGCT":"IW", "GGGA":"IX",
         "GGGC":"I9", "GGGG":"IY", "GGGT":"IZ", "GGTA":"2B", "GGTC":"2D",
         "GGTG":"20", "GGTT":"2E", "GTAA":"2F", "GTAC":"21", "GTAG":"2H",
         "GTAT":"2I", "GTCA":"22", "GTCC":"2J", "GTCG":"2K", "GTCT":"23",
         "GTGA":"2L", "GTGC":"2M", "GTGG":"24", "GTGT":"2N", "GTTA":"2O",
         "GTTC":"25", "GTTG":"2P", "GTTT":"2Q", "TAAA":"26", "TAAC":"2R",
         "TAAG":"2S", "TAAT":"27", "TACA":"2U", "TACC":"2V", "TACG":"28",
         "TACT":"2W", "TAGA":"2X", "TAGC":"29", "TAGG":"2Y", "TAGT":"2Z",
         "TATA":"JB", "TATC":"JD", "TATG":"J0", "TATT":"JE", "TCAA":"JF",
         "TCAC":"J1", "TCAG":"JH", "TCAT":"JI", "TCCA":"J2", "TCCC":"JJ",
         "TCCG":"JK", "TCCT":"J3", "TCGA":"JL", "TCGC":"JM", "TCGG":"J4",
         "TCGT":"JN", "TCTA":"JO", "TCTC":"J5", "TCTG":"JP", "TCTT":"JQ",
         "TGAA":"J6", "TGAC":"JR", "TGAG":"JS", "TGAT":"J7", "TGCA":"JU",
         "TGCC":"JV", "TGCG":"J8", "TGCT":"JW", "TGGA":"JX", "TGGC":"J9",
         "TGGG":"JY", "TGGT":"JZ", "TGTA":"KB", "TGTC":"KD", "TGTG":"K0",
         "TGTT":"KE", "TTAA":"KF", "TTAC":"K1", "TTAG":"KH", "TTAT":"KI",
         "TTCA":"K2", "TTCC":"KJ", "TTCG":"KK", "TTCT":"K3", "TTGA":"KL",
         "TTGC":"KM", "TTGG":"K4", "TTGT":"KN", "TTTA":"KO", "TTTC":"K5",
         "TTTG":"KP", "TTTT":"KQ"}


def get_id(seq):
    result = 'isomiRNA-' + str(len(seq)) + '-'
    while seq:
        result += codes[seq[:5]]
        seq = seq[5:]
    return result


def compare_chars(fastq, consensus):
    return fastq in comparison[consensus]


def compare_strings(fastq, consensus):
    fastq_len = len(fastq)
    if fastq_len != len(consensus):
        return False
    for i in range(fastq_len):
        if not compare_chars(fastq[i], consensus[i]):
            return False
    return True


def calc_trimming(seq_end, consensus_end):
    for i in range(0, min(len(seq_end), len(consensus_end))):
        if not compare_chars(seq_end[i], consensus_end[i]):
            return len(consensus_end) - i
    if len(seq_end) < len(consensus_end):
        return len(consensus_end) - len(seq_end)
    return 0


def calc_trimming_nonamb(seq_end, consensus_end):
    for i in range(0, min(len(seq_end), len(consensus_end))):
        if seq_end[i]!= consensus_end[i]:
            return len(consensus_end) - i
    if len(seq_end) < len(consensus_end):
        return len(consensus_end) - len(seq_end)
    return 0


def calc_trimming_5p(seq_end, consensus_end):
    for i in range(0, min(len(seq_end), len(consensus_end))):
        if not compare_chars(seq_end[len(seq_end) - i - 1], consensus_end[len(
                consensus_end) - i - 1]):
            return len(consensus_end) - i
    if len(seq_end) < len(consensus_end):
        return len(consensus_end) - len(seq_end)
    return 0


def calc_trimming_5p_nonamb(seq_end, consensus_end):
    for i in range(0, min(len(seq_end), len(consensus_end))):
        if seq_end[len(seq_end) - i - 1] != consensus_end[len(
                consensus_end) - i - 1]:
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


def motif_consensus_to_dict(file):
    ordered_dict = co.OrderedDict()
    handle = open(file)
    seq_records = list(SeqIO.parse(handle, 'fasta'))
    for seq_record in seq_records:
        mirna = seq_record.description.split()[0]
        motif = seq_record.description.split()[1]
        consensus = str(seq_record.seq)
        if motif in ordered_dict:
            logging.info('Motif ' + motif + ' encountered more than once in ' +
                         file + ': only first occurrence is kept')
            continue
        ordered_dict[motif] = [mirna, consensus]
    return ordered_dict


def other_motifs_pulled_seq(dict_mirna_consensus, line, motif, motif_list):
    annotation = ""
    matching_motifs = [x for x in motif_list if contains_motif(
        line, x) if x != motif]
    for matched_motif in matching_motifs:
        annotation += dict_mirna_consensus[matched_motif][0]
        annotation += " "
    return annotation


def input_name(prefix, sufix):
    answer = []
    for sample in SAMPLES:
        answer.append(prefix + sample + sufix)
    return answer


def sample_name(input, prefix, sufix):
    return input[len(prefix):-len(sufix)]


configfile:
    'config.yaml'

if config['ambiguous_letters']:
    def find_in_string(str1, str2):
        len2 = len(str2)
        for i in range(len(str1)-len2+1):
            if compare_strings(str1[i:i+len2], str2):
                return i
        return -1
else:
    def find_in_string(str1, str2):
        return str1.find(str2)


def contains_motif(line, motif):
    if find_in_string(line, motif) == -1:
        return False
    return True


def contains_ambiguous_letters(str):
    for letter in ambiguous_letters:
        if str.find(letter) != -1:
            return True
    return False


def variant(row):
    result = ''
    reference = row['REFERENCE']
    sequence = row['SEQUENCE']
    var_5p = int(row['MIRNA_START'] - row['START_G'])
    var_3p = int(row['END_G'] - row['MIRNA_END'])
    comp_start = 0
    comp_end = len(reference)
    if var_5p > 0:
        result += 'iso_5p:+' + str(var_5p) + ','
    elif var_5p < 0:
        result += 'iso_5p:' + str(var_5p) + ','
        comp_start = -var_5p
    if var_3p > 0:
        result += 'iso_3p:+' + str(var_3p) + ','
    elif var_3p < 0:
        result += 'iso_3p:' + str(var_3p) + ','
        comp_end += var_3p
    int_1_s = max(comp_start, 2)
    int_1_e = min(comp_end, 8)
    int_2_s = max(comp_start, 9)
    int_2_e = min(comp_end, 13)
    int_3_s = max(comp_start, 13)
    int_3_e = min(comp_end, 18)
    if sequence[int_1_s + var_5p:int_1_e + var_5p] != reference[
                                                      int_1_s:int_1_e]:
        result += 'iso_snp_seed,'
    if var_5p >= -8 and 8 + var_5p < len(sequence) and sequence[
        8 + var_5p] != reference[8]:
        result += 'iso_snp_central_offset,'
    if sequence[int_2_s + var_5p:int_2_e + var_5p] != reference[
                                                      int_2_s:int_2_e]:
        result += 'iso_snp_central,'
    if sequence[int_3_s + var_5p:int_3_e + var_5p] != reference[
                                                      int_3_s:int_3_e]:
        result += 'iso_snp_central_supp,'
    if comp_start <= 2 and sequence[
                           comp_start + var_5p:2 + var_5p] != reference[
                                                              comp_start:2]:
        result += 'iso_snp,'
    elif sequence[18 + var_5p:comp_end + var_5p] != reference[18:comp_end]:
        result += 'iso_snp,'
    if result:
        return result[:-1]
    return 'NA'


if config['edit_distance_3p'] == -1:
    if config['edit_distance_5p'] == -1:
        def calc_edit_distance(seq_3p, seq_5p, consensus_3p, consensus_5p,
                               delete_costs, substitute_costs, insert_costs):
            lev_3p = lev(str1=seq_3p, str2=consensus_3p,
                         delete_costs=delete_costs,
                         substitute_costs=substitute_costs,
                         insert_costs=insert_costs)
            lev_5p = lev(str1=seq_5p, str2=consensus_5p,
                         delete_costs=delete_costs,
                         substitute_costs=substitute_costs,
                         insert_costs=insert_costs)
            return lev_3p + lev_5p
    else:
        def calc_edit_distance(seq_3p, seq_5p, consensus_3p, consensus_5p,
                               delete_costs, substitute_costs, insert_costs):
            lev_5p = lev(str1=seq_5p, str2=consensus_5p,
                         delete_costs=delete_costs,
                         substitute_costs=substitute_costs,
                         insert_costs=insert_costs)
            if lev_5p > config['edit_distance_5p']:
                return -5
            lev_3p = lev(str1=seq_3p, str2=consensus_3p,
                         delete_costs=delete_costs,
                         substitute_costs=substitute_costs,
                         insert_costs=insert_costs)
            return lev_3p + lev_5p
else:
    if config['edit_distance_5p'] == -1:
        def calc_edit_distance(seq_3p, seq_5p, consensus_3p, consensus_5p,
                               delete_costs, substitute_costs, insert_costs):
            lev_3p = lev(str1=seq_3p, str2=consensus_3p,
                         delete_costs=delete_costs,
                         substitute_costs=substitute_costs,
                         insert_costs=insert_costs)
            if lev_3p > config['edit_distance_3p']:
                return -3
            lev_5p = lev(str1=seq_5p, str2=consensus_5p,
                         delete_costs=delete_costs,
                         substitute_costs=substitute_costs,
                         insert_costs=insert_costs)
            return lev_3p + lev_5p
    else:
        def calc_edit_distance(seq_3p, seq_5p, consensus_3p, consensus_5p,
                               delete_costs, substitute_costs, insert_costs):
            lev_3p = lev(str1=seq_3p, str2=consensus_3p,
                         delete_costs=delete_costs,
                         substitute_costs=substitute_costs,
                         insert_costs=insert_costs)
            if lev_3p > config['edit_distance_3p']:
                return -3
            lev_5p = lev(str1=seq_5p, str2=consensus_5p,
                         delete_costs=delete_costs,
                         substitute_costs=substitute_costs,
                         insert_costs=insert_costs)
            if lev_5p > config['edit_distance_5p']:
                return -5
            return lev_3p + lev_5p


def chr_sort(x):
    if x.lower() == 'chrx':
        return 998
    if x.lower() == 'chry':
        return 999
    return int(x[3:])


#setting edit distance from config file
insert_costs = np.ones(128, dtype=np.float64)
delete_costs = np.ones(128, dtype=np.float64)
substitute_costs = np.ones((128, 128), dtype=np.float64)
if config['insertion_score']:
    for base_ord in base_ords:
        insert_costs[base_ord] = config['insertion_score']
if config['deletion_score']:
    for base_ord in base_ords:
        delete_costs[base_ord] = config['deletion_score']
if config['substitution_AG']:
    substitute_costs[A, G] = config['substitution_AG']
if config['substitution_GA']:
    substitute_costs[G, A] = config['substitution_GA']
if config['substitution_CT']:
    substitute_costs[C, T] = config['substitution_CT']
if config['substitution_TC']:
    substitute_costs[T, C] = config['substitution_TC']
if config['substitution_AT']:
    substitute_costs[A, T] = config['substitution_AT']
if config['substitution_TA']:
    substitute_costs[T, A] = config['substitution_TA']
if config['substitution_AC']:
    substitute_costs[A, C] = config['substitution_AC']
if config['substitution_CA']:
    substitute_costs[C, A] = config['substitution_CA']
if config['substitution_GC']:
    substitute_costs[G, C] = config['substitution_GC']
if config['substitution_CG']:
    substitute_costs[C, G] = config['substitution_CG']
if config['substitution_GT']:
    substitute_costs[G, T] = config['substitution_GT']
if config['substitution_TG']:
    substitute_costs[T, G] = config['substitution_TG']
for base_ord in comparison_ord.keys():
    substitute_costs[N, base_ord] = 0
    substitute_costs[base_ord, N] = 0
for base_ord in comparison_ord.keys():
    for base_ord1 in comparison_ord.keys():
        if base_ord ^ N and base_ord1 ^ N and not (base_ord in {
            A, C, G, T} and base_ord1 in acgt_ords):
            substitute_costs[base_ord, base_ord1] = min(
                [substitute_costs[x, y] for x in
                 comparison_ord[
                     base_ord] for y in comparison_ord[base_ord1]])

input_folder = config['data']
if input_folder[-1] != '/':
    input_folder += '/'


rule all:
    input:
        expand('results/{A}.isomir.tsv', A=SAMPLES),
        expand('results/{A}.isomir.sequence_info.tsv', A=SAMPLES),
        expand('group_results/' + config['group_output_name'] + '.isomir.tsv'),
        expand('group_results/' + config['group_output_name'] + '.isomir.sequence_info.tsv'),
        expand('group_results/' + config['group_output_name'] + '.isomir.nucleotide_dist.tsv'),
        expand('results/{A}.gff', A=SAMPLES)


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
        input_files=input_folder+'{A}',
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
        # SECTION | SETUP #####################################################
        dict_mirna_consensus = motif_consensus_to_dict(input.motif_consensus)

        # dict_ambiguous = {}
        # dict_nonambiguous = {}
        # for key, value in dict_mirna_consensus.items():
        #     if contains_ambiguous_letters(key) and

        motif_list = list(dict_mirna_consensus.keys())

        dict_mirna_consensus_amb = {}
        dict_mirna_consensus_nonamb = {}
        for motif, value in dict_mirna_consensus.items():
            if contains_ambiguous_letters(value[1]):
                dict_mirna_consensus_amb[motif] = value
            else:
                dict_mirna_consensus_nonamb[motif] = value
        total_reads_in_sample = sum(1 for line in open(input.input_files)) / 4
        first_ds = True
        first_dsi = True
        first_dnd = True
        open(output[0], 'a').close()
        open(output[1], 'a').close()
        open(output[2], 'a').close()

        # SECTION | MIRNA LOOP ################################################

        table_out = {}
        freq_nt_all = {}
        with open(str(input.collapsed_fasta), "rt") as sample:
            # calculate total reads
            for line in sample:

                seq = line.rpartition(' ')[2].rstrip()
                num_reads = int(line.rpartition(' ')[0])
                ratio = float(num_reads)
                len_read = len(seq)

                for motif, value in dict_mirna_consensus_amb.items():
                    mirna = value[0]
                    consensus = value[1]

                    if mirna not in freq_nt_all.keys():
                        freq_nt_all[mirna] = co.defaultdict(lambda: co.defaultdict(int))

                    ##logging.debug("Start motif: " + motif + ' mirna: ' + mirna)

                    if not contains_motif(seq, motif):
                        continue

                    # ascertain sequences pulled in by several miRNA motifs
                    has_other = other_motifs_pulled_seq(
                        dict_mirna_consensus, line, motif, motif_list)

                    # sequence manipulations
                    consensus_index_5p = str.find(consensus, motif)
                    consensus_index_3p = consensus_index_5p + len(motif)
                    seq_index_5p = find_in_string(seq, motif)
                    seq_index_3p = seq_index_5p + len(motif)
                    consensus_end_3p = consensus[consensus_index_3p:]
                    seq_end_3p = seq[seq_index_3p:]
                    consensus_end_5p = consensus[:consensus_index_5p]
                    seq_end_5p = seq[:seq_index_5p]

                    # calculation of single-statistics
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

                    # option to check if same sequence mathches multiple mirna
                    # if destructive pull is TRUE, assign seq to mirna/motif
                    # with best distance metric
                    dist = calc_edit_distance(seq_end_3p, seq_end_5p,
                                              consensus_end_3p, consensus_end_5p,
                                              delete_costs, substitute_costs,
                                              insert_costs)
                    if dist == -3:
                        logging.warning(
                        'Skipped (3p sequencing error) ' + seq + ' ' +
                            str(len_trim) + ' ' + mirna)
                        continue
                    if dist == -5:
                        logging.warning(
                        'Skipped (5p substitution) ' + seq + ' ' + mirna)
                        continue
                    if (config['destructive_motif_pull'] and len(has_other) > 0):
                        best_matching_mirna = ""
                        for k, v in dict_mirna_consensus.items():
                            if v[0] in has_other.split(" "):
                                if (lev(seq, v[1], delete_costs = delete_costs,
                                    substitute_costs = substitute_costs,
                                    insert_costs = insert_costs) < dist):
                                    best_matching_mirna = v[0]
                        if (len(best_matching_mirna) > 0):
                            logging.warning(
                                'Skipped (' + seq + ') better matches mirna: ' +
                                best_matching_mirna)
                            continue
                    # calculation of nt frequencies at each position

                    nt_offset = seq_index_5p - consensus_index_5p
                    for index, nt in enumerate(seq):
                        freq_nt_all[mirna][nt][index - nt_offset] += num_reads
                        for nt2 in ['A', 'C', 'G', 'T', 'N']:
                            if nt2!=nt and (freq_nt_all[mirna][nt] is None or freq_nt_all[
                                mirna][nt2][index - nt_offset] is None):
                                freq_nt_all[mirna][nt2][index - nt_offset] = 0

                    # add to display queue
                    if mirna in table_out:
                        table_out[mirna].append([mirna, seq, len_read, num_reads,
                                ratio, len_trim, len_tail, seq_tail,
                                vari_5p, has_other, dist])
                    else:
                        table_out[mirna] = [[mirna, seq, len_read, num_reads,
                                ratio, len_trim, len_tail, seq_tail,
                                vari_5p, has_other, dist]]

                if 'N' in seq:
                    for motif, value in dict_mirna_consensus_nonamb.items():
                        mirna = value[0]
                        consensus = value[1]

                        if mirna not in freq_nt_all.keys():
                            freq_nt_all[mirna] = co.defaultdict(
                                lambda: co.defaultdict(int))

                        ##logging.debug("Start motif: " + motif + ' mirna: '
                        # + mirna)

                        if not contains_motif(seq, motif):
                            continue

                        # ascertain sequences pulled in by several miRNA motifs
                        has_other = other_motifs_pulled_seq(
                            dict_mirna_consensus, line, motif, motif_list)

                        # sequence manipulations
                        consensus_index_5p = str.find(consensus, motif)
                        consensus_index_3p = consensus_index_5p + len(motif)
                        seq_index_5p = find_in_string(seq, motif)
                        seq_index_3p = seq_index_5p + len(motif)
                        consensus_end_3p = consensus[consensus_index_3p:]
                        seq_end_3p = seq[seq_index_3p:]
                        consensus_end_5p = consensus[:consensus_index_5p]
                        seq_end_5p = seq[:seq_index_5p]

                        # calculation of single-statistics
                        len_trim = calc_trimming(seq_end_3p,
                                                 consensus_end_3p)
                        len_tail = calc_tailing(
                            seq_end_3p, consensus_end_3p, len_trim)
                        seq_tail = get_tailing_seq(seq, len_tail)
                        len_trim_5p = calc_trimming_5p(
                            seq_end_5p, consensus_end_5p)
                        vari_5p = max(len_trim_5p,
                                      calc_tailing(
                                          seq_end_5p, consensus_end_5p,
                                          len_trim_5p))

                        # option to check if same sequence mathches multiple
                        #  mirna
                        # if destructive pull is TRUE, assign seq to
                        # mirna/motif
                        # with best distance metric
                        dist = calc_edit_distance(seq_end_3p, seq_end_5p,
                                                  consensus_end_3p,
                                                  consensus_end_5p,
                                                  delete_costs,
                                                  substitute_costs,
                                                  insert_costs)
                        if dist == -3:
                            logging.warning(
                                'Skipped (3p sequencing error) ' + seq + ' ' +
                                str(len_trim) + ' ' + mirna)
                            continue
                        if dist == -5:
                            logging.warning(
                                'Skipped (5p substitution) ' + seq + ' ' +
                                mirna)
                            continue
                        if (config['destructive_motif_pull'] and len(
                                has_other) > 0):
                            best_matching_mirna = ""
                            for k, v in dict_mirna_consensus.items():
                                if v[0] in has_other.split(" "):
                                    if (lev(seq, v[1],
                                            delete_costs=delete_costs,
                                            substitute_costs=substitute_costs,
                                            insert_costs=insert_costs) < dist):
                                        best_matching_mirna = v[0]
                            if (len(best_matching_mirna) > 0):
                                logging.warning(
                                    'Skipped (' + seq + ') better matches '
                                                        'mirna: ' +
                                    best_matching_mirna)
                                continue

                        # calculation of nt frequencies at each position
                        nt_offset = seq_index_5p - consensus_index_5p
                        for index, nt in enumerate(seq):
                            freq_nt_all[mirna][nt][index - nt_offset] += num_reads
                            for nt2 in ['A', 'C', 'G', 'T', 'N']:
                                if nt2 != nt and (
                                        freq_nt_all[mirna][nt] is None or freq_nt_all[
                                    mirna][nt2][index - nt_offset] is None):
                                    freq_nt_all[mirna][nt2][index - nt_offset] = 0

                        # add to display queue
                        if mirna in table_out:
                            table_out[mirna].append(
                                [mirna, seq, len_read, num_reads,
                                 ratio, len_trim, len_tail, seq_tail,
                                 vari_5p, has_other, dist])
                        else:
                            table_out[mirna] = [[mirna, seq, len_read, num_reads,
                                 ratio, len_trim, len_tail, seq_tail,
                                 vari_5p, has_other, dist]]

                else:
                    for motif, value in dict_mirna_consensus_nonamb.items():
                        mirna = value[0]
                        consensus = value[1]

                        if mirna not in freq_nt_all.keys():
                            freq_nt_all[mirna] = co.defaultdict(
                                lambda: co.defaultdict(int))

                        ##logging.debug("Start motif: " + motif + ' mirna: '
                        # + mirna)

                        if motif not in seq:
                            continue

                        # ascertain sequences pulled in by several miRNA motifs
                        has_other = other_motifs_pulled_seq(
                            dict_mirna_consensus, line, motif, motif_list)

                        # sequence manipulations
                        consensus_index_5p = str.find(consensus, motif)
                        consensus_index_3p = consensus_index_5p + len(motif)
                        seq_index_5p = str.find(seq, motif)
                        seq_index_3p = seq_index_5p + len(motif)
                        consensus_end_3p = consensus[consensus_index_3p:]
                        seq_end_3p = seq[seq_index_3p:]
                        consensus_end_5p = consensus[:consensus_index_5p]
                        seq_end_5p = seq[:seq_index_5p]

                        # calculation of single-statistics
                        len_trim = calc_trimming_nonamb(seq_end_3p,
                                                 consensus_end_3p)
                        len_tail = calc_tailing(
                            seq_end_3p, consensus_end_3p, len_trim)
                        seq_tail = get_tailing_seq(seq, len_tail)
                        len_trim_5p = calc_trimming_5p_nonamb(
                            seq_end_5p, consensus_end_5p)
                        vari_5p = max(len_trim_5p,
                                      calc_tailing(
                                          seq_end_5p, consensus_end_5p,
                                          len_trim_5p))

                        # option to check if same sequence mathches multiple
                        #  mirna
                        # if destructive pull is TRUE, assign seq to
                        # mirna/motif
                        # with best distance metric
                        dist = calc_edit_distance(seq_end_3p, seq_end_5p,
                                                  consensus_end_3p,
                                                  consensus_end_5p,
                                                  delete_costs,
                                                  substitute_costs,
                                                  insert_costs)
                        if dist == -3:
                            logging.warning(
                                'Skipped (3p sequencing error) ' + seq + ' ' +
                                str(len_trim) + ' ' + mirna)
                            continue
                        if dist == -5:
                            logging.warning(
                                'Skipped (5p substitution) ' + seq + ' ' +
                                mirna)
                            continue
                        if (config['destructive_motif_pull'] and len(
                                has_other) > 0):
                            best_matching_mirna = ""
                            for k, v in dict_mirna_consensus.items():
                                if v[0] in has_other.split(" "):
                                    if (lev(seq, v[1],
                                            delete_costs=delete_costs,
                                            substitute_costs=substitute_costs,
                                            insert_costs=insert_costs) < dist):
                                        best_matching_mirna = v[0]
                            if (len(best_matching_mirna) > 0):
                                logging.warning(
                                    'Skipped (' + seq + ') better matches '
                                                        'mirna: ' +
                                    best_matching_mirna)
                                continue

                        # calculation of nt frequencies at each position
                        nt_offset = seq_index_5p - consensus_index_5p
                        for index, nt in enumerate(seq):
                            freq_nt_all[mirna][nt][index - nt_offset] += num_reads
                            for nt2 in ['A', 'C', 'G', 'T', 'N']:
                                if nt2 != nt and (
                                        freq_nt_all[mirna][nt] is None or freq_nt_all[
                                    mirna][nt2][index - nt_offset] is None):
                                    freq_nt_all[mirna][nt2][index - nt_offset] = 0

                        # add to display queue
                        if mirna in table_out:
                            table_out[mirna].append(
                                [mirna, seq, len_read, num_reads,
                                 ratio, len_trim, len_tail, seq_tail,
                                 vari_5p, has_other, dist])
                        else:
                            table_out[mirna] = [
                                [mirna, seq, len_read, num_reads,
                                 ratio, len_trim, len_tail, seq_tail,
                                 vari_5p, has_other, dist]]

    # SECTION | MOVE STATISTICS INTO DATAFRAME ################################
        for motif, value in dict_mirna_consensus.items():
            mirna = value[0]
            consensus = value[1]
            if mirna not in table_out:
                continue
            table_out_mirna = table_out[mirna]
            sequence_info = pd.DataFrame(table_out_mirna,
                                         columns=["MIRNA",
                                                  "SEQUENCE",
                                                  "LEN_READ",
                                                  "READS",
                                                  "RATIO",
                                                  "LEN_TRIM",
                                                  "LEN_TAIL",
                                                  "SEQ_TAIL",
                                                  "VAR_5P",
                                                  "MATCH",
                                                  "DISTANCE"])
            if len(table_out_mirna) > 0:
                sequence_info.sort_values(by="READS", ascending=0, inplace=1)
            else:
                continue
            # calculate total reads
            total_reads = float(sequence_info['READS'].sum())
            if total_reads == 0:
                logging.info("Mirna " + mirna + " skipped: no supporting/matched reads")
                continue

            # calculate ratio
            sequence_info['RATIO'] = sequence_info['READS'].apply(
                lambda x: round(100 * x / total_reads, 2))
            nucleotide_dist = pd.DataFrame(freq_nt_all[mirna])
            nucleotide_dist['MIRNA'] = mirna
            for base in acgtn:
                if base not in nucleotide_dist:
                    nucleotide_dist[base] = 0.0000
                else:
                    nucleotide_dist[base] = nucleotide_dist[base].fillna(0.0000)
            if len(table_out_mirna) > 0:
                nucleotide_dist['READS'] = nucleotide_dist.sum(axis=1)
                nucleotide_dist.loc[
                :, acgtn] = nucleotide_dist.loc[:, acgtn].div(
                    nucleotide_dist["READS"], axis=0)
                nucleotide_dist = np.round(nucleotide_dist, decimals=4)
                nucleotide_dist.index.name = 'NT_POSITION'
                nucleotide_dist.set_index('MIRNA', append=True, inplace=True)
                nucleotide_dist = nucleotide_dist.swaplevel(0, 1)
                nucleotide_dist = nucleotide_dist[["A", "C", "G", "T", "N",
                                                   "READS"]]

# SECTION | GENERATE SUMMARY STATISTICS ###################################
            # calculate 5' fidelity score
            vals_vari_5p = sequence_info['VAR_5P'].tolist()
            vals_reads = sequence_info['READS'].tolist()
            fidelity = 0
            for vari_5p, read in zip(vals_vari_5p, vals_reads):
                fidelity += (vari_5p * read)
            fidelity = round((fidelity / total_reads), 4)

            # calculate individual nt tailing ratios
            vals_tail_seq = sequence_info['SEQ_TAIL'].tolist()
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
            vals_len_trim = sequence_info['LEN_TRIM'].tolist()
            vals_len_tail = sequence_info['LEN_TAIL'].tolist()
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
            sequence_info = sequence_info[(sequence_info.RATIO > config['min_ratio']) |
                                          (sequence_info.READS > config['min_read'])]
            sequence_info['RATIO'] = sequence_info['RATIO'].apply(lambda x: str(x) + '%')

            # calculate number of isomirs above threshold
            total_isomirs = sequence_info.shape[0]
            if (sequence_info['LEN_TRIM'].iloc[0] == 0 and sequence_info[
                'LEN_TAIL'].iloc[0] == 0):
                total_isomirs -= 1

# SECTION | DISPLAY HEADER AND SUMMARY STATISTICS #########################
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
            summary = pd.DataFrame(summary_out,
                                   columns = ["MIRNA",
                                              "MOTIF",
                                              "CONSENSUS",
                                              "TOTAL_READS",
                                              "TOTAL_ISOMIRS",
                                              "FIDELITY_5P",
                                              "A_TAILING",
                                              "C_TAILING",
                                              "G_TAILING",
                                              "T_TAILING",
                                              "SEQUENCE_TRIMMING_ONLY",
                                              "SEQUENCE_TRIMMING",
                                              "SEQUENCE_TAILING_ONLY",
                                              "SEQUENCE_TAILING",
                                              "SEQUENCE_TRIMMING_AND_TAILING",
                                              "TOTAL_READS_IN_SAMPLE"])
            if first_ds:
                summary_all = summary.copy()
                first_ds = False
            else:
                summary_all = summary_all.append(summary)

# SECTION | DISPLAY SEQUENCES AND SINGLE STATISTICS #######################
            if first_dsi:
                sequence_info_all = sequence_info.copy()
                first_dsi = False
            else:
                sequence_info_all = sequence_info_all.append(sequence_info)
            if first_dnd:
                nucleotide_dist_all = nucleotide_dist.copy()
                first_dnd = False
            else:
                nucleotide_dist_all = nucleotide_dist_all.append(nucleotide_dist)

    # SECTION | CALCULATE CPM AND RPKM OF READS AND OUTPUT RESULTS ############
        if config['display_summary'] and not first_ds:
            with open(output[0], 'a') as out:
                summary_all.to_csv(out, sep='\t', index=False, header=True)
        if config['display_sequence_info'] and not first_dsi:
            with open(output[1], 'a') as out:
                sequence_info_all[
                    'CPM'] = sequence_info_all['READS'] / total_reads_in_sample * float(10 ** 6)
                sequence_info_all[
                    'RPKM'] = sequence_info_all['CPM'] / sequence_info_all['LEN_READ'] * float(10 ** 3)
                sequence_info_all = sequence_info_all[["MIRNA",
                                                       "SEQUENCE",
                                                       "LEN_READ",
                                                       "READS",
                                                       "RATIO",
                                                       "LEN_TRIM",
                                                       "LEN_TAIL",
                                                       "SEQ_TAIL",
                                                       "VAR_5P",
                                                       "MATCH",
                                                       "CPM",
                                                       "RPKM",
                                                       "DISTANCE"]]
                sequence_info_all.to_csv(
                    out, sep='\t', index=False, header=True)
        if config['display_nucleotide_dist'] and not first_dnd:
            with open(output[2], 'a') as out:
                nucleotide_dist_all.to_csv(
                    out, sep='\t', index=True, header=True)


rule gff_file:
    input:
        input_sample=input_folder+'{A}',
        sequence_info='results/{A}.isomir.sequence_info.tsv',
        reference_file=config['reference_file']
    output:
        'results/{A}.gff'
    log:
        os.path.join("logs/", TIMESTAMP)
    run:
        ref = pd.read_csv(input.reference_file, sep='\t')[['MIRNA',
                                                     'SEQUENCE',
                                                     'MOTIF.13',
                                                     'CHROMOSOME',
                                                     'X.COORDINATE',
                                                     'Y.COORDINATE',
                                                     'UNIQUE.MOTIF',
                                                     'PRIMIRNA',
                                                     'PRI.SEQUENCE']]

        ref.columns = ['MIRNA', 'REFERENCE', 'MOTIF.13', 'CHROMOSOME',
                       'X.COORDINATE', 'Y.COORDINATE', 'UNIQUE.MOTIF',
                       'PRIMIRNA', 'PRI.SEQUENCE']

        ref = ref.dropna()
        ref['MOTIF_START'] = ref.apply(lambda row: row['X.COORDINATE'] + row[
            'PRI.SEQUENCE'].index(row['MOTIF.13']), axis=1)
        ref['MIRNA_START'] = ref.apply(lambda row: row['X.COORDINATE'] + row[
            'PRI.SEQUENCE'].index(row['REFERENCE']), axis=1)
        ref['MIRNA_END'] = ref.apply(
            lambda row: row['MIRNA_START'] + len(row['REFERENCE']), axis=1)

        res = pd.read_csv(input.sequence_info, sep='\t')[["MIRNA",
                                                    "SEQUENCE",
                                                    "LEN_READ",
                                                    "READS",
                                                    "DISTANCE"]]

        gff = pd.merge(res, ref, left_on='MIRNA', right_on='MIRNA', how='left')
        gff['START_G'] = gff.apply(lambda row: int(row['MOTIF_START'] - row[
            'SEQUENCE'].index(row['MOTIF.13'])), axis=1)
        gff['END_G'] = gff.apply(lambda row: row['START_G'] + len(
            row['SEQUENCE']), axis=1)
        gff['SOURCE'] = 'miRBase' + str(config['mirbase_version'])
        gff['TYPE'] = gff.apply(lambda row: 'SO:0002166' if row[
            'SEQUENCE'] == row['REFERENCE'] else 'SO:0002167', axis=1)
        gff['STRAND'] = '+'
        gff['PHASE'] = '.'
        gff['CHR_NO'] = gff['CHROMOSOME'].apply(chr_sort)

        gff['START_H'] = gff.apply(lambda row: int(row['START_G'] - row[
            'X.COORDINATE']), axis=1)
        gff['END_H'] = gff.apply(lambda row: int(row['END_G'] - row[
            'X.COORDINATE']), axis=1)
        gff = gff.sort_values(['CHR_NO', 'START_G'])
        gff['POS'] = gff.apply(lambda row: row['CHROMOSOME'] + ':' +
                                           str(row['START_G']) + '-' +
                                           str(row['END_G']), axis=1)
        gff['UID'] = gff['SEQUENCE'].apply(get_id)
        no_paralogs = gff.groupby(['MIRNA', 'SEQUENCE']).apply(
            lambda x: x.POS.nunique())
        gff['NUMBER_OF_PARALOGS'] = gff.apply(lambda row: no_paralogs[
            row['MIRNA'], row['SEQUENCE']], axis=1)
        no_hits = gff.groupby(['SEQUENCE']).apply(lambda x: x.POS.nunique())
        gff['HITS'] = gff.apply(lambda row: no_hits[row['SEQUENCE']], axis=1)
        gff['VARIANT'] = gff.apply(variant, axis=1)
        gff['ATTRIBUTES'] = gff.apply(lambda row: '; '.join([
            'UID=' + row['UID'],
            'Name=' + row['MIRNA'],
            'Parent=' + row['PRIMIRNA'],
            'Variant=' + row['VARIANT'],
            'Hits=' + str(row['HITS']),
            'Genomic=' + row['POS'],
            'Expression=' + str(row['READS']),
            'Filter=Pass',
            'sequence=' + row['SEQUENCE'],
            'number_of_paralogs=' + str(row['NUMBER_OF_PARALOGS'])
        ]), axis=1)

        gff = gff[['PRIMIRNA', 'SOURCE', 'TYPE', 'START_H', 'END_H',
                   'DISTANCE', 'STRAND', 'PHASE', 'ATTRIBUTES']]
        gff.columns = ['seqID',	'source', 'type', 'start', 'end', 'score',
                       'strand', 'phase', 'attributes']

        header = '## VERSION: 1.0\n'
        header += '## source-ontology: ' + config['source_ontology'].rstrip() + '\n'
        header += '## COLDATA: ' + input.input_sample + '\n'

        with open(output[0], 'a') as w:
            w.write(header)
            gff.to_csv(w, index=False, sep='\t', header=False)



rule group_outputs:
    input:
        input_name("results/", ".isomir.tsv"),
        input_name("results/", ".isomir.sequence_info.tsv"),
        input_name("results/", ".isomir.nucleotide_dist.tsv")
    output:
        'group_results/' + config['group_output_name'] + '.isomir.tsv',
        'group_results/' + config['group_output_name'] + '.isomir.sequence_info.tsv',
        'group_results/' + config['group_output_name'] + '.isomir.nucleotide_dist.tsv'
    log:
        os.path.join("logs/", TIMESTAMP)
    run:
        prefix = "results/"
        sufix = [".isomir.tsv",
                 ".isomir.sequence_info.tsv",
                 ".isomir.nucleotide_dict.tsv"]
        condition = [config['display_group_output'] and config['display_summary'],
                     config['display_group_output'] and config['display_sequence_info'],
                     config['display_group_output'] and config['display_nucleotide_dist']]
        logging.basicConfig(
            filename=log[0],
            level=logging.DEBUG,
            format='%(levelname)s: %(message)s')
        for i in range(3):
            open(output[i], 'a').close()
            if condition[i]:
                first = True
                with open(output[i], 'a') as out:
                    for inp in input[i*len(SAMPLES) : (i+1)*len(SAMPLES)]:
                        if first:
                            try:
                                df = pd.read_csv(inp, delimiter="\t", header = 0)
                                first = False
                            except:
                                logging.warning("File " + inp + " is empty")
                                continue
                            df['SAMPLE'] = sample_name(
                            input[i*len(SAMPLES)], prefix, sufix[i])
                            cols = df.columns.tolist()
                            cols = cols[-1:] + cols[:-1]
                            df = df[cols]
                        else:
                            try:
                                df_rest = pd.read_csv(inp, delimiter="\t", header = 0)
                            except:
                                logging.warning("File " + inp + " is empty")
                                continue
                            df_rest['SAMPLE'] = sample_name(inp, prefix, sufix[i])
                            cols = df_rest.columns.tolist()
                            cols = cols[-1:] + cols[:-1]
                            df_rest = df_rest[cols]
                            df = df.append(df_rest)
                    if 'df' in vars():
                        df.to_csv(out, sep='\t', index=False, header=True)
