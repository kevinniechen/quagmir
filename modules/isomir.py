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

