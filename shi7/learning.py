#!/usr/bin/env python
from shi7 import read_fastq
import os

def subsample_fastqs(path_fastqs, num_files=10, num_sequences=100):
    for i, path_fastq in enumerate(path_fastqs):
        with open(path_fastq) as fastq_inf:
            if i >= num_files:
                break
            fastq_gen = read_fastq(fastq_inf)
            yield limit_fastq(fastq_gen, num_sequences=num_sequences)


def limit_fastq(fastq_gen, num_sequences=100):
    for i in range(num_sequences):
        yield next(fastq_gen)


def get_seq_length_qual_scores(path_fastqs, num_files=10, num_sequences=1000):
    subsampled_fastqs = subsample_fastqs(path_fastqs, num_files=num_files, num_sequences=num_sequences)
    sequence_len_sum = 0.
    quality_sum = 0
    num_sequences = 0.
    # TODO: Write the subsampled fastqs
    # sequences = []
    # qualities = []
    for fastq_gen in subsampled_fastqs:
        for header, sequence, quality in fastq_gen:
            # TODO: Write the subsampled fastqs
            # sequences.append(sequences)
            # qualities.append(qualities)
            sequence_len_sum += len(sequence)
            quality_sum += sum([ord(i) for i in quality])
            num_sequences += 1.
    # Return (average length of sequences, average quality score)
    return sequence_len_sum/num_sequences, quality_sum/sequence_len_sum


def count_num_lines(path):
    with open(path) as path_inf:
        return sum(1 for line in path_inf)


def get_file_size(path):
    return os.path.getsize(path)

'''
def check_sequence_name(path):
    ids = []
    with open(path) as path_inf:
        fastq_gen = read_fastq(path_inf)
        for gen in fastq_gen:
            for header in gen: #why not header, sequence, quality?
                ids.append(header[-9])

    return ids #why the output has size of 1497
'''

def detect_paired_end(path_fastqs):
    if len(path_fastqs) % 2 == 1:
        return False
    path_fastqs = sorted(path_fastqs)
    path_R1_fastqs, path_R2_fastqs = path_fastqs[::2], path_fastqs[1::2]

    R1_lines_num = []
    R2_lines_num = []
    R1_files_size = []
    R2_files_size = []
    R1_seqs_name = []
    R2_seqs_name = []
    for path_R1_fastq in path_R1_fastqs:
        R1_lines_num.append(count_num_lines(path_R1_fastq))
        R1_files_size.append(get_file_size(path_R1_fastq))
        #R1_seqs_name.append(check_sequence_name(path_R1_fastq,'1'))
    for path_R2_fastq in path_R2_fastqs:
        R2_lines_num.append(count_num_lines(path_R2_fastq))
        R2_files_size.append(get_file_size(path_R2_fastq))
        #R2_seqs_name.append(check_sequence_name(path_R2_fastq,'2'))

    print(R1_seqs_name)
    print(R2_seqs_name)
    if not R1_lines_num == R2_lines_num or not R1_files_size == R2_files_size: #or not R1_seqs_name == R2_seqs_name:
        return False

    return True


