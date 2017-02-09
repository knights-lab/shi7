#!/usr/bin/env python
from shi7 import read_fastq


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
