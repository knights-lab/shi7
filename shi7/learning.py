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
