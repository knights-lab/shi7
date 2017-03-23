#!/usr/bin/env python
from shi7 import *
import multiprocessing
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


def get_seq_length_qual_scores(path_fastqs, output_path, num_files=10, num_sequences=1000):
    subsampled_fastqs = subsample_fastqs(path_fastqs, num_files=num_files, num_sequences=num_sequences)
    sequence_len_sum = 0.
    quality_sum = 0
    num_sequences = 0.

    # sequences = []
    # qualities = []
    for fastq_path, fastq_gen in zip(path_fastqs, subsampled_fastqs):
        with open(os.path.join(output_path, os.path.basename(fastq_path)), 'w') as outf:
            for header, sequence, quality in fastq_gen:
                outf.write("@%s\n%s\n+\n%s\n" % (header, sequence, quality))
                # sequences.append(sequences)
                # qualities.append(qualities)
                sequence_len_sum += len(sequence)
                quality_sum += sum([ord(i) for i in quality])
                num_sequences += 1.
    # Return (average length of sequences, average quality score)
    return sequence_len_sum/num_sequences, quality_sum/sequence_len_sum

# def get_all_fastq_seqs_qual(path_fastqs):
#     for i, path_fastq in enumerate(path_fastqs):
#         with open(path_fastq) as fastq_inf:
#             fastq_gen = read_fastq(fastq_inf)
#             yield next(fastq_gen)


def count_num_lines(path):
    with open(path) as path_inf:
        return sum(1 for line in path_inf)


def get_file_size(path):
    return os.path.getsize(path)


def check_sequence_name(path_R1, path_R2):
    with open(path_R1) as path_inf_R1, open(path_R2) as path_inf_R2:
        fastq_gen_R1 = read_fastq(path_inf_R1)
        fastq_gen_R2 = read_fastq(path_inf_R2)
        for gen_R1, gen_R2 in zip(fastq_gen_R1,fastq_gen_R2):
            title_R1, title_R2 = gen_R1[0], gen_R2[0]
            if len(title_R1) != len(title_R2):
                return False
            diff_idx = [i for i in range(len(title_R1)) if title_R1[i] != title_R2[i]]
            if len(diff_idx) != 1:
                return False
            if int(title_R2[diff_idx[0]]) - int(title_R1[diff_idx[0]]) != 1:
                return False
    return True


def detect_paired_end(path_fastqs):
    #path_fastqs = [f for f in path_fastqs if f.endswith('.fastq') or f.endswith('.fq')]
    if len(path_fastqs) % 2 == 1:
        return False
    path_fastqs = sorted(path_fastqs)
    path_R1_fastqs, path_R2_fastqs = path_fastqs[::2], path_fastqs[1::2]
    if len(path_R1_fastqs) != len(path_R2_fastqs) or len(path_R1_fastqs) < 1:
        return False
    R1_lines_num = []
    R2_lines_num = []
    R1_files_size = []
    R2_files_size = []
    seqs_name = []
    for path_R1_fastq in path_R1_fastqs:
        R1_lines_num.append(count_num_lines(path_R1_fastq))
        R1_files_size.append(get_file_size(path_R1_fastq))
    for path_R2_fastq in path_R2_fastqs:
        R2_lines_num.append(count_num_lines(path_R2_fastq))
        R2_files_size.append(get_file_size(path_R2_fastq))
    for path_R1_fastq, path_R2_fastq in zip(path_R1_fastqs, path_R2_fastqs):
        seqs_name.append(check_sequence_name(path_R1_fastq, path_R2_fastq))
    if not R1_lines_num == R2_lines_num or not R1_files_size == R2_files_size or False in seqs_name:
        return False
    return True


def get_directory_size(path):
    total_size = 0
    for fastq in os.listdir(path):
        total_size += get_file_size(os.path.join(path, fastq))
    return total_size


def remove_directory_contents(path):
    for f in os.listdir(path):
        os.remove(os.path.join(path, f))


def choose_axe_adaptors(path_subsampled_fastqs, output_path):
    # assume we already have the directory that contains subsample fastqs
    adapters = ['Nextera', 'TruSeq2', 'TruSeq3', 'TruSeq3-2']
    threads = min(multiprocessing.cpu_count(),16)
    print('directory name: ', os.path.dirname(path_subsampled_fastqs[0]))
    original_size = get_directory_size(os.path.dirname(path_subsampled_fastqs[0]))
    print('Original size = ', original_size)
    best_size = original_size
    best_adap = None
    print(path_subsampled_fastqs)
    if detect_paired_end(path_subsampled_fastqs):
        print("Entering paired end ...")
        for adapter in adapters:
            axe_adaptors_paired_end(path_subsampled_fastqs, output_path, adapter, threads, shell=False)
            fastqs_path_size = get_directory_size(output_path)
            print(adapter, fastqs_path_size)
            if fastqs_path_size < best_size:
                best_size = fastqs_path_size
                best_adap = adapter
            remove_directory_contents(output_path)
    else:
        print("Entering single end ...")
        for adapter in adapters:
            axe_adaptors_single_end(path_subsampled_fastqs, output_path, adapter, threads, shell=False)
            fastqs_path_size = get_directory_size(output_path)
            print(adapter, fastqs_path_size)
            if fastqs_path_size < best_size:
                best_size = fastqs_path_size
                best_adap = adapter
            remove_directory_contents(output_path)

    if best_size < 0.995*original_size:
        return best_adap, best_size
    else:
        return None, original_size

