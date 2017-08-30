#!/usr/bin/env python
from __future__ import print_function, division
import multiprocessing
import os
import csv
import datetime
import logging
from datetime import datetime
import argparse

from shi7 import VERSION, read_fastq, axe_adaptors_single_end, axe_adaptors_paired_end, flash_part1, flash_part2, split_fwd_rev

#TODO: Finish the inner array SD and Mean


def make_arg_parser():
    parser = argparse.ArgumentParser(description='This is the commandline interface for shi7_learning',
                                     usage='shi7_learning v{version} -i <input> -o <output> ...'.format(version=VERSION))
    parser.add_argument('-i', '--input', help='Set the directory path of the fastq directory OR oligos.txt if splitting', required=True)
    parser.add_argument('-o', '--output', help='Set the directory path of the output (default: cwd)', default=os.getcwd())
    parser.add_argument('-t', '--threads', help='Set the number of threads (default: %(default)s)',
                        default=min(multiprocessing.cpu_count(), 16))
    parser.set_defaults()
    return parser

# template.format(input=input, adapters=adapters)

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
    path_fastqs = [f for f in path_fastqs if f.endswith('.fastq') or f.endswith('.fq')]
    if len(path_fastqs) % 2 == 1:
        return False
    path_fastqs = sorted(path_fastqs)
    path_R1_fastqs, path_R2_fastqs = path_fastqs[::2], path_fastqs[1::2]
    print(path_R2_fastqs)
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
    return sum([get_file_size(os.path.join(path, fastq)) for fastq in os.listdir(path)])


def remove_directory_contents(path):
    for f in os.listdir(path):
        os.remove(os.path.join(path, f))

def choose_axe_adaptors(path_subsampled_fastqs, output_path):
    adapters = ['Nextera', 'TruSeq2', 'TruSeq3', 'TruSeq3-2']
    threads = min(multiprocessing.cpu_count(),16)
    original_size = get_directory_size(os.path.dirname(path_subsampled_fastqs[0]))
    print('Original size of the subsampled_fastqs = ', original_size)
    best_size = original_size
    best_adap = None
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


def flash_stitchable_and_check_outies(adapter_output_filenames, flash_output_path):
    num_threads = min(multiprocessing.cpu_count(),16)
    flash_output_str = flash_part1(adapter_output_filenames, flash_output_path, max_overlap=700, min_overlap=20, allow_outies=True, threads=num_threads, shell=False)

    allow_outies_count = 0
    for flash_out in flash_output_str:
        flash_str_list = flash_out.strip().split('\n')
        outies_info = flash_str_list[-8]
        outies_percent = float(outies_info[outies_info.find('(')+1:outies_info.find('%')])
        print('outies_percent:', outies_percent)
        if outies_percent >= 15:
            allow_outies_count += 1

    path_flash_fqs = flash_part2(flash_output_str, flash_output_path)
    path_R1_fastqs, _ = split_fwd_rev(adapter_output_filenames)
    print('flash files:', path_flash_fqs)
    print('original files:', path_R1_fastqs)

    matched_count = 0
    for original_fq, flash_fq in zip(path_R1_fastqs, path_flash_fqs):
        print('Flash number of lines =', count_num_lines(flash_fq), '\nOriginal number of lines =', count_num_lines(original_fq))
        if count_num_lines(flash_fq) > count_num_lines(original_fq)*0.3:
            matched_count = matched_count + 1

    return matched_count/len(path_flash_fqs) >= 0.75, allow_outies_count/len(flash_output_str) >= 0.75


def flash_check_cv(flash_output_path):
    hist_files = [os.path.join(flash_output_path, f) for f in os.listdir(flash_output_path) if f.endswith('.hist')]
    total_cv = 0
    for f in hist_files:
        with open(f) as inf:
            csv_inf = csv.reader(inf, delimiter="\t")
        rows = [row for row in csv_inf]
        x2f = 0
        sum = 0
        cnt = 0
        for row in rows:
            cnt = cnt + row[1]
            sum = sum + row[0] * row[1]
            x2f = x2f + row[0] * row[0] * row[1]
        mean = sum/cnt
        std = math.sqrt((x2f - sum*sum/cnt)/cnt-1)
        cv = std/mean
        total_cv = total_cv + cv
    total_files = len(hist_files)
    return total_cv/hist_len


def trimmer_learning(flash_output_filenames):
    filter_q_sum = 0
    trim_q_sum = 0
    length = 0
    for fq_path in flash_output_filenames:
        with open(fq_path) as fq_inf:
            fq_gen = read_fastq(fq_inf)
            for gen in fq_gen:
                qualities = gen[2]
                length = length + len(qualities)
                qualities = [ord(qual) for qual in qualities]
                filter_q_sum = filter_q_sum + sum(qualities)
                trim_q_sum = trim_q_sum + sum(qualities[:10]) + sum(qualities[-10:])

    print('filter_q_sum:', filter_q_sum)
    print('trim_q_sum:', trim_q_sum)
    print('total length:', length)
    print('filter_q:', filter_q_sum/length)
    print('trim_q:', trim_q_sum/length)

    return round(filter_q_sum/length), round(trim_q_sum/length)

def template_input(input):
    input = os.path.abspath(input)
    # input, input_cmd
    return "input\t{}".format(input), "--input {}".format(input)

def template_paired_end(bool):
    # bool, paired_end
    if bool:
        cmd = ""
    else:
        cmd = "-SE"
    return "paired_end\t{}".format(str(bool)), cmd

def template_output(output):
    # output, output_cmd
    output = os.path.abspath(output)
    return "output\t{}".format(output), "--output {}".format(output)

def main():
    start_time = datetime.now()

    parser = make_arg_parser()
    args = parser.parse_args()

    #learning_params = ['input_cmd', 'paired_end_cmd', 'output_cmd']
    learning_params = ["shi7.py"]
    #learning_dict = dict(zip(learning_params, [""]*len(learning_params)))

    input = os.path.abspath(args.input)
    output = os.path.abspath(args.output)

    if not os.path.exists(output):
        os.makedirs(output)

    # Put in the logging file
    logging.basicConfig(filename=os.path.join(output, 'shi7_learning.log'), filemode='w', level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    logger = logging.getLogger()

    # Record the input
    results, addon = template_input(input)
    logger.info(results)
    if addon:
        learning_params.extend(addon)

    path_fastqs = [os.path.join(input, f) for f in os.listdir(args.input) if f.endswith('fastq') or f.endswith('fq')]

    if len(path_fastqs) == 0:
        msg = "No FASTQS found in input folder {}".format(input)
        logger.critical(msg)
        raise IOError(msg)

    # Detect if paired end
    results, addon = template_paired_end(detect_paired_end(path_fastqs))
    logger.info(results)
    if addon:
        learning_params.extend(addon)

    # Detect output folder
    results, addon = template_output(output)
    logger.info(results)
    if addon:
        learning_params.extend(addon)

    with open(os.path.join(args.output, "shi7_cmd.sh"), "w") as output:
        cmd = " ".join(learning_params)
        output.write(cmd)
        print(cmd)

if __name__ == "__main__":
    main()
