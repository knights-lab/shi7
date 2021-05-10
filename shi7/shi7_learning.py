#!/usr/bin/env python
from __future__ import print_function, division
import multiprocessing
import os
import csv
import datetime
import logging
from datetime import datetime
import argparse
import shutil
import math
from glob import glob
import gzip
from shi7 import __version__

from shi7.shi7 import TRUE_FALSE_DICT, read_fastq, axe_adaptors_single_end, axe_adaptors_paired_end, flash_part1, \
    flash_part2, split_fwd_rev, match_pairs, link_manicured_names

def make_arg_parser():
    parser = argparse.ArgumentParser(description='This is the commandline interface for shi7_learning',
                                     usage='shi7_learning v{version}\nshi7_learning.py -i <input> -o <output> ...'.format(version=__version__))
    parser.add_argument('-i', '--input', help='Set the directory path of the fastq directory OR oligos.txt if splitting', required=True)
    parser.add_argument('-o', '--output', help='Set the directory path of the output (default: cwd)', default=os.getcwd())
    parser.add_argument('--debug', help='Retain all intermediate files (default: Disabled)', dest='debug', action='store_true')
    parser.add_argument('-t', '--threads', help='Set the number of threads (default: %(default)s)',
                        default=min(multiprocessing.cpu_count(), 16))
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    parser.set_defaults()
    return parser


def subsample_fastqs(path_fastqs, num_files=10, num_sequences=1000):
    for i, path_fastq in enumerate(path_fastqs):
        if i >= num_files:
            return
        with open(path_fastq) as fastq_inf:
            fastq_gen = read_fastq(fastq_inf)
            yield limit_fastq(fastq_gen, num_sequences=num_sequences)


def limit_fastq(fastq_gen, num_sequences=1000):
    for i in range(num_sequences):
        try:
            yield next(fastq_gen)
        except StopIteration:
            return


def get_seq_length_qual_scores(path_fastqs, output_path, num_files=10, num_sequences=1000):
    subsampled_fastqs = subsample_fastqs(path_fastqs, num_files=num_files, num_sequences=num_sequences)
    sequence_len_sum = 0.
    quality_sum = 0
    num_sequences = 0.

    for fastq_path, fastq_gen in zip(path_fastqs, subsampled_fastqs):
        with open(os.path.join(output_path, os.path.basename(fastq_path)), 'w') as outf:
            for header, sequence, quality in fastq_gen:
                outf.write("@%s\n%s\n+\n%s\n" % (header, sequence, quality))
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
    path_fastqs = [f for f in path_fastqs if f.endswith('.fastq') or f.endswith('.fq') or f.endswith('.fastq.gz') or f.endswith('.fq.gz')]
    if len(path_fastqs) % 2 == 1: return False, [path_fastqs, None, None, None]
    pair_obj = match_pairs(path_fastqs, True)
    path_fastqs = pair_obj[0]
    if pair_obj[1]==None: return False, pair_obj
    return True, pair_obj

def get_directory_size(path):
    return sum([get_file_size(os.path.join(path, fastq)) for fastq in os.listdir(path)])

def remove_directory_contents(path):
    for f in os.listdir(path):
        os.remove(os.path.join(path, f))

def choose_axe_adaptors(path_subsampled_fastqs, paired_end, output_path, threads):
    adapters = ['TruSeq2', 'TruSeq3', 'TruSeq3-2', 'Nextera']
    threads = min(threads, multiprocessing.cpu_count(), 16)
    original_size = get_directory_size(os.path.dirname(path_subsampled_fastqs[0]))
    logging.info('Original size of the subsampled_fastqs = ' + str(original_size))
    best_size = original_size
    best_adap = None
    for adapter in adapters:
        if paired_end:
            axe_adaptors_paired_end(path_subsampled_fastqs, output_path, adapter, threads, shell=False)
        else:
            axe_adaptors_single_end(path_subsampled_fastqs, output_path, adapter, threads, shell=False)
        fastqs_path_size = get_directory_size(output_path)
        logging.info("Adapters: {adapter}\tFile Size: {filesize}".format(adapter=adapter, filesize=fastqs_path_size))
        if fastqs_path_size <= best_size:
            best_size = fastqs_path_size
            best_adap = adapter

    if best_size < 0.995*original_size:
        # Actually write the best files again for use in later steps
        logging.info("Best Adapters: {adapter}\tFile Size: {filesize}".format(adapter=best_adap, filesize=best_size))

        if paired_end:
            files = axe_adaptors_paired_end(path_subsampled_fastqs, output_path, best_adap, threads, shell=False)
        else:
            files = axe_adaptors_single_end(path_subsampled_fastqs, output_path, best_adap, threads, shell=False)
        return best_adap, best_size, files
    else:
        return None, original_size, path_subsampled_fastqs


def flash_stitchable_and_check_outies(adapter_output_filenames, flash_output_path, threads):
    flash_output_str = flash_part1(adapter_output_filenames, flash_output_path, max_overlap=700, \
        min_overlap=10, allow_outies=True, threads=threads, shell=False)

    allow_outies_count = 0
    for flash_out in flash_output_str:
        flash_str_list = flash_out.strip().split('\n')
        outies_info = flash_str_list[-8]
        outies_percent = float(outies_info[outies_info.find('(')+1:outies_info.find('%')])
        if outies_percent >= 15:
            allow_outies_count += 1

    path_flash_fqs = flash_part2(flash_output_str, flash_output_path)
    path_R1_fastqs, _ = split_fwd_rev(adapter_output_filenames)

    matched_count = 0
    for original_fq, flash_fq in zip(path_R1_fastqs, path_flash_fqs):
        if count_num_lines(flash_fq) > count_num_lines(original_fq)*0.3:
            matched_count = matched_count + 1

    return matched_count/len(path_flash_fqs) >= 0.75, allow_outies_count/len(flash_output_str) >= 0.75, path_flash_fqs


def flash_check_cv(flash_output_path):
    hist_files = [os.path.join(flash_output_path, f) for f in os.listdir(flash_output_path) if f.endswith('.hist')]
    total_cv = total_mean = 0
    for f in hist_files:
        with open(f) as inf:
            csv_inf = csv.reader(inf, delimiter="\t")
            x2f = 0
            sum = 0
            cnt = 0
            for row in csv_inf:
                row = [int(r) for r in row]
                cnt = cnt + row[1]
                sum = sum + row[0] * row[1]
                x2f = x2f + row[0] * row[0] * row[1]
            mean = sum/cnt
            std = math.sqrt((x2f - sum*sum/cnt)/(cnt-1))
            cv = std/mean
            total_cv = total_cv + cv
            total_mean = total_mean + mean
    total_files = len(hist_files)
    return total_cv/total_files, total_mean/total_files


def trimmer_learning(flash_output_filenames):
    filter_q_sum = 0
    trim_q_sum = 0
    totbases = 0
    tottrim = 0
    num = 0
    for fq_path in flash_output_filenames:
        with open(fq_path) as fq_inf:
            fq_gen = read_fastq(fq_inf)
            for gen in fq_gen:
                num = num + 1
                qualities = gen[2]
                totbases = totbases + len(qualities)
                qualities = [ord(qual)-33 for qual in qualities]
                filter_q_sum = filter_q_sum + sum(qualities)
                if (len(qualities) >= 20):
                    trim_q_sum = trim_q_sum + sum(qualities[:10]) + sum(qualities[-10:])
                    tottrim = tottrim + 20
    logging.info('num seqs: %d' % num)
    logging.info('filter_q_sum: %d' % filter_q_sum)
    logging.info('trim_q_sum: %d' % trim_q_sum)
    logging.info('total bases considered: %d (trim: %d)' % (totbases, tottrim))
    logging.info('filter_q: %d' % (filter_q_sum/totbases))
    logging.info('trim_q: %d' % (trim_q_sum/tottrim))

    filter_q = math.floor(filter_q_sum/totbases)
    trim_q = math.floor(trim_q_sum/tottrim)-1
    trim_q = trim_q if trim_q > filter_q - 3 else filter_q - 3

    return filter_q, trim_q

def template_input(input):
    input = os.path.abspath(input)
    # input, input_cmd
    return "input\t{}".format(input), ["--input", input]

def template_paired_end(bool):
    # bool, paired_end
    if bool:
        return "paired_end\t{}".format(str(bool)), None
    else:
        return "paired_end\t{}".format(str(bool)), ["-SE"]

def template_trim(filt_q, trim_q):
    return "filt_q: %d, trim_q: %d" % (filt_q, trim_q), ["--filter_qual", str(filt_q), "--trim_qual", str(trim_q)]

def template_cv(minstitch, maxstitch):
    return "minstitch: %d, maxstitch: %d" % (minstitch, maxstitch), ["--min_overlap", str(minstitch), "--max_overlap", str(maxstitch)]

def template_output(output):
    # output, output_cmd
    output = os.path.abspath(output)
    return "output\t{}".format(output), ["--output", output]

def template_choose_axe_adaptors(best_adapt, best_size):
   if best_adapt:
       return "axe_adaptors\t" + best_adapt, ["--adaptor", best_adapt]
   else:
       return "axe_adaptors\tNA", ["--adaptor", "None"]

def template_flash(stitches, do_outies):
    return "stitches: %s, outies: %s" % (stitches, do_outies), ["--flash", str(stitches), "--allow_outies", str(do_outies)]

def main():
    start_time = datetime.now()

    parser = make_arg_parser()
    args = parser.parse_args()

    learning_params = ["shi7.py"]
    learning_pretty = ["SHI7 version", __version__]

    input = os.path.abspath(args.input)
    output = os.path.abspath(args.output)

    # Make output folder
    if not os.path.exists(output):
        os.makedirs(output)

    # Put in the logging file
    logging.basicConfig(filename=os.path.join(output, 'shi7_learning.log'), filemode='w', level=logging.DEBUG, \
        format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    # Make temp outfolder
    if os.path.exists(os.path.join(args.output, 'temp')):
        shutil.rmtree(os.path.join(args.output, 'temp'))
        logging.info('Existing temp directory deleted.')
        os.makedirs(os.path.join(args.output, 'temp'))
    else:
        os.makedirs(os.path.join(args.output, 'temp'))


    path_fastqs = [os.path.join(input, f) for f in os.listdir(input) if f.endswith('fastq') or f.endswith('fq') or f.endswith('fq.gz') or f.endswith('fastq.gz')]

    if len(path_fastqs) == 0:
        msg = "No FASTQS found in input folder {}".format(input)
        logging.critical(msg)
        raise IOError(msg)

    # Record the input
    results, addon = template_input(input)
    logging.info(results)
    if addon:
        learning_params.extend(addon)

    # Write temp subsampled fastqs
    subsampled_fastq_path = os.path.join(output, 'temp', 'subsampled')
    os.makedirs(subsampled_fastq_path)
    totbases = totseqs = 0
    for file in path_fastqs:
        basename = os.path.basename(file)
        if(file.endswith('.fastq') or file.endswith('.fq')):
            fastq_inf = open(file)
        else:
            fastq_inf = gzip.open(file, 'rt')
        fastq_gen = read_fastq(fastq_inf)
        if(basename.endswith('.gz')):
            basename = basename[:-3]
        with open(os.path.join(subsampled_fastq_path, basename), 'w') as outf:
            for header, seq, quality in limit_fastq(fastq_gen):
                outf.write("@{header}\n{seq}\n+\n{quality}\n".format(header=header, seq=seq, quality=quality))
                totbases += len(seq)
                totseqs += 1
    avlen = totbases/totseqs
    path_fastqs = glob(os.path.join(subsampled_fastq_path , "*"))



    # Detect if paired end
    paired_end, pair_obj = detect_paired_end(path_fastqs)
    path_fastqs = pair_obj[0]
    link_outdir = os.path.join(output, 'temp', 'link')
    os.makedirs(link_outdir)
    snames = [os.path.basename(n) for n in path_fastqs]
    path_fastqs = link_manicured_names(path_fastqs, snames, link_outdir, not paired_end, pair_obj[1:])

    results, addon = template_paired_end(paired_end)
    logging.info(results)
    if addon: learning_params.extend(addon)
    learning_pretty += ["Paired end",paired_end]

    # Detect adapters
    axe_adaptors_path = os.path.join(output, 'temp', 'axe_adaptors')
    os.makedirs(axe_adaptors_path)
    best_adap, best_size, fastq_paths = choose_axe_adaptors(path_fastqs, paired_end, axe_adaptors_path, int(args.threads))
    results, addon = template_choose_axe_adaptors(best_adap, best_size)
    logging.info(results)
    if addon: learning_params.extend(addon)
    learning_pretty += ["Detected adaptors",best_adap]

    # Detect output folder
    results, addon = template_output(output)
    logging.info(results)
    if addon: learning_params.extend(addon)

    # Detect stitching
    stitched_path = os.path.join(output, 'temp', 'flash')
    os.makedirs(stitched_path)
    if paired_end:
        stitches, do_outies, fastq_paths = flash_stitchable_and_check_outies(fastq_paths, stitched_path, int(args.threads))
    else: stitches, do_outies = False, False
    results, addon = template_flash(stitches, do_outies)
    logging.info(results)
    if addon: learning_params.extend(addon)
    if paired_end:
        learning_pretty += ["Stitching",stitches]
        if stitches: learning_pretty += ["Outies allowed",do_outies]

    filt_q, trim_q = trimmer_learning(fastq_paths)
    results, addon = template_trim(int(filt_q), int(trim_q))
    logging.info(results)
    if addon: learning_params.extend(addon)
    learning_pretty += ["Filter quality",filt_q,"Trimming quality",trim_q]

    # Check whether to implement stitching bounds
    if stitches:
        cv, mean = flash_check_cv(stitched_path)
        if cv < 0.1:
            learning_pretty += ["Amplicon mode",True]
            logging.info("CV: %f, Mean: %f, Avlen: %f" % (cv, mean, avlen))
            if avlen > mean: avlen = mean
            mr = math.ceil(cv*mean)
            logging.info("SD was: %d" % mr)
            minstitch, maxstitch = int(2*avlen - mean-mr), int(2*avlen - mean+mr)
            if minstitch < 8: minstitch = 8
            logging.info("Amplicon mode: stitch range [%d, %d]" % (minstitch, maxstitch))
            results, addon = template_cv(minstitch, maxstitch)
            logging.info(results)
            if addon: learning_params.extend(addon)
            learning_pretty += ["Amplicon stitch minimum",minstitch]
            learning_pretty += ["Amplicon stitch maximum",maxstitch]
        else: learning_pretty += ["Amplicon mode",False]

    #print(str(learning_params))
    with open(os.path.join(args.output, "shi7_cmd.sh"), "w") as output:
        cmd = " ".join(learning_params)
        output.write(cmd)
        print(cmd)

    with open(os.path.join(args.output, "learning_params.txt"),"w") as output:
        for ix in range(0,len(learning_pretty),2):
            output.write(str(learning_pretty[ix]) + "\t" + str(learning_pretty[ix+1]) + "\n")

    if not args.debug:
        shutil.rmtree(os.path.join(args.output, 'temp'))
    logging.info('Execution time: %s' % (datetime.now() - start_time))

if __name__ == "__main__":
    main()
