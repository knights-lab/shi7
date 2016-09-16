#!/usr/bin/env python
from __future__ import print_function, division
import argparse
import subprocess
import os
import re
import shutil
import multiprocessing
from datetime import datetime

import logging

def make_arg_parser():
    parser = argparse.ArgumentParser(description='This is the commandline interface for NINJA-SHI7',
                                     usage='python ninja_shi7.py -i <input> -o <output> -t_trim <threads>...')

    # parser.add_argument('-adaptor', '--adaptor_type', help='Set the type of the adaptor (default: None)', choices=[None, 'Nextera', 'TruSeq3', 'TruSeq2'], default=None)
    parser.add_argument('--axe_adaptors', help='Path to the adaptor', default=None)
    # parser.add_argument('--no_axe_adaptors', help='Disable trimmomatic axe adaptors (default: Enabled)', dest='axe_adaptors', action='store_false')
    parser.add_argument('--no_flash', help='Disable FLASH stiching (default: Enabled)', dest='flash', action='store_false')
    parser.add_argument('--no_trim', help='Disable the TRIMMER (default: Enabled)', dest='trim', action='store_false')
    parser.add_argument('--no_allow_outies', help='Disable "outie" orientation (default: Enabled)', dest='allow_outies', action='store_false')
    parser.add_argument('--no_convert_fasta', help='Disable convert FASTQS to FASTA (default: Enabled)', dest='convert_fasta', action='store_false')
    parser.add_argument('--no_combine_fasta', help='Disable the FASTA append mode (default: Enabled)', dest='combine_fasta', action='store_false')
    parser.add_argument('--shell', help='Use shell in Python system calls, NOT RECOMMENDED (default: Disabled)', dest='shell', action='store_true')
    parser.add_argument('-i', '--input', help='Set the directory path of the fastq directory', required=True)
    parser.add_argument('-o', '--output', help='Set the directory path of the output (default: cwd)', default=os.getcwd())
    parser.add_argument('-t', '--threads', help='Set the number of threads (default: %(default)s)',
                        default=multiprocessing.cpu_count())
    parser.add_argument('-m', '--min_overlap',
                        help='Set the minimum overlap length between two reads (default: %(default)s)', default=20, type=int)
    parser.add_argument('-M', '--max_overlap',
                        help='Set the maximum overlap length between two reads (default: %(default)s)', default=700, type=int)
    parser.add_argument('-trim_l', '--trim_length', help='Set the trim length (default: %(default)s)', default=150, type=int)
    parser.add_argument('-trim_q', '--trim_qual', help='Set the trim qual (default: %(default)s)', default=20, type=int)
    parser.set_defaults(flash=True, trim=True, allow_outies=True, convert_fasta=True, combine_fasta=True, shell=False)
    return parser


def run_command(cmd, shell=False):
    """
    Run prepared behave command in shell and return its output.
    :param cmd: Well-formed behave command to run.
    :param shell: Force subprocess to use shell, not recommended
    :return: Command output as string.
    """

    try:
        cmd = [str(i) for i in cmd]
        logging.debug(' '.join(cmd))
        output = subprocess.check_output(
            cmd,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            shell=shell,
            cwd=os.getcwd(),
        )

    except subprocess.CalledProcessError as e:
        output = e.output

    return output


def read_fastq(fh):
    line = next(fh)
    while line:
        title = line[1:].strip()
        data = ''
        qualities = ''
        flag = True
        line = next(fh)
        while not line[0] == '@' and line:
            if line[0] == '+':
                flag = False
            elif flag:
                data += line.strip()
            else:
                qualities += line.strip()
            line = next(fh)
        yield title, data, qualities


def split_fwd_rev(paths):
    path_R1_fastqs = [f for f in paths if 'R1' in os.path.basename(f)]
    path_R2_fastqs = [f for f in paths if 'R2' in os.path.basename(f)]
    if len(path_R1_fastqs) != len(path_R2_fastqs) or len(path_R1_fastqs) < 1:
        raise ValueError('Error: The input directory %s must contain at least one pair of R1 & R2 fastq file!' % os.path.dirname(paths[0]))
    return path_R1_fastqs, path_R2_fastqs


def axe_adaptors(input_fastqs, output_path, adapters, threads=1, shell=False):
    path_R1_fastqs, path_R2_fastqs = split_fwd_rev(input_fastqs)
    output_filenames = []
    for input_path_R1, input_path_R2 in zip(path_R1_fastqs, path_R2_fastqs):
        output_path_R1 = os.path.join(output_path, os.path.basename(input_path_R1))
        unpaired_R1 = os.path.join(output_path, 'unpaired.%s' % os.path.basename(input_path_R1))
        output_path_R2 = os.path.join(output_path, os.path.basename(input_path_R2))
        unpaired_R2 = os.path.join(output_path, 'unpaired.%s' % os.path.basename(input_path_R1))
        trim_cmd = ['trimmomatic', 'PE', input_path_R1, input_path_R2, output_path_R1, unpaired_R1,  output_path_R2, unpaired_R2, 'ILLUMINACLIP:%s:2:30:10:2:true' % adapters, '-threads', threads]
        logging.info(run_command(trim_cmd, shell=shell))
        output_filenames.append(output_path_R1)
        output_filenames.append(output_path_R2)
    return output_filenames


def flash(input_fastqs, output_path, max_overlap, min_overlap, allow_outies, threads=1, shell=False):
    path_R1_fastqs, path_R2_fastqs = split_fwd_rev(input_fastqs)
    for input_path_R1, input_path_R2 in zip(path_R1_fastqs, path_R2_fastqs):
        flash_cmd = ['flash', input_path_R1, input_path_R2, '-o', re.sub('_R1_+.fastq', '', os.path.basename(input_path_R1)), '-d', output_path, '-M', max_overlap, '-m', min_overlap, '-t', threads]
        if allow_outies:
            flash_cmd.append('-O')
        logging.info(run_command(flash_cmd, shell=shell))
    return [f for f in os.listdir(output_path) if f.endswith('extendedFrags.fastq')]


def trimmer(input_fastqs, output_path, trim_length, trim_qual, threads=1, shell=False):
    [logging.info(filename) for filename in input_fastqs]
    output_filenames = []
    for path_input_fastq in input_fastqs:
        path_output_fastq = os.path.join(output_path, re.sub('.extendedFrags.fastq', '.trimmed.fastq', os.path.basename(path_input_fastq)))
        ninja_shi7_cmd = ['ninja_shi7', path_input_fastq, path_output_fastq, trim_length, trim_qual, 'FLOOR', 5, 'ASS_QUALITY', 30]
        logging.info(run_command(ninja_shi7_cmd, shell=shell))
        output_filenames.append(path_output_fastq)
    return output_filenames


def convert_fastqs(input_fastqs, output_path):
    output_filenames = []
    for path_input_fastq in input_fastqs:
        with open(path_input_fastq) as inf_fastq:
            gen_fastq = read_fastq(inf_fastq)
            output_filename = os.path.join(output_path, re.sub('fastq', 'fna', os.path.basename(path_input_fastq)))
            with open(output_filename, 'w') as outf_fasta:
                for title, seq, quals in gen_fastq:
                    outf_fasta.write('>%s\n%s\n' % (title, seq))
        output_filenames.append(output_filename)
    return output_filenames


def convert_combine_fastqs(input_fastqs, output_path, basenames):
    output_filenames = []
    for i, path_input_fastq, basename in enumerate(zip(input_fastqs, basenames)):
        output_filename = os.path.join(output_path, 'combined_seqs.fna')
        with open(output_filename, 'w') as outf_fasta:
            with open(path_input_fastq) as inf_fastq:
                gen_fastq = read_fastq(inf_fastq)
                for title, seq, quals in gen_fastq:
                    outf_fasta.write('>%s_%i %s\n%s\n' % (basename, i, title, seq))
        output_filenames.append(output_filename)
    return output_filenames


def main():
    # gcc -m64 -O3 ninja_shi7.c -o ninja_shi7
    # gcc -m64 -O3 -msse4.1 ninja_shi7.c -o ninja_shi7
    start_time = datetime.now()

    parser = make_arg_parser()
    args = parser.parse_args()

    logging.basicConfig(filename=os.path.join(args.output, 'ninja_shi7.log'), filemode='w', level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    # FIRST CHECK IF THE INPUT AND OUTPUT PATH EXIST. IF DO NOT, RAISE EXCEPTION AND EXIT
    if not os.path.exists(args.input):
        raise ValueError('Error: Input directory %s doesn\'t exist!' % args.input)

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    if os.path.exists(os.path.join(args.output, 'temp')):
        shutil.rmtree(os.path.join(args.output, 'temp'))
        print('Existing temp directory deleted.')
        os.makedirs(os.path.join(args.output, 'temp'))
    else:
        os.makedirs(os.path.join(args.output, 'temp'))

    # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # AXE_ADAPTORS

    path_fastqs = [f for f in os.listdir(args.input) if f.endswith('fastq')]
    logging.debug(path_fastqs)

    if args.axe_adaptors:
        path_fastqs = axe_adaptors(path_fastqs, os.path.join(args.output, 'temp'), args.axe_adaptors, threads=args.threads, shell=args.shell)
        logging.info('AXE_ADAPTORS done!')

    # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # FLASH
    if args.flash:
        path_fastqs = flash(path_fastqs, os.path.join(args.output, 'temp'), args.max_overlap, args.min_overlap, args.allow_outies, threads=args.threads, shell=args.shell)
        logging.info('FLASH done!')

    # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # CREATE_TRIMMER_GENERAL

    if args.trim:
        path_fastqs = trimmer(path_fastqs, os.path.join(args.output, 'temp'), args.trim_length, args.trim_qual, threads=args.threads, shell=args.shell)
        logging.info('CREATE_TRIMMER_GENERAL done!')

    # ">[SAMPLENAME]_[INDX starting at 0] HEADER"
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # CONVERT FASTA TO FASTQ
    if args.convert_fasta:
        if args.combine_fasta:
            basenames = ['.'.join(
                re.sub('[^0-9a-zA-Z]+', '.', re.sub('_L001', '', re.sub('_001', '', os.path.basename(f)))).split('.')[:-1])
                                for f in path_fastqs]
            path_fastqs = convert_combine_fastqs(path_fastqs, os.path.join(args.output, 'temp'), basenames=basenames)
            logging.info('Convert FASTQs to FASTAs done!')
            logging.info('Combine FASTAs done!')
        else:
            path_fastqs = convert_fastqs(path_fastqs, os.path.join(args.output, 'temp'))
            logging.info('Convert FASTQs to FASTAs done!')

    [shutil.move(file, args.output) for file in path_fastqs]
    shutil.rmtree(os.path.join(args.output, 'temp'))
    logging.info('Execution time: %s' % (datetime.now() - start_time))


if __name__ == '__main__':
    main()
