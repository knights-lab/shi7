#!/usr/bin/env python
from __future__ import print_function, division
import argparse
import subprocess
import os
import re
import shutil
import multiprocessing
import pkg_resources
from datetime import datetime

import logging

t_f_values = {"True": True, "False": False}

def convert_t_or_f(value):
    return t_f_values[value]


def make_arg_parser():

    parser = argparse.ArgumentParser(description='This is the commandline interface for shi7en',
                                     usage='shi7en -i <input> -o <output> -t_trim <threads>...')
    parser.add_argument('--debug', help='Enable debug (default: Disabled)', dest='debug', action='store_true')
    parser.add_argument('--adaptor', help='Set the type of the adaptor (default: None)', choices=[None, 'Nextera', 'TruSeq3', 'TruSeq2'], default=None)
    parser.add_argument('-SE', help='Run in Single End mode (default: Disabled)', dest='single_end', action='store_true')
    parser.add_argument('--flash', help='Enable (True) or Disable (False) FLASH stiching (default: True)', choices=[True,False].keys(), default='True', type=convert_t_or_f)
    parser.add_argument('--trim', help='Enable (True) or Disable (False) the TRIMMER (default: True)', choices=[True,False].keys(), default='True', type=convert_t_or_f)
    parser.add_argument('--allow_outies', help='Enable (True) or Disable (False) the "outie" orientation (default: True)', choices=[True,False].keys(), default='True', type=convert_t_or_f)
    parser.add_argument('--convert_fasta', help='Enable (True) or Disable (False) the conversion of FASTQS to FASTA (default: True)', choices=[True,False].keys(), default='True', type=convert_t_or_f)
    parser.add_argument('--combine_fasta', help='Enable (True) or Disable (False) the FASTA append mode (default: True)', choices=[True,False].keys(), default='True', type=convert_t_or_f)
    parser.add_argument('--shell', help='Use shell in Python system calls, NOT RECOMMENDED (default: Disabled)', dest='shell', action='store_true')
    parser.add_argument('-i', '--input', help='Set the directory path of the fastq directory', required=True)
    parser.add_argument('-o', '--output', help='Set the directory path of the output (default: cwd)', default=os.getcwd())
    parser.add_argument('-t', '--threads', help='Set the number of threads (default: %(default)s)',
                        default=multiprocessing.cpu_count())
    # TODO: Table of presets for different variable regions
    parser.add_argument('-m', '--min_overlap',
                        help='Set the minimum overlap length between two reads. If V4 set to 285 (default: %(default)s)', default=20, type=int)
    parser.add_argument('-M', '--max_overlap',
                        help='Set the maximum overlap length between two reads. If V4 set to 300 (default: %(default)s)', default=700, type=int)
    parser.add_argument('-trim_l', '--trim_length', help='Set the trim length (default: %(default)s)', default=150, type=int)
    parser.add_argument('-trim_q', '--trim_qual', help='Set the trim qual (default: %(default)s)', default=20, type=int)
    parser.set_defaults(shell=False, single_end=False)
    
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


def format_basename(filename):
    return '.'.join(re.sub('[^0-9a-zA-Z]+', '.', re.sub('_L001', '', re.sub('_001', '', os.path.basename(filename)))).split('.')[:-1])


def whitelist(dir, whitelist):
    for root, subdirs, files in os.walk(dir):
            for file in files:
                if os.path.join(root, file) not in whitelist:
                    os.remove(os.path.join(root, file))


def read_fastq(fh):
    line = next(fh)
    while line:
        title = line[1:].strip()
        data = ''
        qualities = ''
        flag = True
        line = next(fh)
        while line and (flag or len(data) != len(qualities)):
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

def resolve_adapter_path(adaptor_name, paired_end):
    adaptor_dict_PE = {'Nextera': 'NexteraPE-PE.fa', 'TruSeq2': 'TruSeq2-PE.fa', 'TruSeq3': 'TruSeq3-PE.fa', 'TruSeq32': 'TruSeq3-PE-2.fa'}
    adaptor_dict_SE = {'TruSeq2': 'TruSeq2-SE.fa', 'TruSeq3': 'TruSeq3-SE.fa'}
    if paired_end:
        adap_data = pkg_resources.resource_filename(__name__, os.path.join('adapters',adaptor_dict_PE[adaptor_name]))
    else:
        adap_data = pkg_resources.resource_filename(__name__, os.path.join('adapters',adaptor_dict_SE[adaptor_name]))
    return adap_data   


def axe_adaptors_single_end(input_fastqs, output_path, adapters, threads=1, shell=False):
    adap_data = resolve_adapter_path(adapters, 0)
    output_filenames = []
    for fastq in input_fastqs:
        output_fastq = os.path.join(output_path, format_basename(fastq) + '.fastq')
        trim_cmd = ['trimmomatic', 'SE', fastq, output_fastq, 'ILLUMINACLIP:%s:2:30:10:2:true' % adap_data, '-threads', threads]
        logging.info(run_command(trim_cmd, shell=shell))
        output_filenames.append(output_fastq)
    return output_filenames 


def axe_adaptors_paired_end(input_fastqs, output_path, adapters, threads=1, shell=False):
    adap_data = resolve_adapter_path(adapters, 1)
    path_R1_fastqs, path_R2_fastqs = split_fwd_rev(input_fastqs)
    output_filenames = []
    for input_path_R1, input_path_R2 in zip(path_R1_fastqs, path_R2_fastqs):
        output_path_R1 = os.path.join(output_path, format_basename(input_path_R1) + '.fastq')
        unpaired_R1 = os.path.join(output_path, 'unpaired.%s' % os.path.basename(input_path_R1))
        output_path_R2 = os.path.join(output_path, format_basename(input_path_R2) + '.fastq')
        unpaired_R2 = os.path.join(output_path, 'unpaired.%s' % os.path.basename(input_path_R1))
        trim_cmd = ['trimmomatic', 'PE', input_path_R1, input_path_R2, output_path_R1, unpaired_R1,  output_path_R2, unpaired_R2, 'ILLUMINACLIP:%s:2:30:10:2:true' % adap_data, '-threads', threads]
        logging.info(run_command(trim_cmd, shell=shell))
        output_filenames.append(output_path_R1)
        output_filenames.append(output_path_R2)
    return output_filenames


def flash(input_fastqs, output_path, max_overlap, min_overlap, allow_outies, threads=1, shell=False):
    # TODO: Can we run a two-pass approach?
    # TODO: Wiki table and hard code most common presets 
    path_R1_fastqs, path_R2_fastqs = split_fwd_rev(input_fastqs)
    for input_path_R1, input_path_R2 in zip(path_R1_fastqs, path_R2_fastqs):
        flash_cmd = ['flash', input_path_R1, input_path_R2, '-o', format_basename(re.sub('R1', '', input_path_R1)), '-d', output_path, '-M', max_overlap, '-m', min_overlap, '-t', threads]
        if allow_outies:
            flash_cmd.append('-O')
        logging.info(run_command(flash_cmd, shell=shell))
    output_filenames = []
    for f in os.listdir(output_path):
        if f.endswith('extendedFrags.fastq'):
            dest = re.sub('.extendedFrags', '', os.path.join(output_path, f))
            shutil.move(os.path.join(output_path, f), dest)
            output_filenames.append(dest)
    return output_filenames


def trimmer(input_fastqs, output_path, trim_length, trim_qual, threads=1, shell=False):
    [logging.info(filename) for filename in input_fastqs]
    output_filenames = []
    for path_input_fastq in input_fastqs:
        path_output_fastq = os.path.join(output_path, format_basename(path_input_fastq) + '.fastq')
        shi7en_cmd = ['shi7en_trimmer', path_input_fastq, path_output_fastq, trim_length, trim_qual, 'FLOOR', 5, 'ASS_QUALITY', 30]
        logging.info(run_command(shi7en_cmd, shell=shell))
        output_filenames.append(path_output_fastq)
    return output_filenames


def convert_fastqs(input_fastqs, output_path):
    output_filenames = []
    for path_input_fastq in input_fastqs:
        with open(path_input_fastq) as inf_fastq:
            gen_fastq = read_fastq(inf_fastq)
            output_filename = os.path.join(output_path, format_basename(path_input_fastq) + '.fna')
            with open(output_filename, 'w') as outf_fasta:
                for title, seq, quals in gen_fastq:
                    outf_fasta.write('>%s\n%s\n' % (title, seq))
        output_filenames.append(output_filename)
    return output_filenames


def convert_combine_fastqs(input_fastqs, output_path):
    output_filename = os.path.join(output_path, 'combined_seqs.fna')
    with open(output_filename, 'w') as outf_fasta:
        for path_input_fastq in input_fastqs:
                basename = format_basename(path_input_fastq)
                with open(path_input_fastq) as inf_fastq:
                    gen_fastq = read_fastq(inf_fastq)
                    for i, (title, seq, quals) in enumerate(gen_fastq):
                        outf_fasta.write('>%s_%i %s\n%s\n' % (basename, i, title, seq))
    return [output_filename]


def main():
    # gcc -m64 -O3 ninja_shi7.c -o ninja_shi7
    # gcc -m64 -O3 -msse4.1 ninja_shi7.c -o ninja_shi7
    start_time = datetime.now()

    parser = make_arg_parser()
    args = parser.parse_args()


    # FIRST CHECK IF THE INPUT AND OUTPUT PATH EXIST. IF DO NOT, RAISE EXCEPTION AND EXIT
    if not os.path.exists(args.input):
        raise ValueError('Error: Input directory %s doesn\'t exist!' % args.input)

    if not args.convert_fasta and args.combine_fasta:
        raise ValueError('Error: convert_fasta must be enabled if combine_fasta is enabled!')

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # Put in the logging file
    logging.basicConfig(filename=os.path.join(args.output, 'shi7en.log'), filemode='w', level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    if args.debug:
        logging.info('Debug Mode is Enabled. Retaining intermediate files.')

    if os.path.exists(os.path.join(args.output, 'temp')):
        shutil.rmtree(os.path.join(args.output, 'temp'))
        logging.info('Existing temp directory deleted.')
        os.makedirs(os.path.join(args.output, 'temp'))
    else:
        os.makedirs(os.path.join(args.output, 'temp'))

    # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # AXE_ADAPTORS

    path_fastqs = [os.path.join(args.input, f) for f in os.listdir(args.input) if f.endswith('fastq')]
    # TODO: Filename to samplename map?
    logging.debug(path_fastqs)

    if args.adaptor:
        axe_output = os.path.join(args.output, 'temp', 'axe')
        os.makedirs(axe_output)
        if args.single_end:
            path_fastqs = axe_adaptors_single_end(path_fastqs, axe_output, args.adaptor, threads=args.threads, shell=args.shell)
        else:
            path_fastqs = axe_adaptors_paired_end(path_fastqs, axe_output, args.adaptor, threads=args.threads, shell=args.shell)
        if not args.debug:
            whitelist(os.path.join(args.output, 'temp'), path_fastqs)
        logging.info('AXE_ADAPTORS done!')
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # FLASH
    if args.flash:
        if not args.single_end:
            flash_output = os.path.join(args.output, 'temp', 'flash')
            os.makedirs(flash_output)
            path_fastqs = flash(path_fastqs, flash_output, args.max_overlap, args.min_overlap, args.allow_outies, threads=args.threads, shell=args.shell)
            if not args.debug:
                whitelist(os.path.join(args.output, 'temp'), path_fastqs)
            logging.info('FLASH done!')
        else:
            logging.warning('Single End mode enabled with FLASH. Skipping this step.')

    # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # SHI7EN

    if args.trim:
        trimmer_output = os.path.join(args.output, 'temp', 'trimmer')
        os.makedirs(trimmer_output)
        path_fastqs = trimmer(path_fastqs, trimmer_output, args.trim_length, args.trim_qual, threads=args.threads, shell=args.shell)
        if not args.debug:
            whitelist(os.path.join(args.output, 'temp'), path_fastqs)
        logging.info('CREATE_TRIMMER_GENERAL done!')

    # ">[SAMPLENAME]_[INDX starting at 0] HEADER"
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # CONVERT FASTQ TO FASTA
    if args.convert_fasta:
        convert_output = os.path.join(args.output, 'temp', 'convert')
        os.makedirs(convert_output)
        if args.combine_fasta:
            path_fastqs = convert_combine_fastqs(path_fastqs, convert_output)
            logging.info('Convert FASTQs to FASTAs done!')
            logging.info('Combine FASTAs done!')
        else:
            path_fastqs = convert_fastqs(path_fastqs, convert_output)
            logging.info('Convert FASTQs to FASTAs done!')

    for file in path_fastqs:
        dest = os.path.join(args.output, os.path.basename(file))
        if os.path.exists(dest):
            os.remove(dest)
        shutil.move(file, args.output)
    if not args.debug:
        shutil.rmtree(os.path.join(args.output, 'temp'))
    logging.info('Execution time: %s' % (datetime.now() - start_time))


if __name__ == '__main__':
    main()
