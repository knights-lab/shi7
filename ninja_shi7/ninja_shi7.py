#!/usr/bin/env python
from __future__ import print_function, division
import argparse
import subprocess
import os
import re
import shutil
import multiprocessing
from datetime import datetime


def make_arg_parser():
    parser = argparse.ArgumentParser(description='This is the commandline interface for NINJA-SHI7',
                                     usage='python ninja_shi7.py -i <input> -o <output> -t_trim <threads>...')

    parser.add_argument('-adaptor', '--adaptor_type', help='Set the type of the adaptor (default: %(default)s)',
                        choices=['None', 'Nextera', 'TruSeq3'], default='None')
    parser.add_argument('-f', '--flash', help='Enable or disable the FLASH (default: %(default)s)',
                        choices=['enabled', 'disabled'], default='enabled')
    # path to sequences for trimming
    parser.add_argument('-trim', '--trimmer', help='Enable or disable the TRIMMER (default: %(default)s)',
                        choices=['enabled', 'disabled'], default='enabled')
    parser.add_argument('-fasta', '--make_fastas', help='Enable or disable the MAKE_FASTAS (default: %(default)s)',
                        choices=['enabled', 'disabled'], default='enabled')
    parser.add_argument('-qiime', '--QIIME',
                        help='Enable or disable the QIIME (default: %(default)s) NOTE: MAKE_FASTAS needs to be enabled in order to run QIIME.',
                        choices=['enabled', 'disabled'], default='enabled')
    parser.add_argument('-append', '--append_fna', help='Enable or disable the fna append mode (default: %(default)s)',
                        choices=['enabled', 'disabled'], default='disabled')

    parser.add_argument('-i', '--input', help='Set the directory path of the fastq directory', required=True)
    parser.add_argument('-o', '--output',
                        help='Set the directory path of the destination (default: the directory path where fastq files are located)')
    parser.add_argument('-t', '--threads', help='Set the number of threads (default: %(default)s)',
                        default=str(multiprocessing.cpu_count()))
    parser.add_argument('-O', '--allow_outies',
                        help='Enable "outie" orientation Choose one option: enabled/disabled (deafult: %(default)s)',
                        default='enabled')
    parser.add_argument('-m', '--min_overlap',
                        help='Set the minimum overlap length between two reads (default: %(default)s)', default='20')
    parser.add_argument('-M', '--max_overlap',
                        help='Set the maximum overlap length between two reads (default: %(default)s)', default='700')
    parser.add_argument('-trim_l', '--trim_length', help='Set the trim length (default: %(default)s)', default='150')
    parser.add_argument('-trim_q', '--trim_qual', help='Set the trim qual (default: %(default)s)', default='20')
    # print args.threads_trimmomatic, args.min_overlap, args.max_overlap, args.input, args.allow_outies
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


def main():
    # gcc -m64 -O3 ninja_shi7.c -o ninja_shi7
    # gcc -m64 -O3 -msse4.1 ninja_shi7.c -o ninja_shi7
    startTime = datetime.now()

    parser = make_arg_parser()
    args = parser.parse_args()

    script_path = os.path.dirname(os.path.realpath(__file__))
    print('script_path:', script_path)
    fastq_path = args.input
    print('fastq_path:', fastq_path)

    if args.output:
        fna_path = args.output
    else:
        fna_path = os.path.join(os.path.dirname(os.path.dirname(fastq_path)), 'seqs.fna')
        args.output = fna_path
    print('output_path:', fna_path)

    # FIRST CHECK IF THE INPUT AND OUTPUT PATH EXIST. IF DO NOT, RAISE EXCEPTION AND EXIT
    if not os.path.exists(fastq_path):
        print('Error:', fastq_path, 'doesn\'t exist!')
        exit()

    if args.output and not os.path.exists(os.path.dirname(args.output)):
        print('Error:', os.path.dirname(args.output), 'doesn\'t exist!')
        exit()

    if args.output and not args.output.endswith('.fna'):
        print('Error: output file must be a .fna file!')
        exit()

    if os.path.exists(os.path.join(os.path.dirname(args.output), 'temp')):
        shutil.rmtree(os.path.join(os.path.dirname(args.output), 'temp'))
        print('Existing temp directory deleted.')

    # check if MAKE_FASTAS is enabled when QIIME is enabled
    if args.QIIME == 'enabled' and args.make_fastas == 'disabled':
        print('MAKE_FASTAS needs to be ENABLED in order to run QIIME!')
        exit()

    # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # AXE_ADAPTORS

    fastqs = [f for f in os.listdir(fastq_path) if f.endswith('fastq')]
    print(fastqs)
    R1_fqs = [f for f in fastqs if re.search('R1', f)]
    R2_fqs = [f for f in fastqs if re.search('R2', f)]
    print(R1_fqs, '\n')
    print(R2_fqs, '\n')
    if len(R1_fqs) != len(R2_fqs) or len(R1_fqs) < 1:
        print('Error:', fastq_path + ':', 'The input directory must contain at least one pair of R1 & R2 fastq file!')
        exit()

    #
    if args.adaptor_type == 'None':
        fwdp_fnas = R1_fqs
    else:
        for R1, R2 in zip(R1_fqs, R2_fqs):
            trim_command = ['trimmomatic', 'PE', R1, R2, os.path.join('temp', re.sub('R1', 'fwdp', os.path.basename(R1))), os.path.join('temp', re.sub('R1', 'revp', os.path.basename(R1))), 'ILLUMINACLIP:' + args.adapters + ':2:30:10:2:true', '--threads', args.threads]
            run_command(trim_command, shell=args.shell)
        fwdp_fnas = [f for f in os.listdir(os.path.join(os.path.dirname(args.output), 'temp')) if re.search('_fwdp.+.fastq', f)]

    print('AXE_ADAPTORS done!\n')

    # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # FLASH
    print(fwdp_fnas, '\n')
    if args.flash == 'enabled':
        for f in fwdp_fnas:
            flash_command = ['flash', re.sub("fwdp", "revp", f), '-o ', os.path.join('temp', re.sub('_fwdp.+.fastq', '', os.path.basename(f))), '-M', args.max_overlap, '-m', args.min_overlap]
            if args.allow_outies == 'enabled':
                flash_command.append('-O')
            subprocess.call(flash_command, shell=True)
            trimmer_fnas = [f for f in os.listdir(os.path.join(os.path.dirname(fna_path), 'temp')) if f.endswith('.extendedFrags.fastq')]
    else:
        trimmer_fnas = fwdp_fnas

    print('FLASH done!\n')

    # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # CREATE_TRIMMER_GENERAL
    print(trimmer_fnas, '\n')
    if args.trimmer == 'enabled':
        if os.path.exists(os.path.join(fastq_path, 'temp', 'ninja_shi7_report.log')):
            os.remove(os.path.join(fastq_path, 'temp', 'ninja_shi7_report.log'))
        for f in trimmer_fnas:
            s1 = 'cd ' + os.path.join(os.path.dirname(fna_path), 'temp') + '; echo ' + f + ' >> ninja_shi7_report.log'
            s2 = os.path.join(script_path, 'ninja_shi7_mac ') + f + ' ' + re.sub('.extendedFrags.fastq',
                                                                                 '.trimmed.fastq',
                                                                                 f) + ' ' + args.trim_length + ' ' + args.trim_qual + ' ' + args.trim_qual + 'FLOOR 5 ASS_QUALITY 30 >> ninja_shi7_report.log'
            subprocess.call(s1 + ' && ' + s2, shell=True)
    else:
        for f in trimmer_fnas:
            shutil.copy(os.path.join(os.path.dirname(fna_path), 'temp', f),
                        os.path.join(os.path.dirname(fna_path), 'temp',
                                     re.sub('.extendedFrags.fastq', '.trimmed.fastq', f)))

    print('CREATE_TRIMMER_GENERAL done!\n')

    # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # MAKE_FASTAS

    if os.path.exists(os.path.join(os.path.dirname(fna_path), 'fasta')):
        shutil.rmtree(os.path.join(os.path.dirname(fna_path), 'fasta'))
    os.mkdir(os.path.join(os.path.dirname(fna_path), 'fasta'))
    FASTA = [f for f in os.listdir(os.path.join(os.path.dirname(fna_path), 'temp')) if f.endswith('.trimmed.fastq')]
    print(FASTA, '\n')
    if args.make_fastas == 'enabled':
        iter = 0
        while iter < len(FASTA):
            with open(os.path.join(os.path.dirname(fna_path), 'temp', FASTA[iter])) as fo:
                record = fo.readlines()
                trim = []
                i = 0
                while i < len(record):
                    trim.append(record[i])
                    trim.append(record[i + 1])
                    i = i + 4
                trim_record = [re.sub('^@', '>', line) for line in trim]
                trim_string = ''.join(trim_record)
                # fo = open(fastq_path + 'temp/fasta/'+re.sub('fastq','fasta',FASTA[iter]),'w')
            with open(os.path.join(os.path.dirname(fna_path), 'fasta', re.sub('fastq', 'fasta', FASTA[iter])),
                      'w') as fo2:
                # fo = open(os.path.join(fastq_path, 'temp', 'fasta', re.sub('fastq','fasta',FASTA[iter]),'w'))
                fo2.write(trim_string)
            iter += 1

    print('MAKE_FASTAS done!\n')

    # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # QIIME to .fna

    File = [f for f in os.listdir(os.path.join(os.path.dirname(fna_path), 'fasta')) if f.endswith('.fasta')]
    print(File, '\n')

    if args.QIIME == 'enabled':
        if args.append_fna == 'enabled' and os.path.exists(fna_path):
            # open to read the last sample index
            with open(fna_path) as fo:
                # with open(fna_path,'a') as fo_app: #open the fna
                last_sample_name = fo.readlines()[-2]
                print('last sample name: ', last_sample_name)
                s = re.split(' ', last_sample_name)
                s = re.split('_', s[0])
                last_idx = int(s[1])
                print('The last index of the existing fna file:', last_idx)
            sample_idx = last_idx + 1
        else:
            # fo_app = open(fna_path,'w') #open the fna
            sample_idx = 0
            print('idx:', sample_idx)

        iter2 = 0
        fastq = []
        while iter2 < len(File):
            with open(os.path.join(os.path.dirname(fna_path), 'fasta', File[iter2])) as fo_fasta:
                record = fo_fasta.readlines()
                record = [re.sub('^>', '', line) for line in record]
                sample_name = re.sub('.trimmed.fasta', '', File[iter2])
                sample_name = '>' + re.sub('[-_]', '.', sample_name) + '_'
                for idx, line in enumerate(record):
                    if idx % 2 == 0:
                        line = sample_name + str(sample_idx) + ' ' + line
                        sample_idx += 1
                    fastq.append(line)
            iter2 += 1

        fastq_string = ''.join(fastq)
        if args.QIIME == 'enabled':
            with open(fna_path, 'a') as fo_app:
                fo_app.write(fastq_string)
        else:
            with open(fna_path, 'w') as fo_app:
                fo_app.write(fastq_string)

    print('QIIME done!')
    print('Execution time:', datetime.now() - startTime)


if __name__ == '__main__':
    main()
