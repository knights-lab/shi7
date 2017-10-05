#!/usr/bin/env python
from __future__ import print_function, division
import argparse
import subprocess
import os
import re
import sys
import shutil
import multiprocessing
import pkg_resources
from datetime import datetime

import logging

from shi7 import __version__

TRUE_FALSE_DICT = {
              "true": True,
              "false": False,
              "t": True,
              "f": False,
              "1": True,
              "0": False
              }


def convert_t_or_f(value):
    value = value.lower()
    return TRUE_FALSE_DICT[value]

def make_arg_parser():
    # TODO: Preset modes will get precedence over default values, but lose to explicit settings from user
    parser = argparse.ArgumentParser(description='This is the commandline interface for shi7',
                                     usage='shi7 %s -i <input> -o <output> ...' % __version__)
    parser.add_argument('--gotta_split', help='Split one giant fastq (or one pair of R1/R2) into 1 fastq per sample', dest='split', choices=[True,False], default='False', type=convert_t_or_f)
    parser.add_argument('--gotta_split_output', help='output directory for the newly-split fastqs')
    parser.add_argument('--gotta_split_r1', help='r1 to split by sample names in oligos.txt')
    parser.add_argument('--gotta_split_r2', help='r2 to split by sample names in oligos.txt')
    parser.add_argument('--debug', help='Retain all intermediate files (default: Disabled)', dest='debug', action='store_true')
    parser.add_argument('--adaptor', help='Set the type of the adaptors to remove (default: None)', choices=['None', 'Nextera', 'TruSeq3', 'TruSeq2', 'TruSeq3-2'], default=None)
    parser.add_argument('-SE', help='Run in Single End mode (default: Disabled)', dest='single_end', action='store_true')
    parser.add_argument('--flash', help='Enable (True) or Disable (False) FLASH stitching (default: True)', choices=[True, False], default='True', type=convert_t_or_f)
    parser.add_argument('--trim', help='Enable (True) or Disable (False) the TRIMMER (default: True)', choices=[True, False], default='True', type=convert_t_or_f)
    parser.add_argument('-outies','--allow_outies', help='Enable (True) or Disable (False) the "outie" orientation (default: True)', choices=[True, False], default='True', type=convert_t_or_f)
    parser.add_argument('--convert_fasta', help='Enable (True) or Disable (False) the conversion of FASTQS to FASTA (default: True)', choices=[True, False], default='True', type=convert_t_or_f)
    parser.add_argument('--combine_fasta', help='Enable (True) or Disable (False) the FASTA append mode (default: True)', choices=[True, False], default='True', type=convert_t_or_f)
    #parser.add_argument('--shell', help='Use shell in Python system calls, NOT RECOMMENDED (default: Disabled)', dest='shell', action='store_true')
    parser.add_argument('-i', '--input', help='Set the directory path of the fastq directory OR oligos.txt if splitting', required=True)
    parser.add_argument('-o', '--output', help='Set the directory path of the output (default: cwd)', default=os.getcwd())
    parser.add_argument('-t', '--threads', help='Set the number of threads (default: %(default)s)',
                        default=min(multiprocessing.cpu_count(),16))
    parser.add_argument('-s', '--strip_delim', help='Prune names at delim d\'s nth occurance (usage -s d,n)',default="",type=str)
    # TODO: Table of presets for different variable regions
    parser.add_argument('-m', '--min_overlap',
                        help='Set the minimum overlap length between two reads. If V4 w/primers, try 285 (default: %(default)s)', default=10, type=int)
    parser.add_argument('-M', '--max_overlap',
                        help='Set the maximum overlap length between two reads. If V4 w/primers, try 300 (default: %(default)s)', default=700, type=int)
    parser.add_argument('-filter_l', '--filter_length', help='Set the length of reads to retain (default: %(default)s)', default=80, type=int)
    parser.add_argument('-filter_q', '--filter_qual', help='Set the avg quality of the reads to retain (default: %(default)s)', default=35, type=int)
    parser.add_argument('-trim_q', '--trim_qual', help='Trim read ends until they reach trim_q (default: %(default)s)', default=32, type=int)
    parser.add_argument('--drop_r2', help='When combining FASTAs, drop R2 reads and remove "R1" from read name (default: False)', choices=[True, False], default='False', type=convert_t_or_f)
    parser.set_defaults(shell=False, single_end=False)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
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
    return re.sub('[^0-9a-zA-Z]+', '.', '.'.join(os.path.basename(filename).split('.')[:-1]))


def whitelist(dir, whitelist):
    for root, subdirs, files in os.walk(dir):
            for file in files:
                if os.path.join(root, file) not in whitelist:
                    os.remove(os.path.join(root, file))


def read_fastq(fh):
    # Assume linear FASTQS
    while True:
        title = next(fh)
        while title[0] != '@':
            title = next(fh)
        # Record begins
        if title[0] != '@':
            raise IOError('Malformed FASTQ files; verify they are linear and contain complete records.')
        title = title[1:].strip()
        sequence = next(fh).strip()
        garbage = next(fh).strip()
        if garbage[0] != '+':
            raise IOError('Malformed FASTQ files; verify they are linear and contain complete records.')
        qualities = next(fh).strip()
        if len(qualities) != len(sequence):
            raise IOError('Malformed FASTQ files; verify they are linear and contain complete records.')
        yield title, sequence, qualities

def split_fwd_rev(paths):
    if len(paths) % 2: 
        raise ValueError('Sequence number is uneven; these aren\'t pairs!')
    return paths[::2], paths[1::2]


def resolve_adapter_path(adaptor_name, paired_end):
    adaptor_dict_PE = {'Nextera': 'NexteraPE-PE.fa', 'TruSeq2': 'TruSeq2-PE.fa', 'TruSeq3': 'TruSeq3-PE.fa', 'TruSeq3-2': 'TruSeq3-PE-2.fa'}
    adaptor_dict_SE = {'TruSeq2': 'TruSeq2-SE.fa', 'TruSeq3': 'TruSeq3-SE.fa', 'Nextera': 'NexteraPE-PE.fa', 'TruSeq3-2': 'TruSeq3-PE-2.fa'}
    if paired_end:
        adap_data = pkg_resources.resource_filename(__name__, os.path.join('adapters', adaptor_dict_PE[adaptor_name]))
    else:
        adap_data = pkg_resources.resource_filename(__name__, os.path.join('adapters', adaptor_dict_SE[adaptor_name]))
    return adap_data


def axe_adaptors_single_end(input_fastqs, output_path, adapters, threads=1, shell=False):
    # Can we automate the adaptors?
    adap_data = resolve_adapter_path(adapters, 0)
    output_filenames = []
    for fastq in input_fastqs:
        output_fastq = os.path.join(output_path, os.path.basename(fastq))
        trim_cmd = ['trimmomatic', 'SE', '-threads', threads, fastq, output_fastq, 'ILLUMINACLIP:%s:2:30:10:2:true' % adap_data]
        logging.info(run_command(trim_cmd, shell=shell))
        output_filenames.append(output_fastq)
    return output_filenames


def axe_adaptors_paired_end(input_fastqs, output_path, adapters, threads=1, shell=False):
    adap_data = resolve_adapter_path(adapters, 1)
    path_R1_fastqs, path_R2_fastqs = split_fwd_rev(input_fastqs)
    output_filenames = []
    for input_path_R1, input_path_R2 in zip(path_R1_fastqs, path_R2_fastqs):
        name_R1 = os.path.basename(input_path_R1)
        name_R2 = os.path.basename(input_path_R2)
        output_path_R1 = os.path.join(output_path, name_R1)
        unpaired_R1 = os.path.join(output_path, 'unpaired.%s' % name_R1)
        output_path_R2 = os.path.join(output_path, name_R2)
        unpaired_R2 = os.path.join(output_path, 'unpaired.%s' % name_R2)
        trim_cmd = ['trimmomatic', 'PE', '-threads', threads, input_path_R1, input_path_R2, output_path_R1, unpaired_R1, output_path_R2, unpaired_R2, 'ILLUMINACLIP:%s:2:20:5:1:true' % adap_data]
        logging.info(run_command(trim_cmd, shell=shell))
        output_filenames.append(output_path_R1)
        output_filenames.append(output_path_R2)
    return output_filenames


def flash_part1(input_fastqs, output_path, max_overlap, min_overlap, allow_outies, threads=1, shell=False):
    # TODO: Can we run a two-pass approach?
    # TODO: Wiki table and hard code most common presets
    flash_output_str = []
    path_R1_fastqs, path_R2_fastqs = split_fwd_rev(input_fastqs)
    for input_path_R1, input_path_R2 in zip(path_R1_fastqs, path_R2_fastqs):
        name_R1 = os.path.basename(input_path_R1)
        name_R2 = os.path.basename(input_path_R2)
        k = name_R1.rfind(".R1")
        new = name_R1[:k] + name_R2[k+3:]
        flash_cmd = ['flash', input_path_R1, input_path_R2, '-o', new, '-d', output_path, '-M', max_overlap, '-m', min_overlap, '-t', threads]
        if allow_outies:
            flash_cmd.append('-O')
        flash_output_str.append(run_command(flash_cmd, shell=shell))

    return flash_output_str


def flash_part2(flash_output, output_path, log = True):
    if log: [logging.info(s) for s in flash_output]
    output_filenames = []
    for f in os.listdir(output_path):
        f = os.path.join(output_path,f)
        if f.endswith('.extendedFrags.fastq'):
            dest = f[:-20]
            shutil.move(f, dest)
            output_filenames.append(dest)
    return output_filenames


def trimmer(input_fastqs, output_path, filter_length, trim_qual, filter_qual, threads=1, shell=False):
    output_filenames = []
    for path_input_fastq in input_fastqs:
        fq_name = os.path.basename(path_input_fastq)
        path_output_fastq = os.path.join(output_path, fq_name)
        shi7_cmd = ['shi7_trimmer', path_input_fastq, path_output_fastq, filter_length, trim_qual, 'FLOOR', 5, 'ASS_QUALITY', filter_qual]
        logging.info(run_command(shi7_cmd, shell=shell))
        output_filenames.append(path_output_fastq)
    return output_filenames


def convert_fastqs(input_fastqs, output_path):
    output_filenames = []
    for path_input_fastq in input_fastqs:
        with open(path_input_fastq, 'r') as inf_fastq:
            gen_fastq = read_fastq(inf_fastq)
            outfile = os.path.basename(path_input_fastq)[:-3]
            output_filename = os.path.join(output_path, outfile + '.fna')
            with open(output_filename, 'w') as outf_fasta:
                for title, seq, quals in gen_fastq:
                    outf_fasta.write('>%s\n%s\n' % (title, seq))
        output_filenames.append(output_filename)
    return output_filenames


def convert_combine_fastqs(input_fastqs, output_path, drop_r2=False):
    output_filename = os.path.join(output_path, 'combined_seqs.fna')
    with open(output_filename, 'w') as outf_fasta:
        for path_input_fastq in input_fastqs:
                basename = os.path.basename(path_input_fastq)[:-3]
                with open(path_input_fastq) as inf_fastq:
                    gen_fastq = read_fastq(inf_fastq)
                    for i, (title, seq, quals) in enumerate(gen_fastq):
                        if drop_r2: 
                            if basename.endswith('.R2'): continue
                            outf_fasta.write('>%s_%i %s\n%s\n' % (basename[:-3], i, title, seq))
                        else: outf_fasta.write('>%s_%i %s\n%s\n' % (basename, i, title, seq))
    return [output_filename]


def gotta_split(r1_input, r2_input, input_fastqs, output_path):
    gotta_split_cmd = ['gotta_split', r1_input, r2_input, input_fastqs, output_path]
    logging.info(run_command(gotta_split_cmd))


def strip_delim(path_fastqs, token):
    stripped_names = [os.path.basename(n) for n in path_fastqs]
    if token == "": return stripped_names
    sp_token = token.split(",")
    if (len(sp_token) != 2 or (int)(sp_token[1]) < 1):
        raise ValueError("ERROR: strip_delim arg must be 'd,n' (d=delimeter, n=which one).")
    delim = sp_token[0]
    count = int(sp_token[1])
    logging.info("Splitting names on occurance #%d of delimiter: %s" % (count,delim))
    path = os.path.dirname(path_fastqs[0])
    for i in range(0,len(stripped_names)):
        name = stripped_names[i]
        splits = name.split(delim)
        if len(splits) > count:
            stripped_names[i] = delim.join(splits[0:count]) + ".fq"
    return stripped_names


def match_pairs(path_fastqs, doSE):
    orig = set([os.path.basename(p) for p in path_fastqs])
    nfiles = len(path_fastqs)
    if not doSE and nfiles % 2: 
        raise ValueError("ERROR: Odd number of fastq files in paired-end mode.")
    Pat1 = [".R1","_R1","-R1","R1",".1","_1","-1",".0","_0","-0",".F","_F","-F"]
    Pat2 = [".R2","_R2","-R2","R2",".2","_2","-2",".1","_1","-1",".R","_R","-R"]
    path = os.path.dirname(path_fastqs[0])
    for i in range(0,len(Pat1)):
        p1 = Pat1[i]
        p2 = Pat2[i]
        matched = 0
        where = []
        for j in range(0,nfiles):
            n = os.path.basename(path_fastqs[j])
            Ixs = [m.start() for m in re.finditer(p1,n)]
            if len(Ixs) == 0: continue
            for k in range(0,len(Ixs)):
                trans = n[0:Ixs[k]] + p2 + n[(Ixs[k]+len(p1)):]
                #print("Found pattern " + p1 + " in " + n + " at # " + str(Ixs[k]))
                if trans in orig:
                    where.append(Ixs[k])
                    matched = matched + 1
                    break
        if matched == nfiles // 2:
            matched = 0
            minp = min(where) 
            plen = len(p1)
            neworig = []
            nwhere = []
            for k in range(0,nfiles):
                if matched >= nfiles // 2: break
                n = os.path.basename(path_fastqs[k])
                ix = where[matched]
                if n[ix:ix+plen] != p1: continue
                trans = os.path.join(path, n[0:ix] + p2 + n[ix+plen:])
                neworig.append(path_fastqs[k])
                neworig.append(trans)
                nwhere += [ix,ix]
                matched = matched + 1
            if matched == nfiles // 2: 
                if doSE: 
                    logging.warning("WARNING: SE is toggled on apparently paired reads.")
                    if nfiles % 2: neworig = sorted(path_fastqs)
                return [neworig, p1, p2, nwhere]
    if not doSE: raise ValueError("ERROR: No known pattern reliably distinguishes mate pairs.")
    return [path_fastqs, None, None, None]


def link_manicured_names(orig_paths, snames, subdir, doSE, delimCtr):
    nfiles = len(snames)
    ctr = 0
    Names = []
    Rep2 = ["_R1","_R2"]
    which = delimCtr[2]
    oset = set(snames)
    ambig = len(snames) > len(oset)
    plen = len(delimCtr[0])
    nmatches = len(which)
    if ambig or not doSE:
        for i in range(0,nfiles):
            n = snames[i] 
            p = delimCtr[ctr]
            ix = which[i] if i < nmatches else -1
            #print("File %d ('%s') supposedly has delim %s at pos %d" % (i,n,p,ix))
            if ix < 0 or len(n) <= ix + 3: Names.append(n[0:n.rfind('.')]+Rep2[ctr])
            else:
                mate = n[0:ix]+delimCtr[not ctr]+n[ix+plen:]
                if mate not in oset:
                    raise ValueError('No pair found: %s, mate %s' % (n,mate))
                n = n[0:ix]+n[ix+plen:]
                Names.append(n[0:n.rfind('.')]+Rep2[ctr])
            ctr = not ctr
        Names = [re.sub('[^0-9a-zA-Z]+', '.', e)+".fq" for e in Names]
        for x in range(nfiles):
            bn = os.path.basename(orig_paths[x])
            logging.info("%s --> %s" % (bn,Names[x]))
    else: Names = snames

    if (len(Names) > len(set(Names))):
        raise ValueError('ERROR: processed file names are not unique.')
    Names = [os.path.join(subdir,n) for n in Names]
    for i in range(0,nfiles): os.link(orig_paths[i],Names[i])
    return Names


def main():
    # gcc -m64 -O3 ninja_shi7.c -o ninja_shi7
    start_time = datetime.now()

    parser = make_arg_parser()
    args = parser.parse_args()

    # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Perform validation of input parameters

    outdir = os.path.join(args.output)
    tmpdir = os.path.join(outdir,"temp")
    if not os.path.exists(args.input):
        raise ValueError('ERROR: Input directory %s doesn\'t exist!' % args.input)

    if not args.convert_fasta and args.combine_fasta:
        raise ValueError('ERROR: convert_fasta must be enabled if combine_fasta is enabled!')

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if args.split:
        if not args.r1 and not args.r2:
            raise ValueError('ERROR: arguments -r1 and -r2 are required for split!')
        elif not args.r1:
            raise ValueError('ERROR: argument -r1 is required for split!')
        elif not args.r2:
            raise ValueError('ERROR: argument -r2 is required for split!')
        if not os.path.exists(args.r1):
            raise ValueError('ERROR: r1 directory %s doesn\'t exist!' % args.r1)
        if not os.path.exists(args.r2):
            raise ValueError('ERROR: r1 directory %s doesn\'t exist!' % args.r2)
        if args.gotta_split_output:
            if not os.path.exists(args.gotta_split_output):
                raise ValueError('ERROR: Gotta_split output directory %s doesn\'t exist!' % args.gotta_split_output)

    # Initialize the logging file
    logging.basicConfig(filename=os.path.join(outdir, 'shi7.log'), filemode='w', level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    logging.info('SHI7 invocation: %s' % " ".join(sys.argv))

    if args.drop_r2 and args.flash: 
        logging.warning('WARNING: Not dropping R2 reads because flash is enabled')
        args.drop_r2 = False

    if args.debug:
        logging.info('Debug Mode is Enabled. Retaining intermediate files.')

    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)
        logging.info('Existing temp directory deleted.')
        os.makedirs(tmpdir)
    else:
        os.makedirs(tmpdir)
    

    # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Demuliplex reads

    if args.split:
        if args.gotta_split_output:
            splitted_output = args.gotta_split_output
        else:
            splitted_output = os.path.join(os.path.dirname(args.input), 'split_fastqs') #the args.input here is the oligos path

        gotta_split(args.gotta_split_r1, args.gotta_split_r2, args.input, splitted_output)
        args.input = splitted_output
        logging.info('Gotta Split done!')

    # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Prepare file paths/names for subsequent processing

    path_fastqs = [os.path.join(args.input, f) for f in os.listdir(args.input) if f.endswith('fastq') or f.endswith('fq')]

    # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Match pairs appropriately

    pp_paths = match_pairs(path_fastqs, args.single_end)
    if (pp_paths[1]==None): 
        raise ValueError("No pattern found for distinguishing mate pairs. Try -SE")
    logging.info("Detected pairs match on delimiter %s" % pp_paths[1])
    path_fastqs = pp_paths[0]
    snames = strip_delim(path_fastqs, args.strip_delim)
    link_outdir = os.path.join(tmpdir, 'link')
    os.makedirs(link_outdir)
    path_fastqs = link_manicured_names(path_fastqs, snames, link_outdir, args.single_end, pp_paths[1:])
    if not args.debug:
        whitelist(link_outdir, path_fastqs)

    
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # AXE_ADAPTORS
    # TODO: Filename to samplename map?

    if args.adaptor and args.adaptor != str(None):
        axe_output = os.path.join(tmpdir, 'axe')
        os.makedirs(axe_output)
        if args.single_end:
            path_fastqs = axe_adaptors_single_end(path_fastqs, axe_output, args.adaptor, threads=args.threads, shell=args.shell)
        else:
            path_fastqs = axe_adaptors_paired_end(path_fastqs, axe_output, args.adaptor, threads=args.threads, shell=args.shell)
        if not args.debug:
            whitelist(tmpdir, path_fastqs)
        logging.info('Adaptor removal complete!')

    # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # FLASH
    if args.flash:
        if not args.single_end:
            flash_output = os.path.join(tmpdir, 'flash')
            os.makedirs(flash_output)
            flash_output_str = flash_part1(path_fastqs, flash_output, args.max_overlap, args.min_overlap, args.allow_outies, threads=args.threads, shell=args.shell)
            path_fastqs = flash_part2(flash_output_str, flash_output)
            if not args.debug:
                whitelist(tmpdir, path_fastqs)
            logging.info('Stitching of paired ends with FLASH complete!')
        else:
            logging.warning('Single End mode enabled. Skipping FLASH stitching.')

    # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Trimmer

    if args.trim:
        trimmer_output = os.path.join(tmpdir, 'trimmer')
        os.makedirs(trimmer_output)
        path_fastqs = trimmer(path_fastqs, trimmer_output, args.filter_length, args.trim_qual, args.filter_qual, threads=args.threads, shell=args.shell)
        if not args.debug:
            whitelist(tmpdir, path_fastqs)
        logging.info('Final quality trimming complete!')

    # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # CONVERT FASTQ TO FASTA

    if args.convert_fasta:
        drop_r2 = args.drop_r2
        convert_output = os.path.join(tmpdir, 'convert')
        os.makedirs(convert_output)
        if not args.debug and args.combine_fasta:
            path_fastqs = convert_combine_fastqs(path_fastqs, convert_output, drop_r2=drop_r2)
        elif args.combine_fasta:
            convert_fastqs(path_fastqs, convert_output)
            path_fastqs = convert_combine_fastqs(path_fastqs, convert_output, drop_r2=drop_r2)
        else:
            path_fastqs = convert_fastqs(path_fastqs, convert_output)
        logging.info('Convert ' + ('and combine ' if args.combine_fasta else '') + 'FASTQs done!')

    for file in path_fastqs:
        dest = os.path.join(outdir, os.path.basename(file))
        if os.path.exists(dest):
            os.remove(dest)
        shutil.move(file, outdir)
    if not args.debug:
        shutil.rmtree(tmpdir)
    logging.info('Execution time: %s' % (datetime.now() - start_time))


if __name__ == '__main__':
    main()
