#!/usr/bin/env python
from __future__ import print_function, division
import argparse
import os

def make_arg_parser():
    parser = argparse.ArgumentParser(description='This is the commandline interface for stitch_combine',
                                     usage='combine_stitched_se.py --input_se <input_se> --input_pe <input> -o <output> ...')
    parser.add_argument('--input_se', help='Set the directory path of the fasta pe file', required=True)
    parser.add_argument('--input_pe', help='Set the directory path of the fasta se file', required=True)
    parser.add_argument('-o', '--output', help='Set the directory path of the output to place the combined_seqs.all.fna (default: cwd)', default=os.getcwd())
    return parser


def read_fasta(fh):
    """
    :return: tuples of (title, seq)
    """
    title = None
    data = None
    for line in fh:
        if line[0] == ">":
            if title:
                yield (title.strip(), data)
            title = line[1:]
            data = ''
        else:
            data += line.strip()
    if not title:
        yield None
    yield title.strip(), data
    
def combine_stitched_se(output_file, pe_inf, se_inf):
    paired_success = set()

    with open(pe_inf, "r") as inf:
        with open(output_file, "w") as outf_fasta:
            gen_fasta = read_fasta(inf)
            for title, seq in gen_fasta:
                    # Grab the name of the sequence
                    paired_success.add(title.split()[0])
                    outf_fasta.write('>%s\n%s\n' % (title, seq))
    
    num_not_stitched = 0
    with open(se_inf, "r") as inf:
        with open(output_file, "a") as outf_fasta:
            gen_fasta = read_fasta(inf)
            for title, seq in gen_fasta:
                title_check = title.split()[0].replace(".R1", "").replace(".R2", "")
                if not title_check in paired_success:
                    num_not_stitched += 1
                    print(title)
                    outf_fasta.write('>%s\n%s\n' % (title, seq))
    print("Number not stitched %d" % num_not_stitched)

def main():
    parser = make_arg_parser()
    args = parser.parse_args()
    
    combine_stitched_se(args.output, args.input_pe, args.input_se)

if __name__ == '__main__':
    main()
