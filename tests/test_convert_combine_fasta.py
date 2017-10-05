import unittest
from shi7.shi7 import convert_fastqs, format_basename

import os
from tempfile import TemporaryDirectory

from nose.tools import assert_equals


def read_fastq_binary(fh):
    line = next(fh)
    while line:
        title = line[1:].strip()
        data = b''
        qualities = b''
        flag = True
        line = next(fh)
        while line and (flag or len(data) != len(qualities)):
            if line[0] == b'+':
                flag = False
            elif flag:
                data += line.strip()
            else:
                qualities += line.strip()
            line = next(fh)
        yield title, data, qualities


def convert_fastqs_binary(input_fastqs, output_path):
    output_filenames = []
    for path_input_fastq in input_fastqs:
        with open(path_input_fastq, 'br') as inf_fastq:
            gen_fastq = read_fastq_binary(inf_fastq)
            output_filename = os.path.join(output_path, format_basename(path_input_fastq) + '.fna')
            with open(output_filename, 'bw') as outf_fasta:
                for title, seq, quals in gen_fastq:
                    outf_fasta.write(b'>%s\n%s\n' % (title, seq))
        output_filenames.append(output_filename)
    return output_filenames


class ConvertFastqs(unittest.TestCase):
    def test(self):
        path_fastqs = [os.path.join('testfq', f) for f in os.listdir('testfq') if f.endswith('fastq')]
        with TemporaryDirectory() as outdir:
            # convert_fastqs(path_fastqs, outdir)
            convert_fastqs_binary(path_fastqs, outdir)
            #TODO md5checksum on files
            assert_equals(None, None)
