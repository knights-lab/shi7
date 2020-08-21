import unittest
from shi7.shi7 import convert_fastqs, format_basename, read_fastq, convert_fastqs

import os
from tempfile import TemporaryDirectory

from nose.tools import assert_equals

class ConvertFastqs(unittest.TestCase):
    def test(self):
        path_fastqs = [os.path.join('testfq', f) for f in os.listdir('testfq') if f.endswith('fastq')]
        with TemporaryDirectory() as outdir:
            convert_fastqs(path_fastqs, outdir)
            
            # TODO md5checksum on files
            assert_equals(None, None)
