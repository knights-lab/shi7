import unittest
from shi7.shi7 import flash

import os
from tempfile import TemporaryDirectory

from nose.tools import assert_equals


class Flash(unittest.TestCase):
    def test(self):
        path_fastqs = [os.path.join('testfq', f) for f in os.listdir('testfq') if f.endswith('fastq')]
        with TemporaryDirectory() as outdir:
            flash(path_fastqs, outdir, 10, 100, False)
            #TODO md5checksum on files
            assert_equals(None, None)
