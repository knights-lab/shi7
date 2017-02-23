from shi7.learning import read_fastq
import os
import unittest

from nose.tools import assert_equals


class TestFastqGen(unittest.TestCase):
    def test(self):
        path_fastqs = [os.path.join('testfq', f) for f in os.listdir('testfq') if f.endswith('fastq')]
        i = 0
        j = 0
        for path_fastq in path_fastqs:
            with open(path_fastq) as inf:
                fastq_gen = read_fastq(inf)
                for title, data, qualities in fastq_gen:
                    i += 1
            with open(path_fastq) as inf:
                for line in inf:
                    j += 1
        assert i == j/4
