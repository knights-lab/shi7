from shi7.learning import subsample_fastqs
import os
import math

import nose

def test_subsample():
    path_fastqs = [os.path.join('testfq', f) for f in os.listdir('testfq') if f.endswith('fastq')]
    subsampled_fastqs = subsample_fastqs(path_fastqs)
    sequence_len_sum = 0.
    quality_sum = 0
    num_sequences = 0.
    for fastq_gen in subsampled_fastqs:
        for header, sequence, quality in fastq_gen:
            sequence_len_sum += len(sequence)
            quality_sum += sum([ord(i) for i in quality])
            num_sequences += 1.
    print('AVG Sequence Length: %d.3' % (sequence_len_sum/num_sequences))
    print('Quality AVG: %d.3' % (quality_sum/sequence_len_sum))
    assert(None, None)

test_subsample()
