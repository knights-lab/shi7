from shi7.learning import get_seq_length_qual_scores
import os


def test_subsample():
    path_fastqs = [os.path.join('testfq', f) for f in os.listdir('testfq') if f.endswith('fastq')]
    avg_sequence_len, qual_scores = get_seq_length_qual_scores(path_fastqs)
    print(avg_sequence_len)
    print(qual_scores)
    assert(None, None)

test_subsample()
