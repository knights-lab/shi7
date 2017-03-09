from shi7.learning import get_seq_length_qual_scores
import os
import shutil


def test_subsample():
    path_fastqs = [os.path.join('testfq', f) for f in os.listdir('testfq') if f.endswith('fastq')]
    output_path = os.path.join("testfq", "temp", "learning")
    if os.path.exists(output_path):
        shutil.rmtree(output_path)
    os.makedirs(output_path)
    avg_sequence_len, qual_scores = get_seq_length_qual_scores(path_fastqs, output_path, num_sequences=200, num_files=3)
    shutil.rmtree(output_path)
    assert(None, None)

test_subsample()
