from shi7.learning import *
import os
import shutil


def test_subsample():
    # Fetch example fastqs
    path_fastqs = [os.path.join('testfq', f) for f in os.listdir('testfq') if f.endswith('fastq')]
    # Create a fake output directory
    output_path = os.path.join("testfq", "temp", "learning")

    if os.path.exists(output_path):
        shutil.rmtree(output_path)
    os.makedirs(output_path)
    avg_sequence_len, qual_scores = get_seq_length_qual_scores(path_fastqs, output_path, num_sequences=200, num_files=3)
    shutil.rmtree(output_path)
    assert(None, None)


def test_axe_adaptors():
    path_fastqs = [os.path.join('testfq', f) for f in os.listdir('testfq') if f.endswith('fastq')]
    subsampled_fastqs = os.path.join("testfq", "temp", "subsampled_fastqs")
    output_path = os.path.join("testfq", "temp", "axe_adaptors")
    if os.path.exists(subsampled_fastqs):
        shutil.rmtree(subsampled_fastqs)
    if os.path.exists(output_path):
        shutil.rmtree(output_path)
    os.makedirs(subsampled_fastqs)
    os.makedirs(output_path)
    get_seq_length_qual_scores(path_fastqs, subsampled_fastqs, num_sequences=200)
    path_subsampled_fastqs = [os.path.join(subsampled_fastqs, f) for f in os.listdir(subsampled_fastqs) if f.endswith('fastq')]
    best_adap, best_size = choose_axe_adaptors(path_subsampled_fastqs, output_path)
    print("Best adaptor:", best_adap)
    print("Best adaptor size:", best_size)


test_subsample()
test_axe_adaptors()
