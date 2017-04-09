from shi7.learning import *
from shi7.shi7 import *
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


def test_axe_adaptors(subsampled_fastqs):
    path_fastqs = [os.path.join('testfq', f) for f in os.listdir('testfq') if f.endswith('fastq')]
    output_path = os.path.join("testfq", "temp", "axe_adaptors")
    if os.path.exists(output_path):
        shutil.rmtree(output_path)
    os.makedirs(output_path)
    get_seq_length_qual_scores(path_fastqs, subsampled_fastqs, num_sequences=200) #To create the subsampled_fastqs
    path_subsampled_fastqs = [os.path.join(subsampled_fastqs, f) for f in os.listdir(subsampled_fastqs) if f.endswith('fastq')]
    best_adap, best_size = choose_axe_adaptors(path_subsampled_fastqs, output_path)
    print("Best adaptor:", best_adap)
    print("Best adaptor size:", best_size)
    return best_adap


def test_flash(adapter_output_filenames):
    print(adapter_output_filenames)
    flash_output_path = os.path.join("testfq", "temp", "flash")
    if os.path.exists(flash_output_path):
        shutil.rmtree(flash_output_path)
    os.makedirs(flash_output_path)
    is_stitchable = flash_is_stitchable(adapter_output_filenames, flash_output_path)
    print(is_stitchable)
    if is_stitchable:
        flash_check_cv(flash_output_path)
    # incomplete


if __name__ == '__main__':
    test_subsample()
    subsampled_fastqs = os.path.join("testfq", "temp", "subsampled_fastqs")
    if os.path.exists(subsampled_fastqs):
        shutil.rmtree(subsampled_fastqs)
    os.makedirs(subsampled_fastqs)
    best_adap = test_axe_adaptors(subsampled_fastqs) # This will also generate the subsampled fastqs
    path_subsampled_fastqs = [os.path.join(subsampled_fastqs, f) for f in os.listdir(subsampled_fastqs) if f.endswith('fastq')]
    adapter_output_path = os.path.join("testfq", "temp", "axe_adaptors")
    if os.path.exists(adapter_output_path):
        shutil.rmtree(adapter_output_path)
    os.makedirs(adapter_output_path)
    if best_adap:
        threads = min(multiprocessing.cpu_count(),16)
        adapter_output_filenames = axe_adaptors_paired_end(path_subsampled_fastqs, adapter_output_path, best_adap, threads, shell=False)
        if detect_paired_end(path_subsampled_fastqs):
            test_flash(adapter_output_filenames)
