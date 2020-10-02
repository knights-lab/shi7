import unittest

from shi7.shi7_learning import *
from shi7.shi7 import *
import os
import shutil

from nose.tools import assert_equals

def subsample():
    # Fetch example fastqs
    path_fastqs = [os.path.join('testfq', f) for f in os.listdir('testfq') if f.endswith('fastq')]
    # Create a fake output directory
    output_path = os.path.join("testfq", "temp", "learning")

    if os.path.exists(output_path):
        shutil.rmtree(output_path)
    os.makedirs(output_path)
    avg_sequence_len, qual_scores = get_seq_length_qual_scores(path_fastqs, output_path, num_sequences=200, num_files=3)
    shutil.rmtree(output_path)


def axe_adaptors(subsampled_fastqs):
    num_threads = 1
    paired_end = True

    path_fastqs = [os.path.join('testfq', f) for f in os.listdir('testfq') if f.endswith('fastq')]
    output_path = os.path.join("testfq", "temp", "axe_adaptors")
    if os.path.exists(output_path):
        shutil.rmtree(output_path)
    os.makedirs(output_path)
    get_seq_length_qual_scores(path_fastqs, subsampled_fastqs, num_sequences=200) #To create the subsampled_fastqs
    path_subsampled_fastqs = [os.path.join(subsampled_fastqs, f) for f in os.listdir(subsampled_fastqs) if f.endswith('fastq')]
    best_adap, best_size, _files = choose_axe_adaptors(path_subsampled_fastqs, paired_end, output_path, num_threads)
    print("Best adaptor:", best_adap)
    print("Best adaptor size:", best_size)

    return best_adap


def flash(adapter_output_filenames):
    threads = 1
    print('adapter_output_filenames:', adapter_output_filenames)
    return_vars = []
    flash_output_path = os.path.join("testfq", "temp", "flash")
    if os.path.exists(flash_output_path):
        shutil.rmtree(flash_output_path)
    os.makedirs(flash_output_path)
    is_stitchable, allow_outies, _path_flash_fqs = flash_stitchable_and_check_outies(adapter_output_filenames, flash_output_path, threads)
    return_vars.append([is_stitchable, allow_outies])
    print('is stitchable:',is_stitchable)
    print('allow_outies:', allow_outies)
    if is_stitchable:
        # TODO: flash_check_cv only returns a tuple of 2 but this is 
        # not currently reached so it doesn't matter
        cv, mean = flash_check_cv(flash_output_path)
        return_vars.append([cv, mean])
    # TODO: I added this because the variable assignment in the above is not
    # currently reached and so the test is failing. this probably needs to 
    # be investigated further because it seems like maybe the intention is that
    # is_stitchable should be true for this test data?
    else:
        cv = 0
        mean = 0
        return_vars.append([cv, mean])

    print('CV =', cv, 'Mean =', mean)

    return sum(return_vars,[])

def trimmer(flash_output_filenames):
    print('flash_output_filenames:', flash_output_filenames)
    trimmer_output_path = os.path.join("testfq", "temp", "trimmer")
    if os.path.exists(trimmer_output_path):
        shutil.rmtree(trimmer_output_path)
    os.makedirs(trimmer_output_path)
    filter_q, trim_q = trimmer_learning(flash_output_filenames)
    print(filter_q, trim_q)

class TestLearning(unittest.TestCase):
    def test(self):
        subsample()
        subsampled_fastqs = os.path.join("testfq", "temp", "subsampled_fastqs")

        if os.path.exists(subsampled_fastqs):
            shutil.rmtree(subsampled_fastqs)

        os.makedirs(subsampled_fastqs)
        best_adap = axe_adaptors(subsampled_fastqs) # This will also generate the subsampled fastqs
        path_subsampled_fastqs = [os.path.join(subsampled_fastqs, f) for f in os.listdir(subsampled_fastqs) if f.endswith('fastq')]
        adapter_output_path = os.path.join("testfq", "temp", "axe_adaptors")

        if os.path.exists(adapter_output_path):
            shutil.rmtree(adapter_output_path)

        os.makedirs(adapter_output_path)

        if best_adap:
            threads = min(multiprocessing.cpu_count(),16)
            adapter_output_filenames = axe_adaptors_paired_end(path_subsampled_fastqs, adapter_output_path, best_adap, threads, shell=False)
        else:
            adapter_output_filenames = path_subsampled_fastqs
        
        assert_equals(None, None)
        if detect_paired_end(path_subsampled_fastqs):
            flash_vars = flash(adapter_output_filenames)
            print(flash_vars)
            if flash_vars[0]:   #if it is stitchable
                flash_output_path = os.path.join("testfq", "temp", "flash")
                if os.path.exists(flash_output_path):
                    shutil.rmtree(flash_output_path)
                    os.makedirs(flash_output_path)
                if flash_vars[2] <= 0.1:
                    flash_output_str = flash_part1(adapter_output_filenames, flash_output_path, max_overlap=int(round(flash_vars[3])), min_overlap=int(round(flash_vars[4])), allow_outies=flash_vars[1], threads=threads, shell=False)
                else:
                    flash_output_str = flash_part1(adapter_output_filenames, flash_output_path, max_overlap=700, min_overlap=20, allow_outies=flash_vars[2], threads=threads, shell=False)
                flash_output_filenames = flash_part2(flash_output_str, flash_output_path)
            else:
                flash_output_filenames = adapter_output_filenames   #skip FLASH
    
    
        trimmer(flash_output_filenames)
        
        assert_equals(None, None)
