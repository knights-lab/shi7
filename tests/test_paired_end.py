from shi7.shi7_learning import detect_paired_end
import os


def test_paired_end():
    path_fastqs = [os.path.join('testfq', f) for f in os.listdir('testfq') if f.endswith('fastq')]
    is_paired_end = detect_paired_end(path_fastqs)
    if is_paired_end:
        print('Fastq files are paired end.')
    else:
        print('Fastq files are not paired end.')

test_paired_end()
