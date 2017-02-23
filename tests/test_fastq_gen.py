from shi7.learning import read_fastq
import hashlib
import os


def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def test_fastq_gen():
    path_fastqs = [os.path.join('testfq', f) for f in os.listdir('testfq') if f.endswith('fastq')]
    i = 0
    j = 0
    for path_fastq in path_fastqs:
        with open(path_fastq) as inf:
            fastq_gen = read_fastq(inf)
            for title, data, qualities in fastq_gen:
                i += 1
                print(title)
        with open(path_fastq) as inf:
            for line in inf:
                j += 1
    print(i)
    print(j/4)

test_fastq_gen()
