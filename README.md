# Installation

These installation instructions are streamlined for Linux. The tool SHI7EN is installable on OSX/Windows with a few minor tweaks to this tutorial. This package requires anaconda, which is a system agnostic package and virtual environment manager. Follow the installation instructions for your system at <http://conda.pydata.org/miniconda.html>.

Once anaconda is installed, create a new virtual environment with python3.

```
conda create -n shi7en python=3
```

Now activate the environment.

```
# OSX, Linux
source activate shi7en
```

With the shogun environment activated, install the developmental SHI7EN toolchain.

```
# Put condas and pip setuptools in sync
pip install -I --upgrade setuptools

# If you want to use flash
conda install -c bioconda flash

# If you want to use trimmomatic
conda install -c bioconda trimmomatic

# Install shi7en
pip install git+https://github.com/knights-lab/shi7en --upgrade --no-cache-dir
```

With the flags provided to pip, copying and pasting any of these commands will redo the installation if a failure happened.

The final step of the procedure is to add the binary shi7en_trimmer to your path. That binary is available on the release page.

Now that everything is installed, the command 'shi7en' will be on your path when the conda environment is active. Here is the helpfile for the command:

```
$ shi7en --help
shi7en --help
usage: shi7en -i <input> -o <output> -t_trim <threads>...

This is the commandline interface for shi7en

optional arguments:
  -h, --help            show this help message and exit
  --gotta_split {True,False}
                        Split one giant fastq (well, one pair -- an R1 and R2)
                        into samples
  --debug               Enable debug (default: Disabled)
  --adaptor {None,Nextera,TruSeq3,TruSeq2,TruSeq3-2}
                        Set the type of the adaptor (default: None)
  -SE                   Run in Single End mode (default: Disabled)
  --flash {True,False}  Enable (True) or Disable (False) FLASH stiching
                        (default: True)
  --trim {True,False}   Enable (True) or Disable (False) the TRIMMER (default:
                        True)
  --allow_outies {True,False}
                        Enable (True) or Disable (False) the "outie"
                        orientation (default: True)
  --convert_fasta {True,False}
                        Enable (True) or Disable (False) the conversion of
                        FASTQS to FASTA (default: True)
  --combine_fasta {True,False}
                        Enable (True) or Disable (False) the FASTA append mode
                        (default: True)
  --shell               Use shell in Python system calls, NOT RECOMMENDED
                        (default: Disabled)
  -i INPUT, --input INPUT
                        Set the directory path of the fastq directory
  -o OUTPUT, --output OUTPUT
                        Set the directory path of the output (default: cwd)
  -t THREADS, --threads THREADS
                        Set the number of threads (default: 4)
  -m MIN_OVERLAP, --min_overlap MIN_OVERLAP
                        Set the minimum overlap length between two reads. If
                        V4 set to 285 (default: 20)
  -M MAX_OVERLAP, --max_overlap MAX_OVERLAP
                        Set the maximum overlap length between two reads. If
                        V4 set to 300 (default: 700)
  -trim_l TRIM_LENGTH, --trim_length TRIM_LENGTH
                        Set the trim length (default: 150)
  -trim_q TRIM_QUAL, --trim_qual TRIM_QUAL
                        Set the trim qual (default: 20)

```
