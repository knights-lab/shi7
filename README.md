# Prerequisites
1. Python 2.7+
2. Java

# Installation

## The CONDA way (personal install)
1. Follow steps 1 and 2 of https://bioconda.github.io/ (including installing MiniConda 3.6 if you don't have miniconda)
2. Follow the install instructions here: https://anaconda.org/knights-lab/shi7


## Portable/server install:
Grab the latest [release](https://github.com/knights-lab/shi7/releases) package (not source), extract, then add to PATH. You should be able to execute shi7.py on the commandline. 

How to add to PATH? Well, like this:

```
echo 'PATH=$PATH: <path_to_binary>' >> ~/.bashrc
. ~/.bashrc
```

### Here are some more specific installation instructions that we use on our local supercomputer:
Installation on MSI (only do this once!):

1. Log into MSI (ssh user@login.msi.umn.edu), then gain an interactive shell: 
`isub -n nodes=1:ppn=16 -m 22GB -w 12:00:00`
(Or replace this with logging into your own machine...)
1. (optional; advanced). Change directory into where you'd like to install shi7 like "cd ~/bin"
1. Download and unpack the latest release:
```
wget https://github.com/knights-lab/shi7/releases/download/v0.92/shi7_0.92a_linux_release.zip
unzip shi7_0.92a_linux_release.zip
chmod +x shi7_0.92_linux_release/*
```
1. Add SHI7 binaries to your PATH so they can be called on the commandline anywhere:
`echo "PATH=$PWD/shi7_0.92_linux_release:$PATH" >> ~/.bashrc`
1. Reload your terminal environment and test shi7.py:
```
. ~/.bashrc
shi7.py -h
```
At this point you should see the help screen printed out.

# Now, to use SHI7:
1. Get on a system with about 1GB/ram per core you want to run SHI7 with:
(our machine): 
```
isub -n nodes=1:ppn=16 -m 22GB -w 12:00:00`
module load python
```
2. Learn the appropriate shi7 parameters from the data:
`shi7_learning.py -i myFastqFolder -o learnt`
3. Run shi7.py with the output of 2b (if run; otherwise fresh: `shi7.py -i myFastqFolder -o myOutput ...other commands`)

## Other usage examples:

Assuming you have a bunch of fastq files, of forward and reverse reads, split up by sample, that have Nextera adaptors: 

`shi7.py -i MyFastQFolder -o MyOutputFolder --adaptor Nextera`

Assuming you only have R1 reads (no paired end):

`shi7.py -i MyFastQFolder -o MyOutputFolder --adaptor Nextera -SE`

This sets the minimum read length to 285 and the maximum to 300 when stitching, which is the canonical HMP V4 16S primer coverage region. This can be a powerful QC step in and of itself. Note: if using the [EMP V4 protocol](http://press.igsb.anl.gov/earthmicrobiome/protocols-and-standards/16s/), omit these arguments.

If you have shotgun sequences, you might want to try not stitching (we recommend trying first and seeing how many stitch -- "percent combined" in the `shi7.log` file):

`--flash False`

Including `--drop_r2 True` in this case returns only R1 reads.

We recommend the following format for sequence file names:
```
sampleID_other_information_R1.fastq
sampleID_other_information_R2.fastq
```
Then, using `strip_underscore True` will return processed reads with just the sampleID, simplifying downstream processing. For example, an efficient command for non-stitching shotgun sequences:

`shi7.py -i MyFastQFolder -o MyOutputFolder --adaptor Nextera --flash False --strip_underscore True --drop_r2 True`

# Cite

To cite SHI7:
`Al-Ghalith GA, Ang K, Hillmann B, Shields-Cutler R, Knights D. (2017). SHI7: A Streamlined short-read iterative trimming pipeline. DOI:10.5281/zenodo.808832`  [![DOI](https://zenodo.org/badge/66102758.svg)](https://zenodo.org/badge/latestdoi/66102758)



# Installation (old way)

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
# If you want to use flash
conda install -c bioconda flash

# If you want to use trimmomatic
conda install -c bioconda trimmomatic

# Install shi7en
pip install git+https://github.com/knights-lab/shi7en --upgrade --no-cache-dir
```

With the flags provided to pip, copying and pasting any of these commands will redo the installation if a failure happened.

The final step of the procedure is to add the binary shi7en_trimmer to your path. The binary is available on the [release page](https://github.com/knights-lab/shi7en/releases). It is either ```ninja_shi7_linux``` or ```ninja_shi7_mac```, depending on your machine. Please rename it to ```shi7en_trimmer```. The tutorial for adding the binary to your path is shown as following:

```
echo 'PATH=$PATH: <path_to_binary>' >> ~/.bashrc
. ~/.bashrc
```

### Example
If your binary is in your ```/home/username/Downloads/shi7en_trimmer```, you can add the binary to your path in this way:
```
echo 'PATH=$PATH:/home/username/Downloads/shi7en_trimmer' >> ~/.bashrc
.~/.bashrc
```

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
