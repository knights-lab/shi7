![](https://github.com/knights-lab/shi7/workflows/build/badge.svg)

# New!
Now the *shi7_trimmer binary* in bin/ is capable of standalone paired-end QC on a per-sample basis without any dependencies, prerequisites, or installation. 

It has the following features:
1. Native gzip support built in; just feed it fastq.gz file(s)
2. Adaptor autodetection and trimming
3. New automatic adapter support including poly-G repeats on newer Illumina 2-color instruments
4. Detection and removal of any novel adapters or suspected technical sequence via PALINCUT mode (PE only)
5. Support for stripping/anonymizing all headers (sanitizing)
6. Support for outputting a set maximum depth (rarefaction) of QC'd reads
7. Support for fasta as well as fastq output formats directly
8. Support for discarding reads upon encountering an Illumina machine discard signal (an ultra-low QC code embedded in the read)
9. A variety of base-pair quality trimming modes and windows, including FLOOR and blind cut as well as different rolling average varieties
10. Support for different trimming quality thresholds at front and end of reads
11. (Almost) all other features in python-wrapped shi7
12. Over an order of magnitude faster in many cases, with no intermediates!

# Prerequisites (classic shi7 and shi7_learning pipeline)
1. Python 2.7+
2. Java

# Installation

## The CONDA way (personal install), recommended for Linux or Windows WSL
1. Follow steps 1 and 2 of https://bioconda.github.io/ (including installing MiniConda 3.6 if you don't have miniconda)
2. Do this in a terminal:
```
conda create -n shi7 -c knights-lab shi7
source activate shi7
```

## Alternative portable/server install:
Grab the latest [portable release](https://github.com/knights-lab/shi7/releases/tag/v1.0.1) package (not source), extract, then add to PATH. You should be able to execute shi7.py on the commandline. 

Confused or new to UNIX? Or something not working right? Try following the super-specific directions below:

### Here are some super specific installation instructions for the portable install:
For our purposes, we will install and use SHI7 on an interactive shell on our supercomputer, MSI, like so: `isub -n nodes=1:ppn=16 -m 22GB -w 12:00:00` (skip this if installing on your own computer, just open a (bash) terminal and optionally enter a directory where you want to install SHI7 with `cd`). You'll want to update the URLs below with the latest version on the release page above!

0. Have python (you will most likely have this already!) and importantly (at least on MacOS) the Java SDK:  http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html
(If you are installing java for the first time, you MUST REBOOT to have it be recognized).

1. Download and unpack the latest release: 
(if on MacOS, replace the word "linux" with "mac" near the end of the wget link!)
 ```
wget https://github.com/knights-lab/shi7/releases/download/v0.9.9/shi7_0.9.9_linux_release.zip
unzip shi7_*_release.zip
chmod +x shi7_*_linux_release/*
 ```
2. Add SHI7 binaries to your PATH so they can be called on the commandline anywhere:
On Linux:
```
echo "PATH=$PWD/shi7_0.9.9_linux_release:$PATH" >> ~/.bashrc
```
On Mac:
```
echo "PATH=$PWD/shi7_0.9.9_mac_release:$PATH" >> ~/.bash_profile
```
3. Reload your terminal environment and test shi7.py:
```
. ~/.bash_profile
. ~/.bashrc
shi7.py -h
```
At this point you should see the help screen printed out and SHI7 should be installed.

# Using SHI7 (simplest method):
*Note: if you installed via the CONDA method, omit the ".py" extension at the end of shi7 and shi7_learning*

1. To run interactively on a supercomputer like MSI (skip this step otherwise): 
```
isub -n nodes=1:ppn=16 -m 22GB -w 12:00:00`
module load python
```
2. Learn the appropriate shi7 parameters from the raw (unzipped) fastq data (replace 'myFastqFolder' with your actual data directory):
```
shi7_learning.py -i myFastqFolder -o learnt
```
3. Run shi7.py with the output of 2b:
`chmod +x learnt/shi7_cmd.sh && ./learnt/shi7_cmd.sh`
(or just copy the command that shi7_learning prints out when it finishes and run that)

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
Then, using `strip_underscore True` will return processed reads with just the sampleID, simplifying downstream processing. For example, an efficient command for non-stitching shotgun sequences with the sample names in the filenames before the first underscore character:

`shi7.py -i MyFastQFolder -o MyOutputFolder --adaptor Nextera --flash False --strip_delim _,1 --drop_r2 True`

# Cite

To cite SHI7:
Please cite the article published here: http://msystems.asm.org/content/3/3/e00202-17 

Al-Ghalith, Gabriel A., Benjamin Hillmann, Kaiwei Ang, Robin Shields-Cutler, and Dan Knights. 2018. “SHI7 Is a Self-Learning Pipeline for Multipurpose Short-Read DNA Quality Control.” mSystems 3 (3). American Society for Microbiology Journals: e00202–17. 

