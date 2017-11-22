# Prerequisites
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
Grab the latest [release](https://github.com/knights-lab/shi7/releases) package (not source), extract, then add to PATH. You should be able to execute shi7.py on the commandline. 

Confused or new to UNIX? Or something not working right? Try following the super-specific directions below:

### Here are some super specific installation instructions for the portable install:
For our purposes, we will install and use SHI7 on an interactive shell on our supercomputer, MSI, like so: `isub -n nodes=1:ppn=16 -m 22GB -w 12:00:00` (skip this if installing on your own computer, just open a (bash) terminal and optionally enter a directory where you want to install SHI7 with `cd`).

0. Have python (you will most likely have this already!) and importantly (at least on MacOS) the Java SDK:  http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html
(If you are installing java for the first time, you MUST REBOOT to have it be recognized).

1. Download and unpack the latest release: 
(if on MacOS, replace the word "linux" with "mac" near the end of the wget link!)
 ```
wget https://github.com/knights-lab/shi7/releases/download/v0.9.5/shi7_0.9.3c_linux_release.zip
unzip shi7_0.9.3_*_release.zip
chmod +x shi7_0.9.3_*_release/*
 ```
2. Add SHI7 binaries to your PATH so they can be called on the commandline anywhere:
On Linux:
```
echo "PATH=$PWD/shi7_0.9.3_linux_release:$PATH" >> ~/.bashrc
```
On Mac:
```
echo "PATH=$PWD/shi7_0.9.3_mac_release:$PATH" >> ~/.bash_profile
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
`Al-Ghalith GA, Ang K, Hillmann B, Shields-Cutler R, Knights D. (2017). SHI7: A Streamlined short-read iterative trimming pipeline. DOI:10.5281/zenodo.808832`  [![DOI](https://zenodo.org/badge/66102758.svg)](https://zenodo.org/badge/latestdoi/66102758)

