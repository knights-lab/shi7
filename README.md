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
