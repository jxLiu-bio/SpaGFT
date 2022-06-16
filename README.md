<p align="center">

   <img src="https://github.com/jxLiu-bio/SpaGFT/blob/master/source/images/SpaGFT_logo.svg?raw=True" width="200">

</p>

# SpaGFT : Graph fourier transformer for representation, analysis, and interpretation of spatially variable genes

<img src="https://img.shields.io/badge/Platform-Linux-green"> <img src="https://img.shields.io/badge/Language-python3-green"> <img src="https://img.shields.io/badge/License-MIT-green"><img src="https://img.shields.io/badge/notebooks-passing-green"><img src="https://img.shields.io/badge/docs-passing-green">

Given a gene expression matrix that consists of $n$ spots as well as their spatial coordinates and $n$ genes, ```SpaGFT``` can detect spatially variable genes (SVG) and identify tissue modules that are determined by a group of SVGs with similar spatial patterns.

<p align="center">

   <img src="https://github.com/jxLiu-bio/SpaGFT/blob/master/source/images/SpaGFT_workflow.svg?raw=True" width="1000">

</p>

## System Requirments

### Hardware Requirements

```SpaGFT``` is friendly to hardware. All functions in _SpaGFT_ need the minimum
requirements of a CPU with 4 cores and 4G RAM. For large datasets, a large RAM is
required to avoid memory overflow.

### OS requirements

_SpaGFT_ can run on Windows, Linux, Mac os. The package has been tested on 
the following systems:

- Linux: Ubuntu 20.04
- Windows: Windows 10

### Python Dependencies

```SpaGFT``` requires python version >= 3.7.

```{txt}

kneed==0.7.0
louvain==0.7.1
matplotlib==3.5.2
networkx==2.8
numba==0.55.1
numpy =1.21.5
pandas==1.4.2
plotnine==0.8.0
scanpy==1.9.1
scikit-learn==1.0.2
scipy==1.8.0
gseapy==0.10.8
```

## Installation Guide

### Create a virtual environment

The virtual environment is recommended before installing ```SpaGFT```. Users can
install ```anaconda``` by following this tutorial. [https://www.anaconda.com/]

If users do not have conda please install Miniconda first:

```bash
cd /path/to/software
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Create a separated virtual environment

```shell
conda create -n spagft_env python==3.8.0
conda activate spagft_env
```

If users want to quit this virtual environment, just run ``` conda deactivate ```

## Install ```SpaGFT```

### Approach 1: install ```SpaGFT``` by `pip`

Users can install ```SpaGFT``` easily in this way by

```bash
pip install SpaGFT
```

Users can also install ```SpaGFT``` via Github if there exists any problems.

### Approach 2: install ```SpaGFT``` via Github

Before installing ```SpaGFT``` formally, the dependency packages should be installed.

Users can install all dependencies by:

```bash
git clone https://github.com/OSU-BMBL/SpaGFT
cd SpaGFT
pip install -r requirements.txt
```

Next, run

```bash
python setup.py install
```

Note that we recommend [jupyter](https://jupyter.org/) for interactive usage. It can be installed and configured by

```bash
conda install jupyter
python -m ipykernel install --user --name=spagft_env --display-name=spagft_env
```

## Usage and Tutorials

The tutorial of ```SpaGFT``` could be found [here](https://spagft.readthedocs.io/en/latest/).
