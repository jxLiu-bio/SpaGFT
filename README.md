# SpaGFT : Graph fourier transformer for representation, analysis, and interpretation of spatially variable genes
<img src="https://img.shields.io/badge/SpaGFT-v1.0.0-green"><img src="https://img.shields.io/badge/Platform-Linux-green"> <img src="https://img.shields.io/badge/Platform-Windows-green"> <img src="https://img.shields.io/badge/Language-python3-green"> <img src="https://img.shields.io/badge/License-MIT-green">

Given gene expression matrix that consists $n$ spots as well as their spatial coordinates and $n$ genes, ```SpaGFT``` can detect spatially variable genes (SVG) and identify tissue modules which determined by a group of SVGs with similar spatial patterns.
## System requirments
### Hardware Requirements
```SpaGFT``` is friendly to hardware. All functions in _SpaGFT_ need the minimum
requriments of a CPU with 4 cores and 4G RAM. For large datasets, a large RAM is
required to avoid memory overflow.

### OS requirements
_SpaGFT_ can run on Windows, Linux, Mac os. The package has been tested on 
following systems:

- Linux: Ubuntu 20.04
- Windows: Windows 10, Windows 11

### Python Dependencies
```SpaGFT``` requries python version >= 3.7.

``` 
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
```
## Installation Guide
### Create a vitual environment 

The virtual environment is recommended before installing ```SpaGFT```. You can
achieve this step easily by ```annconda```. [https://www.anaconda.com/]

If you do not have conda please install Miniconda first:
cd /path/to/software
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

```shell
conda create -n spagft python==3.8
conda activate spagft
```
If you want to quit this virtual environment, just run ``` conda deacticate ```

### Install dependency packages
Before intalling ```SpaGFT``` formally, the dependency packages should be installed.

You can install all dependencies by:
```shell
git clone https://github.com/jxLiu-bio/SpaGFT
cd SpaGFT
pip install -r requirements.txt
```
### Install SpaGFT from GitHub
```shell
python setup.py install
```

## Usage and Tutorials
The tutorial of using ```SpaGFT``` could be found here https://spagft.readthedocs.io/en/latest/.

Please report bugs via https://github.com/jxLiu-bio/SpaGFT/issues and ask any usage questions in https://github.com/jxLiu-bio/SpaGFT/discussions.

