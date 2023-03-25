<p align="center">

   <img src="https://user-images.githubusercontent.com/26455910/227746120-9bf87226-7320-49c3-8108-8ed09763cdf8.png" width="700">

</p>

# Spatial omics representation and functional tissue module inference using graph Fourier transform


<img src="https://img.shields.io/badge/Platform-Linux-green"> <img src="https://img.shields.io/badge/Language-python3-green"> <img src="https://img.shields.io/badge/License-MIT-green"><img src="https://img.shields.io/badge/notebooks-passing-green"><img src="https://img.shields.io/badge/docs-passing-green">

Tissue module (TM) is a spatially organized tissue region and executes specialized biological functions, recurring and varying at different tissue sites. However, the computational identification of TMs poses challenges due to their convoluted biological functions, poorly-defined molecular features, and varying spatially organized patterns. Here, we present a hypothesis-free graph Fourier transform model, ```SpaGFT```, to represent spatially organized features using the Fourier coefficients, leading to an accurate representation of spatially variable genes and proteins and the characterization of TM at a fast computational speed. We implemented sequencing-based and imaging-based spatial transcriptomics, spatial-CITE-seq, and spatial proteomics to identify spatially variable genes and proteins, define TM identities, and infer convoluted functions among TMs in mouse brains and human lymph nodes.  We collected a human tonsil sample and performed CODEX to accurately demonstrate molecular and cellular variability within the secondary follicle structure. The superior accuracy, scalability, and interpretability of SpaGFT indicate that it is an effective representation of spatially-resolved omics data and an essential tool for bringing new insights into molecular tissue biology.

<p align="center">

   <img src="https://user-images.githubusercontent.com/26455910/226082582-cd77af6b-13b2-4f8c-9003-d96f0a28bd4d.svg" width="1000">

</p>

## System Requirments

### Hardware Requirements

```SpaGFT``` is friendly to hardware. All functions in _SpaGFT_ need the minimum
requirements of a CPU with 4 cores and 4G RAM. For large datasets, a large RAM is
required to avoid memory overflow.

### OS requirements

_SpaGFT_ can run on Windows, Linux, Mac os. The package has been tested on 
the following systems:

- Linux: Ubuntu 20.04 (recommend)
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
igraph==0.9.10
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
