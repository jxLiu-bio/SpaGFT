## Installation


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

### Install ```SpaGFT```

#### Approach 1: Install from PyPI
```bash
pip install SpaGFT
```

#### Approach 2: Install from source

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

### Install jupyter (optional)
Note that we recommend [jupyter](https://jupyter.org/) for interactive usage. It can be installed and configured by

```bash
conda install jupyter
python -m ipykernel install --user --name=spagft_env --display-name=spagft_env
```

