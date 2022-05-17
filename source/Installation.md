## Installation

### Configure the virtual environment 

The virtual environment is recommended before installing ```SpaGFT```. User can
configure the virtual environment by [annconda](https://www.anaconda.com/).

#### 1. Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
If the user does not install conda, please install a free minimal installer for conda (named Miniconda):

```
cd /path/to/software
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

#### 2. Create the virtual envrionment
```shell
conda create -n spagft_env python==3.8.0
conda activate spagft_env
```
If you want to quit this virtual environment, just run ``` conda deactivate ```

#### 3. Install dependency packages
Before installing ```SpaGFT``` formally, the dependency packages should be installed.

You can install all dependencies by:
```shell
git clone https://github.com/OSU-BMBL/SpaGFT
cd SpaGFT
pip install -r requirements.txt
```
#### 4. Install SpaGFT from GitHub 
```shell
python setup.py install
```

Note: We recommend [jupyter](https://jupyter.org/) for interactive usage. It can be installed and configured by
```
pip install jupyter
python -m ipykernel install --user --name=spagft_env --display-name=spagft_env
```

