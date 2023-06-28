## Installation

### Configure the virtual environment 

The virtual environment is recommended before installing ```SpaGFT```. Users can
configure the virtual environment by [anaconda](https://www.anaconda.com/).

#### 1. Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html)[optional]
If the user does not install conda, please install a free minimal installer for conda (named Miniconda):

```shell
cd /path/to/software
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

#### 2. Create the virtual environment
```shell
conda create -n spagft_env python==3.8.0
conda activate spagft_env
```
If you want to quit this virtual environment, just run ``` conda deactivate ```

#### 3. Install ```SpaGFT```
##### 3.1 Approach 1: install ```SpaGFT``` by Pipy.
You can install ```SpaGFT``` easily in this way by
```shell
pip install SpaGFT
```
You can also install ```SpaGFT``` via Github if there exists any problems.
##### 3.2 Approach 2: install ```SpaGFT``` via Github
Before installing ```SpaGFT``` formally, the dependency packages should be installed.

You can install all dependencies by:
```shell
git clone https://github.com/OSU-BMBL/SpaGFT
cd SpaGFT
pip install -r requirements.txt
```
Next, run
```shell
python setup.py install
```

We recommend [jupyter](https://jupyter.org/) for interactive usage. It can be installed and configured by
```shell
conda install jupyter
python -m ipykernel install --user --name=spagft_env --display-name=spagft_env
```

