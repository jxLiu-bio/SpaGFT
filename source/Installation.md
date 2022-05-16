## Installation

### Create a vitual environment 

The virtual environment is recommended before installing ```SpaGFT```. You can
achieve this step easily by ```annconda```. [https://www.anaconda.com/]

If you do not have conda please install Miniconda first:

```
cd /path/to/software
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Create a separated virtual envrionment
```shell
conda create -n spagft_env python==3.8.0
conda activate spagft_env
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

For using tutorial easily, you can install jupyter and add it this environment by
```
pip install jupyter
python -m ipykernel install --user --name=spagft_env --display-name=spagft_env
```

