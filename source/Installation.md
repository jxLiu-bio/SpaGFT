# Installation

## Create a vitual environment 

The virtual environment is recommended before installing ```SpaGFT```. You can
achieve this step easily by ```annconda```. [https://www.anaconda.com/]
```shell
conda create -n spagft python==3.8
conda activate spagft
```
If you want to quit this virtual environment, just run ``` conda deacticate ```

## Install dependency packages
Before installing ```SpaGFT``` formally, the following dependency packages are 
required:

``` 
numpy==1.20.3 
pandas==1.3.4
scipy==1.7.1
networkx==2.6.3
scanpy==1.8.2
statsmodels==0.13.1 
scikit-learn==1.0.1
```
You can install all dependencies by:
```shell
git clone https://github.com/jxLiu-bio/SpaGFT
cd SpaGFT/SpaGFT_packages
pip install -r requirments.txt
```
## Install SpaGFT from GitHub
```shell
python install setup.py
```

