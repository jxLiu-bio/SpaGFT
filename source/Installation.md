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
kneed==0.7.0
louvain==0.7.1
matplotlib==3.5.2
networkx==2.8
numba==0.55.1
numpy ==1.21.5
pandas==1.4.2
plotnine==0.8.0
scanpy==1.9.1
scikit-learn==1.0.2
scipy==1.8.0
```
You can install all dependencies by:
```shell
git clone https://github.com/jxLiu-bio/SpaGFT
cd SpaGFT
pip install -r requirements.txt
```
## Install SpaGFT from GitHub
```shell
python setup.py install
```

