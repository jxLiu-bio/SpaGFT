## System Requirements

### Hardware Requirements
```SpaGFT``` is friendly to hardware. All functions in _SpaGFT_ need the minimum
requirements of a CPU with 4 cores and 4G RAM. For large datasets, a large RAM is
required to avoid memory overflow.

### OS requirements
_SpaGFT_ can run on Windows, Linux, and Mac OS. We suggested analyzing users' data on ***Linux: Ubuntu 20.04*** to reproduce results from our paper.

### Python Dependencies
```SpaGFT``` requries python version >= 3.7.

```python
kneed==0.7.0
louvain==0.7.1
leidenalg==0.8.10
matplotlib==3.5.2
networkx==2.8
numba==0.55.1
numpy==1.21.5
pandas==1.4.2
plotnine==0.8.0
scanpy==1.9.1
scikit-learn==1.0.2
scipy==1.8.0
gseapy==0.10.8
igraph==0.9.10
# Note that in the current requirements there are two incorrectly
# configured dependencies. So we need to manually specify them.
chardet==5.1.0
charset-normalizer==3.1.0
```
