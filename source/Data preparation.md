## Data preparation

### Spatial transcriptomics data
_SpaGFT_ adopts _scanpy_ input methods, and the following two files are needed:
- Raw count matrix, which indicates gene expression information on all spots/pixels. 
It might be ```.h5ad```,  ```.h5``` , ```.csv```, which can be loaded by 
_scanpy_  normally for creating the AnnData object. 
- Tissue_positions_list file: A ```.csv``` or  ```.txt``` file contains coordinate
information of spots. More details could be found [here (_scanpy_)](https://scanpy.readthedocs.io/en/stable/api.html#reading).

Optionally, stained image information could be used to visualize as the background of SVGs or TMs. The method of adding them into 
AnnData object can be found in [scanpy](https://scanpy.readthedocs.io/en/stable/generated/scanpy.read_visium.html) or  [stLearn](https://stlearn.readthedocs.io/en/latest/stlearn.add.image.html).

We recommend load Visium data by:
```
import scanpy as sc
adata = sc.read_visium(path_to_visium_dataset)
```
For all spatial transcriptomics datasets, it should be pointed out that raw count matrix needs to be found at _adata.X_ and the spatial coordinate information needs to be found at _adata.obs_ or _adata.obsm_


