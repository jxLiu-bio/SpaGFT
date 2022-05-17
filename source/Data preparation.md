## Data preparation

### Spatial transcriptomics data
_SpaGFT_ adopts _scanpy_ input methods, and the following two files are needed:
- Raw count matrix, which indicates gene expression information on all spots/pixels. 
It might be ```.h5ad```,  ```.h5``` , ```.csv```, which could be load by 
scanpy for creating the AnnData object. 
- Tissue_positions_list file: A ```.csv``` or  ```.txt``` file contains coordinate
information of spots. More details could be found [here (_scanpy_)](https://scanpy.readthedocs.io/en/stable/api.html#reading).

Optionally, stained image information could be used to visualize as the background of SVGs or TMs. The method of adding them into 
AnnData object can be found in [scanpy](https://scanpy.readthedocs.io/en/stable/generated/scanpy.read_visium.html) or  [stLearn](https://stlearn.readthedocs.io/en/latest/stlearn.add.image.html).


