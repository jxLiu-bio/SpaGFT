# Data preparation

### Spatial transcriptomics data
The following two files are needed at least:
- Raw count matrix which indicates gene expresion information on all spots/pixels. 
It could be ```.h5ad```,  ```.h5``` , ```.csv``` etc. once it could be load by 
scanpy correctly to used as AnnData object. 
- Tissue_positions_list file: A ```.csv``` or  ```.txt``` file contains coordinate
information of spots. More details could be found [here](https://scanpy.readthedocs.io/en/stable/api.html#reading).

Additionally, stained image information could be used to enhance the 
visualization of genes or tissue modules. The method of adding them into 
AnnData object can found in [scanpy](https://scanpy.readthedocs.io/en/stable/generated/scanpy.read_visium.html) or  [stLearn](https://stlearn.readthedocs.io/en/latest/stlearn.add.image.html).

