## About
```SpaGFT``` is a python package to process spatial transcriptomics data by
graph Fourier transform. 

Given the count matrix of spatial transcriptomics data and corresponding spatial
coordinate information, _SpaGFT_ transfer gene expression from the spatial domain
to frequency domain (also called spectral main), i.e., finding a new representation
method for genes which combined gene expression and spatial information 
simultaneously. _SpaGFT_ can transform complex gene expression patterns into simple but informative signals (frequency signals), identify spatially variable genes (SVG), detect tissue modules (TM), and enhance gene expression.

![](./images/SpaGFT_workflow.svg)

### SpaGFT's analyses objects
- Transform gene expression to frequency signals.

- Identify SVG.

- Characterize TM.

- Enhance gene expression (optional).

### SpaGFT's advantages
- Hypothesis-free. _SpaGFT_ is a hypothesis-free graph Fourier transform framework (GFT) for SVG and TM identification from spatial transcriptomics data without assuming any spatial patterns. 

- Computational effectiveness and efficiency. _SpaGFT_ can run on PC and still keep high accuracy and efficiency (time complexity is O(n^2)). For Visium data, it will take ~20 seconds to identify spatially variable genes. 

- Characterizing TM from a gene-centric perspective. By selecting and grouping the low-frequency FM signals, SVGs with similar spatial patterns will be grouped into one cluster and form TM. Due to TM-associated SVGs, the gene signatures and biological functions of characterized TM can be further studied.
