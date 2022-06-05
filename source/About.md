## About
```SpaGFT``` is a python package to process spatial transcriptomics data by graph Fourier transform. 

Given the raw count matrix of spatial transcriptomics data and corresponding spatial coordinate information, _SpaGFT_ transforms gene expression information from the spatial domain
to frequency domain (also called spectral domain), i.e., finding a new representation
method of genes which combined gene expression and spatial information of spots
simultaneously. _SpaGFT_ can transform complex gene expression patterns into one-dimensional, simple but informative signals (frequency signals). In this way, several tasks can be achieved by analyzing such the frequency signals, include identifying spatially variable genes (SVG), detecting tissue modules (TM), and enhancing gene expression.

![](./images/SpaGFT_workflow.svg)

### SpaGFT's functions
- Transform gene expression to frequency signals.

- Identify SVG.

- Characterize TM.

- Enhance gene expression.

### SpaGFT's advantages
- Hypothesis-free. _SpaGFT_ is a hypothesis-free graph Fourier transform framework (GFT) for SVG and TM identification from spatial transcriptomics data without assuming any spatial distribution patterns. 

- Computational effectiveness and efficiency. _SpaGFT_ can run on PC directly and still keep high accuracy and efficiency (time complexity is O(n^2)). For Visium data, it will take ~20 seconds to identify spatially variable genes. 

- Characterizing TM from a gene-centric perspective. By selecting and grouping the low-frequency FM signals, SVGs with similar spatial patterns will be grouped into one cluster and form TM. Due to TM-associated SVGs, the gene signatures and biological functions of characterized TM can be further studied.
