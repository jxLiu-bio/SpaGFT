## Spg script
This script can be used to run examples shown in the jypyter notebooks.

Below is some example script usage. Run from project root:



### Usage ###
```shell
# usage: SpaGFT [-h] [-env] [-v] [-data PATH] [-o PATH] [-omcs {gene,protein}] [-org {Mouse,Human}] [-flp] [-fcnt] [-fnorm] [-flog] [-cell2loc PATH] [-clustering_alg ALG] [-fident_svos] [-fident_tms] [-fcmp_tm_clustering] [-save FNAME] [-gsvosum] [-gsvos PATH] [-gtmsum] [-gtms PATH]
#
# SpaGFT is a python package to analyze spatial transcriptomics data. It uses a hypothesis-free graph Fourier transform model to identify spatially variable genes/proteins, tissue modules (TM) and functional enrichment analysis to explore the underlying biological processes in these TMs.
#
# optional arguments:
#   -h, --help            show this help message and exit
#   -env                  python environment
#   -v                    program version
#
# data:
#   -data PATH            input file/folder in .h5ad/visium format
#   -o PATH               output directory
#
# analyze:
#   -omcs {gene,protein}  specify the omics kind
#   -org {Mouse,Human}    specify the organism kind
#   -flp                  prepare dataset: low pass filter in the graph Fourier domain
#   -fcnt                 prepare dataset: kNN filter
#   -fnorm                prepare dataset: normalize
#   -flog                 prepare dataset: logarithmize
#   -cell2loc PATH        augment dataset: add cell2location data specified in .csv file
#   -clustering_alg ALG   specify the algorithm used for clustering tms
#   -fident_svos          identify spatially variable omics
#   -fident_tms           identify tissue modules
#   -fcmp_tm_clustering   compare tm clustering with ground truth specified in the dataset annotations
#   -save FNAME           save analyzed dataset to output directory with the given .h5ad filename
#
# graph:
#   -gsvosum              graph the svo summary
#   -gsvos PATH           graph selected svos specified in .csv file
#   -gtmsum               graph the tissue module summary
#   -gtms PATH            graph selected tissue modules specified in .csv file

python "./source/spatial/spg.py" -v
python "./source/spatial/spg.py" -env
```



### Example 1 - mouse_brain_coronal ###

Dataset:
- [mouse-he-coronal/](https://drive.google.com/drive/folders/1-B9ipgx5gdWU-go38-acCdS9xHCK9D3E)
- [mouse-he-coronal_cell2location_results.csv](https://drive.google.com/drive/folders/1-B9ipgx5gdWU-go38-acCdS9xHCK9D3E)

Folders:
```python
# source/
# +- spatial/
#   +- spg.py   # script
#   +- data/
#     +- mouse-he-coronal/   # (initial) input dataset
#     +- mouse-he-coronal_cell2location_results.csv   # input
#     +- spg/
#       +- mouse_brain_coronal/   # input
#         +- svo_selection.csv
#         +- tm_selection.csv
#   +- results/
#     +- spg
#       +- mouse_brain_coronal/   # output + (modified) input dataset
```

svo_selection.csv
```csv
SVO_NAME
Ahi1
Nsmf
Tbr1
Ano2
Cartpt
Cbln2
Ttr
Pmch
```

tm_selection.csv
```csv
TM_NAME,TM_INDEX,TM_SVOS,OVERLAP_TMS
tm_6,6,Ahi1 Ap1s2 Arxes2 Bex1 Ctxn2 Hap1,8 9
tm_8,8,Ahi1 Ap1s2 Arxes2 Bex1 Ctxn2 Hap1,6 9
tm_9,9,Ahi1 Ap1s2 Arxes2 Bex1 Ctxn2 Hap1,6 8
```

Investigate svos:
```shell
# identify svos and get svo summary
python "./source/spatial/spg.py" -data "./source/spatial/data/mouse-he-coronal" -o "./source/spatial/results/spg/mouse_brain_coronal" -omcs "gene" -fcnt -fnorm -flog -fident_svos -save "dataset.h5ad" -gsvosum

# graph specific svos
python "./source/spatial/spg.py" -data "./source/spatial/results/spg/mouse_brain_coronal/dataset.h5ad" -o "./source/spatial/results/spg/mouse_brain_coronal" -gsvos "./source/spatial/data/spg/mouse_brain_coronal/svo_selection.csv"
```

Investigate tms:
```shell
# identify tms and get tm summary
python "./source/spatial/spg.py" -data "./source/spatial/results/spg/mouse_brain_coronal/dataset.h5ad" -o "./source/spatial/results/spg/mouse_brain_coronal" -flp -org "Mouse" -cell2loc "./source/spatial/data/mouse-he-coronal_cell2location_results.csv" -fident_tms -save "dataset.h5ad" -gtmsum

# graph specific tms
python "./source/spatial/spg.py" -data "./source/spatial/results/spg/mouse_brain_coronal/dataset.h5ad" -o "./source/spatial/results/spg/mouse_brain_coronal" -gtms "./source/spatial/data/spg/mouse_brain_coronal/tm_selection.csv"
```



### Example 2 - human_lymph_node ###

Dataset:
- [V1_Human_Lymph_Node/](https://support.10xgenomics.com/spatial-gene-expression/datasets/1.0.0/V1_Human_Lymph_Node)
- [lymphnode_cell2location_results.csv](https://drive.google.com/drive/folders/1-B9ipgx5gdWU-go38-acCdS9xHCK9D3E)

Folders:
```python
# source/
# +- spatial/
#   +- spg.py   # script
#   +- data/
#     +- V1_Human_Lymph_Node/   # (initial) input dataset
#     +- lymphnode_cell2location_results.csv   # input
#     +- spg/
#       +- human_lymph_node/   # input
#         +- svo_selection.csv
#         +- tm_selection.csv
#   +- results/
#     +- spg
#       +- human_lymph_node/   # output + (modified) input dataset
```

svo_selection.csv
```csv
SVO_NAME
CD3E
IL7R
CCR7
PCNA
CDK1
CDC20
CD19
CD79B
```

tm_selection.csv
```csv
TM_NAME,TM_INDEX,TM_SVOS,OVERLAP_TMS
tm_3,3,CD3E IL7R CCR7,5 7
tm_5,5,CD3E IL7R CCR7,3 7
tm_7,7,CD3E IL7R CCR7,3 5
```

Investigate svos:
```shell
# identify svos and get svo summary
python "./source/spatial/spg.py" -data "./source/spatial/data/V1_Human_Lymph_Node" -o "./source/spatial/results/spg/human_lymph_node" -omcs "gene" -fcnt -fnorm -flog -fident_svos -save "dataset.h5ad" -gsvosum

# graph specific svos
python "./source/spatial/spg.py" -data "./source/spatial/results/spg/human_lymph_node/dataset.h5ad" -o "./source/spatial/results/spg/human_lymph_node" -gsvos "./source/spatial/data/spg/human_lymph_node/svo_selection.csv"
```

Investigate tms:
```shell
# identify tms and get tm summary
python "./source/spatial/spg.py" -data "./source/spatial/results/spg/human_lymph_node/dataset.h5ad" -o "./source/spatial/results/spg/human_lymph_node" -flp -org "Human" -cell2loc "./source/spatial/data/lymphnode_cell2location_results.csv" -fident_tms -save "dataset.h5ad" -gtmsum

# graph specific tms
python "./source/spatial/spg.py" -data "./source/spatial/results/spg/human_lymph_node/dataset.h5ad" -o "./source/spatial/results/spg/human_lymph_node" -gtms "./source/spatial/data/spg/human_lymph_node/tm_selection.csv"
```



### Example 3 - codex_lymph_node_A6 ###

Dataset:
- [A6_compress_codex.h5ad](https://drive.google.com/drive/folders/1-B9ipgx5gdWU-go38-acCdS9xHCK9D3E)

Folders:
```python
# source/
# +- spatial/
#   +- spg.py   # script
#   +- data/
#     +- A6_compress_codex.h5ad   # (initial) input dataset
#     +- spg/
#       +- codex_lymph_node_A6/   # input
#         +- svo_selection.csv
#         +- tm_selection.csv
#   +- results/
#     +- spg
#       +- codex_lymph_node_A6/   # output + (modified) input dataset
```

svo_selection.csv
```csv
SVO_NAME
reg001_cyc002_ch004_CD56
reg001_cyc008_ch003_CD5
reg001_cyc016_ch003_CD57
reg001_cyc019_ch004_CD11b
reg001_cyc009_ch004_PD-1
reg001_cyc023_ch004_Podoplanin
```

tm_selection.csv
```csv
TM_NAME,TM_INDEX,TM_SVOS,OVERLAP_TMS
tm_2,2,reg001_cyc019_ch004_CD11b reg001_cyc011_ch003_CD11c reg001_cyc023_ch004_Podoplanin,1
```

Investigate svos:
```shell
# identify svos and get svo summary
python "./source/spatial/spg.py" -data "./source/spatial/data/A6_compress_codex.h5ad" -o "./source/spatial/results/spg/codex_lymph_node_A6" -omcs "protein" -fcnt -fident_svos -save "dataset.h5ad" -gsvosum

# graph specific svos
python "./source/spatial/spg.py" -data "./source/spatial/results/spg/codex_lymph_node_A6/dataset.h5ad" -o "./source/spatial/results/spg/codex_lymph_node_A6" -gsvos "./source/spatial/data/spg/codex_lymph_node_A6/svo_selection.csv"
```

Investigate tms:
```shell
# identify tms and get tm summary
python "./source/spatial/spg.py" -data "./source/spatial/results/spg/codex_lymph_node_A6/dataset.h5ad" -o "./source/spatial/results/spg/codex_lymph_node_A6" -org "Human" -fident_tms -save "dataset.h5ad" -gtmsum

# graph specific tms
python "./source/spatial/spg.py" -data "./source/spatial/results/spg/codex_lymph_node_A6/dataset.h5ad" -o "./source/spatial/results/spg/codex_lymph_node_A6" -gtms "./source/spatial/data/spg/codex_lymph_node_A6/tm_selection.csv"
```



### Example 4 - e95_e1s1_mosta ###

[Dataset](https://db.cngb.org/stomics/mosta/download/):
- [E9.5_E1S1.MOSTA.h5ad](https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/E9.5_E1S1.MOSTA.h5ad)

Folders:
```python
# source/
# +- spatial/
#   +- spg.py   # script
#   +- data/
#     +- E9.5_E1S1.MOSTA.h5ad   # (initial) input dataset
#     +- E9.5_E1S1.MOSTA.h5ad:Zone.Identifier
#   +- results/
#     +- spg
#       +- e95_e1s1_mosta/   # output + (modified) input dataset
```

Investigate svos + tms:
```shell
# identify svos and get svo summary
python "./source/spatial/spg.py" -data "./source/spatial/data/E9.5_E1S1.MOSTA.h5ad" -o "./source/spatial/results/spg/e95_e1s1_mosta" -omcs "gene" -fcnt -fnorm -flog -fident_svos -save "dataset.h5ad" -gsvosum

# identify tms with louvain algorithm and get tm summary with clustering comparison
python "./source/spatial/spg.py" -data "./source/spatial/results/spg/e95_e1s1_mosta/dataset.h5ad" -o "./source/spatial/results/spg/e95_e1s1_mosta" -org "Mouse" -clustering_alg "louvain" -fident_tms -fcmp_tm_clustering -save "dataset.h5ad" -gtmsum

# identify tms with leiden algorithm and get tm summary with clustering comparison
python "./source/spatial/spg.py" -data "./source/spatial/results/spg/e95_e1s1_mosta/dataset.h5ad" -o "./source/spatial/results/spg/e95_e1s1_mosta" -org "Mouse" -clustering_alg "leiden" -fident_tms -fcmp_tm_clustering -save "dataset.h5ad" -gtmsum
```



### Example 5 - dorsal_midbrain_cell ###

[Dataset](https://db.cngb.org/stomics/mosta/download/):
- [Dorsal_midbrain_cell_bin.h5ad](https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/Dorsal_midbrain_cell_bin.h5ad)

Folders:
```python
# source/
# +- spatial/
#   +- spg.py   # script
#   +- data/
#     +- Dorsal_midbrain_cell_bin.h5ad   # (initial) input dataset
#     +- Dorsal_midbrain_cell_bin.h5ad:Zone.Identifier
#   +- results/
#     +- spg
#       +- dorsal_midbrain_cell/   # output + (modified) input dataset
```

Investigate svos + tms:
```shell
# identify svos and get svo summary
python "./source/spatial/spg.py" -data "./source/spatial/data/Dorsal_midbrain_cell_bin.h5ad" -o "./source/spatial/results/spg/dorsal_midbrain_cell" -omcs "gene" -fcnt -fnorm -flog -fident_svos -save "dataset.h5ad" -gsvosum

# identify tms with louvain algorithm and get tm summary with clustering comparison
python "./source/spatial/spg.py" -data "./source/spatial/results/spg/dorsal_midbrain_cell/dataset.h5ad" -o "./source/spatial/results/spg/dorsal_midbrain_cell" -org "Mouse" -clustering_alg "louvain" -fident_tms -fcmp_tm_clustering -save "dataset.h5ad" -gtmsum

# identify tms with leiden algorithm and get tm summary with clustering comparison
python "./source/spatial/spg.py" -data "./source/spatial/results/spg/dorsal_midbrain_cell/dataset.h5ad" -o "./source/spatial/results/spg/dorsal_midbrain_cell" -org "Mouse" -clustering_alg "leiden" -fident_tms -fcmp_tm_clustering -save "dataset.h5ad" -gtmsum
```




