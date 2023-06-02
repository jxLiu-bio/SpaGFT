import os
import platform
import unicodedata
import argparse as ar
import numpy as np
import pandas as pd
import scanpy as sc
import sklearn as sk
import matplotlib.pyplot as plt
import SpaGFT as spg
from typing import List, Union, Tuple, Any, Callable
from enum import Enum
from pathlib import Path
from anndata import AnnData
from pandas import DataFrame

spg_version = "SpaGFT 0.1.2"
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi_save=600, facecolor="white")


### ____________________________________________________________________________
### Utilities
class OmicsType(Enum):
    GENE = 1
    PROTEIN = 2


def _omics_type(type: str) -> OmicsType:
    if type == "gene":
        return OmicsType.GENE
    elif type == "protein":
        return OmicsType.PROTEIN
    else:
        raise NotImplementedError()


def _norm_part(path_part: str) -> str:
    path_part = (
        unicodedata.normalize("NFKD", path_part)
        .encode("ascii", "replace")
        .decode("ascii")
        .lower()
    )
    return path_part


def _ensure_ext(path: str, ext: str) -> str:
    if path.endswith(ext):
        return path
    return f"{path}{ext}"


def _read_visium(path: Union[str, Path]) -> AnnData:
    adata = sc.read_visium(path)
    adata.var_names_make_unique()
    adata.raw = adata
    return adata


def _read_h5ad(filepath: Union[str, Path]) -> AnnData:
    adata = sc.read_h5ad(filepath)
    adata.var_names_make_unique()
    adata.raw = adata
    return adata


def _read_csv(
    filepath: Union[str, Path], index_col: Union[int, None] = None
) -> DataFrame:
    data = pd.read_csv(filepath, index_col=index_col)
    return data


def _write_h5ad(adata: AnnData, filepath: Union[str, Path]) -> None:
    filepath = _ensure_ext(filepath, ".h5ad")
    adata.write_h5ad(filepath)


def _write_csv(data: DataFrame, filepath: Union[str, Path], index: bool = True) -> None:
    filepath = _ensure_ext(filepath, ".csv")
    data.to_csv(filepath, index=index)


def _write_json(data: DataFrame, filepath: Union[str, Path]) -> None:
    filepath = _ensure_ext(filepath, ".json")
    data.to_json(filepath)


### ____________________________________________________________________________
### SpaGFT functionalities - insights


def read_dataset(path: Union[str, Path]) -> AnnData:
    if os.path.isfile(path):
        return _read_h5ad(path)
    elif os.path.isdir(path):
        return _read_visium(path)
    else:
        raise NotImplementedError()


def read_cell2location(adata: AnnData, filepath: Union[str, Path]) -> None:
    # load deconvolution results from cell2location.
    deconv_df = _read_csv(filepath, index_col=0)
    # add deconvolution results to anndata object
    adata.obsm["cell_type_proportion"] = deconv_df


def write_dataset(adata: AnnData, filepath: Union[str, Path]) -> None:
    filepath = _ensure_ext(filepath, ".h5ad")
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    _write_h5ad(adata, filepath)


def preprocess(
    adata: AnnData,
    filter_low_pass: bool,
    filter_by_count: bool,
    normalize: bool,
    logarithmize: bool,
) -> None:
    if filter_low_pass:
        spg.low_pass_enhancement(adata, c=0.01, ratio_low_freq=15, inplace=True)
    if filter_by_count:
        sc.pp.filter_genes(adata, min_cells=10, inplace=True)
    if normalize:
        sc.pp.normalize_total(adata, inplace=True)
    if logarithmize:
        sc.pp.log1p(adata, copy=False)


def identify_svo_list(
    adata: AnnData, omics_type: OmicsType, spatial_info: Union[str, Tuple[str, str]]
) -> None:
    # determine the number of low-frequency FMs and high-frequency FMs
    low_cutoff, high_cutoff = spg.gft.determine_frequency_ratio(
        adata, ratio_neighbors=1, spatial_info=spatial_info
    )
    adata.uns["low_cutoff"] = low_cutoff
    adata.uns["high_cutoff"] = high_cutoff
    print(f"Cutoff frequencies: {low_cutoff}/{high_cutoff}")

    # calculation
    score_df = spg.detect_svg(
        adata,
        spatial_info=spatial_info,
        ratio_low_freq=low_cutoff,
        ratio_high_freq=high_cutoff,
        ratio_neighbors=1,
        filter_peaks=True,
        S=6,  # kneedle sensitivity
    )
    adata.uns["score_df"] = score_df

    # extract spatially variable omics
    svo_list = None
    if omics_type == OmicsType.GENE:
        svo_list = score_df[score_df.cutoff_gft_score][
            score_df.qvalue < 0.05
        ].index.to_list()
    elif omics_type == OmicsType.PROTEIN:
        # for non-whole-transcriptomics datasets, use qvalue as cutoff to obtain spatial variable features
        svo_list = score_df[score_df.qvalue < 0.05].index.to_list()
    else:
        raise NotImplementedError()
    print(f"The number of SVOs: {len(svo_list)}")
    adata.uns["svo_list"] = svo_list


def identify_tm_list(
    adata: AnnData,
    svo_list: List[str],
    spatial_info: Union[str, Tuple[str, str]],
    low_cutoff: Union[float, Any],
    resolution: Union[float, Tuple[float, ...]],
    ratio_neighbors: int,
    n_neighbors: int,
    tm_clustering_alg: Union[str, None],
) -> None:
    if tm_clustering_alg is None:
        tm_clustering_alg = "louvain"

    tm_svo_df, _ = spg.identify_tissue_module(
        adata,
        svg_list=svo_list,
        ratio_fms=low_cutoff,
        spatial_info=spatial_info,
        resolution=resolution,
        ratio_neighbors=ratio_neighbors,
        n_neighbors=n_neighbors,
        algorithm=tm_clustering_alg,
    )
    adata.uns["tm_svo_df"] = tm_svo_df
    adata.uns["tm_clustering_alg"] = tm_clustering_alg
    print(f"Algorithm used for svo clustering: {tm_clustering_alg}")


def compare_tm_clustering(
    adata: AnnData,
    svo_list: List[str],
) -> None:
    pass
    current_genes = adata.uns["detect_TM_data"]["gft_umap_tm"].index.tolist()
    if set(svo_list) <= set(current_genes):
        svo_list = np.intersect1d(svo_list, current_genes)

    clustering_df = pd.concat(
        (
            adata.uns["detect_TM_data"]["gft_umap_tm"].loc[svo_list, :],
            adata.var.loc[svo_list, :].tissue_module,
        ),
        axis=1,
    )

    categories = [eval(i) for i in np.unique(clustering_df.tissue_module)]
    categories = np.sort(np.array(categories))
    categories = categories.astype(str)
    clustering_df.tissue_module = pd.Categorical(
        clustering_df.tissue_module, categories=categories
    )

    tm_clustering_score = sk.metrics.adjusted_rand_score(
        adata.obs["annotation"], clustering_df
    )

    adata.uns["tm_clustering_score"] = tm_clustering_score
    print(f"Svo clustering score (compared to reference in annotations): {tm_clustering_score}")


### ____________________________________________________________________________
### SpaGFT functionalities - plotting


def write_svo_summary(
    dirpath: Union[str, Path],
    adata: AnnData,
    svo_list: List[str],
) -> None:
    os.makedirs(dirpath, exist_ok=True)

    # svo list
    svo_df = pd.DataFrame(data={"SVO_NAME": svo_list})
    _write_csv(svo_df, f"{dirpath}/svo_all.csv", index=False)

    # svo umap
    spg.plot.gene_signal_umap(
        adata,
        svg_list=svo_list,
        size=30,
        return_fig=True,
        show_fig=False,
        save_path=f"{dirpath}/svo_all.umap.png",
    )
    plt.close()


def read_svo_selection(filepath: Union[str, Path]) -> List[str]:
    svo_df = _read_csv(filepath)
    svo_list = svo_df["SVO_NAME"].values.tolist()
    return svo_list


def write_svo_selection(
    dirpath: Union[str, Path],
    adata: AnnData,
    svo_selection: List[str],
    omics_type: OmicsType,
) -> None:
    os.makedirs(dirpath, exist_ok=True)

    # svo list
    svo_df = pd.DataFrame(data={"SVO_NAME": svo_selection})
    _write_csv(svo_df, f"{dirpath}/svo_selection.csv", index=False)

    # individual svos
    svo_dirpath = f"{dirpath}/svo"
    os.makedirs(svo_dirpath, exist_ok=True)
    for svo in svo_selection:
        write_svo(svo_dirpath, adata, svo, omics_type)


def write_svo(
    dirpath: Union[str, Path],
    adata: AnnData,
    svo: str,
    omics_type: OmicsType,
) -> None:
    os.makedirs(dirpath, exist_ok=True)

    if omics_type == OmicsType.GENE:
        # svgs
        sc.pl.spatial(
            adata,
            color=svo,
            size=1.6,
            cmap="magma",
            use_raw=False,
            return_fig=True,
            show=False,
            save=False,
        )
        plt.savefig(f"{dirpath}/{_norm_part(svo)}.svo.png")
        plt.close()
    elif omics_type == OmicsType.PROTEIN:
        # svps
        spg.plot.scatter_gene(
            adata,
            gene=svo,
            spatial_info=["x", "y"],
            cmap="magma",
            return_fig=True,
            show_fig=False,
            save_path=f"{dirpath}/{_norm_part(svo)}.svo.png",
        )
        plt.close()
    else:
        raise NotImplementedError()

    # freq signals
    spg.plot.gene_freq_signal(
        adata,
        gene=svo,
        return_fig=True,
        show_fig=False,
        save_path=f"{dirpath}/{_norm_part(svo)}.freq.png",
    )
    plt.close()


def write_tm_summary(
    dirpath: Union[str, Path],
    adata: AnnData,
    svo_list: List[str],
    tm_info: Union[dict, None] = None,
) -> None:
    os.makedirs(dirpath, exist_ok=True)

    # tissue module general info
    if tm_info:
        tm_info_df = pd.DataFrame.from_dict(tm_info, orient="index")
        _write_json(tm_info_df, f"{dirpath}/tm_info.json")

    # the genes which support corresponding tissue modules
    tm_gene_df = adata.var.loc[svo_list, :]
    _write_csv(tm_gene_df, f"{dirpath}/tm_svo.csv")

    # tissue module supporting svos information
    tm_binary_df = adata.obsm["tm_binary"]
    _write_csv(tm_binary_df, f"{dirpath}/tm_supporting_svos.csv")

    # svo clustering
    spg.plot.scatter_umap_clustering(
        adata,
        svo_list,
        return_fig=True,
        show_fig=False,
        save_path=f"{dirpath}/tm_svo.umap.png",
    )
    plt.close()


def read_tm_selection(filepath: Union[str, Path]) -> DataFrame:
    tm_selection = _read_csv(filepath)
    return tm_selection


def write_tm_selection(
    dirpath: Union[str, Path],
    adata: AnnData,
    tm_selection: DataFrame,
    spatial_info: Union[str, Tuple[str, str]],
    organism: str,
    cell2loc: bool,
) -> None:
    os.makedirs(dirpath, exist_ok=True)

    # tm list
    _write_csv(tm_selection, f"{dirpath}/tm_selection.csv", index=False)

    # plot tms
    tm_name_list = tm_selection["TM_NAME"].tolist()
    tm_colors = np.array(["#9aabe1", "#e2b97c", "#5ac096"], dtype=str)
    tm_colors = np.tile(tm_colors, np.ceil(len(tm_name_list) / len(tm_colors))).tolist()

    spg.plot.scatter_tm(
        adata,
        tm=tm_name_list,
        tm_color=tm_colors,
        spatial_info=spatial_info,
        size=15,
        return_fig=True,
        show_fig=False,
        save_path=f"{dirpath}/tm_selection.png",
    )
    plt.close()

    # individual tms
    tm_dirpath = f"{dirpath}/tm"
    os.makedirs(tm_dirpath, exist_ok=True)
    for i in range(len(tm_name_list)):
        tm_name = str(tm_selection["TM_NAME"][i])
        tm_idx = str(tm_selection["TM_INDEX"][i])
        tm_svo_list = str(tm_selection["TM_SVOS"][i]).split(" ")
        tm_overlap_show = str(tm_selection["OVERLAP_TMS"][i]).split(" ")

        write_tm(
            tm_dirpath,
            adata,
            tm_name,
            tm_idx,
            tm_svo_list,
            tm_overlap_show,
            spatial_info,
            organism,
            cell2loc,
        )


def write_tm(
    dirpath: Union[str, Path],
    adata: AnnData,
    tm_name: str,
    tm_idx: str,
    tm_svo_list: List[str],
    tm_overlap_show: Union[List[str], None],
    spatial_info: Union[str, Tuple[str, str]],
    organism: str,
    cell2loc: bool,
) -> None:
    os.makedirs(dirpath, exist_ok=True)

    # plot tm and corresponding svos
    spg.plot.scatter_tm_gene(
        adata,
        tm=tm_name,
        tm_color="#9aabe1",
        gene=tm_svo_list,
        spatial_info=spatial_info,
        size=13,
        return_fig=True,
        show_fig=False,
        save_path=f"{dirpath}/{_norm_part(tm_name)}.svo.png",
    )
    plt.close()

    # plot tissue module id card
    if cell2loc:
        # TODO: Fix corresponding svgs not showing.
        tm_colors = np.array(
            ["#9aabe1", "#e2b97c", "#5ac096", "#c17a71", "#c265b6"], dtype=str
        )
        tm_colors = np.tile(
            tm_colors,
            np.ceil(
                len(tm_overlap_show if tm_overlap_show is not None else [])
                / len(tm_colors)
            ),
        ).tolist()
        spg.plot.draw_tissue_module_id_card(
            adata,
            svg_list=tm_svo_list,
            tm=tm_idx,
            tm_overlap_show=tm_overlap_show,
            spatial_info=spatial_info,
            tm_color=tm_colors,
            organism=organism,
            deconvolution_key="cell_type_proportion",
            return_fig=True,
            show_fig=False,
            save_path=None,
        )
        plt.savefig(f"{dirpath}/{_norm_part(tm_name)}.idcard.png")
        plt.close()


### ____________________________________________________________________________
### SpaGFT functionalities - python environment


def env() -> None:
    print(f"python {platform.python_version()}")
    print("### Numerical dependencies:")
    sc.logging.print_header()


def version() -> None:
    print(spg_version)


### ____________________________________________________________________________
### Parser


def _v_dir_path(path: str) -> str:
    normpath = os.path.abspath(path)
    if not os.path.exists(normpath):
        raise ar.ArgumentTypeError(f'"{path}" does not exist')

    if not os.path.isdir(normpath):
        raise ar.ArgumentTypeError(f'"{path}" is not a directory')

    return normpath


def _v_file_path(ext: Union[str, None] = None) -> Callable[[str], str]:
    def validate(path: str):
        normpath = os.path.abspath(path)
        if not os.path.exists(normpath):
            raise ar.ArgumentTypeError(f'"{path}" does not exist')

        if not os.path.isfile(normpath):
            raise ar.ArgumentTypeError(f'"{path}" is not a file')

        if ext and not Path(normpath).suffix == f".{ext}":
            raise ar.ArgumentTypeError(f'"{path}" is not a {ext} file')

        return normpath

    return validate


def _v_dataset_path(path: str) -> str:
    normpath = os.path.abspath(path)
    if not os.path.exists(normpath):
        raise ar.ArgumentTypeError(f'"{path}" does not exist')

    if os.path.isfile(normpath):
        if not Path(normpath).suffix == f".h5ad":
            raise ar.ArgumentTypeError(f'"{path}" is not a .h5ad file')
    elif os.path.isdir(normpath):
        pass
    else:
        raise ar.ArgumentTypeError(f'"{path}" is not a file/directory')

    return normpath


def parser() -> ar.ArgumentParser:
    # usage: SpaGFT [-h] [-env] [-v] [-data PATH] [-o PATH] [-omcs {gene,protein}] [-org {Mouse,Human}] [-flp] [-fcnt] [-fnorm] [-flog] [-cell2loc PATH] [-clustering_alg ALG] [-fident_svos]
    #               [-fident_tms] [-fcmp_tm_clustering] [-save FNAME] [-gsvosum] [-gsvos PATH] [-gtmsum] [-gtms PATH]
    #
    # SpaGFT is a python package to analyze spatial transcriptomics data. It uses a hypothesis-free graph Fourier transform model to identify spatially variable genes/proteins, tissue modules
    # (TM) and functional enrichment analysis to explore the underlying biological processes in these TMs.
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
    parser = ar.ArgumentParser(
        prog="SpaGFT",
        description="""
    SpaGFT is a python package to analyze spatial transcriptomics data. It uses a
    hypothesis-free graph Fourier transform model to identify spatially variable
    genes/proteins, tissue modules (TM) and functional enrichment analysis to
    explore the underlying biological processes in these TMs.
    """,
    )

    parser.add_argument(
        "-env",
        action="store_true",
        help="python environment",
    )
    parser.add_argument(
        "-v",
        action="store_true",
        help="program version",
    )

    group_data = parser.add_argument_group("data")

    group_data.add_argument(
        "-data",
        metavar="PATH",
        type=_v_dataset_path,
        help="input file/folder in .h5ad/visium format",
    )
    group_data.add_argument("-o", metavar="PATH", type=str, help="output directory")

    group_analyze = parser.add_argument_group("analyze")

    group_analyze.add_argument(
        "-omcs",
        type=str,
        choices=["gene", "protein"],
        help="specify the omics kind",
    )
    group_analyze.add_argument(
        "-org",
        type=str,
        choices=["Mouse", "Human"],
        help="specify the organism kind",
    )

    group_analyze.add_argument(
        "-flp",
        action="store_true",
        help="prepare dataset: low pass filter in the graph Fourier domain",
    )
    group_analyze.add_argument(
        "-fcnt",
        action="store_true",
        help="prepare dataset: kNN filter",
    )
    group_analyze.add_argument(
        "-fnorm",
        action="store_true",
        help="prepare dataset: normalize",
    )
    group_analyze.add_argument(
        "-flog",
        action="store_true",
        help="prepare dataset: logarithmize",
    )

    group_analyze.add_argument(
        "-cell2loc",
        metavar="PATH",
        type=_v_file_path(ext="csv"),
        help="augment dataset: add cell2location data specified in .csv file",
    )
    group_analyze.add_argument(
        "-clustering_alg",
        metavar="ALG",
        choices=["louvain", "leiden"],
        default=None,
        help="specify the algorithm used for clustering tms",
    )

    group_analyze.add_argument(
        "-fident_svos",
        action="store_true",
        help="identify spatially variable omics",
    )
    group_analyze.add_argument(
        "-fident_tms",
        action="store_true",
        help="identify tissue modules",
    )
    group_analyze.add_argument(
        "-fcmp_tm_clustering",
        action="store_true",
        help="compare tm clustering with ground truth specified in the dataset annotations",
    )

    group_analyze.add_argument(
        "-save",
        metavar="FNAME",
        type=str,
        help="save analyzed dataset to output directory with the given .h5ad filename",
    )

    group_graph = parser.add_argument_group("graph")

    group_graph.add_argument(
        "-gsvosum",
        action="store_true",
        help="graph the svo summary",
    )
    group_graph.add_argument(
        "-gsvos",
        metavar="PATH",
        type=_v_file_path(ext="csv"),
        help="graph selected svos specified in .csv file",
    )
    group_graph.add_argument(
        "-gtmsum",
        action="store_true",
        help="graph the tissue module summary",
    )
    group_graph.add_argument(
        "-gtms",
        metavar="PATH",
        type=_v_file_path(ext="csv"),
        help="graph selected tissue modules specified in .csv file",
    )

    return parser


def process_actions(parser: ar.ArgumentParser, args: ar.Namespace) -> None:
    ### info ###
    if not vars(args):
        parser.print_help()
        exit(0)

    if args.v:
        version()
        exit(0)

    if args.env:
        env()
        exit(0)

    ### data ###
    if args.data is None:
        raise ar.ArgumentTypeError("dataset path not specified")
    adata: AnnData = read_dataset(args.data)

    output_dir: Union[str, None] = args.o
    output_dataset_fname: Union[str, None] = args.save

    ### analyze ###
    omics_type: Union[OmicsType, None] = None
    if True:
        omics_type_str: Union[str, None] = args.omcs
        if omics_type_str is not None:
            adata.uns["omics_type"] = omics_type_str
        else:
            if "omics_type" in adata.uns:
                omics_type_str = adata.uns["omics_type"]

        if omics_type_str:
            omics_type = _omics_type(omics_type_str)

    organism: Union[str, None] = args.org
    if organism is not None:
        adata.uns["organism"] = organism
    else:
        if "organism" in adata.uns:
            organism = adata.uns["organism"]

    spatial_info = None
    spg_key_gene = ["array_row", "array_col"]
    spg_key_protein = ["x", "y"]
    key_default = "spatial"
    # Required for compatibility with SpaGFT examples.
    if (
        spg_key_gene[0] in adata.obs_keys() and spg_key_gene[1] in adata.obs_keys()
    ) or (
        spg_key_gene[0] in adata.obsm_keys() and spg_key_gene[1] in adata.obsm_keys()
    ):
        spatial_info = spg_key_gene
    elif (
        spg_key_protein[0] in adata.obs_keys()
        and spg_key_protein[1] in adata.obs_keys()
    ) or (
        spg_key_protein[0] in adata.obsm_keys()
        and spg_key_protein[1] in adata.obsm_keys()
    ):
        spatial_info = spg_key_protein
    elif key_default in adata.obs_keys() or key_default in adata.obsm_keys():
        spatial_info = key_default
    else:
        raise ar.ArgumentTypeError("could not find spatial key in dataset")

    cell2loc: bool = False
    if args.cell2loc is not None:
        cell2loc = True
        read_cell2location(adata, args.cell2loc)
    elif "cell_type_proportion" in adata.obsm:
        cell2loc = True

    filter_low_pass = args.flp or False
    filter_by_count = args.fcnt or False
    normalize = args.fnorm or False
    logarithmize = args.flog or False
    if filter_low_pass or filter_by_count or normalize or logarithmize:
        preprocess(adata, filter_low_pass, filter_by_count, normalize, logarithmize)

    low_cutoff: Union[float, int, None] = None
    high_cutoff: Union[float, int, None] = None
    score_df: Union[DataFrame, None] = None
    svo_list: Union[List[str], None] = None
    if args.fident_svos:
        if omics_type is None:
            raise ar.ArgumentTypeError("omics type not specified")
        identify_svo_list(adata, omics_type, spatial_info)
    if "low_cutoff" in adata.uns:
        low_cutoff = adata.uns["low_cutoff"]
    if "high_cutoff" in adata.uns:
        high_cutoff = adata.uns["high_cutoff"]
    if "score_df" in adata.uns:
        score_df = adata.uns["score_df"]
    if "svo_list" in adata.uns:
        svo_list = adata.uns["svo_list"]

    tm_svo_df: Union[DataFrame, None] = None
    tm_clustering_alg: Union[str, None] = args.clustering_alg
    if args.fident_tms:
        if organism is None:
            raise ar.ArgumentTypeError("organism not specified")
        if omics_type is None:
            raise ar.ArgumentTypeError("omics type not specified")
        if svo_list is None or low_cutoff is None:
            raise ar.ArgumentTypeError(
                "need to identify spatially variable omics first"
            )

        # TODO: This should be read as parameters.
        resolution: Union[Tuple[float, float, float], None] = None
        ratio_neighbors: Union[int, None] = None
        n_neighbors: Union[int, None] = None
        # mouse_brain_coronal
        if organism == "Mouse" and omics_type == OmicsType.GENE:
            resolution = (0.7, 1.8, 0.1)
            ratio_neighbors = 2
            n_neighbors = 15
        # lymphnode_tutorial
        elif organism == "Human" and omics_type == OmicsType.GENE:
            resolution = (0.5, 1.5, 0.1)
            ratio_neighbors = 2
            n_neighbors = 15
        # codex_A6
        elif organism == "Human" and omics_type == OmicsType.PROTEIN:
            resolution = (1.5, 2.0, 0.1)
            ratio_neighbors = 1
            n_neighbors = 5
        else:
            # TODO: Remove these not particularly meaningful defaults.
            resolution = (0.5, 1.5, 0.1)
            ratio_neighbors = 2
            n_neighbors = 15

        identify_tm_list(
            adata,
            svo_list,
            spatial_info,
            low_cutoff,
            resolution,
            ratio_neighbors,
            n_neighbors,
            tm_clustering_alg,
        )
    if "tm_svo_df" in adata.uns:
        tm_svo_df = adata.uns["tm_svo_df"]
    if "tm_clustering_alg" in adata.uns:
        tm_clustering_alg = adata.uns["tm_clustering_alg"]

    tm_clustering_score = None
    if args.fcmp_tm_clustering:
        if svo_list is None:
            raise ar.ArgumentTypeError(
                "need to identify spatially variable omics first"
            )
        if tm_svo_df is None:
            raise ar.ArgumentTypeError("need to identify tissue modules first")

        compare_tm_clustering(adata, svo_list)
    if "tm_clustering_score" in adata.uns:
        tm_clustering_score = adata.uns["tm_clustering_score"]

    if args.save:
        if output_dir is None:
            raise ar.ArgumentTypeError("output directory for dataset not specified")
        if output_dataset_fname is None:
            raise ar.ArgumentTypeError("output dataset filename not specified")

        write_dataset(adata, f"{output_dir}/{output_dataset_fname}")

    ### graph ###
    if args.gsvosum:
        if output_dir is None:
            raise ar.ArgumentTypeError("output directory for svo summary not specified")
        if svo_list is None:
            raise ar.ArgumentTypeError(
                "need to identify spatially variable omics first"
            )

        write_svo_summary(
            output_dir,
            adata,
            svo_list,
        )

    if args.gsvos:
        if output_dir is None:
            raise ar.ArgumentTypeError(
                "output directory for svo selection not specified"
            )
        if omics_type is None:
            raise ar.ArgumentTypeError("omics type not specified")

        svo_selection: List[str] = read_svo_selection(args.gsvos)
        write_svo_selection(
            output_dir,
            adata,
            svo_selection,
            omics_type,
        )

    if args.gtmsum:
        if output_dir is None:
            raise ar.ArgumentTypeError("output directory for tm summary not specified")
        if svo_list is None:
            raise ar.ArgumentTypeError(
                "need to identify spatially variable omics first"
            )

        tm_info = dict()
        if tm_clustering_alg is not None:
            tm_info["tm_clustering_alg"] = tm_clustering_alg
        if tm_clustering_score is not None:
            tm_info["tm_clustering_score"] = tm_clustering_score

        write_tm_summary(
            output_dir,
            adata,
            svo_list,
            tm_info,
        )

    if args.gtms:
        if output_dir is None:
            raise ar.ArgumentTypeError(
                "output directory for tm selection not specified"
            )
        if omics_type is None:
            raise ar.ArgumentTypeError("omics type not specified")
        if organism is None:
            raise ar.ArgumentTypeError("organism not specified")

        tm_selection = read_tm_selection(args.gtms)
        write_tm_selection(
            output_dir,
            adata,
            tm_selection,
            spatial_info,
            organism,
            cell2loc,
        )


def main():
    par = parser()
    args = par.parse_args()
    process_actions(par, args)


if __name__ == "__main__":
    main()
