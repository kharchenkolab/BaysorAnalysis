from collections import defaultdict
import argparse
import os
import json

import numpy as np
import pandas as pd
from tqdm import tqdm

import scanpy as sc


def octagon(poly):
    poly_as_int = poly.astype('int')
    min_x = int(np.min(poly_as_int[:, [0]]))
    max_x = int(np.max(poly_as_int[:, [0]]))
    min_y = int(np.min(poly_as_int[:, [1]]))
    max_y = int(np.max(poly_as_int[:, [1]]))

    summed = np.sum(poly_as_int, axis=1)
    diffed = np.diff(poly_as_int, axis=1)

    min_sum = int(np.min(summed))
    max_sum = int(np.max(summed))
    min_diff = int(np.min(diffed))
    max_diff = int(np.max(diffed))

    return [
        [min_x, min_sum - min_x],
        [min_x, max_diff + min_x],  # ^ Left
        [max_y - max_diff, max_y],
        [max_sum - max_y, max_y],  # ^ Botton
        [max_x, max_sum - max_x],
        [max_x, min_diff + max_x],  # ^ Right
        [min_y - min_diff, min_y],
        [min_sum - min_y, min_y]  # ^ Top
    ]


def cells_dict(df:pd.DataFrame, cell_col:str='cell', gene_col:str='gene', x_col:str='x', y_col:str='y'):
    res = {}

    df_genes = df.groupby([cell_col, gene_col]).size().reset_index(name='count')
    genes_pair_list = {k: {s: int(i) for s,i in df[[gene_col, 'count']].values} for k,df in tqdm(df_genes.groupby(cell_col))}
    xy_dict = {k: list(df[[x_col, y_col]].values.astype(float)) for k,df in tqdm(df[[cell_col, x_col, y_col]].groupby(cell_col))}
    gene_dict_def = dict.fromkeys(df.gene.unique(), 0)

    for cell in xy_dict:
        gd = gene_dict_def.copy()
        gd.update(dict(genes_pair_list[cell]))
        res[cell] = {
            "mappings": {},
            "genes": gd,
            "xy": list(map(np.mean, zip(*xy_dict[cell]))),
            "factors": {},
            "poly": octagon(np.array(xy_dict[cell]))
        }

    return res


def molecules_dict(df, gene_col='gene', x_col='x', y_col='y'):
    molecules_dict = {x[0]: [[int(v) for v in a] for a in x[1][[x_col, y_col]].values.astype(np.int32)]
                      for x in tqdm(df.groupby(gene_col))}
    return molecules_dict


def embedding_dict(emb: pd.DataFrame):
    return {i: tuple([float(v) for v in c[:2]]) for i,c in zip(emb.index, emb.values)}


def get_genes(metadata):
    genes = defaultdict(lambda: {'max': 0, 'cells': {}})
    for cell_id, cell_data in metadata.items():
        for gene_id, expression_level in cell_data['genes'].items():
            gene_data = genes[gene_id]
            gene_data['cells'][cell_id] = expression_level
            if gene_data['max'] < expression_level:
                gene_data['max'] = expression_level
    return genes


# def image_dict():
#     cloud_target = open('cloud_target.txt').read().strip()
#     url = 'https://s3.amazonaws.com/{}/wang.images/info.json'.format(
#         cloud_target
#     )
#     image_dict = {
#         'MERrfish': {
#             'sample': 1,
#             'tileSource': url
#         }
#     }

#     return image_dict


def generate_config(name:str, url:str, description:str=None, uid:str=None, version:str = "1.0.0",
                    spatial_center:tuple=None, pca_center:tuple=None, umap_center:tuple=None):
    if uid is None:
        uid = name

    if description is None:
        description = name

    if spatial_center is None:
        spatial_center = (0, 0)

    if pca_center is None:
        pca_center = (0, 0)

    if umap_center is None:
        umap_center = (0, 0)

    config = {
        "name": name,
        "version": version,
        "description": description,
        "public": True,
        "datasets": [{
            "uid": uid,
            "name": name,
            "files": [
                {
                    "type": "cells",
                    "fileType": "cells.json",
                    "url": f"{url}/cells.json"
                },
                {
                    "type": "cell-sets",
                    "fileType": "cell-sets.json",
                    "url": f"{url}/cell-sets.json"
                },
                {
                    "type": "molecules",
                    "fileType": "molecules.json",
                    "url": f"{url}/molecules.json"
                },
                {
                    "type": "expression-matrix",
                    "fileType": "genes.json",
                    "url": f"{url}/genes.json"
                }
            ]
        }],
        "initStrategy": "auto",
        "coordinationSpace": {
            "embeddingZoom": { "PCA": -5, "UMAP": 4 },
            "embeddingType": { "PCA": "PCA", "UMAP": "UMAP" },
            "embeddingTargetX": { "PCA": pca_center[0], "UMAP": umap_center[0] },
            "embeddingTargetY": { "PCA": pca_center[1], "UMAP": -umap_center[1] },
            "spatialZoom": { "A": -2 },
            "spatialTargetX": { "A": spatial_center[0] },
            "spatialTargetY": { "A": spatial_center[1] },
            "spatialLayers": {
                "A": [
                    {"type": "molecules", "radius": 2, "opacity": 1, "visible": True},
                    {"type": "cells", "opacity": 1, "radius": 50, "visible": True, "stroked": False}
                ]
            }
        },
        "layout": [
            {
                "component": "description",
                "props": { "description": description },
                "x": 0, "y": 0, "w": 2, "h": 1
            },
            { "component": "layerController", "coordinationScopes": { "spatialLayers": "A" }, "x": 0, "y": 1, "w": 2, "h": 4 },
            { "component": "status", "x": 0, "y": 5, "w": 2, "h": 1 },
            {
                "component": "spatial",
                "coordinationScopes": { "spatialZoom": "A", "spatialLayers": "A", "spatialTargetX": "A", "spatialTargetY": "A" },
                "x": 2, "y": 0, "w": 4, "h": 4
            },
            { "component": "genes", "x": 9, "y": 0, "w": 3, "h": 2 },
            { "component": "cellSets", "x": 9, "y": 3, "w": 3, "h": 2 },
            {
                "component": "heatmap",
                "props": { "transpose": True },
                "x": 2, "y": 4, "w": 10, "h": 2
            },
            {
                "component": "scatterplot",
                "coordinationScopes": { "embeddingType": "PCA", "embeddingZoom": "PCA", "embeddingTargetY": "PCA", "embeddingTargetX": "PCA" },
                "x": 6, "y": 0, "w": 3, "h": 2
            },
            {
                "component": "scatterplot",
                "coordinationScopes": { "embeddingType": "UMAP", "embeddingZoom": "UMAP", "embeddingTargetY": "UMAP", "embeddingTargetX": "UMAP" },
                "x": 6, "y": 2, "w": 3, "h": 2
            }
        ]
    }

    return config


def preprocess_data(df: pd.DataFrame, cell_col="cell", gene_col="gene", k:int=30, n_pcs:int=20):
    df["v"] = np.ones(df.shape[0], int)
    adata = sc.AnnData(df[[cell_col, gene_col, "v"]].groupby([cell_col, gene_col]).sum().pivot_table("v", index=cell_col, columns=gene_col).fillna(0).iloc[1:])
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=k, n_pcs=n_pcs)

    sc.tl.umap(adata)
    sc.tl.leiden(adata)

    return adata


def add_to_metadata(metadata, adata):
    for (k,v) in embedding_dict(pd.DataFrame(adata.obsm['X_umap'], index=adata.obs.index)).items():
        metadata[k]['mappings']['UMAP'] = v

    for (k,v) in embedding_dict(pd.DataFrame(adata.obsm['X_pca'], index=adata.obs.index)).items():
        metadata[k]['mappings']['PCA'] = v

    for (k,v) in dict(adata.obs.leiden).items():
        metadata[k]['factors']['cluster'] = v

    return metadata


def parse_arguments():
    parser = argparse.ArgumentParser(description='Convert a segmented spot table to the vitessce format.')
    parser.add_argument('-c', '--cell-col', default='cell', type=str, help='cell column name')
    parser.add_argument('-g', '--gene-col', default='gene', type=str, help='gene column name')
    parser.add_argument('-x', '--x-col', default='x', type=str, help='x column name')
    parser.add_argument('-y', '--y-col', default='y', type=str, help='y column name')
    parser.add_argument('-o', '--output', required=True, type=str, help='output folder')
    parser.add_argument('-u', '--url', required=True, type=str, help='dataset url prefix')
    parser.add_argument('--n-pcs', default=10, type=int, help='number of principal components for scanpy')
    parser.add_argument('-n', '--name', required=True, type=str, help='dataset name')
    parser.add_argument('--description', type=str, help='dataset description')
    parser.add_argument('--uid', type=str, help='dataset uid')

    parser.add_argument('spots', type=str, help='path to the spot data frame')

    args = parser.parse_args()

    return args


def main():
    args = parse_arguments()
    df = pd.read_csv(args.spots)
    df[args.cell_col] = df[args.cell_col].astype(str)

    out_path = args.output
    os.makedirs(out_path, exist_ok=True)

    metadata = cells_dict(df, cell_col=args.cell_col, gene_col=args.gene_col, x_col=args.x_col, y_col=args.y_col)
    mol_dict = molecules_dict(df, gene_col=args.gene_col, x_col=args.x_col, y_col=args.y_col)
    adata = preprocess_data(df, cell_col=args.cell_col, gene_col=args.gene_col, n_pcs=args.n_pcs)

    cell_sets = {
        "version": "0.1.2",
        "datatype": "cell",
        "tree": [
            {
                "name": "leiden",
                "children":[{
                    "name": n,
                    "children": [{"name": n, "set": list(v.values)}]
                } for n,v in adata.obs.leiden.index.groupby(adata.obs.leiden.values).items()]
            }
        ]
    }

    metadata = add_to_metadata(metadata, adata)

    del metadata['0']

    config = generate_config(args.name, url=args.url, description=args.description, uid=args.uid,
                             spatial_center=(df.x.median(), df.y.median()),
                             pca_center=np.mean(adata.obsm['X_pca'], axis=0)[:2].astype(float),
                             umap_center=np.mean(adata.obsm['X_umap'], axis=0).astype(float))

    with open(f"{out_path}/cells.json", "w") as f:
        json.dump(metadata, f, indent=1)

    with open(f"{out_path}/molecules.json", "w") as f:
        json.dump(mol_dict, f, indent=1)

    with open(f"{out_path}/genes.json", "w") as f:
        json.dump(get_genes(metadata), f, indent=1)

    with open(f"{out_path}/cell-sets.json", "w") as f:
        json.dump(cell_sets, f, indent=1)

    with open(f"{out_path}/config.json", "w") as f:
        json.dump(config, f, indent=1)


if __name__ == '__main__':
    main()