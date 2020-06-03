---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.4.1
  kernelspec:
    display_name: Julia 1.3.1
    language: julia
    name: julia-1.3
---

```julia
import Baysor
import Colors
import CSV
import Clustering
import Images
import MAT
import MultivariateStats
import Plots
import PlotThemes

using DataFrames
using DataFramesMeta
using NearestNeighbors
using ProgressMeter
using RCall
using SparseArrays
using Statistics
using StatsBase
using StatsPlots

B = Baysor;
```

```julia
R"""
library(ggplot2)
library(ggrastr)
library(ggridges)
theme_set(theme_bw())
""";
```

```julia
@time df_spatial, gene_names = B.load_df("/home/vpetukhov/spatial/Benchmarking/run_results/iss_hippo/ca1_no_prior/segmentation.csv");
@time dapi_arr = Float16.(Images.load("/home/vpetukhov/data/spatal/iss/hippocampus/CA1/Viktor/CA1DapiBoundaries_4-3_right_med_filt.tif"));
@time dapi_seg_labels = Matrix(B.load_segmentation_mask("/home/vpetukhov/data/spatal/iss/hippocampus/CA1/Viktor/CA1DapiBoundaries_4-3_right_watershed.tif"));
df_spatial[!, :cell_dapi] = denserank(B.staining_value_per_transcript(df_spatial, dapi_seg_labels)) .- 1;
df_spatial[!, :cell_prior] = B.load_df("/home/vpetukhov/spatial/Benchmarking/run_results/iss_hippo/ca1_watershed_prior/segmentation.csv")[1].cell;
```

```julia
df_spatial[!, :cell_prior] = B.load_df("/home/vpetukhov/spatial/Benchmarking/run_results/iss_hippo/ca1_watershed_prior/segmentation.csv")[1].cell;
```

```julia
@time neighb_cm = B.neighborhood_count_matrix(df_spatial, 30);
@time color_transformation = B.gene_composition_transformation(neighb_cm);
@time gene_colors = B.gene_composition_colors(neighb_cm, color_transformation; color_range=400);
df_spatial[!, :color] = gene_colors;
```

```julia
@time paper_polys = B.boundary_polygons(df_spatial, df_spatial.cell_dapi, bandwidth=2.0, grid_step=1.0, min_border_length=50);
```

```julia
# # B.plot_comparison_for_cell(df_spatial, (5625, 5950), (2000, 2250), dapi_seg_labels, dapi_arr; ms=4.0, alpha=0.75,
# B.plot_comparison_for_cell(df_spatial, (5650, 5800), (2100, 2250), dapi_seg_labels, dapi_arr; ms=4.0, alpha=0.75,
#     bandwidth=2.0, grid_step=1.0, ticks=true, size_mult=3.0, polygon_line_width=3.0, polygon_alpha=0.5)
```

```julia
# B.plot_comparison_for_cell(df_spatial, (5650, 5800), (2100, 2250), nothing, dapi_arr; ms=4.0, alpha=0.75, paper_polys=paper_polys,
#     bandwidth=2.0, grid_step=1.0, ticks=true, size_mult=3.0, polygon_line_width=3.0, polygon_alpha=0.5)
```

```julia
# B.plot_comparison_for_cell(df_spatial, (5650, 5800), (2100, 2250), nothing, dapi_arr; ms=4.0, alpha=0.75, paper_polys=paper_polys,
#     bandwidth=2.0, grid_step=1.0, ticks=true, size_mult=3.0, polygon_line_width=3.0, polygon_alpha=0.5, cell_col=:cell_prior)
```

```julia
# xs, ys = (5650, 5800), (2100, 2250)
xs, ys = (5700, 5800), (2100, 2200)
size_mult = 4.0
ms = 6.0
plts = [B.plot_comparison_for_cell(df_spatial, xs, ys, dapi_seg_labels, dapi_arr; ms=ms, alpha=0.75, dapi_alpha=0.5, cell_col=cs,
    bandwidth=2.0, grid_step=1.0, ticks=true, size_mult=size_mult, polygon_line_width=3.0, polygon_alpha=0.65, plot_raw_dapi=false) for cs in [:cell ,:cell_prior]];

plt = Plots.plot(plts..., size=(800, 400))
Plots.savefig(plt, "./plots/iss_dense.png")
plt
```
