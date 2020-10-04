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
import Images
import Plots

using DataFrames
using DataFramesMeta
using NearestNeighbors
using ProgressMeter
using Statistics
using StatsBase

B = Baysor;
```

```julia
@time df_spatial, gene_names = B.load_df("/home/vpetukhov/spatial/Benchmarking/run_results/spacejam2/allen_sm_fish/mask_prior_1/segmentation.csv");
@time dapi_arr = Float16.(Images.load("/home/vpetukhov/data/spatal/SpaceJam2Full/allen_sm_fish/dapi_merged.tiff"));
@time dapi_seg_labels = B.load_segmentation_mask("/home/vpetukhov/data/spatal/SpaceJam2Full/allen_sm_fish/dapi_merged_watershed.tif");

df_spatial[!, :cell_watershed] = denserank(B.staining_value_per_transcript(df_spatial, dapi_seg_labels)) .- 1;
df_spatial[!, :cell_watershed_prior] = deepcopy(df_spatial.cell);
df_spatial.cell .= B.load_df("/home/vpetukhov/spatial/Benchmarking/run_results/spacejam2/allen_sm_fish/no_dapi/segmentation.csv")[1].cell;
```

```julia
@time neighb_cm = B.neighborhood_count_matrix(df_spatial, 50);
@time color_transformation = B.gene_composition_transformation(neighb_cm);
@time gene_colors = B.gene_composition_colors(neighb_cm, color_transformation);
df_spatial[!, :color] = gene_colors;
```

```julia
Plots.pyplot()
```

```julia
B.plot_cell_borders_polygons(@where(df_spatial, :x .> 12000, :x .< 14000, :y .> 12000, :y .< 14000), color=:color, size=(1000, 1000), ticks=true)
```

```julia
xs, ys = (12300, 13300), (12500, 13500)
plts = [B.plot_subset(df_spatial, dapi_arr, xs, ys; ms=3.0, alpha=0.75, dapi_alpha=0.5, cell_col=cs, noise=true, #color_col=cs,
        bandwidth=10.0, grid_step=3.0, min_border_length=20, dens_threshold=1e-3, ticks=false, size_mult=4.0, 
        polygon_line_width=3.0, polygon_alpha=0.65, plot_raw_dapi=false, xlabel="x", ylabel="y", guidefontsize=20) 
    for (cs, bw) in zip([:cell, :cell_watershed_prior], [3.0, 1.75])];

plt = Plots.plot(plts..., size=(400, 800), layout=(2,1))
Plots.savefig(plt, "./plots/allen_smfish.png")
Plots.closeall()
plt
```

```julia
plt = B.plot_subset(df_spatial, dapi_arr, xs, ys; ms=1.25 * ms, alpha=0.75, dapi_alpha=0.5, cell_col=:cell, noise=true, #color_col=:cell,
        bandwidth=2.0, grid_step=0.5, min_border_length=50, dens_threshold=1e-5, ticks=false, size_mult=size_mult, 
        polygon_line_width=3.0, polygon_alpha=0.65, plot_raw_dapi=false, xlabel="x", ylabel="y", guidefontsize=20)

plt = Plots.plot(plt, size=(500, 500))
Plots.savefig(plt, "./plots/iss_no_prior.png")
Plots.closeall()
plt
```

```julia
# plts = [B.plot_comparison_for_cell(df_spatial, xs, ys, nothing, dapi_arr; ms=ms, alpha=0.75, dapi_alpha=0.5, cell_col=cs, paper_polys=polys_paper,
#     bandwidth=2.0, grid_step=1.0, ticks=false, size_mult=size_mult, polygon_line_width=3.0, polygon_alpha=0.65, plot_raw_dapi=false,
#     xlabel="x", ylabel="y", guidefontsize=20) for cs in [:cell ,:cell_paper_prior]];

# plt = Plots.plot(plts..., size=(800, 400))
# Plots.closeall()
# plt
```
