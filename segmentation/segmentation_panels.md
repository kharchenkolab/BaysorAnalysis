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
import MAT
import MultivariateStats
import Plots
import CSV

using DataFrames
using DataFramesMeta
using NearestNeighbors
using ProgressMeter
using Statistics
using StatsBase
```

<!-- #region toc-hr-collapsed=true toc-nb-collapsed=true -->
## MERFISH MTG
<!-- #endregion -->

```julia
@time df_spatial, gene_names = Baysor.load_df("../run_results/merfish_mtg_103/segmentation.csv");
```

```julia
grid_step = 5.0
@time polygons = Baysor.boundary_polygons(df_spatial, df_spatial.cell, grid_step=grid_step, verbose=true, method=:knn);
```

```julia
@time neighb_cm = Baysor.neighborhood_count_matrix(df_spatial, 50);
@time color_transformation = Baysor.gene_composition_transformation(neighb_cm[:, df_spatial.confidence .> 0.95]);
```

```julia
# @time mtx_trans = Baysor.transform(color_transformation, neighb_cm);
# @views mtx_trans .-= minimum(mtx_trans[:, df_spatial.confidence .> 0.9], dims=2)
# @views mtx_trans ./= maximum(mtx_trans[:, df_spatial.confidence .> 0.9], dims=2);
# mtx_colors = deepcopy(mtx_trans)
# mtx_colors[1,:] .*= 100
# mtx_colors[2:3,:] .-= 0.5
# mtx_colors[2:3,:] .*= 500
# gene_colors = vec(mapslices(col -> Colors.Lab(col...), mtx_colors, dims=1));

@time gene_colors = Baysor.gene_composition_colors(neighb_cm, color_transformation; confidences=df_spatial.confidence, color_range=500);
```

```julia
Baysor.plot_dataset_colors(df_spatial, gene_colors, polygons=polygons, min_molecules_per_cell=30,
    min_pixels_per_cell=10, polygon_alpha=0.5)[1]
```

```julia
# n_mols_per_cell = Baysor.count_array(df_spatial.cell .+ 1)[2:end]
# plt_n_mols = Plots.histogram(n_mols_per_cell[(n_mols_per_cell .> 1) .& (n_mols_per_cell .< quantile(n_mols_per_cell, 0.99) / 0.99)], xlabel="Num. molecules per cell", ylabel="Num. cells")
```

<!-- #region toc-hr-collapsed=true toc-nb-collapsed=true -->
### Insets
<!-- #endregion -->

```julia
module T

using DataFrames
using DataFramesMeta

import Baysor
import Plots

function get_plot_inset(df_spatial::DataFrame, min_x::T, max_x::T, min_y::T, max_y::T; target_x::T, target_y::T, grid_step::Float64, upscale::Float64=4.0) where T <: Real
    subs_df = @where(df_spatial, :x .> min_x, :x .< max_x, :y .> min_y, :y .< max_y)
    subs_polygons = Baysor.boundary_polygons(subs_df, subs_df.cell, grid_step=grid_step);
    (px_min, px_max), (py_min, py_max) = [Baysor.val_range(df_spatial[!,s]) for s in (:x, :y)];
    pwidth, pheight = px_max - px_min, py_max - py_min

    nw, nh = min.([upscale * (max_x - min_x) / pwidth, upscale * (max_y - min_y) / pheight], 1.0);
    inset = Plots.bbox((target_x - px_min) / pwidth, (target_y - py_min) / pheight, nw, nh, :bottom, :left);
    return subs_df, subs_polygons, (1, inset)
end

end
```

```julia

```

```julia
subs_df1, subs_polygons1, inset1 = T.get_plot_inset(df_spatial, -3500, -3000, -500, -200; target_x=-2500, target_y=-3000, grid_step=grid_step);
subs_df2, subs_polygons2, inset2 = T.get_plot_inset(df_spatial, 400, 900, 0, 200; target_x=-200, target_y=-1400, grid_step=grid_step);

insets = [inset1, inset2]

plt = Baysor.plot_dataset_colors(df_spatial, gene_colors, min_molecules_per_cell=30, inset=insets)[1]

Baysor.plot_cell_borders_polygons!(subs_df1, subs_polygons1, color=subs_df1.color,
    xaxis=false, yaxis=false, subplot=2)

Baysor.plot_cell_borders_polygons!(subs_df2, subs_polygons2, color=subs_df2.color,
    xaxis=false, yaxis=false, subplot=3)
```

## MERFISH Hippocampus

```julia
@time df_spatial, gene_names = Baysor.load_df("../run_results/merfish_moffit/segmentation.csv");
```

```julia
@time df_adj, gn_adj = Baysor.load_df("/home/vpetukhov/data/spatal/merfish_moffit/merfish_coords_adj.csv");
all(gn_adj[df_adj.gene] .== gene_names[df_spatial.gene])
```

```julia
@time dapi = Images.load("/home/vpetukhov/data/spatal/merfish_moffit/dapi_merged.tiff");
```

```julia
grid_step = 1.0
@time polygons = Baysor.boundary_polygons(df_spatial, df_spatial.cell, grid_step=grid_step);
length(polygons)
```

```julia
@time neighb_cm = Baysor.neighborhood_count_matrix(df_spatial, 70);
@time color_transformation = Baysor.gene_composition_transformation(neighb_cm, df_spatial.confidence);
```

```julia
@time gene_colors = Baysor.gene_composition_colors(neighb_cm, color_transformation; color_range=1500);
```

```julia
@time plt = Baysor.plot_dataset_colors(df_spatial, gene_colors, min_molecules_per_cell=50,
    min_pixels_per_cell=20, ms=3.0);

@time Plots.savefig(plt, "./plots/merfish_hippocampus/gene_coloring.png");
```

```julia
@time plt = Baysor.plot_dataset_colors(df_spatial, gene_colors, polygons=polygons, min_molecules_per_cell=50,
    min_pixels_per_cell=10, polygon_alpha=0.5, ms=1.5, polygon_line_width=1.5)[1];

@time Plots.savefig(plt, "./plots/merfish_hippocampus/polygons.png");
```

```julia
dapi_arr = Float64.(dapi);
```

```julia
df_adj[!, :cell_paper] = ifelse.(ismissing.(df_adj.cell), 0, denserank(df_adj.cell));
df_adj[!, :cell] = df_spatial.cell;
df_adj[!, :cluster] = df_spatial.cluster;
df_adj[!, :color] = gene_colors;
df_adj[!, :xr] = df_spatial.x;
df_adj[!, :yr] = df_spatial.y;
```

```julia
module T

using DataFrames
using DataFramesMeta

import Plots
import Baysor

function plot_subset(df_adj::DataFrame, dapi_arr::Matrix{Float64}, (xs, xe), ((ys, ye)); plot_polygons::Bool=true)
    df_subs = @where(df_adj, :xr .> xs, :xr .< xe, :yr .> ys, :yr .< ye);
    pol_subs_kde = plot_polygons ? Baysor.boundary_polygons(df_subs, df_subs.cell, grid_step=5.0, min_molecules_per_cell=10, method=:kde, verbose=true, bandwidth=8.0) : Matrix{Float64}[]

    ysa, yea = Int(minimum(df_subs.y)), Int(maximum(df_subs.y));
    xsa, xea = Int(minimum(df_subs.x)), Int(maximum(df_subs.x));
    
    xticks = range(0, xea-xsa, length=5)
    yticks = range(0, yea-ysa, length=5)

    plt1 = Plots.heatmap(1.0 .- dapi_arr[ysa:yea, xsa:xea], color=:grayscale, colorbar=:none, size=((xea-xsa), yea-ysa) ./ 3,
                alpha=0.9, format=:png, ticks=false, legend=:none)
    Baysor.plot_cell_borders_polygons!(df_subs, pol_subs_kde, color=df_subs.color, ms=2.0, alpha=0.2, offset=(-xsa, -ysa),
        polygon_line_width=2, polygon_alpha=0.75, is_noise=df_subs.cell .== 0, noise_kwargs=Dict(:ms => 1.0))

    Plots.vline!(xticks, color="black", alpha=0.5)
    Plots.hline!(yticks, color="black", alpha=0.5)

    # Possible colorschemes: tarn, diff, lime_grad, thermal
    plt2 = Plots.heatmap(dapi_arr[ysa:yea, xsa:xea], color=:diff, colorbar=:none,
            alpha=0.9, format=:png, ticks=false, legend=:none)

    Plots.plot!([Plots.Shape(pg[:,1] .- xsa, pg[:,2] .- ysa) for pg in pol_subs_kde], 
        fill=(0, 0.0), linewidth=2.0, linecolor="black", alpha=0.4, label="", xlims=(0, (xea-xsa)), ylims=(0, (yea-ysa)));

    Plots.vline!(xticks, color="black", alpha=0.5)
    Plots.hline!(yticks, color="black", alpha=0.5)

    return Plots.plot(plt1, plt2, layout=2, size=(2 * (xea-xsa), yea-ysa) ./ 3)
end

rectangle((xs, xe), (ys, ye)) = Plots.Shape([xs, xe, xe, xs], [ys, ys, ye, ye])

end
```

```julia
subset_coords = [((-2475, -2235), (-3300, -3050)), ((-2872, -2822), (-2900, -2800)), ((-2825, -2685), (-3000, -2900)), ((-2090, -1950), (-2750, -2600)), ((-2928, -2825), (-3850, -3700)), ((-2857, -2807), (-2800, -2700))];
```

```julia
@time plt = Baysor.plot_dataset_colors(df_spatial, gene_colors, polygons=polygons, min_molecules_per_cell=50,
    min_pixels_per_cell=10, polygon_alpha=0.5, ms=1.5, polygon_line_width=1.5)[1];

plt = Plots.plot!([T.rectangle(sc...) for sc in subset_coords], fill=(0, 0.0), color=colorant"#5603fc", lw=3.0)

@time Plots.savefig(plt, "./plots/merfish_hippocampus/polygons.png");
```

```julia
plt = T.plot_subset(df_adj, dapi_arr, subset_coords[1]...)
Plots.savefig(plt, "./plots/merfish_hippocampus/dapi_subs1.png");
plt
```

```julia
plt = T.plot_subset(df_adj, dapi_arr, subset_coords[2]...)
Plots.savefig(plt, "./plots/merfish_hippocampus/dapi_subs2.png");
plt
```

```julia
plt = T.plot_subset(df_adj, dapi_arr, subset_coords[6]...)
Plots.savefig(plt, "./plots/merfish_hippocampus/dapi_subs2_2.png");
plt
```

```julia
plt = T.plot_subset(df_adj, dapi_arr, subset_coords[3]...)
Plots.savefig(plt, "./plots/merfish_hippocampus/dapi_subs3.png");
plt
```

```julia
plt = T.plot_subset(df_adj, dapi_arr, subset_coords[4]...)
Plots.savefig(plt, "./plots/merfish_hippocampus/dapi_subs4.png");
plt
```

```julia
plt = T.plot_subset(df_adj, dapi_arr, subset_coords[5]...)
Plots.savefig(plt, "./plots/merfish_hippocampus/dapi_subs5.png");
plt
```

```julia
brightness_per_mol = dapi_arr[CartesianIndex.(Int.(df_adj.y), Int.(df_adj.x))];
```

```julia
plt = Plots.histogram(brightness_per_mol[df_adj.cell .> 0], bins=80, label="Assigned", xlabel="DAPI brightness", ylabel="Num. of molecules")
plt = Plots.histogram!(brightness_per_mol[df_adj.cell .== 0], bins=80, label="Noise")
# Plots.savefig(plt, "./plots/merfish_hippocampus/dapi_brightness.png");
plt
```

<!-- #region toc-hr-collapsed=true toc-nb-collapsed=true -->
## ISS Mouse
<!-- #endregion -->

```julia
# @time dapi = MAT.matread("/home/vpetukhov/data/spatal/SpaceJam2/iss_1_HCA_11_1_hm_133_c1-1.tif.mat")["I"];
```

```julia
# @time plt = Baysor.plot_dataset_colors(df_spatial, gene_colors, polygons=polygons, min_molecules_per_cell=5,
#     min_pixels_per_cell=30, polygon_alpha=0.5, ms=3.0)[1];

# @time Plots.savefig(plt, "test2.png");
# close()
```

```julia
# @time df_spatial, gene_names = Baysor.load_df("../run_results/iss_1_s10/segmentation.csv");
@time df_spatial, gene_names = Baysor.load_df("../run_results/iss_mouse_3/segmentation.csv");
```

```julia
# # @time dapi = Images.load("/home/vpetukhov/spatial/JamboreeV2/iss_mouse_3_dapi.tif");
# @time dapi = Images.load("/home/vpetukhov/spatial/JamboreeV2/iss_mouse_1_HCA_11_3_m_133_c1-1.tif");
# dapi_arr = Float64.(dapi);
```

```julia
df_spatial[!, :x] = round.(Int, df_spatial.x / 0.1625);
df_spatial[!, :y] = round.(Int, df_spatial.y / 0.1625);
```

```julia
maximum(df_spatial.y), minimum(df_spatial.y)
```

```julia
grid_step = 10.0
@time polygons = Baysor.boundary_polygons(df_spatial, df_spatial.cell, grid_step=grid_step, verbose=true, bandwidth=10.0);
length(polygons)
```

```julia
@time neighb_cm = Baysor.neighborhood_count_matrix(df_spatial, 50);
@time color_transformation = Baysor.gene_composition_transformation(neighb_cm, df_spatial.confidence);
```

```julia
# @time mtx_trans = Baysor.transform(color_transformation, neighb_cm);
# min_vals = minimum(color_transformation.embedding, dims=2)
# mtx_trans .-= min_vals
# mtx_trans ./= (maximum(color_transformation.embedding, dims=2) - min_vals);

# mtx_colors = deepcopy(mtx_trans)
# mtx_colors[1,:] .*= 100
# mtx_colors[2:3,:] .-= 0.5
# mtx_colors[2:3,:] .*= 200
# gene_colors = vec(mapslices(col -> Colors.Lab(col...), mtx_colors, dims=1));

@time gene_colors = Baysor.gene_composition_colors(neighb_cm, color_transformation; color_range=200);
```

```julia
@time plt = Baysor.plot_dataset_colors(df_spatial, gene_colors, polygons=polygons, min_molecules_per_cell=5,
    min_pixels_per_cell=15, polygon_alpha=0.5, ms=2.5, alpha=0.3, polygon_line_width=2.0)[1];

@time Plots.savefig(plt, "./plots/iss_mouse_3/polygons_m3.png");
```

<!-- #region toc-hr-collapsed=true toc-nb-collapsed=true -->
## STARmap
<!-- #endregion -->

```julia
@time df_spatial, gene_names = Baysor.load_df("../run_results/star_map_vis_1020/segmentation.csv");
size(df_spatial)
```

```julia
grid_step = 10.0
@time polygons = Baysor.boundary_polygons(df_spatial, df_spatial.cell, grid_step=grid_step, method=:kde_opt, verbose=true, bandwidth=10.0);
length(polygons)
```

```julia
@time neighb_cm = Baysor.neighborhood_count_matrix(df_spatial, 50);
@time color_transformation = Baysor.gene_composition_transformation(neighb_cm, df_spatial.confidence);
```

```julia
@time mtx_trans = Baysor.transform(color_transformation, neighb_cm);
min_vals = minimum(color_transformation.embedding, dims=2)
mtx_trans .-= min_vals
mtx_trans ./= (maximum(color_transformation.embedding, dims=2) - min_vals);

mtx_colors = deepcopy(mtx_trans)
mtx_colors[1,:] .*= 100
mtx_colors[2:3,:] .-= 0.5
mtx_colors[2:3,:] .*= 500
gene_colors = vec(mapslices(col -> Colors.Lab(col...), mtx_colors, dims=1));

# @time gene_colors = Baysor.gene_composition_colors(neighb_cm, color_transformation; color_range=500);
```

```julia
@time plt = Baysor.plot_dataset_colors(df_spatial, gene_colors, polygons=polygons, min_molecules_per_cell=50,
    min_pixels_per_cell=15, polygon_alpha=0.5, ms=3.0, alpha=0.075, polygon_line_width=2.0)[1];

@time Plots.savefig(plt, "./plots/star_map_vis_1020/polygons.png");
```

## osm-FISH

```julia
@time df_spatial, gene_names = Baysor.load_df("../run_results/osm_fish/segmentation.csv");
size(df_spatial)
```

```julia
grid_step = 10.0
@time polygons = Baysor.boundary_polygons(df_spatial, df_spatial.cell, grid_step=grid_step, verbose=true, bandwidth=15.0, min_molecules_per_cell=10);
length(polygons)
```

```julia
@time neighb_cm = Baysor.neighborhood_count_matrix(df_spatial, 50);
@time color_transformation = Baysor.gene_composition_transformation(neighb_cm, df_spatial.confidence);
```

```julia
@time mtx_trans = Baysor.transform(color_transformation, neighb_cm);
min_vals = minimum(color_transformation.embedding, dims=2)
mtx_trans .-= min_vals
mtx_trans ./= (maximum(color_transformation.embedding, dims=2) - min_vals);
```

```julia
mtx_colors = deepcopy(mtx_trans)
mtx_colors[1,:] .*= 100
mtx_colors[2:3,:] .-= 0.5
mtx_colors[2:3,:] .*= 400
gene_colors = vec(mapslices(col -> Colors.Lab(col...), mtx_colors, dims=1));

# @time gene_colors = Baysor.gene_composition_colors(neighb_cm, color_transformation; color_range=500);
```

```julia
@time plt = Baysor.plot_dataset_colors(df_spatial, gene_colors, polygons=polygons, min_molecules_per_cell=30,
    min_pixels_per_cell=15, polygon_alpha=0.5, ms=3.0, alpha=0.075, polygon_line_width=2.0)[1];

@time Plots.savefig(plt, "./plots/osm_fish/polygons.png");
```

<!-- #region toc-hr-collapsed=true toc-nb-collapsed=true -->
## Test plots
<!-- #endregion -->

```julia
import GR
import Plots
```

```julia
# x = 0:pi/100:2*pi;
# y = sin.(x);
# 
# GR.beginprint("test_gr.png")
# GR.setlinewidth(0.01)
# GR.setmarkersize(0.0)
# GR.polyline(x, y)
# GR.polymarker(x, sqrt.(abs.(y)))
# GR.show()
# GR.endprint()
```

```julia
Plots.pyplot()
```

```julia
plt = Baysor.plot_dataset_colors(sdf, sdf.color, polygons=spols, min_molecules_per_cell=5,
    min_pixels_per_cell=30, polygon_alpha=0.5, ms=3.0)[1]
Plots.plot!(3000:3500, 4000:4500, lw=0.1)

@time Plots.savefig(plt, "test_idr.png")
```

```julia
plt = Baysor.plot_dataset_colors(sdf, sdf.color, polygons=spols, min_molecules_per_cell=5,
    min_pixels_per_cell=60, polygon_alpha=0.5, ms=3.0)[1]
Plots.plot!(3000:3500, 4000:4500, lw=0.001)

@time Plots.savefig(plt, "test.png")
```

```julia
# Baysor.plot_dataset_colors(sdf, sdf.color, polygons=spols, min_molecules_per_cell=5,
#     min_pixels_per_cell=30, polygon_alpha=0.5, ms=3.0)[1]
```

```julia
plt = Baysor.plot_dataset_colors(df_spatial, gene_colors, polygons=polygons[1:5000], min_molecules_per_cell=5,
    min_pixels_per_cell=30, polygon_alpha=0.5, ms=3.0)[1];

@time Plots.savefig(plt, "test.png")
```

```julia
plt = Baysor.plot_dataset_colors(df_spatial, gene_colors, polygons=polygons, min_molecules_per_cell=5,
    min_pixels_per_cell=30, polygon_alpha=0.5, ms=3.0)[1];

@time Plots.savefig(plt, "test.png")
```

## Allen sm-FISH

```julia
@time df_spatial, gene_names = Baysor.load_df("../run_results/spacejam2/allen_sm_fish/no_dapi/segmentation.csv");
df_spatial[!, :x] = round.(Int, 10 .* (df_spatial.x .- minimum(df_spatial.x)));
df_spatial[!, :y] = round.(Int, 10 .* (df_spatial.y .- minimum(df_spatial.y)));
```

```julia
@time dapi = Images.load("/home/vpetukhov/data/spatal/SpaceJam2Full/allen_sm_fish/dapi_merged.tiff");
dapi_arr = Float64.(dapi);
```

```julia
@time seg_labels = Images.load("/home/vpetukhov/data/spatal/SpaceJam2Full/allen_sm_fish/segmentation_labels_from_json.tiff");
seg_labels = round.(Int, Float64.(seg_labels) ./ 1.5259021896696422e-5);
```

```julia
grid_step = 10.0
@time polygons = Baysor.boundary_polygons(df_spatial, df_spatial.cell, grid_step=grid_step, bandwidth=10.0);
length(polygons)
```

```julia
@time neighb_cm = Baysor.neighborhood_count_matrix(df_spatial, 40);
@time color_transformation = Baysor.gene_composition_transformation(neighb_cm, df_spatial.confidence; sample_size=20000, spread=2.0, min_dist=0.1);
```

```julia
# pc2 = MultivariateStats.transform(MultivariateStats.fit(MultivariateStats.PCA, neighb_cm, maxoutdim=2), neighb_cm);
# t_ids = Baysor.select_ids_uniformly(pc2[1,:], pc2[2,:], df_spatial.confidence, 20000)

# mtx_trans = deepcopy(color_transformation.embedding)
# @views mtx_trans .-= minimum(mtx_trans, dims=2)
# @views mtx_trans ./= maximum(mtx_trans, dims=2);

# mtx_trans[1,:] .*= 100
# mtx_trans[2:3,:] .-= 0.5
# mtx_trans[2:3,:] .*= 500
# gene_colors = vec(mapslices(col -> Colors.Lab(col...), mtx_trans, dims=1));
```

```julia
# @time color_transformation = Baysor.gene_composition_transformation(neighb_cm, df_spatial; sample_size=20000, spread=2.0, min_dist=0.1);
@time gene_colors = Baysor.gene_composition_colors(neighb_cm, color_transformation; color_range=500);
```

```julia
@time plt = Baysor.plot_dataset_colors(df_spatial, gene_colors, polygons=polygons, min_molecules_per_cell=50,
    polygon_line_width=2, min_pixels_per_cell=30, ms=1.25, alpha=1.0)[1];

@time Plots.savefig(plt, "./plots/allen_sm_fish/polygons.png");
```

```julia
borders_per_label = Baysor.grid_borders_per_label(seg_labels);
borders_per_label = borders_per_label[map(length, borders_per_label) .>= 3]
paths = Baysor.longest_paths.(Baysor.border_mst.(borders_per_label));

polygons_paper = vcat([[borders[p] for p in cur_paths] for (borders, cur_paths) in zip(borders_per_label, paths)]...);
polygons_paper = [Float64.(hcat(p...)') for p in polygons_paper];
```

```julia
@time plt = Baysor.plot_dataset_colors(df_spatial, gene_colors, polygons=polygons_paper, min_molecules_per_cell=50,
    min_pixels_per_cell=50, ms=5.0)[1];

@time Plots.savefig(plt, "./plots/allen_sm_fish/paper_polygons.png");
close()
```

```julia
# @time begin
#     plt = Plots.heatmap(1. .- dapi, size=(2300, 2100), color=:grayscale, format=:png)
#     plt = Baysor.plot_cell_borders_polygons!(df_spatial, polygons, color=gene_colors, ms=5.0, alpha=0.2,
#         polygon_line_width=3, polygon_alpha=0.75, is_noise=df_spatial.cell .== 0, noise_kwargs=Dict(:markersize => 3.0))

#     Plots.savefig(plt, "./plots/merfish_hippocampus/dapi_polygons.png");
#     close()
# end
```

```julia
@time begin
    plt = Plots.plot(dapi_arr_inv, size=(2300, 2100), color=:grayscale)
    plt = Baysor.plot_cell_borders_polygons!(df_spatial, polygons, color=gene_colors, ms=5.0, alpha=0.2,
        polygon_line_width=3, polygon_alpha=0.75, is_noise=df_spatial.cell .== 0, noise_kwargs=Dict(:markersize => 3.0))

    Plots.savefig(plt, "./plots/merfish_hippocampus/dapi_polygons.png");
    close()
end
```

```julia
xs, xe = -2500, -2100;
ys, ye = -3300, -3050;

df_subs = @where(df_spatial, :x .> xs, :x .< xe, :y .> ys, :y .< ye);
# @time pol_subs = Baysor.boundary_polygons(df_subs, df_subs.cell, grid_step=10.0);
@time pol_subs_kde = Baysor.boundary_polygons(df_subs, df_subs.cell, grid_step=5.0, min_molecules_per_cell=10, method=:kde, verbose=true, bandwidth=5.0);
@time pol_subs_paper = Baysor.boundary_polygons(df_subs, df_subs.cell_paper, grid_step=5.0, min_molecules_per_cell=10, method=:kde, verbose=true, bandwidth=5.0);

ysa, yea = Int(minimum(df_subs.y)), Int(maximum(df_subs.y));
xsa, xea = Int(minimum(df_subs.x)), Int(maximum(df_subs.x));
```

```julia
@time begin
    plt = Plots.heatmap(1.0 .- dapi_arr[ysa:yea, xsa:xea], size=(xea-xsa, yea-ysa), color=:grayscale)
    plt = Baysor.plot_cell_borders_polygons!(df_subs, pol_subs_kde, color=df_subs.color, ms=5.0, alpha=0.2, offset=(-xsa, -ysa),
        polygon_line_width=3, polygon_alpha=0.75, is_noise=df_subs.cell .== 0, noise_kwargs=Dict(:markersize => 3.0))

    Plots.savefig(plt, "./plots/merfish_hippocampus/dapi_subs.png");
    close()
end
```

```julia
@time begin
    plt = Plots.heatmap(1.0 .- dapi_arr[ysa:yea, xsa:xea], size=(xea-xsa, yea-ysa), color=:grayscale)
    plt = Baysor.plot_cell_borders_polygons!(df_subs, pol_subs_paper, color=df_subs.color, ms=5.0, alpha=0.2, offset=(-xsa, -ysa),
        polygon_line_width=3, polygon_alpha=0.75, is_noise=df_subs.cell_paper .== 0, noise_kwargs=Dict(:markersize => 3.0))

    Plots.savefig(plt, "./plots/merfish_hippocampus/dapi_subs_paper.png");
    close()
end
```

```julia
brightness_per_mol = dapi_arr[CartesianIndex.(Int.(df_adj.y), Int.(df_adj.x))];
```

```julia
Plots.histogram(brightness_per_mol[df_adj.cell .> 0], bins=80, label="Assigned", xlabel="DAPI brightness", ylabel="Num. of molecules")
Plots.histogram!(brightness_per_mol[df_adj.cell .== 0], bins=80, label="Noise")
Plots.savefig(plt, "./plots/merfish_hippocampus/dapi_brightness.png");
```

## Ex-Seq

```julia
@time df_spatial, gene_names = Baysor.load_df("/home/vpetukhov/data/spatal/SpaceJam2Full/ex_seq/spottable_exseq.csv";
    x_col=:PositionsGlobalPix_1, y_col=:PositionsGlobalPix_2, gene_col=:GeneName);
size(df_spatial)
```

```julia
Baysor.append_confidence!(df_spatial, nn_id=30);
```

```julia
@time neighb_cm = Baysor.neighborhood_count_matrix(df_spatial, 40);
@time color_transformation = Baysor.gene_composition_transformation(neighb_cm, df_spatial.confidence; sample_size=20000, spread=2.0, min_dist=0.1);

@time gene_colors = Baysor.gene_composition_colors(neighb_cm, color_transformation; color_range=500);
```

```julia
Baysor.val_range(df_spatial.PositionsUm_3)
```

```julia
Baysor.val_range(df_spatial.PositionsGlobalPix_3)
```

```julia
@time plt = Baysor.plot_dataset_colors(df_spatial, gene_colors, min_molecules_per_cell=50,
    min_pixels_per_cell=20, ms=3.0);

@time Plots.savefig(plt, "/home/vpetukhov/tmp/ex_seq_gene_coloring.png");
```

```julia
df_spatial_z1 = @where(df_spatial, :PositionsGlobalPix_3 .< 50);

@time ncm1 = Baysor.neighborhood_count_matrix(df_spatial_z1, 70);
@time ct1 = Baysor.gene_composition_transformation(ncm1, df_spatial_z1.confidence; sample_size=10000, spread=2.0, min_dist=0.1);
@time gc1 = Baysor.gene_composition_colors(ncm1, ct1; color_range=500);
```

```julia
plt = Baysor.plot_cell_borders_polygons(df_spatial_z1, color=gc1, ms=2.0, size=(1000, 1000), xlims=(0, 17500), ylims=(0, 18000))
@time Plots.savefig(plt, "/home/vpetukhov/tmp/ex_seq_gene_coloring_z1.png");
plt
```

```julia
df_spatial_z2 = @where(df_spatial, :PositionsGlobalPix_3 .>= 50, :PositionsGlobalPix_3 .< 100);

@time ncm2 = Baysor.neighborhood_count_matrix(df_spatial_z2, 70);
@time gc2 = Baysor.gene_composition_colors(ncm2, ct1; color_range=500);
```

```julia
plt = Baysor.plot_cell_borders_polygons(df_spatial_z2, color=gc2, ms=2.0, size=(1000, 1000), xlims=(0, 17500), ylims=(0, 18000))
@time Plots.savefig(plt, "/home/vpetukhov/tmp/ex_seq_gene_coloring_z2.png");
plt
```

```julia
df_spatial_z3 = @where(df_spatial, :PositionsGlobalPix_3 .>= 100);

@time ncm3 = Baysor.neighborhood_count_matrix(df_spatial_z3, 70);
@time gc3 = Baysor.gene_composition_colors(ncm2, ct1; color_range=500);
```

```julia
plt = Baysor.plot_cell_borders_polygons(df_spatial_z3, color=gc3, ms=2.0, size=(1000, 1000), xlims=(0, 17500), ylims=(0, 18000))
@time Plots.savefig(plt, "/home/vpetukhov/tmp/ex_seq_gene_coloring_z3.png");
plt
```

```julia
Baysor.plot_expression_vectors(Baysor.count_array(df_spatial.gene), gene_names=gene_names)
```

```julia
@time plt = Baysor.plot_dataset_colors(df_spatial_z1, gc1, min_molecules_per_cell=30,
    min_pixels_per_cell=20, ms=2.0);

@time Plots.savefig(plt, "/home/vpetukhov/tmp/ex_seq_gene_coloring_z1.png");
```
