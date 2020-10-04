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

```julia execution={"iopub.status.busy": "2020-09-22T09:56:20.011000+02:00", "iopub.execute_input": "2020-09-22T09:56:20.011000+02:00", "iopub.status.idle": "2020-09-22T09:57:02.990000+02:00"}
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

```julia execution={"iopub.status.busy": "2020-09-22T09:57:02.990000+02:00", "iopub.execute_input": "2020-09-22T09:57:05.721000+02:00", "iopub.status.idle": "2020-09-22T09:57:07.350000+02:00"}
R"""
library(ggplot2)
library(ggrastr)
theme_set(theme_bw())
""";
```

```julia execution={"iopub.status.busy": "2020-09-22T09:57:07.350000+02:00", "iopub.execute_input": "2020-09-22T09:57:07.350000+02:00", "iopub.status.idle": "2020-09-22T09:57:07.841000+02:00"}
module T

using DataFrames
using StatsBase
using Statistics

import Baysor
B = Baysor

function estimate_local_complexity(df_spatial::DataFrame, k::Int; n_samples::Int=100, cell_col::Symbol=:cell)
    neighb_cm = B.neighborhood_count_matrix(df_spatial, k; normalize_by_dist=false);

    n_genes = maximum(df_spatial.gene)
    est_ent_per_cell = [mean([entropy(B.prob_array(rand(gids, k), max_value=n_genes)) for i in 1:n_samples]) 
        for gids in B.split(df_spatial.gene, df_spatial[!, cell_col] .+ 1)[2:end]];

    obs_ents_per_mol = vec(mapslices(entropy, neighb_cm, dims=1));
    
    est_n_genes_per_cell = [mean([sum(length(unique(rand(gids, k)))) for i in 1:n_samples]) 
        for gids in B.split(df_spatial.gene, df_spatial[!, cell_col] .+ 1)[2:end]];

    obs_n_genes_per_mol = vec(sum(neighb_cm .> 1e-10, dims=1));
    
    return (est_ent_per_cell=est_ent_per_cell, obs_ents_per_mol=obs_ents_per_mol, est_n_genes_per_cell=est_n_genes_per_cell, 
        obs_n_genes_per_mol=obs_n_genes_per_mol, k=k)
end

end
```

```julia execution={"iopub.status.busy": "2020-09-22T09:57:07.841000+02:00", "iopub.execute_input": "2020-09-22T09:57:07.841000+02:00", "iopub.status.idle": "2020-09-22T09:57:08.649000+02:00"}
k_vals = [5, 10, 50]
```

## STARMAP VISp 160 20171120_BF4_light

```julia execution={"iopub.status.busy": "2020-09-22T09:57:08.649000+02:00", "iopub.execute_input": "2020-09-22T09:57:08.649000+02:00", "iopub.status.idle": "2020-09-22T09:57:29.787000+02:00"}
@time df_spatial, gene_names = B.load_df("/home/vpetukhov/spatial/Benchmarking/run_results/star_map/visp160/cl0/segmentation.csv");
# @time df_spatial, gene_names = B.load_df("$PROJECT_DIR/run_results/test_small_cell_penalty/star_map_visp1020/segmentation.csv");
df_spatial.x ./= 2
df_spatial.y ./= 2

@time dapi_seg_labels = B.load_segmentation_mask("/home/vpetukhov/data/spatal/star_map/vis160_20171120_BF4_light/segmentation.tiff");
df_spatial[!, :cell_dapi] = denserank(B.staining_value_per_transcript(df_spatial, dapi_seg_labels)) .- 1;
df_spatial.cell .= denserank(df_spatial.cell) .- 1;

# dapi_seg_labels = nothing;
GC.gc()
```

```julia
comp_per_k = @showprogress map(k -> T.estimate_local_complexity(df_spatial, k; cell_col=:cell_dapi), k_vals);
```

```julia execution={"iopub.status.busy": "2020-09-22T09:57:29.787000+02:00", "iopub.execute_input": "2020-09-22T09:57:29.788000+02:00", "iopub.status.idle": "2020-09-22T09:57:36.991000+02:00"}
Plots.pyplot()
```

```julia
real_mol_mask = (df_spatial.cell_dapi .> 0);
plts = Plots.Plot[]

for (k, cv) in zip(k_vals, comp_per_k)
    t_bins = B.hist_bins(cv.est_ent_per_cell, cv.obs_ents_per_mol, min_val=-0.05, m_quantile=1.0, n_bins=30)
    plt1 = Plots.barhist(cv.est_ent_per_cell[df_spatial.cell_dapi[real_mol_mask]], bins=t_bins, legend=(k == 5 ? :top : :none), label="Expected",
        xlabel="Entropy", ylabel="Num. of\nmolecules", bg_legend=Plots.RGBA(1, 1, 1, 0.1), legendfontsize=12, guidefontsize=13, tickfontsize=10)
    Plots.barhist!(cv.obs_ents_per_mol[real_mol_mask], bins=t_bins, alpha=0.7, label="Observed")
#     Plots.barhist!(cv.obs_ents_per_mol[.!real_mol_mask], bins=t_bins, alpha=0.5, label="Observed, not assigned")
    
    t_bins = B.hist_bins(cv.est_n_genes_per_cell, cv.obs_n_genes_per_mol, min_val=0.99, m_quantile=1.0, n_bins=30)
    plt2 = Plots.barhist(cv.est_n_genes_per_cell[df_spatial.cell_dapi[real_mol_mask]], bins=t_bins, legend=:none, label="Expected", 
        xlabel="Number of unique genes", guidefontsize=13, tickfontsize=10)
    Plots.barhist!(cv.obs_n_genes_per_mol[real_mol_mask], bins=t_bins, alpha=0.7, label="Observed")
#     Plots.barhist!(cv.obs_n_genes_per_mol[.!real_mol_mask], bins=t_bins, alpha=0.5, label="Observed, not assigned")

    y_max = max(Plots.ylims(plt1)[2], Plots.ylims(plt1)[2])
    Plots.ylims!(plt1, 0, y_max)
    Plots.ylims!(plt2, 0, y_max)
    Plots.annotate!(plt1, 0.1, 0.9 * y_max, "K=$k")

    push!(plts, plt1)
    push!(plts, plt2)
end

plt = Plots.plot(plts..., layout=(3,2), size=(800, 550))
Plots.savefig(plt, "./plots/entropy_starmap.pdf")
Plots.closeall()
plt
```

```julia
# p_df = vcat([vcat(
#     DataFrame(
#         :Entropy => cv.est_ent_per_cell[df_spatial.cell_dapi[real_mol_mask]],
#         :NGenes => cv.est_n_genes_per_cell[df_spatial.cell_dapi[real_mol_mask]],
#         :Type => "Expected", :K => cv.k
#     ),
#     DataFrame(
#         :Entropy => cv.obs_ents_per_mol[real_mol_mask],
#         :NGenes => cv.obs_n_genes_per_mol[real_mol_mask],
#         :Type => "Observed, assigned", :K => cv.k
#     ),
#     DataFrame(
#         :Entropy => cv.obs_ents_per_mol[.!real_mol_mask],
#         :NGenes => cv.obs_n_genes_per_mol[.!real_mol_mask],
#         :Type => "Observed, not assigned", :K => cv.k
#     )
# ) for cv in comp_per_k]...);

# RCall.ijulia_setdevice(MIME("image/svg+xml"), width=8, height=4)

# R"""
# ggplot($p_df) +
#     geom_density_ridges(aes(x=Entropy, y=as.factor(K), height=stat(density), fill=Type), stat='binline', bins=50, alpha=0.5, scale=0.95) +
#     theme(legend.pos=c(1, 0), legend.justification=c(1, 0)) + ylab('')
# """
```

```julia
# @time neighb_cm = B.neighborhood_count_matrix(df_spatial, 30);
# @time color_transformation = B.gene_composition_transformation(neighb_cm);
# @time gene_colors = B.gene_composition_colors(neighb_cm, color_transformation; color_range=400);
# df_spatial[!, :color] = gene_colors;
```

```julia execution={"iopub.status.busy": "2020-09-22T10:06:51.680000+02:00", "iopub.execute_input": "2020-09-22T10:06:51.680000+02:00", "iopub.status.idle": "2020-09-22T10:07:02.729000+02:00"}
# polygons_dapi = B.boundary_polygons(p_df, p_df.cell_dapi, grid_step=1.0, bandwidth=4.0, min_border_length=75);
grid_step = 1
@time polygons_dapi = B.extract_polygons_from_label_grid(Matrix(dapi_seg_labels[1:grid_step:end, 1:grid_step:end]'), grid_step=grid_step);
```

```julia execution={"iopub.status.busy": "2020-09-22T11:58:52.107000+02:00", "iopub.execute_input": "2020-09-22T11:58:52.107000+02:00", "iopub.status.idle": "2020-09-22T11:58:54.333000+02:00"}
p_df = @where(df_spatial, :x .>= 2620, :x .<= 2720, :y .>= 950, :y .<= 1130)

t_genes = sortperm(B.prob_array(p_df.gene), rev=true)[1:15]
p_df = @where(p_df, in.(:gene, Ref(t_genes)));

ps = (1.75 * diff(collect(B.val_range(p_df.x)))[1], diff(collect(B.val_range(p_df.y)))[1]) .* 3.0
plt = B.plot_cell_borders_polygons(p_df, polygons_dapi, annotation=gene_names[p_df.gene], ticks=false,
    size=ps, ms=7.0, xlabel="x", ylabel="y", legend=:outerright, legendfontsize=14, guidefontsize=20)

Plots.savefig(plt, "./plots/starmap_grouping.pdf")
Plots.closeall()
plt
```

```julia
# # p_df = @where(df_spatial, :x .< 3000, :x .> 1000, :y .< 2000, :y .> 0)
# # p_df = @where(df_spatial, :x .< 2700, :x .> 2550, :y .< 1000, :y .> 900)
# # p_df = @where(df_spatial, :x .< 2650, :x .> 2550, :y .< 1000, :y .> 900)
# p_df = @where(df_spatial, :x .< 2620, :x .> 2535, :y .< 1056, :y .> 895)

# grid_step, bandwidth = 1.0, 5.0
# polygons_dapi = B.boundary_polygons(p_df, p_df.cell_dapi, grid_step=grid_step, bandwidth=bandwidth);
# polygons = B.boundary_polygons(p_df, p_df.cell, grid_step=grid_step, bandwidth=bandwidth);

# size_mult = 3.0
# ms = 5.0
# ps = (diff(collect(B.val_range(p_df.x)))[1], diff(collect(B.val_range(p_df.y)))[1]) .* size_mult
# plt1 = B.plot_cell_borders_polygons(p_df, polygons_dapi, color=:gene, size=ps, ms=ms, alpha=0.3,
#     ticks=false, xlabel="X", ylabel="Y")

# t_genes = sortperm(B.prob_array(p_df.gene), rev=true)[1:7]
# p_df2 = @where(p_df, in.(:gene, Ref(t_genes)));
# plt2 = B.plot_cell_borders_polygons(p_df2, polygons_dapi, annotation=gene_names[p_df2.gene], size=ps, ms=ms, bg_legend=Plots.RGBA(1,1,1,0.9), 
#     ticks=false, xlabel="X", ylabel="Y")
# Plots.xlims!(plt2, Plots.xlims(plt1))
# Plots.ylims!(plt2, Plots.ylims(plt1))

# # tw1 = 0.38
# # plt = Plots.plot(plt1, plt2, size=(ps[1] * 2.5, ps[2]), layout=Plots.grid(1, 2, widths=[tw1, 1 - tw1]))
# plt = Plots.plot(plt1, plt2, size=(ps[1] * 2, ps[2]), layout=(1, 2))
# Plots.savefig(plt, "./plots/starmap_grouping.pdf")
# Plots.closeall()
# plt
```

```julia
t_mol = df_spatial[findall(comp_per_k[3].obs_n_genes_per_mol .== 1)[200],:];
p_df = @where(df_spatial, :x .<= t_mol.x .+ 100, :x .>= t_mol.x .- 100, :y .<= t_mol.y .+ 100, :y .>= t_mol.y .- 100)

grid_step, bandwidth = 1.0, 5.0
polygons_dapi = B.boundary_polygons(p_df, p_df.cell_dapi, grid_step=grid_step, bandwidth=bandwidth);
polygons = B.boundary_polygons(p_df, p_df.cell, grid_step=grid_step, bandwidth=bandwidth);

size_mult = 1.5
ms = 5.0
ps = (diff(collect(B.val_range(p_df.x)))[1], diff(collect(B.val_range(p_df.y)))[1]) .* size_mult
plt1 = B.plot_cell_borders_polygons(p_df, polygons_dapi, color=:gene, size=ps, ms=ms, alpha=0.3)

t_genes = sortperm(B.prob_array(p_df.gene), rev=true)[1:4]
p_df2 = @where(p_df, in.(:gene, Ref(t_genes)));
plt2 = B.plot_cell_borders_polygons(p_df2, polygons_dapi, annotation=gene_names[p_df2.gene], size=ps, ms=ms, bg_legend=Plots.RGBA(1,1,1,0.9), legend=:bottomright)
Plots.xlims!(plt2, Plots.xlims(plt1))
Plots.ylims!(plt2, Plots.ylims(plt1))

Plots.plot(plt1, plt2, size=(ps[1] * 2, ps[2]), layout=(1, 2))
```

```julia
plt1 = B.plot_cell_borders_polygons(p_df, polygons, color=:gene, size=ps, ms=ms, alpha=0.3)
plt2 = B.plot_cell_borders_polygons(p_df2, polygons, annotation=gene_names[p_df2.gene], size=ps, ms=ms, bg_legend=Plots.RGBA(1,1,1,0.9), legend=:bottomright)
Plots.xlims!(plt2, Plots.xlims(plt1))
Plots.ylims!(plt2, Plots.ylims(plt1))

Plots.plot(plt1, plt2, size=(ps[1] * 2, ps[2]), layout=(1, 2))
```

<!-- #region toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true -->
## Other STARmap
<!-- #endregion -->

<!-- #region toc-hr-collapsed=true toc-nb-collapsed=true -->
### mpfc_160_20180424_BZ14_control
<!-- #endregion -->

```julia
@time df_spatial, gene_names = B.load_df("/home/vpetukhov/data/spatal/star_map/mpfc_160_20180424_BZ14_control/molecules.csv");
@time dapi_seg_labels = B.load_segmentation_mask("/home/vpetukhov/data/spatal/star_map/mpfc_160_20180424_BZ14_control/segmentation.tiff");
df_spatial[!, :cell_dapi] = denserank(B.staining_value_per_transcript(df_spatial, dapi_seg_labels)) .- 1;

dapi_seg_labels = nothing;
GC.gc()
```

```julia
# p_df = @where(df_spatial, :x .<= 3000, :x .>= 1000, :y .<= 2000, :y .>= 0)
p_df = @where(df_spatial, :x .<= 1250, :x .>= 1000, :y .<= 1930, :y .>= 1650)

@time polygons_dapi = B.boundary_polygons(p_df, p_df.cell_dapi, grid_step=1.0, bandwidth=5.0);

size_mult = 3.0
ms = 5.0
ps = (diff(collect(B.val_range(p_df.x)))[1], diff(collect(B.val_range(p_df.y)))[1]) .* size_mult
plt1 = B.plot_cell_borders_polygons(p_df, polygons_dapi, color=:gene, size=ps, ms=ms, alpha=0.3)

t_genes = sortperm(B.prob_array(p_df.gene), rev=true)[1:7]
p_df2 = @where(p_df, in.(:gene, Ref(t_genes)));
plt2 = B.plot_cell_borders_polygons(p_df2, polygons_dapi, annotation=gene_names[p_df2.gene], size=ps, ms=ms, bg_legend=Plots.RGBA(1,1,1,0.9))
Plots.xlims!(plt2, Plots.xlims(plt1))
Plots.ylims!(plt2, Plots.ylims(plt1))

Plots.plot(plt1, plt2, size=(ps[1] * 2, ps[2]), layout=(1, 2))
```

### visual_160_20170915_BD2_dark

```julia
@time df_spatial, gene_names = B.load_df("/home/vpetukhov/data/spatal/star_map/visual_160_20170915_BD2_dark/molecules.csv");
df_spatial.x ./= 2;
df_spatial.y ./= 2;

@time dapi_seg_labels = B.load_segmentation_mask("/home/vpetukhov/data/spatal/star_map/visual_160_20170915_BD2_dark/segmentation.tiff");
df_spatial[!, :cell_dapi] = denserank(B.staining_value_per_transcript(df_spatial, dapi_seg_labels)) .- 1;

dapi_seg_labels = nothing;
GC.gc()
```

```julia
# p_df = @where(df_spatial, :x .<= 3000, :x .>= 1000, :y .<= 2000, :y .>= 0)
# # p_df = @where(df_spatial, :x .<= 1250, :x .>= 1000, :y .<= 1930, :y .>= 1650)

# @time polygons_dapi = B.boundary_polygons(p_df, p_df.cell_dapi, grid_step=1.0, bandwidth=5.0);

# size_mult = 3.0
# ms = 5.0
# ps = (diff(collect(B.val_range(p_df.x)))[1], diff(collect(B.val_range(p_df.y)))[1]) .* size_mult
# plt1 = B.plot_cell_borders_polygons(p_df, polygons_dapi, color=:gene, size=ps, ms=ms, alpha=0.3)

# t_genes = sortperm(B.prob_array(p_df.gene), rev=true)[1:7]
# p_df2 = @where(p_df, in.(:gene, Ref(t_genes)));
# plt2 = B.plot_cell_borders_polygons(p_df2, polygons_dapi, annotation=gene_names[p_df2.gene], size=ps, ms=ms, bg_legend=Plots.RGBA(1,1,1,0.9))
# Plots.xlims!(plt2, Plots.xlims(plt1))
# Plots.ylims!(plt2, Plots.ylims(plt1))

# Plots.plot(plt1, plt2, size=(ps[1] * 2, ps[2]), layout=(1, 2))
```

### visual_1020_20180410-BY3_1kgenes

```julia
@time df_spatial, gene_names = B.load_df("/home/vpetukhov/data/spatal/star_map/visual_1020_20180410-BY3_1kgenes/molecules.csv");
@time dapi_seg_labels = B.load_segmentation_mask("/home/vpetukhov/data/spatal/star_map/visual_1020_20180410-BY3_1kgenes/segmentation.tiff");
df_spatial[!, :cell_dapi] = denserank(B.staining_value_per_transcript(df_spatial, dapi_seg_labels)) .- 1;

dapi_seg_labels = nothing;
GC.gc()
```

```julia
comp_per_k = @showprogress map(k -> T.estimate_local_complexity(df_spatial, k; cell_col=:cell_dapi), k_vals);
```

```julia
real_mol_mask = (df_spatial.cell_dapi .> 0);
plts = Plots.Plot[]

for (k, cv) in zip(k_vals, comp_per_k)
    t_bins = B.hist_bins(cv.est_ent_per_cell, cv.obs_ents_per_mol, min_val=-0.05, m_quantile=1.0, n_bins=30)
    plt1 = Plots.barhist(cv.est_ent_per_cell[df_spatial.cell_dapi[real_mol_mask]], bins=t_bins, legend=(k == 5 ? :top : :none), label="Expected",
        xlabel="Entropy", ylabel="Num. of transcripts", bg_legend=Plots.RGBA(1, 1, 1, 0.1))
    Plots.barhist!(cv.obs_ents_per_mol[real_mol_mask], bins=t_bins, alpha=0.7, label="Observed")
#     Plots.barhist!(cv.obs_ents_per_mol[.!real_mol_mask], bins=t_bins, alpha=0.5, label="Observed, not assigned")
    
    t_bins = B.hist_bins(cv.est_n_genes_per_cell, cv.obs_n_genes_per_mol, min_val=0.99, m_quantile=1.0, n_bins=30)
    plt2 = Plots.barhist(cv.est_n_genes_per_cell[df_spatial.cell_dapi[real_mol_mask]], bins=t_bins, legend=:none, label="Expected", 
        xlabel="Number of unique genes")
    Plots.barhist!(cv.obs_n_genes_per_mol[real_mol_mask], bins=t_bins, alpha=0.7, label="Observed")
#     Plots.barhist!(cv.obs_n_genes_per_mol[.!real_mol_mask], bins=t_bins, alpha=0.5, label="Observed, not assigned")

    y_max = max(Plots.ylims(plt1)[2], Plots.ylims(plt1)[2])
    Plots.ylims!(plt1, 0, y_max)
    Plots.ylims!(plt2, 0, y_max)
    Plots.annotate!(plt1, 0.1, 0.9 * y_max, "K=$k")

    push!(plts, plt1)
    push!(plts, plt2)
end

Plots.plot(plts..., layout=(3,2), size=(800, 550))
```

```julia
p_df = @where(df_spatial, :x .<= 3000, :x .>= 2500, :y .<= 1750, :y .>= 1000)
# p_df = @where(df_spatial, :x .<= 1250, :x .>= 1000, :y .<= 1930, :y .>= 1650)

@time polygons_dapi = B.boundary_polygons(p_df, p_df.cell_dapi, grid_step=1.0, bandwidth=5.0);

size_mult = 3.0
ms = 5.0
ps = (diff(collect(B.val_range(p_df.x)))[1], diff(collect(B.val_range(p_df.y)))[1]) .* size_mult
plt1 = B.plot_cell_borders_polygons(p_df, polygons_dapi, color=:gene, size=ps, ms=ms, alpha=0.3)

t_genes = sortperm(B.prob_array(p_df.gene), rev=true)[1:7]
p_df2 = @where(p_df, in.(:gene, Ref(t_genes)));
plt2 = B.plot_cell_borders_polygons(p_df2, polygons_dapi, annotation=gene_names[p_df2.gene], size=ps, ms=ms, bg_legend=Plots.RGBA(1,1,1,0.9))
Plots.xlims!(plt2, Plots.xlims(plt1))
Plots.ylims!(plt2, Plots.ylims(plt1))

Plots.plot(plt1, plt2, size=(ps[1] * 2, ps[2]), layout=(1, 2))
```

<!-- #region toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true -->
## MERFISH
<!-- #endregion -->

```julia
@time df_spatial_merfish = B.load_df("/home/vpetukhov/spatial/Benchmarking/run_results/merfish_moffit/segmentation.csv")[1];
```

```julia
comp_per_k_merfish = @showprogress map(k -> T.estimate_local_complexity(df_spatial_merfish, k), k_vals);
```

```julia
real_mol_mask = (df_spatial_merfish.cell .> 0);
plts = Plots.Plot[]

for (k, cv) in zip(k_vals, comp_per_k_merfish)
    t_bins = B.hist_bins(cv.est_ent_per_cell, cv.obs_ents_per_mol, min_val=-0.05, m_quantile=1.0, n_bins=30)
    plt1 = Plots.barhist(cv.est_ent_per_cell[df_spatial_merfish.cell[real_mol_mask]], bins=t_bins, legend=(k == 5 ? :topleft : :none), label="Expected",
        xlabel="Entropy", ylabel="Num. of transcripts", bg_legend=Plots.RGBA(1, 1, 1, 0.1))
    Plots.barhist!(cv.obs_ents_per_mol[real_mol_mask], bins=t_bins, alpha=0.7, label="Observed")
#     Plots.barhist!(cv.obs_ents_per_mol[.!real_mol_mask], bins=t_bins, alpha=0.5, label="Observed, not assigned")
    
    t_bins = B.hist_bins(cv.est_n_genes_per_cell, cv.obs_n_genes_per_mol, min_val=0.99, m_quantile=1.0, n_bins=30)
    plt2 = Plots.barhist(cv.est_n_genes_per_cell[df_spatial_merfish.cell[real_mol_mask]], bins=t_bins, legend=:none, label="Expected", 
        xlabel="Number of unique genes")
    Plots.barhist!(cv.obs_n_genes_per_mol[real_mol_mask], bins=t_bins, alpha=0.7, label="Observed")
#     Plots.barhist!(cv.obs_n_genes_per_mol[.!real_mol_mask], bins=t_bins, alpha=0.5, label="Observed, not assigned")

    y_max = max(Plots.ylims(plt1)[2], Plots.ylims(plt1)[2])
    Plots.ylims!(plt1, 0, y_max)
    Plots.ylims!(plt2, 0, y_max)
    Plots.annotate!(plt1, 0.1, 0.9 * y_max, "K=$k")

    push!(plts, plt1)
    push!(plts, plt2)
end

plt = Plots.plot(plts..., layout=(3,2), size=(800, 550))
Plots.savefig(plt, "./plots/entropy_merfish.pdf")
plt
```

<!-- #region toc-hr-collapsed=true toc-nb-collapsed=true -->
## Allen smFISH
<!-- #endregion -->

```julia execution={"iopub.status.busy": "2020-09-22T10:09:54.866000+02:00", "iopub.execute_input": "2020-09-22T10:09:54.866000+02:00", "iopub.status.idle": "2020-09-22T10:10:12.014000+02:00"}
@time df_spatial_allen, gene_names = B.load_df("/home/vpetukhov/spatial/Benchmarking/run_results/spacejam2/allen_sm_fish/allen_prior/segmentation.csv");

@time dapi_seg_labels = B.load_segmentation_mask("/home/vpetukhov/data/spatal/SpaceJam2Full/allen_sm_fish/segmentation_labels_from_json_transposed.tiff");
df_spatial_allen[!, :cell_dapi] = denserank(B.staining_value_per_transcript(df_spatial_allen, dapi_seg_labels)) .- 1;

dapi_seg_labels = nothing;
GC.gc()
```

```julia execution={"iopub.status.busy": "2020-09-21T15:07:26.649000+02:00", "iopub.execute_input": "2020-09-21T15:07:26.649000+02:00", "iopub.status.idle": "2020-09-21T15:07:46.414000+02:00"}
comp_per_k_allen = @showprogress map(k -> T.estimate_local_complexity(df_spatial_allen, k, cell_col=:cell_dapi), k_vals);
```

```julia
real_mol_mask = (df_spatial_allen.cell_dapi .> 0);
plts = Plots.Plot[]

for (k, cv) in zip(k_vals, comp_per_k_allen)
    t_bins = B.hist_bins(cv.est_ent_per_cell, cv.obs_ents_per_mol, min_val=-0.05, m_quantile=1.0, n_bins=30)
    plt1 = Plots.barhist(cv.est_ent_per_cell[df_spatial_allen.cell_dapi[real_mol_mask]], bins=t_bins, legend=(k == 5 ? :topleft : :none), label="Expected",
        xlabel="Entropy", ylabel="Num. of transcripts")
    Plots.barhist!(cv.obs_ents_per_mol[real_mol_mask], bins=t_bins, alpha=0.7, label="Observed")
#     Plots.barhist!(cv.obs_ents_per_mol[.!real_mol_mask], bins=t_bins, alpha=0.5, label="Observed, not assigned")
    
    t_bins = B.hist_bins(cv.est_n_genes_per_cell, cv.obs_n_genes_per_mol, min_val=0.99, m_quantile=1.0, n_bins=30)
    plt2 = Plots.barhist(cv.est_n_genes_per_cell[df_spatial_allen.cell_dapi[real_mol_mask]], bins=t_bins, legend=:none, label="Expected", 
        xlabel="Number of unique genes")
    Plots.barhist!(cv.obs_n_genes_per_mol[real_mol_mask], bins=t_bins, alpha=0.7, label="Observed")
#     Plots.barhist!(cv.obs_n_genes_per_mol[.!real_mol_mask], bins=t_bins, alpha=0.5, label="Observed, not assigned")

    y_max = max(Plots.ylims(plt1)[2], Plots.ylims(plt1)[2])
    Plots.ylims!(plt1, 0, y_max)
    Plots.ylims!(plt2, 0, y_max)
    
#     y_ticks = round.(range(0, y_max, length=5), sigdigits=2)
    y_ticks = range(0, round(y_max, sigdigits=1), length=5)
    Plots.yticks!(plt1, y_ticks)
    Plots.yticks!(plt2, y_ticks)
    
    Plots.annotate!(plt1, 0.1, 0.9 * y_max, "K=$k")

    push!(plts, plt1)
    push!(plts, plt2)
end

plt = Plots.plot(plts..., layout=(3,2), size=(800, 550))
Plots.savefig(plt, "./plots/entropy_allen_smfish.pdf")
plt
```

```julia execution={"iopub.status.busy": "2020-09-21T15:07:46.414000+02:00", "iopub.execute_input": "2020-09-21T15:07:46.414000+02:00", "iopub.status.idle": "2020-09-21T15:09:08.577000+02:00"}
@time neighb_cm = B.neighborhood_count_matrix(df_spatial_allen, 30);
@time color_transformation = B.gene_composition_transformation(neighb_cm);
@time gene_colors = B.gene_composition_colors(neighb_cm, color_transformation);
df_spatial_allen[!, :color] = gene_colors;
```

```julia execution={"iopub.status.busy": "2020-09-22T11:54:52.207000+02:00", "iopub.execute_input": "2020-09-22T11:54:52.207000+02:00", "iopub.status.idle": "2020-09-22T11:54:52.722000+02:00"}
import PyPlot

t_genes = ["Sv2c", "Satb2", "Pvalb", "Gad2"];
t_colors = Dict(Pair.(t_genes, ["#601A4A", "#EE442F", "#63ACBE", "#CCBE9F"]))

p_df = @where(df_spatial_allen, :x .< 14650, :x .> 14480, :y .> 11200, :y .< 11400) |> deepcopy
p_df = @where(p_df, in.(gene_names[:gene], Ref(t_genes)))
p_df.z .*= 10

scales = [diff(collect(B.val_range(p_df[!, s])))[1] for s in [:x, :y, :z]]

fig, axes = PyPlot.subplots(2, 1, figsize=(scales[1], sum(scales[2:3])) ./ 30, gridspec_kw=Dict(:height_ratios=>(scales[2] * 0.8, scales[3])))

for g in unique(p_df.gene)
    c_df = @where(p_df, :gene .== g)
    gene = gene_names[g]
    axes[1].scatter(c_df.x, c_df.y, label=gene, color=t_colors[gene])
    axes[2].scatter(c_df.x, c_df.z, label=gene, color=t_colors[gene])
end


axes[1].legend(ncol=4, loc="bottom right", bbox_to_anchor=(1.03, 1.3), handletextpad=0.0, handlelength=1.0, 
    columnspacing=0.5, fontsize=20, markerscale=1.5)
axes[1].set_ylim(11220, 11360)

cur_xticks = axes[1].get_xticks()
cur_xticks = cur_xticks[1:3:length(cur_xticks)]

label_font_size = 24

for ax in axes
    cur_yticks = ax.get_yticks()
    cur_yticks = cur_yticks[1:3:length(cur_yticks)]
    ax.grid(true, alpha=0.5)
    ax.set_xticks(cur_xticks)
    ax.set_yticks(cur_yticks)
    ax.set_xlabel("x", fontsize=label_font_size)
    
    ax.set_xticklabels(Int.(cur_xticks), fontsize=14)
    ax.set_yticklabels(Int.(cur_yticks), fontsize=14)
#     ax.set_xticklabels([])
#     ax.set_yticklabels([])
end

axes[1].set_ylabel("y", fontsize=label_font_size)
axes[2].set_ylabel("z", fontsize=label_font_size)

PyPlot.tight_layout()
PyPlot.savefig("./plots/allen_smfish_3d.pdf")
```

```julia execution={"iopub.status.busy": "2020-09-21T15:20:17.013000+02:00", "iopub.execute_input": "2020-09-21T15:20:17.013000+02:00", "iopub.status.idle": "2020-09-21T15:20:24.305000+02:00"}
Plots.pyplot()
```

```julia execution={"iopub.status.busy": "2020-09-21T17:41:32.730000+02:00", "iopub.execute_input": "2020-09-21T17:41:32.731000+02:00", "iopub.status.idle": "2020-09-21T17:41:34.843000+02:00"}
p_df = @where(df_spatial_allen, :x .< 14650, :x .> 14450, :y .> 11200, :y .< 11400) |> deepcopy
t_genes = ["Sv2c", "Satb2", "Pvalb", "Gad2"];
p_df = @where(p_df, in.(gene_names[:gene], Ref(t_genes)))
p_df.z .*= 10
t_colors = Dict(Pair.(t_genes, ["#601A4A", "#EE442F", "#63ACBE", "#CCBE9F"]))

guidefontsize = 16
size_mult = 2.0
ms = 5.0
x_lims = (14480, 14650)

s1 = (diff(collect(B.val_range(p_df.x)))[1], 80) .* size_mult
p1 = B.plot_cell_borders_polygons(DataFrame(:x => p_df.x, :y => p_df.z), color=p_df.gene, legend=:none,
    ann_colors=t_colors, annotation=gene_names[p_df.gene], ms=ms, size=s1, xlabel="X", ylabel="Z", alpha=1.0,
    xticks=x_lims[1]:50:x_lims[2], xlims=x_lims, guidefontsize=guidefontsize)

s2 = (diff(collect(B.val_range(p_df.x)))[1], diff(collect(B.val_range(p_df.y)))[1]) .* size_mult
p2 = B.plot_cell_borders_polygons(p_df, legend=:topleft, annotation=gene_names[p_df.gene], 
    xticks=x_lims[1]:50:x_lims[2], xlims=x_lims, legendfontsize=12, guidefontsize=guidefontsize,
    ann_colors=t_colors, bg_legend=Plots.RGBA(1,1,1,0.1), ms=ms, size=s2, xlabel="X", ylabel="y", alpha=1.0)

plt = Plots.plot(p2, p1, size=(s2[1], s2[2] + s1[2]), layout=Plots.grid(2, 1, heights=[s2[2], s1[2]] ./ (s2[2] + s1[2])))
# Plots.savefig(plt, "./plots/allen_smfish_3d.pdf")
Plots.closeall()

plt
```

```julia
# import PlotlyJS
# PLY = PlotlyJS;

# PLY.plot(
#     PLY.scatter3d(x=p_df.x, y=p_df.y, z=p_df.z, text=gene_names[p_df.gene], 
# #         marker_color=p_df.color, 
#         marker_color=B.distinguishable_colors(p_df.gene)[:colors], 
#         marker_size=4.0, mode="markers"),
#     PLY.Layout(width=1000, height=1000)
# )
```

```julia
# p_df = @where(df_spatial_allen, :x .< 15000, :x .> 10000, :y .< 15000, :y .> 10000)
# B.plot_cell_borders_polygons(p_df, color=:color)
```
