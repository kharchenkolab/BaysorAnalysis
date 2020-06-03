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
import Images
import Plots
import CSV
import MAT

using DataFrames
using DataFramesMeta
using ProgressMeter
using Statistics
using StatsBase

import Baysor
B = Baysor;

COORD_PATH = "/home/vpetukhov/data/spatal/seq_fish/NIH3T3/seqFISH+_NIH3T3_point_locations/";
```

```julia
gene_names = vec(String.(MAT.matread(COORD_PATH * "all_gene_Names.mat")["allNames"]));
@time rna_locations_run_1 = Array{Float64, 2}.(MAT.matread(COORD_PATH * "RNA_locations_run_1.mat")["tot"]);
# @time rna_locations_run_2 = Array{Float64, 2}.(MAT.matread(COORD_PATH * "RNA_locations_run_2.mat")["tot"]);

size(rna_locations_run_1), size(rna_locations_run_2)
```

```julia
module T

using DataFrames
using ProgressMeter

function convert_seq_fish_cell(cell::Matrix{Float64}, cell_id::Int, gene_name::String)
    df = DataFrame(cell, [:x, :y, :z])
    df[!, :cell] .= cell_id
    df[!, :gene] .= gene_name
    
    return df
end

function seq_fish_matrix_to_df(mat::Array{Array{Float64, 2}, 2}, gene_names::Array{String, 1})
    dfs = [convert_seq_fish_cell(mat[ci, gi], ci, gene_names[gi]) 
        for gi in 1:size(mat, 2) for ci in 1:size(mat, 1) if length(mat[ci, gi]) > 0];

    return vcat(dfs...);
end

function seq_fish_data_to_dfs(arr::Array{Matrix{Float64}, 3}, gene_names::Array{String, 1})
    dfs_spatial = @showprogress map(i -> seq_fish_matrix_to_df(arr[i,:,:], gene_names), 1:size(arr, 1))
    n_cells = 0
    for (i,df) in enumerate(dfs_spatial)
        df.cell .+= n_cells
        df[!, :fov] .= i
        n_cells = maximum(df.cell)
    end

    return dfs_spatial
end

end
```

```julia
dfs_spatial = T.seq_fish_data_to_dfs(rna_locations_run_1, gene_names);
```

```julia
B.plot_cell_borders_polygons(dfs_spatial[5], color=:cell, alpha=0.1, size=(600, 600))
```

## Investigate sub-cellular structure

```julia
df_spatial = deepcopy(dfs_spatial[5]);
gene_names = sort(unique(df_spatial.gene))
df_spatial[!, :gene] = denserank(df_spatial.gene);
size(df_spatial)
```

```julia
expressed_genes = findall(B.count_array(df_spatial.gene) .> 50)
gene_names = gene_names[expressed_genes]
df_spatial = @where(df_spatial, in.(:gene, Ref(Set(expressed_genes))));
df_spatial.gene .= denserank(df_spatial.gene);
size(df_spatial)
```

```julia
@time neighb_cm = B.neighborhood_count_matrix(df_spatial, 500);
@time color_transformation = B.gene_composition_transformation(neighb_cm, n_pcs=50);
@time gene_colors = B.gene_composition_colors(neighb_cm, color_transformation; color_range=200);
```

```julia
B.plot_cell_borders_polygons(df_spatial, color=gene_colors, size=(600, 600))
```

```julia
# adjacent_points, adjacent_weights = B.build_molecule_graph(df_spatial, filter=false);

# clust_res = B.cluster_molecules_on_mrf(df_spatial.gene, adjacent_points, adjacent_weights, ones(size(df_spatial, 1)), n_clusters=2);

# B.plot_cell_borders_polygons(df_spatial, annotation=clust_res.assignment, alpha=0.1)
```

## Run segmentation

```julia
df_spatial[!, :confidence] .= 1.0;
```

```julia
# @time bm_data_init = B.initial_distribution_arr(df_spatial; n_frames=1, scale=200, scale_std="25%", min_molecules_per_cell=1000, use_local_gene_similarities=false, 
@time bm_data = B.initial_distribution_arr(df_spatial; n_frames=1, scale=175, scale_std="50%", min_molecules_per_cell=1000, use_local_gene_similarities=false, 
    n_cells_init=50, confidence_nn_id=0)[1];
B.bmm!(bm_data, n_iters=500, new_component_frac=0.3, min_molecules_per_cell=1000, assignment_history_depth=100, verbose=true, log_step=50);
```

```julia
# B.plot_num_of_cells_per_iterarion(bm_data.tracer)
```

```julia
assignment = B.estimate_assignment_by_history(bm_data)[1];
assignment[B.count_array(assignment .+ 1)[assignment .+ 1] .< 500] .= 0;
B.plot_cell_borders_polygons(df_spatial, annotation=assignment, alpha=0.5, noise_ann=0)
```

```julia
polygons = B.boundary_polygons(df_spatial, assignment, grid_step=3.0, bandwidth=1.0, dens_threshold=1e-7, min_border_length=100);
polygons_paper = B.boundary_polygons(df_spatial, df_spatial.cell, grid_step=3.0, bandwidth=1.0, dens_threshold=1e-7, min_border_length=100);

poly_width = 2.0
plt1 = B.plot_cell_borders_polygons(df_spatial, polygons_paper, color=gene_colors, polygon_line_width=poly_width);
plt2 = B.plot_cell_borders_polygons(df_spatial, polygons, annotation=assignment, legend=:none, alpha=0.5, noise_ann=0, size=(400, 400), polygon_line_width=poly_width)
plt = Plots.plot(plt1, plt2, size=(800, 400), format=:png)
Plots.savefig(plt, "./plots/seq_fish_fibroblasts.png")
plt
```
