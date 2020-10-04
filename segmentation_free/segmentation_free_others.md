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
using Statistics
using StatsBase

B = Baysor;
```

<!-- #region toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true -->
## MERFISH Hippocampus
<!-- #endregion -->

```julia
@time df_spatial, gene_names = Baysor.load_df("../run_results/merfish_moffit/segmentation.csv");
```

```julia
@time neighb_cm = Baysor.neighborhood_count_matrix(df_spatial, 70);
@time color_transformation = Baysor.gene_composition_transformation(neighb_cm, df_spatial.confidence; sample_size=20000, spread=2.0, min_dist=0.1);
```

```julia
@time gene_colors = Baysor.gene_composition_colors(neighb_cm, color_transformation; color_range=1000);
```

```julia
@time plt = Baysor.plot_dataset_colors(df_spatial, gene_colors, min_molecules_per_cell=50, min_pixels_per_cell=15, ms=3.0)[1]
```

### Noise estimation

```julia
using Colors
```

```julia
@time mean_dists, confidences, d1, d2 = Baysor.append_confidence!(df_spatial, nn_id=50);
max_dist = quantile(mean_dists, 0.99);
mean_dists[mean_dists .> max_dist] .= max_dist;
```

```julia
border_dist = mean_dists[findmin(abs.(confidences .- 0.5))[2]];

Baysor.plot_noise_estimation_diagnostics(mean_dists, confidences, d1, d2; confidence_nn_id=50, linewidth=2.0, bins=50)
Plots.vline!([border_dist], label="", color="black")
```

```julia
n1 = round(Int, (border_dist - minimum(mean_dists)) / (max_dist - minimum(mean_dists)) * 100);
offset = 10

# t_palette = vcat(reverse(Colors.sequential_palette(10, n1 + offset)[offset:end]), Colors.sequential_palette(250, 100 + offset - n1)[offset:end]);
t_palette = reverse(Colors.sequential_palette(10, 100)[50:end]);
dist_colors = Baysor.map_to_colors(mean_dists, palette=t_palette);
Baysor.plot_colorbar(dist_colors)
```

```julia
df_spatial[!, :color] = dist_colors[:colors];
df_spatial[!, :gene_color] = gene_colors;
```

```julia
cur_df = @where(df_spatial, :x .< -3250, :y .> -2700);
size(cur_df)
```

```julia
@time polygons = Baysor.boundary_polygons(cur_df, (cur_df.confidence .>= 0.25) .+ 1, grid_step=0.5, method=:kde,
    bandwidth=3.0, exclude_labels=[1], dens_threshold=1e-10);

# length(polygons)
Baysor.plot_cell_borders_polygons(cur_df, polygons, color=cur_df.color, ms=1.0, alpha=0.4)
```

```julia
Baysor.plot_cell_borders_polygons(cur_df, polygons, color=cur_df.gene_color, ms=0.5, alpha=0.4)
```

### Clustering

```julia
@time adjacent_points, adjacent_weights = Baysor.build_molecule_graph(cur_df, filter=false);
```

```julia
mol_clust8_t = @async Baysor.cluster_molecules_on_mrf(cur_df.gene, adjacent_points, adjacent_weights, cur_df.confidence; n_clusters=8, max_iters=200, 
    n_iters_without_update=100, min_mols_per_cell=50);
```

```julia
mol_clust8_t.result
```

```julia
@time mol_clust8 = Baysor.cluster_molecules_on_mrf(cur_df.gene, adjacent_points, adjacent_weights, cur_df.confidence; n_clusters=8, max_iters=2000, 
    n_iters_without_update=100, min_mols_per_cell=50);

Baysor.plot_cell_borders_polygons(cur_df, annotation=mol_clust8.assignment, ms=1.0)
```

<!-- #region toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true -->
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
@time df_spatial, gene_names = Baysor.load_df("../run_results/iss_1_s10/segmentation.csv");
```

```julia
# @time dapi = Images.load("/home/vpetukhov/spatial/JamboreeV2/iss_mouse_3_dapi.tif");
@time dapi = Images.load("/home/vpetukhov/spatial/JamboreeV2/iss_mouse_1_HCA_11_3_m_133_c1-1.tif");
dapi_arr = Float64.(dapi);
```

```julia
df_spatial[!, :x] = round.(Int, df_spatial.x / 0.1625);
df_spatial[!, :y] = round.(Int, df_spatial.y / 0.1625);
```

```julia
maximum(df_spatial.y), minimum(df_spatial.y)
```

```julia
size(dapi)
```

```julia
grid_step = 10.0
@time polygons = Baysor.boundary_polygons(df_spatial, df_spatial.cell, grid_step=grid_step);
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

# @time gene_colors = Baysor.gene_composition_colors(neighb_cm, color_transformation; confidences=df_spatial.confidence, color_range=500);
```

```julia
@time plt = Baysor.plot_dataset_colors(df_spatial, gene_colors, polygons=polygons, min_molecules_per_cell=5,
    min_pixels_per_cell=30, polygon_alpha=0.5, ms=3.0, alpha=0.1)[1];

@time Plots.savefig(plt, "./plots/iss_mouse_3/polygons.png");
close()
```

```julia
df_spatial[!, :color] = gene_colors;
```

```julia
xs, xe = 12500, 15500;
ys, ye = 14000, 17000;

df_subs = @where(df_spatial, :x .> xs, :x .< xe, :y .> ys, :y .< ye);
# @time pol_subs = Baysor.boundary_polygons(df_subs, df_subs.cell, grid_step=10.0);
@time pol_subs_kde = Baysor.boundary_polygons(df_subs, df_subs.cell, grid_step=10.0, min_molecules_per_cell=10, method=:kde, verbose=true, bandwidth=5.0);

ysa, yea = Int(minimum(df_subs.y)), Int(maximum(df_subs.y));
xsa, xea = Int(minimum(df_subs.x)), Int(maximum(df_subs.x));
```

```julia
size(dapi)
```

```julia
maximum(df_spatial.y)
```

```julia
@time begin
    plt = Plots.heatmap(1.0 .- dapi_arr[ysa:yea, xsa:xea], size=(xea-xsa, yea-ysa), color=:grayscale)
    plt = Baysor.plot_cell_borders_polygons!(df_subs, pol_subs_kde, color=df_subs.color, ms=5.0, alpha=0.2, offset=(-xsa, -ysa),
        polygon_line_width=3, polygon_alpha=0.75, is_noise=df_subs.cell .== 0, noise_kwargs=Dict(:markersize => 3.0))

    Plots.savefig(plt, "./plots/iss_mouse_3/dapi_subs.png");
    close()
end
```

## Allen sm-FISH

```julia
using Clustering
using Distances
using SparseArrays
using ProgressMeter
using MatrixMarket

import Gadfly
GF = Gadfly;
```

```julia
@time df_spatial, gene_names = Baysor.load_df("../run_results/spacejam2/allen_sm_fish/no_dapi/segmentation.csv");
df_spatial[!, :x] = round.(Int, 10 .* (df_spatial.x .- minimum(df_spatial.x)));
df_spatial[!, :y] = round.(Int, 10 .* (df_spatial.y .- minimum(df_spatial.y)));
length(gene_names)
```

```julia
@time neighb_cm = Baysor.neighborhood_count_matrix(df_spatial, 40);
@time color_transformation = Baysor.gene_composition_transformation(neighb_cm, df_spatial.confidence; sample_size=20000, spread=2.0, min_dist=0.1);
```

```julia
@time gene_colors = Baysor.gene_composition_colors(neighb_cm, color_transformation; color_range=300);
```

```julia
df_spatial[!, :color] = gene_colors;
```

```julia
# @time Baysor.append_confidence!(df_spatial, nn_id=16);
# @time adjacent_points, adjacent_weights = Baysor.build_molecule_graph(df_spatial, filter=false);
```

### Coloring examples

```julia
sub_lims = [((15000, 17000), (7000, 9000)), ((11000, 14000), (12000, 15000)), ((8000, 9000), (12000, 13000)), ((11900, 12400), (9400, 10000)), ((12100, 12270), (9700, 9900))];
```

```julia
plt = Baysor.plot_cell_borders_polygons(df_spatial, color=gene_colors, ms=1.0, xlims=(5000, 18000), ylims=(2500, 15500), size=(1500, 1500), xlabel="X", ylabel="Y", ticks=false)
for ((xs, xe), (ys, ye)) in sub_lims
    Plots.plot!([xs, xs, xe, xe, xs], [ys, ye, ye, ys, ys], color="black", label="", lw=3.0, alpha=0.75)
end

Plots.savefig("./plots/allen_sm_fish/gene_coloring.png")
plt
```

```julia
for (i,((xs, xe), (ys, ye))) in enumerate(sub_lims[1:(end-1)])
    plt = Baysor.plot_cell_borders_polygons(@where(df_spatial, :x .>= xs, :x .<= xe, :y .>= ys, :y .<= ye), color=:color, ms=2.5, 
        size=(xe - xs, ye - ys) ./ 3, ticks=false) # , xlabel="X", ylabel="Y"
    Plots.savefig("./plots/allen_sm_fish/gene_coloring_e$(i).png")
    display(plt)
end
```

### Method scheme

```julia
sub_lims = [((15000, 17000), (7000, 9000)), ((11000, 14000), (12000, 15000)), ((8000, 9000), (12000, 13000)), ((11900, 12400), (9400, 10000)), 
    ((12150, 12250), (9720, 9900))];

(xs,xe), (ys,ye) = sub_lims[end]

p_df = @where(df_spatial, :x .>= xs, :x .<= xe, :y .>= ys, :y .<= ye);
p_df = @where(p_df, (B.count_array(:gene)[:gene] .> 16) .& (B.count_array(:gene)[:gene] .< 100), .!in.(gene_names[:gene], Ref(["Mpped1"])))
x_ticks = (xs:50:xe, "")
y_ticks = (ys:50:ye, "")

p_size = (xe - xs, ye - ys) .* 2
plt1 = B.plot_cell_borders_polygons(p_df, annotation=gene_names[p_df.gene], ms=5, size=p_size, legend=:bottomright, xticks=x_ticks, yticks=y_ticks)
plt2 = B.plot_cell_borders_polygons(p_df, color=:color, ms=5, size=p_size, xticks=x_ticks, yticks=y_ticks)

# B.plot_cell_borders_polygons(p_df, annotation=gene_names[p_df.gene], ms=5, size=(500, 1000), legend=:bottomright, xticks=x_ticks, yticks=y_ticks)
# Plots.annotate!(p_df.x, p_df.y, 1:size(p_df, 1))

t_pos_data = B.position_data(p_df);
t_id = 103
t_pd = t_pos_data[:, knn(KDTree(t_pos_data), t_pos_data[:,[t_id]], 16, true)[1][1][2:end]];
for pl in [plt1, plt2]
    for col in eachcol(t_pd)
        Plots.plot!(pl, [col[1], t_pos_data[1,t_id]], [col[2], t_pos_data[2,t_id]], label="", color="black", lw=1.5)
    end
end

Plots.savefig(plt1, "./plots/allen_sm_fish/knns_genes.pdf")
Plots.savefig(plt2, "./plots/allen_sm_fish/knns_expr.pdf")

Plots.plot(plt1, plt2, size=(p_size[1] * 2, p_size[2]))
```

<!-- #region toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true -->
### Local Expression Vectors
<!-- #endregion -->

```julia
module T

import Baysor
using DataFrames

function build_coexpression_mtx(expr_mat::Matrix{Int}, g1::Int, g2::Int, max_expr::Int=maximum(expr_mat[[g1, g2],:]))
    coexpr_mat = zeros(Int, max_expr, max_expr);
    for i in 1:size(expr_mat, 2)
        coexpr_mat[expr_mat[g1, i] + 1, expr_mat[g2, i] + 1] += 1
    end

    coexpr_mat[1, 1] = 0;
    return coexpr_mat ./ sum(coexpr_mat)
end

end
```

```julia
t_nm, t_ids = T.select_filtered_local_vectors(df_spatial, adjacent_points, 40; confidence=nothing);

all(t_ids .== 1:size(df_spatial, 1))
all(abs.(t_nm .- neighb_cm) .< 1e-5)
```

```julia
g1, g2 = [findfirst(gene_names .== g) for g in ["Gad2", "Satb2"]];
color_lims = (0., 5.);
```

```julia
@time mol_clust4 = Baysor.cluster_molecules_on_mrf(df_spatial.gene, adjacent_points, adjacent_weights, df_spatial.confidence; 
    n_clusters=4, max_iters=2000, n_iters_without_update=100);
```

```julia
# Baysor.plot_cell_borders_polygons(df_spatial, annotation=mol_clust4.assignment, ms=1.0, xlims=(1800, 23300), ylims=(1500, 21000), size=(700, 700))
```

```julia
t_nm, t_ids = Baysor.extract_filtered_local_vectors(df_spatial, adjacent_points, 40; confidence=nothing, normalize=false);
t_nm2, t_ids2 = Baysor.extract_filtered_local_vectors(df_spatial, adjacent_points, 40; normalize=false);
t_nm3, t_ids3 = Baysor.extract_filtered_local_vectors(df_spatial, adjacent_points, clust_per_mol=mol_clust4.assignment, 40, normalize=false);
```

```julia
plot_df = vcat([DataFrame(Dict(:expr => m[g2, m[g1,:] .> m[g2,:]], :type => t)) for (t, m) in zip([:t1, :t2, :t3], (t_nm, t_nm2, t_nm3))]...);
t_probs = [Baysor.prob_array(m[g2, m[g1,:] .> m[g2,:]] .+ 1) for m in (t_nm, t_nm2, t_nm3)]

plt = Plots.bar((1:length(t_probs[3])) .- 0.25, t_probs[3], bar_width=0.2, xticks=(1:10, ["$(x / 40)" for x in 0:9]), 
    label="Clustering", xlabel="Expression of a wrong marker", ylabel="Fraction of molecules", size=(600, 300))
Plots.bar!(t_probs[2], bar_width=0.2, label="Noise filtration")
Plots.bar!((1:length(t_probs[1])) .+ 0.25, t_probs[1], bar_width=0.2, label="Basic")

Plots.savefig("./plots/allen_sm_fish/patches_wrong_frac.pdf")
plt
```

```julia
plt = Plots.heatmap(log10.(T.build_coexpression_mtx(t_nm3, g1, g2, 32) .+ offset)[:,1:25], clims=color_lims, size=(500, 200), colorbar_title="log10(molecule frac)")
Plots.savefig("./plots/allen_sm_fish/patches_heatmap_legend.pdf")
plt
```

```julia
offset = 1e-5
color_lims = (-5., 0.)
t_plots = []
for (tit,m) in zip(["Noise molecules", "Real molecules, no adjustment", "Real molecules, adjustment by confidence", "Real molecules, adjustment by cluster"], 
        (t_nm[:, df_spatial.confidence .< 0.05], t_nm[:, df_spatial.confidence .> 0.75], t_nm2, t_nm3))
    plt = Plots.heatmap(log10.(T.build_coexpression_mtx(m, g1, g2, 32) .+ offset)[:,1:25], clims=color_lims, legend=:none, 
        xlims=(0.5, 25.5), ylims=(0.5, 32.5), xlabel="Satb2 expression fraction", ylabel="Gad2 expression fraction", title=tit,
        xticks=((0:5:25),  ["$(f/40)" for f in 0:5:25]), yticks=((0:5:30),  ["$(f/40)" for f in 0:5:30]))
    plt = Plots.plot!(range(0, 18.75, length=100), range(25, 0, length=100), color="#6be4ff", lw=2.0, alpha=0.75)
    plt = Plots.plot!(range(0, 15, length=100), range(20, 0, length=100), color="#6be4ff", lw=2.0, alpha=0.75)
    push!(t_plots, plt)
end

plt = Plots.plot(t_plots..., layout=(2, 2), size=(1000, 1000))
Plots.savefig("./plots/allen_sm_fish/patches_heatmaps.pdf")

plt
```

<!-- #region toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true -->
#### Attempts to show border effect
<!-- #endregion -->

```julia
mix_frac_per_vec = [min.(m[g1, :], m[g2, :]) ./ (m[g1, :] .+ m[g2, :]) for m in (t_nm, t_nm2, t_nm3)];
```

```julia
plot_ids = findall(.!isnan.(mix_frac_per_vec[1]));
# plt = Baysor.plot_cell_borders_polygons(df_spatial[t_ids,:], color=gene_colors[t_ids], ms=1.0, size=(1500, 1500), xlims=(1800, 23300), ylims=(1500, 21000), xlabel="X", ylabel="Y")
# plt = Baysor.plot_cell_borders_polygons(df_spatial, color=gene_colors, ms=1.0, xlims=(5000, 18000), ylims=(2500, 15500), size=(1500, 1500), xlabel="X", ylabel="Y")
```

```julia
# t_alphas = mix_frac_per_vec[1][plot_ids] .* 2;
# plt = Baysor.plot_cell_borders_polygons(df_spatial[plot_ids,:], color=gene_colors[plot_ids], alpha=max.(t_alphas, 0.05),
# #     ms=1.0, size=(1500, 1500), xlims=(1800, 23300), ylims=(1500, 21000), xlabel="X", ylabel="Y")
#     ms=1.0, xlims=(5000, 18000), ylims=(2500, 15500), size=(1500, 1500), xlabel="X", ylabel="Y")
```

```julia
# expr_ids = findall(.!isnan.(mix_frac_per_vec[3]))
# plot_ids = t_ids3[expr_ids];
# plot_alphas = mix_frac_per_vec[3][expr_ids] .* 2;
# plot_ann = (t_nm3[g1, :] .> t_nm3[g2, :])[expr_ids]
# # plt = Baysor.plot_cell_borders_polygons(df_spatial[plot_ids,:], color=gene_colors[plot_ids], alpha=max.(plot_alphas, 0.05),
# plt = Baysor.plot_cell_borders_polygons(df_spatial[plot_ids,:], annotation=Int.(plot_ann), alpha=max.(plot_alphas, 0.05),
#     ms=1.0, xlims=(5000, 18000), ylims=(2500, 15500), size=(1500, 1500), xlabel="X", ylabel="Y")
```

```julia
# plt = Baysor.plot_cell_borders_polygons(df_spatial[plot_ids,:], annotation=mol_clust4.assignment[plot_ids], alpha=max.(plot_alphas, 0.1),
#     ms=1.0, xlims=(5000, 18000), ylims=(2500, 15500), size=(1500, 1500), xlabel="X", ylabel="Y", noise_ann="other", noise_kwargs=Dict(:ms => 0.1, :alpha => 0.01))
```

```julia
# c_mask = in.(df_spatial.gene, Ref([g1, g2]))
# plt = Baysor.plot_cell_borders_polygons(df_spatial[c_mask,:], annotation=gene_names[df_spatial.gene][c_mask],
#     ms=1.0, xlims=(5000, 18000), ylims=(2500, 15500), size=(1500, 1500), xlabel="X", ylabel="Y")
```

```julia
# c_gene_annot = ifelse.(in.(df_spatial.gene, Ref([g1, g2])), gene_names[df_spatial.gene], "other");
# plt = Baysor.plot_cell_borders_polygons(df_spatial[plot_ids,:], annotation=c_gene_annot[plot_ids], alpha=max.(plot_alphas, 0.1),
#     ms=1.0, xlims=(5000, 18000), ylims=(2500, 15500), size=(1500, 1500), xlabel="X", ylabel="Y", noise_ann="other", noise_kwargs=Dict(:ms => 0.1, :alpha => 0.01))
```

```julia
# c_df = deepcopy(df_spatial[plot_ids,:])
# c_df[!, :cluster] = mol_clust4.assignment[plot_ids];
# c_df[!, :mix_fracs] = mix_frac_per_vec[3][expr_ids];
# c_df[!, :color] = gene_colors[plot_ids];
# c_df[!, :gene_filt] = c_gene_annot[plot_ids];

# # cur_df = @where(c_df, :x .< 13800, :x .> 13400, :y .< 4450, :y .> 4100);
# cur_df = @where(c_df, :x .< 13490, :x .> 13400, :y .< 4275, :y .> 4200);
# Baysor.plot_cell_borders_polygons(cur_df, 
#     annotation=cur_df.gene_filt, 
# #     annotation=cur_df.cluster, 
# #     annotation=gene_names[cur_df.gene], 
#     alpha=max.(cur_df.mix_fracs, 0.05),
#     ms=15.0, size=(1000, 1000), xlabel="X", ylabel="Y")
# #     ms=5.0, size=(1000, 1000), xlabel="X", ylabel="Y")

# # Plots.annotate!(collect(zip(cur_df.x, cur_df.y, Plots.text.(["$(round(f, digits=3))" for f in cur_df.mix_fracs], 5))))
# Plots.annotate!(collect(zip(cur_df.x, cur_df.y, Plots.text.(["$(round(f, digits=3))" for f in cur_df.mix_fracs], 10))))
```

```julia
# expr_ids = findall(.!isnan.(mix_frac_per_vec[2]))
# plot_ids = t_ids2[expr_ids];
# plot_alphas = mix_frac_per_vec[2][expr_ids] .* 2;
# plt = Baysor.plot_cell_borders_polygons(df_spatial[plot_ids,:], color=gene_colors[plot_ids], alpha=max.(plot_alphas, 0.05),
#     ms=1.0, xlims=(5000, 18000), ylims=(2500, 15500), size=(1500, 1500), xlabel="X", ylabel="Y")
```

```julia
# plt = Baysor.plot_cell_borders_polygons(df_spatial[plot_ids,:], annotation=mol_clust4.assignment[plot_ids], alpha=max.(plot_alphas, 0.05),
#     ms=1.0, xlims=(5000, 18000), ylims=(2500, 15500), size=(1500, 1500), xlabel="X", ylabel="Y")
```

```julia
# c_gene_annot = ifelse.(in.(df_spatial.gene, Ref([g1, g2])), gene_names[df_spatial.gene], "other");
# plt = Baysor.plot_cell_borders_polygons(df_spatial[plot_ids,:], annotation=c_gene_annot[plot_ids], alpha=max.(plot_alphas, 0.1),
#     ms=1.0, xlims=(5000, 18000), ylims=(2500, 15500), size=(1500, 1500), xlabel="X", ylabel="Y", noise_ann="other", noise_kwargs=Dict(:ms => 0.1, :alpha => 0.01))
```

<!-- #region toc-hr-collapsed=true toc-nb-collapsed=true -->
#### Different number of NNs
<!-- #endregion -->

```julia
sample_ids = Baysor.select_ids_uniformly(df_spatial.x, df_spatial.y, df_spatial.confidence, n=20000);

coord_df = deepcopy(df_spatial[sample_ids, 2:end]);
coord_df[!, :gene] = gene_names[coord_df.gene];
CSV.write("../cache/allen_coord_df.csv", coord_df);
p
@showprogress for nn in [5, 50, 300, 500, 2000]
    neighb_cm = Baysor.neighborhood_count_matrix(df_spatial, nn, normalize=false);
    mmwrite("../cache/allen_neighb_$nn.mm", SparseMatrixCSC(neighb_cm))
    mmwrite("../cache/allen_neighb_sample_$nn.mm", SparseMatrixCSC(neighb_cm[:, sample_ids]))
end

CSV.write("../cache/allen_gene_names.csv", DataFrame([gene_names], [:gene]));
```

<!-- #region toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true -->
#### Visualizations
<!-- #endregion -->

```julia
@time neighb_cm = Baysor.neighborhood_count_matrix(df_spatial, 5);

mol_ids = Baysor.select_ids_uniformly(df_spatial.x, df_spatial.y, df_spatial.confidence, 5000);

gene_dists = 1 .- cor(neighb_cm[:, mol_ids]');
patch_dists = 1 .- cor(neighb_cm[:, mol_ids]);

gene_ord = hclust(gene_dists, linkage=:ward).order;
patch_ord = hclust(patch_dists, linkage=:ward).order;

Plots.heatmap(neighb_cm[gene_ord, mol_ids[patch_ord]], yticks=(1:length(gene_names), gene_names[gene_ord]))
```

```julia
@time neighb_cm = Baysor.neighborhood_count_matrix(df_spatial, 10);

mol_ids = Baysor.select_ids_uniformly(df_spatial.x, df_spatial.y, df_spatial.confidence, 5000);

gene_dists = 1 .- cor(neighb_cm[:, mol_ids]');
patch_dists = 1 .- cor(neighb_cm[:, mol_ids]);

gene_ord = hclust(gene_dists, linkage=:ward).order;
patch_ord = hclust(patch_dists, linkage=:ward).order;

Plots.heatmap(neighb_cm[gene_ord, mol_ids[patch_ord]], yticks=(1:length(gene_names), gene_names[gene_ord]))
```

```julia
@time neighb_cm = Baysor.neighborhood_count_matrix(df_spatial, 40);

mol_ids = Baysor.select_ids_uniformly(df_spatial.x, df_spatial.y, df_spatial.confidence, 5000);

gene_dists = 1 .- cor(neighb_cm[:, mol_ids]');
patch_dists = 1 .- cor(neighb_cm[:, mol_ids]);

gene_ord = hclust(gene_dists, linkage=:ward).order;
patch_ord = hclust(patch_dists, linkage=:ward).order;

Plots.heatmap(neighb_cm[gene_ord, mol_ids[patch_ord]], yticks=(1:length(gene_names), gene_names[gene_ord]))
```

```julia
@time neighb_cm = Baysor.neighborhood_count_matrix(df_spatial, 100);

mol_ids = Baysor.select_ids_uniformly(df_spatial.x, df_spatial.y, df_spatial.confidence, 5000);

gene_dists = 1 .- cor(neighb_cm[:, mol_ids]');
patch_dists = 1 .- cor(neighb_cm[:, mol_ids]);

gene_ord = hclust(gene_dists, linkage=:ward).order;
patch_ord = hclust(patch_dists, linkage=:ward).order;

Plots.heatmap(neighb_cm[gene_ord, mol_ids[patch_ord]], yticks=(1:length(gene_names), gene_names[gene_ord]))
```

```julia
@time neighb_cm = Baysor.neighborhood_count_matrix(df_spatial, 200);

mol_ids = Baysor.select_ids_uniformly(df_spatial.x, df_spatial.y, df_spatial.confidence, 5000);

gene_dists = 1 .- cor(neighb_cm[:, mol_ids]');
patch_dists = 1 .- cor(neighb_cm[:, mol_ids]);

gene_ord = hclust(gene_dists, linkage=:ward).order;
patch_ord = hclust(patch_dists, linkage=:ward).order;

Plots.heatmap(neighb_cm[gene_ord, mol_ids[patch_ord]], yticks=(1:length(gene_names), gene_names[gene_ord]))
```

```julia
@time neighb_cm = Baysor.neighborhood_count_matrix(df_spatial, 500);

mol_ids = Baysor.select_ids_uniformly(df_spatial.x, df_spatial.y, df_spatial.confidence, 5000);

gene_dists = 1 .- cor(neighb_cm[:, mol_ids]');
patch_dists = 1 .- cor(neighb_cm[:, mol_ids]);

gene_ord = hclust(gene_dists, linkage=:ward).order;
patch_ord = hclust(patch_dists, linkage=:ward).order;

Plots.heatmap(neighb_cm[gene_ord, mol_ids[patch_ord]], yticks=(1:length(gene_names), gene_names[gene_ord]))
```

<!-- #region toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true -->
### Clustering
<!-- #endregion -->

```julia
confidence_nn_id = Baysor.default_param_value(:confidence_nn_id, 30);
Baysor.append_confidence!(df_spatial, nn_id=confidence_nn_id);
```

```julia
adjacent_points, adjacent_weights = Baysor.build_molecule_graph(df_spatial, filter=false);
```

```julia
@time polygons = Baysor.boundary_polygons(df_spatial, df_spatial.cell; grid_step=2.0, bandwidth=6.5);
```

<!-- #region toc-hr-collapsed=true toc-nb-collapsed=true -->
#### 6 clusters
<!-- #endregion -->

```julia
@time mol_clust6_3 = Baysor.cluster_molecules_on_mrf(df_spatial.gene, adjacent_points, adjacent_weights, df_spatial.confidence;
    n_clusters=6, min_mols_per_cell=30);
```

```julia
Baysor.plot_cell_borders_polygons(df_spatial, annotation=mol_clust6_3.assignment, ms=1.0)
```

```julia
@time mol_clust6_2 = Baysor.cluster_molecules_on_mrf(df_spatial, adjacent_points, adjacent_weights; 
    n_clusters=6, min_mols_per_cell=30);
```

```julia
Baysor.plot_cell_borders_polygons(df_spatial, annotation=mol_clust6_2.assignment, ms=1.0)
```

```julia
Baysor.plot_cell_borders_polygons(df_spatial, annotation=mol_clust6_2.assignment, ms=1.0)
```

```julia
@time mol_clust6 = Baysor.cluster_molecules_on_mrf(df_spatial, adjacent_points, adjacent_weights; 
    n_clusters=6, min_mols_per_cell=30);
```

```julia
@time mol_clust6 = Baysor.cluster_molecules_on_mrf(df_spatial.gene, adjacent_points, adjacent_weights, df_spatial.confidence; 
    n_clusters=6, max_iters=10000, n_iters_without_update=20, min_mols_per_cell=30);
```

```julia
Baysor.plot_cell_borders_polygons(df_spatial, annotation=mol_clust6.assignment, ms=1.0)
```

```julia
@time mol_clust6 = Baysor.cluster_molecules_on_mrf(df_spatial.gene, adjacent_points, adjacent_weights, df_spatial.confidence; 
    n_clusters=6, max_iters=10000, n_iters_without_update=20);
```

```julia
Baysor.plot_cell_borders_polygons(df_spatial, annotation=mol_clust6.assignment, ms=1.0)
```

<!-- #region toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true -->
#### 4 clusters
<!-- #endregion -->

```julia
@time mol_clust4 = Baysor.cluster_molecules_on_mrf(df_spatial.gene, adjacent_points, adjacent_weights; n_clusters=4, max_iters=10000, n_iters_without_update=100);
```

```julia
Baysor.plot_cell_borders_polygons(df_spatial, annotation=mol_clust4.assignment, ms=1.0)
```

<!-- #region toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true -->
#### 8 clusters
<!-- #endregion -->

```julia
@time mol_clust8 = Baysor.cluster_molecules_on_mrf(df_spatial.gene, adjacent_points, adjacent_weights, df_spatial.confidence; n_clusters=8, max_iters=2000, n_iters_without_update=100);
```

```julia
@time mol_clust8_v2 = Baysor.cluster_molecules_on_mrf(df_spatial.gene, adjacent_points, adjacent_weights, df_spatial.confidence; n_clusters=8, max_iters=2000, 
    n_iters_without_update=100, min_mols_per_cell=30);
```

```julia
Baysor.plot_cell_borders_polygons(df_spatial, annotation=mol_clust8.assignment, ms=1.0, xlims=(1800, 23300), ylims=(1500, 21000), size=(1000, 1000))
```

```julia
Baysor.plot_cell_borders_polygons(df_spatial, annotation=mol_clust8_v2.assignment, ms=1.0, xlims=(1800, 23300), ylims=(1500, 21000), size=(1000, 1000))
```

```julia
Baysor.plot_cell_borders_polygons(df_spatial, annotation=mol_clust8_v2.assignment, ms=1.0, xlims=(1800, 23300), ylims=(1500, 21000), size=(1000, 1000))
```

```julia
gene_ord = sortperm(gene_names, rev=true)
Plots.heatmap((mol_clust8_v2.exprs')[gene_ord, :], yticks=(1:length(gene_names), gene_names[gene_ord]), xticks=1:size(mol_clust8_v2.exprs, 1))
```

```julia
gene_ord = sortperm(gene_names, rev=true)
Plots.heatmap((mol_clust8_v2.exprs')[gene_ord, :], yticks=(1:length(gene_names), gene_names[gene_ord]), xticks=1:size(mol_clust8_v2.exprs, 1))
```

<!-- #region toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true -->
#### 10 clusters
<!-- #endregion -->

```julia
@time mol_clust10 = Baysor.cluster_molecules_on_mrf(df_spatial.gene, adjacent_points, adjacent_weights; n_clusters=10, max_iters=10000, n_iters_without_update=100);
```

```julia
Baysor.plot_cell_borders_polygons(df_spatial, annotation=mol_clust10.assignment, ms=1.0)
```

#### 20 clusters

```julia
@time mol_clust20 = Baysor.cluster_molecules_on_mrf(df_spatial.gene, adjacent_points, adjacent_weights, df_spatial.confidence; n_clusters=20, max_iters=10000, 
    n_iters_without_update=100, min_mols_per_cell=30);
```

```julia
clust_order = sortperm(vec(mol_clust20.exprs[:, gene_names .== "Gad2"] .> 0.01));
```

```julia
Baysor.plot_cell_borders_polygons(df_spatial, annotation=sortperm(clust_order)[mol_clust20.assignment], ms=1.0, xlims=(1800, 23300), ylims=(1500, 21000), size=(1000, 1000))
```

```julia
gene_ord = sortperm(gene_names, rev=true)
Plots.heatmap((mol_clust20.exprs')[gene_ord, clust_order], yticks=(1:length(gene_names), gene_names[gene_ord]), xticks=1:size(mol_clust20.exprs, 1))
```

```julia
@time mol_clust20_v2 = Baysor.cluster_molecules_on_mrf(df_spatial.gene, adjacent_points, adjacent_weights, df_spatial.confidence; n_clusters=20, max_iters=10000, 
    n_iters_without_update=100, min_mols_per_cell=30);
```

```julia
clust_order = sortperm(vec(mol_clust20_v2.exprs[:, gene_names .== "Gad2"] .> 0.01) .+ vec(mol_clust20_v2.exprs[:, gene_names .== "Satb2"]));
```

```julia
plt = Baysor.plot_cell_borders_polygons(df_spatial, annotation=sortperm(clust_order)[mol_clust20_v2.assignment], ms=1.0, 
    xlims=(5000, 18000), ylims=(2500, 15500), size=(1500, 1500), xlabel="X", ylabel="Y")
#     xlims=(1800, 23300), ylims=(1500, 20000), size=(1500, 1500), xlabel="X", ylabel="Y")
Plots.savefig("./plots/allen_sm_fish/molecule_clustering.png")
plt
```

```julia
gene_ord = sortperm(gene_names, rev=true)
plt = Plots.heatmap((mol_clust20_v2.exprs')[gene_ord, clust_order], yticks=(1:length(gene_names), gene_names[gene_ord]), xticks=1:size(mol_clust20_v2.exprs, 1), 
    xlabel="Cluster", size=(500, 400))
Plots.savefig("./plots/allen_sm_fish/cluster_centers.png")
plt
```

<!-- #region toc-hr-collapsed=true toc-nb-collapsed=true -->
### Noise estimation
<!-- #endregion -->

```julia
using Colors
```

```julia
@time mean_dists, confidences, d1, d2 = Baysor.append_confidence!(df_spatial, nn_id=16);
max_dist = quantile(mean_dists, 0.99);
mean_dists[mean_dists .> max_dist] .= max_dist;
```

```julia
border_dist = mean_dists[findmin(abs.(confidences .- 0.5))[2]];

plt = Baysor.plot_noise_estimation_diagnostics(mean_dists, confidences, d1, d2; confidence_nn_id=16, linewidth=2.0, bins=50, widen=false, size=(500, 300))
plt = Plots.vline!([border_dist], label="", color="black")

Plots.savefig("./plots/allen_sm_fish/noise_separation.pdf")
plt
```

```julia
# n1 = round(Int, (border_dist - minimum(mean_dists)) / (max_dist - minimum(mean_dists)) * 100);
# offset = 10
# t_palette = vcat(reverse(Colors.sequential_palette(10, n1 + offset)[offset:end]), Colors.sequential_palette(250, 100 + offset - n1)[offset:end]);

# t_palette = reverse(Colors.sequential_palette(10, 100)[50:end]);
# dist_colors = Baysor.map_to_colors(mean_dists, palette=t_palette);
dist_colors = B.map_to_colors(mean_dists, palette=Colors.diverging_palette(10, 250, s=0.75, w=1.0))
plt = Baysor.plot_colorbar(dist_colors, lw=0.0, xlims=(0, maximum(mean_dists)), xlabel="Distance to 16'th nearest neighbor", size=(500, 200))
Plots.savefig("./plots/allen_sm_fish/distance_colorbar.pdf")
plt
```

```julia
df_spatial[!, :color] = dist_colors[:colors];
df_spatial[!, :gene_color] = gene_colors;
```

```julia
cur_df = @where(df_spatial, :x .>= 10000, :x .<= 12500, :y .>= 12500, :y .<= 15000);

@time polygons = Baysor.boundary_polygons(cur_df, (cur_df.confidence .>= 0.25) .+ 1, grid_step=2.0, shape_method=:order, max_dev=30.0,
    bandwidth=8.0, exclude_labels=[1], dens_threshold=1e-8, min_border_length=30);

plt = Baysor.plot_cell_borders_polygons(cur_df, polygons, color=cur_df.color, ms=2.5, alpha=0.4, size=(500, 500), ticks=false)
Plots.savefig("./plots/allen_sm_fish/distance_per_transcript.png")
plt
```

```julia
plt = Baysor.plot_cell_borders_polygons(cur_df, polygons, color=cur_df.gene_color, ms=1.25, alpha=0.4)
Plots.savefig("./plots/allen_sm_fish/noise_borders.png")
plt
```

<!-- #region toc-hr-collapsed=true toc-nb-collapsed=true -->
### Annotation transfer
<!-- #endregion -->

```julia
Baysor.append_confidence!(df_spatial, nn_id=16);
```

```julia
@time adjacent_points, adjacent_weights = Baysor.build_molecule_graph(df_spatial, filter=false, adjacency_type=:both, k_adj=30);
```

```julia
centroid_df = CSV.read("../metadata/allen_sm_fish_visp_subclass_centroids.csv");
centroid_df = @where(centroid_df, :class .!= "Non-Neuronal", .!in.(:subclass, Ref(["CR", "NP", "Meis2", "Serpinf1"]))); # Sncg
centroid_exprs = Matrix(centroid_df[:, Symbol.(gene_names)]);

size(centroid_exprs)
```

```julia
@time type_transfer_init = Baysor.cluster_molecules_on_mrf(df_spatial.gene, adjacent_points, adjacent_weights, df_spatial.confidence; 
    do_maximize=false, max_iters=1000, cell_type_exprs=centroid_exprs, n_iters_without_update=20);
```

```julia
expr_obs = deepcopy(type_transfer_init.exprs)
Baysor.maximize_molecule_clusters!(expr_obs, expr_obs, df_spatial.gene, df_spatial.confidence, type_transfer_init.assignment_probs);
```

```julia
exprs_norm = copy((type_transfer_init.exprs ./ sum(type_transfer_init.exprs, dims=1))');
gene_ord = sortperm(vec(mapslices(x -> findmax(x)[2], exprs_norm, dims=2)), rev=true);
B.clustermap(exprs_norm, gene_names, 
    gene_ord=gene_ord, cell_ord=1:size(type_transfer_init.exprs, 1),
    xticks=(1:size(type_transfer_init.exprs, 1), centroid_df.subclass), color=:OrRd_9
)[1]
```

```julia
exprs_norm = copy((type_transfer_init.exprs ./ sum(type_transfer_init.exprs, dims=1))');
B.clustermap(exprs_norm, gene_names, 
    gene_ord=1:length(gene_names), cell_ord=1:size(type_transfer_init.exprs, 1),
    xticks=(1:size(type_transfer_init.exprs, 1), centroid_df.subclass), color=:OrRd_9
)[1]
```

```julia
B.clustermap((type_transfer_init.exprs ./ sum(type_transfer_init.exprs, dims=1))', gene_names, cell_ord=1:size(type_transfer_init.exprs, 1),
    xticks=(1:size(type_transfer_init.exprs, 1), centroid_df.subclass), color=:OrRd_9
)[1]
```

```julia
gene_ord = sortperm(gene_names, rev=true)
Plots.heatmap((type_transfer_init.exprs')[gene_ord, :], yticks=(1:length(gene_names), gene_names[gene_ord]), 
    xticks=(1:size(type_transfer_init.exprs, 1), centroid_df.subclass), size=(550, 400), color=:OrRd_9)
```

```julia
gene_ord = sortperm(gene_names, rev=true)
plt = Plots.heatmap((expr_obs')[gene_ord, :], yticks=(1:length(gene_names), gene_names[gene_ord]), 
    xticks=(1:size(expr_obs, 1), centroid_df.subclass), size=(550, 400), color=:OrRd_9)

Plots.savefig("./plots/allen_sm_fish/cell_type_centers.png")
plt
```

```julia
plt = Plots.bar(Baysor.count_array(type_transfer_init.assignment), label="Total", ylabel="Num. of molecules", size=(550, 200),
    xticks=(1:length(centroid_df.subclass), centroid_df.subclass), grid=:y)
Plots.bar!(Baysor.count_array(type_transfer_init.assignment[df_spatial.confidence .> 0.95]), label="Confidence > 0.95")

Plots.savefig("./plots/allen_sm_fish/cell_type_nums.pdf")
plt
```

```julia
plt = Baysor.plot_cell_borders_polygons(df_spatial, annotation=centroid_df.subclass[type_transfer_init.assignment], ms=1.0, 
#     xlims=(1800, 23300), ylims=(1500, 21000), size=(1000, 1000))
    xlims=(5000, 18000), ylims=(2500, 15500), size=(1500, 1500), xlabel="X", ylabel="Y", legend=:topright)

Plots.savefig("./plots/allen_sm_fish/molecule_annotation.png")
plt
```

```julia
@time type_transfer_init2 = Baysor.cluster_molecules_on_mrf(df_spatial.gene, adjacent_points, adjacent_weights, df_spatial.confidence; 
    do_maximize=true, max_iters=1000, cell_type_exprs=type_transfer_init.exprs, assignment_probs=type_transfer_init.assignment_probs, n_iters_without_update=20);
```

```julia
Baysor.plot_cell_borders_polygons(df_spatial, annotation=centroid_df.subclass[type_transfer_init2.assignment], ms=1.0, 
    xlims=(1800, 23300), ylims=(1500, 21000), size=(1000, 1000))
```

```julia

```

```julia
class_df = CSV.read("../metadata/allen_sm_fish_visp_class_centroids.csv");
class_exprs = Matrix(class_df[:, Symbol.(gene_names)]);
```
