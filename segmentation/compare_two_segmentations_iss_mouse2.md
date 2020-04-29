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

B = Baysor;
```

### Load data

```julia
@time df_spatial, gene_names = B.load_df("../run_results/spacejam2/iss/mouse_visp2_ISS_2_spot_table/segmentation.csv");
df_spatial[!, :x] = round.(Int, df_spatial.x ./ 0.1625);
df_spatial[!, :y] = round.(Int, df_spatial.y ./ 0.1625);

gene_names_all = gene_names;
```

```julia
@time paper_seg_labels = Float64.(Images.load("/home/vpetukhov/data/spatal/SpaceJam2Full/iss/mouse_visp2_HCA_11_2_m_133_c1-1_segmentation_oleh.tif")');
paper_seg_labels = round.(Int, paper_seg_labels' .* 2^16);

df_spatial[!, :cell_paper] = Baysor.staining_value_per_transcript(df_spatial, paper_seg_labels);
```

```julia
df_spatial[!, :cell_paper] = denserank(df_spatial.cell_paper) .- 1;
df_spatial.cell .= denserank(df_spatial.cell) .- 1;
```

#### Visualization

```julia
@time dapi_arr = Float64.(Images.load("/home/vpetukhov/data/spatal/SpaceJam2Full/iss/mouse_visp2_HCA_11_2_m_133_c1-1.tif"));
```

```julia
@time polygons = Baysor.boundary_polygons(df_spatial, df_spatial.cell, grid_step=10.0, bandwidth=10.0);
```

```julia
@time neighb_cm = Baysor.neighborhood_count_matrix(df_spatial, 20);
@time color_transformation = Baysor.gene_composition_transformation(neighb_cm, df_spatial.confidence);
@time gene_colors = Baysor.gene_composition_colors(neighb_cm, color_transformation; color_range=500);
```

```julia
df_spatial[!, :color] = gene_colors;
```

### QC

```julia
qc_per_cell = Baysor.get_cell_qc_df(df_spatial);
qc_per_cell_paper = Baysor.get_cell_qc_df(df_spatial, df_spatial.cell_paper);

for df in (qc_per_cell, qc_per_cell_paper)
    df[!, :cell_id] = 1:size(df, 1)
    deleterows!(df, findall(df.n_transcripts .< 5));
end

size(qc_per_cell, 1), size(qc_per_cell_paper, 1)
```

```julia
t_bins = 1:1:100
Plots.histogram(qc_per_cell.n_transcripts, label="Baysor", bins=t_bins, widen=false,
    xlabel="Num. of transcripts", ylabel="Num. of cells", xlims=Baysor.val_range(t_bins))
Plots.histogram!(qc_per_cell_paper.n_transcripts, label="DAPI", bins=t_bins, alpha=0.6)
```

```julia
t_bins = 0:0.001:0.1
Plots.histogram(qc_per_cell.density, label="Baysor", bins=t_bins,
    xlabel="Density", ylabel="Num. of cells", xlims=Baysor.val_range(t_bins), widen=false)
Plots.histogram!(qc_per_cell_paper.density, label="DAPI", bins=t_bins, alpha=0.6)
```

```julia
t_bins = 1:0.1:20.0
Plots.histogram(qc_per_cell.elongation, label="Baysor", bins=t_bins,
    xlabel="Elongation", ylabel="Num. of cells", xlims=Baysor.val_range(t_bins), widen=false)
Plots.histogram!(qc_per_cell_paper.elongation, label="DAPI", bins=t_bins, alpha=0.6)
```

```julia
t_bins = 1:3.0:100.0
Plots.histogram(sqrt.(qc_per_cell.area), label="Baysor", bins=t_bins,
    xlabel="sqrt(Area)", ylabel="Num. of cells", xlims=Baysor.val_range(t_bins), widen=false)
Plots.histogram!(sqrt.(qc_per_cell_paper.area), label="DAPI", bins=t_bins, alpha=0.6)
```

### Matching

```julia
assignment_filt = denserank(ifelse.(in.(df_spatial.cell, Ref(qc_per_cell.cell_id)), df_spatial.cell, 0)) .-1 ;
assignment_filt_paper = denserank(ifelse.(in.(df_spatial.cell_paper, Ref(qc_per_cell_paper.cell_id)), df_spatial.cell_paper, 0)) .- 1;

contingency_table = counts(assignment_filt .+ 1, assignment_filt_paper .+ 1);
```

```julia
round.((mean(assignment_filt .== 0), mean(assignment_filt_paper .== 0)), digits=3)
```

```julia
frac_of_noise1 = (contingency_table[1,:] ./ vec(sum(contingency_table, dims=1)))[2:end];
frac_of_noise2 = (contingency_table[:, 1] ./ vec(sum(contingency_table, dims=2)))[2:end];

t_bins = 0.0:0.03:1.0
Plots.histogram(frac_of_noise2, label="Baysor", bins=t_bins, xlims=Baysor.val_range(t_bins), widen=false,
    xlabel="Fraction of molecules, not assigned in other segmentation", ylabel="Num. of cells")
Plots.histogram!(frac_of_noise1, label="DAPI", bins=t_bins, alpha=0.6)
```

```julia
max_overlap1 = maximum(contingency_table[:, 2:end], dims=1)[:] ./ sum(contingency_table, dims=1)[2:end];
max_overlap2 = maximum(contingency_table[2:end, :], dims=2)[:] ./ sum(contingency_table, dims=2)[2:end];

t_bins = 0.0:0.03:1.0
Plots.histogram(max_overlap2, label="Baysor", bins=t_bins, xlims=Baysor.val_range(t_bins), widen=false,
    xlabel="Fraction of molecules, matching to a single cell", ylabel="Num. of cells")
Plots.histogram!(max_overlap1, label="DAPI", bins=t_bins, alpha=0.6)
```

```julia
bin_match = (contingency_table[2:end, 2:end] ./ sum(contingency_table[2:end, 2:end], dims=1) .> 0.05);

t_bins = 1:7
Plots.histogram(sum(bin_match, dims=2)[:], label="Baysor", bins=t_bins, xlims=Baysor.val_range(t_bins), widen=false,
    xlabel="Number of overlapping cells", ylabel="Num. of cells")
Plots.histogram!(sum(bin_match, dims=1)[:], label="DAPI", bins=t_bins, alpha=0.6)
```

```julia
ctn1 = contingency_table ./ sum(contingency_table, dims=1)
ctn2 = contingency_table ./ sum(contingency_table, dims=2)
match_thresholds = 0.02:0.05:1.0
match_nums = [sum(any((ctn1 .> mt) .& (ctn2 .> mt), dims=1)) for mt in match_thresholds];
Plots.plot(match_thresholds, match_nums, xlabel="Minimal fraction of matching molecules", ylabel="Num. of cells", 
    xlims=(match_thresholds[1], match_thresholds[end]), widen=false, legend=:none, lw=2.0)
```

### Overview of matching results

```julia
match_noise = frac_of_noise2 .> 0.5;
match_noise_paper = frac_of_noise1 .> 0.5;
multiple_overlap = (max_overlap2 .< 0.7) .& (frac_of_noise2 .< 0.25);
multiple_overlap_paper = (max_overlap1 .< 0.7) .& (frac_of_noise1 .< 0.25);

@assert !any(match_noise .& multiple_overlap)
```

```julia
stat_df = DataFrame(
    :Metric => ["Num. cells", "Total num. molecules", "Fraction of noise", "Fraction of matching cells", "Fraction of matching to noise", "Fraction of multiple overlap"],
    :Baysor => [size(qc_per_cell, 1), sum(qc_per_cell.n_transcripts), 1 .- sum(qc_per_cell.n_transcripts) / size(df_spatial, 1), 
        mean((.!match_noise) .& (.!multiple_overlap)), mean(match_noise), mean(multiple_overlap)],
    :DAPI => [size(qc_per_cell_paper, 1), sum(qc_per_cell_paper.n_transcripts), 1 .- sum(qc_per_cell_paper.n_transcripts) / size(df_spatial, 1), 
        mean((.!match_noise_paper) .& (.!multiple_overlap_paper)), mean(match_noise_paper), mean(multiple_overlap_paper)]
)

stat_df.Baysor .= round.(stat_df.Baysor, digits=3)
stat_df.DAPI .= round.(stat_df.DAPI, digits=3)
stat_df
```

```julia
t_bins = 4:45
lw = 2.0
Plots.stephist(qc_per_cell.n_transcripts[(.!match_noise) .& (.!multiple_overlap)], label="Baysor, full matching", bins=t_bins, widen=false,
    xlabel="Num. of transcripts", ylabel="Num. of cells", xlims=Baysor.val_range(t_bins), lw=lw)
Plots.stephist!(qc_per_cell.n_transcripts[match_noise], label="Baysor, match noise", bins=t_bins, alpha=0.75, lw=lw)
Plots.stephist!(qc_per_cell.n_transcripts[multiple_overlap], label="Baysor, multiple overlap", bins=t_bins, alpha=0.75, lw=lw)
Plots.stephist!(qc_per_cell_paper.n_transcripts[.!multiple_overlap_paper], label="DAPI, full matching", bins=t_bins, alpha=0.75, lw=lw)
Plots.stephist!(qc_per_cell_paper.n_transcripts[multiple_overlap_paper], label="DAPI, multiple overlap", bins=t_bins, alpha=0.75, lw=lw)
```

```julia
t_bins = 0:0.001:0.06
lw = 2.0
Plots.stephist(qc_per_cell.density[(.!match_noise) .& (.!multiple_overlap)], label="Baysor, full matching", bins=t_bins, widen=false,
    xlabel="Density", ylabel="Num. of cells", xlims=Baysor.val_range(t_bins), lw=lw)
Plots.stephist!(qc_per_cell.density[match_noise], label="Baysor, match noise", bins=t_bins, alpha=0.75, lw=lw)
Plots.stephist!(qc_per_cell.density[multiple_overlap], label="Baysor, partial matching", bins=t_bins, alpha=0.75, lw=lw)
Plots.stephist!(qc_per_cell_paper.density[.!multiple_overlap_paper], label="DAPI, full matching", bins=t_bins, alpha=0.75, lw=lw)
Plots.stephist!(qc_per_cell_paper.density[multiple_overlap_paper], label="DAPI, partial matching", bins=t_bins, alpha=0.75, lw=lw)
```

```julia
# t_bins = 0:0.001:0.1
# Plots.histogram(qc_per_cell.density, label="Baysor, matching", bins=t_bins,
#     xlabel="Density", ylabel="Num. of cells", xlims=Baysor.val_range(t_bins), widen=false)
# Plots.histogram!(qc_per_cell.density[match_noise], label="Baysor, not-matching", bins=t_bins, alpha=0.75)
# Plots.histogram!(qc_per_cell_paper.density, label="DAPI", bins=t_bins, alpha=0.6)
```

```julia
# t_bins = 1:0.1:20.0
# Plots.histogram(qc_per_cell.elongation, label="Baysor, matching", bins=t_bins,
#     xlabel="Elongation", ylabel="Num. of cells", xlims=Baysor.val_range(t_bins), widen=false)
# Plots.histogram!(qc_per_cell.elongation[match_noise], label="Baysor, not-matching", bins=t_bins, alpha=0.75)
# Plots.histogram!(qc_per_cell_paper.elongation, label="DAPI", bins=t_bins, alpha=0.6)
```

### Difference in noise-matching cells

```julia
tr = 0.05;
genes_filt = findall(vec((mean(cm .> 0, dims=2) .> tr) .| (mean(cm_paper .> 0, dims=2) .> tr)));
length(genes_filt)
```

```julia
cm = Baysor.convert_segmentation_to_counts(df_spatial.gene, df_spatial.cell)[:, qc_per_cell.cell_id];
cm_paper = Baysor.convert_segmentation_to_counts(df_spatial.gene, df_spatial.cell_paper)[:, qc_per_cell_paper.cell_id];

cm = cm ./ sum(cm, dims=1)
cm_paper = cm_paper ./ sum(cm_paper, dims=1);

cm = cm[genes_filt, :];
cm_paper = cm_paper[genes_filt, :];
gene_names = gene_names_all[genes_filt];
```

```julia
@time cell_ords, gene_ord = B.joint_ordering(cm[:, .!match_noise], cm[:, match_noise], cm_paper);

Plots.plot(
    Baysor.clustermap(cm[:, .!match_noise], gene_names; gene_ord=gene_ord, cell_ord=cell_ords[1], cbar=false, clims=(0,1), title="Baysor, not noise")[1],
    Baysor.clustermap(cm[:, match_noise], gene_names; gene_ord=gene_ord, cell_ord=cell_ords[2], cbar=false, clims=(0,1), title="Baysor, matching to noise")[1],
    Baysor.clustermap(cm_paper, gene_names; gene_ord=gene_ord, cell_ord=cell_ords[3], cbar=false, clims=(0,1), title="DAPI")[1],
    layout=(1,3), size=(1400, 800)
)
```

```julia
for ci in qc_per_cell.cell_id[match_noise][1:20:200]
    display(B.plot_comparison_for_cell(df_spatial, ci, paper_seg_labels, dapi_arr; ms=7.0, grid_step=2.0, bandwidth=4.0, size_mult=2.0, center_mult=2.0))
end
```

### Differences in partially-matching cells

```julia
cell_ords, gene_ord = B.joint_ordering(cm_paper[:, .!multiple_overlap_paper], cm_paper[:, multiple_overlap_paper], cm);

Plots.plot(
    Baysor.clustermap(cm_paper[:, .!multiple_overlap_paper], gene_names; gene_ord=gene_ord, cell_ord=cell_ords[1], cbar=false, clims=(0,1), title="DAPI, matching")[1],
    Baysor.clustermap(cm_paper[:, multiple_overlap_paper], gene_names; gene_ord=gene_ord, cell_ord=cell_ords[2], cbar=false, clims=(0,1), title="DAPI, multiple overlap")[1],
    Baysor.clustermap(cm, gene_names; gene_ord=gene_ord, cell_ord=cell_ords[3], cbar=false, clims=(0,1), title="Baysor")[1],
    layout=(1,3), size=(1400, 500)
)
```

```julia
for ci in qc_per_cell_paper.cell_id[multiple_overlap_paper][1:10:100]
    display(B.plot_comparison_for_cell(df_spatial, ci, paper_seg_labels, dapi_arr; cell1_col=:cell_paper, cell2_col=:cell, 
            ms=7.0, grid_step=2.0, bandwidth=4.0, size_mult=2.0))
end
```

```julia
cell_ords, gene_ord = B.joint_ordering(cm[:, .!multiple_overlap], cm[:, multiple_overlap], cm_paper);

Plots.plot(
    Baysor.clustermap(cm[:, .!multiple_overlap], gene_names; gene_ord=gene_ord, cell_ord=cell_ords[1], cbar=false, clims=(0,1), title="Baysor, matching")[1],
    Baysor.clustermap(cm[:, multiple_overlap], gene_names; gene_ord=gene_ord, cell_ord=cell_ords[2], cbar=false, clims=(0,1), title="Baysor, multiple overlap")[1],
    Baysor.clustermap(cm_paper, gene_names; gene_ord=gene_ord, cell_ord=cell_ords[3], cbar=false, clims=(0,1), title="DAPI")[1],
    layout=(1,3), size=(1400, 500)
)
```

```julia
for ci in qc_per_cell.cell_id[multiple_overlap][1:20:100]
    display(B.plot_comparison_for_cell(df_spatial, ci, paper_seg_labels, dapi_arr; ms=7.0, grid_step=2.0, bandwidth=4.0, size_mult=2.0))
end
```

### Correlation between parts

```julia
part_overlap_mask = (max_overlap2 .> 0.25) .& (max_overlap2 .< 0.75);
cur_overlaps = max_overlap2[part_overlap_mask]
partially_matching_ids = qc_per_cell.cell_id[part_overlap_mask];
noise_expr = Baysor.prob_array(df_spatial.gene[df_spatial.cell .== 0], max_value=length(gene_names_all));

match_cors = Float64[]
noise_cors = Float64[]

for t_cell_id in partially_matching_ids
    t_df = df_spatial[df_spatial.cell .== t_cell_id,:];
    t_mcp = mode(t_df.cell_paper[t_df.cell_paper .!= 0]);
    
    match_expr = Baysor.prob_array(t_df.gene[t_df.cell_paper .== t_mcp], max_value=length(gene_names_all));
    non_match_expr = Baysor.prob_array(t_df.gene[t_df.cell_paper .!= t_mcp], max_value=length(gene_names_all));
    
    push!(match_cors, cor(non_match_expr, match_expr))
    push!(noise_cors, cor(non_match_expr, noise_expr))
end

t_bins = -0.2:0.02:1.0
Plots.histogram(match_cors, bins=t_bins, widen=false, label="To matching cells", legend=:topleft, xlabel="Correlation", ylabel="Num. of cells")
Plots.histogram!(noise_cors, alpha=0.5, bins=t_bins, label="To noise component")
```

```julia
part_overlap_mask = (max_overlap1 .> 0.25) .& (max_overlap1 .< 0.75);
cur_overlaps_paper = max_overlap1[part_overlap_mask]
partially_matching_ids_paper = qc_per_cell_paper.cell_id[part_overlap_mask];
noise_expr = Baysor.prob_array(df_spatial.gene[df_spatial.cell_paper .== 0], max_value=length(gene_names_all));

match_cors_paper = Float64[]
noise_cors_paper = Float64[]

for t_cell_id in partially_matching_ids_paper
    t_df = df_spatial[df_spatial.cell_paper .== t_cell_id,:];
    t_mcp = mode(t_df.cell[t_df.cell .!= 0]);
    
    match_expr = Baysor.prob_array(t_df.gene[t_df.cell .== t_mcp], max_value=length(gene_names_all));
    non_match_expr = Baysor.prob_array(t_df.gene[t_df.cell .!= t_mcp], max_value=length(gene_names_all));
    
    push!(match_cors_paper, cor(non_match_expr, match_expr))
    push!(noise_cors_paper, cor(non_match_expr, noise_expr))
end

t_bins = -0.3:0.02:1.0
Plots.histogram(match_cors_paper, bins=t_bins, widen=false, label="To matching cells", legend=:topleft, xlabel="Correlation", ylabel="Num. of cells")
Plots.histogram!(noise_cors_paper, alpha=0.5, bins=t_bins, label="To noise component")
```

```julia
# # for (ci,mc,of) in collect(zip(partially_matching_ids, match_cors, cur_overlaps))[match_cors .> 0.75][1:10]
# for (ci,mc,of) in collect(zip(partially_matching_ids, match_cors, cur_overlaps))[match_cors .< 0.25][1:10]
#     display(B.plot_comparison_for_cell(df_spatial, ci, paper_seg_labels, dapi_arr; ms=4.0, bandwidth=10.0, 
#             title="C: $(round(mc, digits=2)), O: $(round(of, digits=2))"))
# end
```

```julia
# for (ci,mc,of) in collect(zip(partially_matching_ids_paper, match_cors_paper, cur_overlaps_paper))[match_cors_paper .< 0.25]
#     display(B.plot_comparison_for_cell(df_spatial, ci, paper_seg_labels, dapi_arr; cell1_col=:cell_paper, cell2_col=:cell, ms=4.0, bandwidth=10.0, 
#             title="C: $(round(mc, digits=2)), O: $(round(of, digits=2))"))
# end
```

### Examples

```julia
plts = []
for (xs, ys) in [((16000, 17000), (16000, 17000)), ((10000, 11000), (10000, 11000)), ((8500, 9500), (13500, 14500))]
    plt = B.plot_comparison_for_cell(df_spatial, xs, ys, paper_seg_labels, dapi_arr; ms=4.0, grid_step=2.0, bandwidth=4.0, size_mult=1.0, 
        alpha=0.6, polygons=false, grid_alpha=0.1, ticks=true, noise=false)
    display(plt)
    push!(plts, plt)
end;
```

```julia
# for i in 1:3
#     Plots.savefig(plts[i], "./plots/iss_mouse_2/example$i.png")
# end
```
