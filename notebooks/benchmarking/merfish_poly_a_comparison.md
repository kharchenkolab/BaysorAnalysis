---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.9.1
  kernelspec:
    display_name: Julia 1.6.0
    language: julia
    name: julia-1.6
---

```julia
using DrWatson
quickactivate(@__DIR__)

import Baysor as B
import BaysorAnalysis as BA
import Colors
import MultivariateStats
import Plots
import Clustering
import ColorSchemes
import CSV
import CairoMakie as MK
import JSON
import Images
import ImageMorphology
import PyPlot as Plt
import Seaborn as Sns

using DataFrames
using DataFramesMeta
using NearestNeighbors
using ProgressMeter
using OrderedCollections
using SparseArrays
using Statistics
using StatsBase
using StatsPlots
using VegaLite

ProgressMeter.ijulia_behavior(:clear);
MK.activate!(type = "png");
BA.set_pyplot_defaults!()

Plots.theme(:wong);
Plots.default(legendfontsize=12, guidefontsize=16, tickfontsize=12, gridlinewidth=1, fg_legend=false)

cell_cols_aliases = Dict(:cell => "Baysor", :cell_watershed => "Watershed", :cell_polya => "poly-A", :cell_pciseq => "pciSeq")
cell_cols = collect(keys(cell_cols_aliases))
color_per_label = BA.method_palette(collect(values(cell_cols_aliases)));

cplotsdir(args...) = plotsdir("benchmarking/merfish_polya", args...);
```

## Load data

```julia
@time data = BA.load_merfish(store_labels=true, watershed=true, min_mols_per_cell=40);
rename!(data[:df], :cell_paper => :cell_polya);

@time data[:polya_labels] = Matrix(B.load_segmentation_mask(datadir("exp_raw/merfish_moffit/moffit_poly_a_mask.tif")));
@time poly_a_arr = Float16.(Images.load(datadir("exp_raw/merfish_moffit/polyA_merged.tiff")));
data[:df][!, :polya_brightness] = B.staining_value_per_transcript(data[:df], poly_a_arr);
data[:watershed_labels] = Matrix(data[:watershed_labels]);
df_spatial = data[:df];
```

```julia
# data[:df] = data[:df][.!in.(data[:df].gene_name, Ref(Set(["Blank-1", "Blank-2", "Blank-3", "Blank-4", "Blank-5", "Krt90"]))),:]
# data[:df] = @orderby(data[:df], :x, :y, :gene_name);

# df_pciseq = DataFrame!(CSV.File(datadir("exp_pro/merfish_moffit/pci_seq_aggr/spots.csv")))[:,2:end];
# df_pciseq = @orderby(df_pciseq, :x, :y, :Gene);
# data[:df][!,:cell_pciseq] = df_pciseq.neighbour;

# df_spatial = data[:df];
```

```julia
@time BA.append_matching_statistics!(data, cell_cols, target_cell_col_names=[:cell_polya]);
```

```julia
@time baysor_labels = B.find_grid_point_labels_kde(B.position_data(df_spatial), data[:assignment_filt][:cell], 
    [0., 0.], Float64.(size(data[:polya_labels])[[2,1]]); grid_step=5.0, bandwidth=5.0, dens_threshold=1e-2);

@time baysor_labels = ImageMorphology.closing(baysor_labels);
@time baysor_labels = B.upscale(baysor_labels, 5)[1:size(data[:polya_labels], 1), 1:size(data[:polya_labels], 2)];
data[:baysor_labels] = baysor_labels;
```

```julia
@time pciseq_labels = B.find_grid_point_labels_kde(B.position_data(df_spatial), data[:assignment_filt][:cell_pciseq],
    [0., 0.], Float64.(size(data[:polya_labels])[[2,1]]); grid_step=5.0, bandwidth=5.0, dens_threshold=1e-2);

@time pciseq_labels = ImageMorphology.closing(pciseq_labels);
@time pciseq_labels = B.upscale(pciseq_labels, 5)[1:size(data[:polya_labels], 1), 1:size(data[:polya_labels], 2)];
data[:pciseq_labels] = pciseq_labels;
```


### Filter segments

```julia tags=[]
lab_keys = [:baysor_labels, :polya_labels, :watershed_labels, :pciseq_labels];
n_pix_per_comp = Dict(k => B.count_array(vec(data[k]), drop_zero=true) for k in lab_keys);

@showprogress for n in lab_keys
    filt_segs = Set(findall(n_pix_per_comp[n] .< 250))
    data[n][in.(data[n], Ref(filt_segs))] .= 0
    vec(data[n]) .= denserank(vec(data[n])) .- 1
    n_pix_per_comp[n] = B.count_array(vec(data[n]), drop_zero=true)
end
```

```julia tags=[]
bins = range(0, 250, length=100)
Plots.histogram(sqrt.(n_pix_per_comp[:baysor_labels]), bins=bins, label="Baysor", xlabel="Num. pixels", ylabel="Num. cells", color=color_per_label["Baysor"])
Plots.histogram!(sqrt.(n_pix_per_comp[:polya_labels]), bins=bins, label="poly-A", alpha=0.75, color=color_per_label["poly-A"])
Plots.histogram!(sqrt.(n_pix_per_comp[:watershed_labels]), bins=bins, label="Watershed", alpha=0.5, color=color_per_label["Watershed"])
Plots.histogram!(sqrt.(n_pix_per_comp[:pciseq_labels]), bins=bins, label="pciSeq", alpha=0.5, color=color_per_label["pciSeq"])
Plots.vline!([sqrt(250)], color="black", label="")
```

## Watershed vs Baysor comparison on transcripts

```julia
n_cells_poly_a = size(data[:qc_per_cell_dfs][:cell_polya], 1)
ovr_counts = [B.count_array(min.(data[s].n_overlaps[1], 3) .+ 1) for s in [:match_res_cell_polya, :match_res_pciseq_polya, :match_res_watershed_polya]];
p_labels = ["Baysor", "pciSeq", "Watershed"]
group_var = repeat(p_labels, inner=length(ovr_counts[1]));

plt = StatsPlots.groupedbar(hcat(ovr_counts...), palette=[color_per_label[s] for s in p_labels], group=group_var, 
    xlabel="Number of overlapping cells", ylabel="Num. of cells", ylim=(0, 6600), yticks=0:1500:6000, xgrid=false, lw=0.5, size=(500, 400))
Plots.hline!([n_cells_poly_a], label="Num. poly-A cells", lw=3, color="black", ylims=(0, n_cells_poly_a + 100))
Plots.xticks!(1:4, ["0", "1", "2", "3+"])

Plots.savefig(cplotsdir("num_overlaps.pdf"))
plt
```

```julia
Plt.figure(figsize=(2, 4))
Plt.bar(1:length(cell_cols), [size(data[:qc_per_cell_dfs][k], 1) for k in cell_cols], color=[color_per_label[cell_cols_aliases[c]] for c in cell_cols])
Plt.xticks([]); Plt.yticks(0:2500:10000);
Plt.xlabel("Method"); Plt.ylabel("Num. of cells");
Plt.gca().set_axisbelow(true)
Plt.grid(true, axis=:y, alpha=0.3)
Plt.tight_layout()

# Plots.xticks!(1:3, ["Baysor", "Watershed", "poly-A"])
Plt.savefig(cplotsdir("num_cells.pdf"))
```

```julia
fig, axes = Plt.subplots(1, 2, figsize=(8, 3))
p_keys = [:match_res_cell_polya => "Baysor", :match_res_watershed_polya => "Watershed", :match_res_pciseq_polya => "pciSeq"]
bins = 0.0:0.04:1.0
alpha = 0.1
lw = 2
for (k, l) in p_keys
    df = BA.precision_recall_df(data[k].contingency)
    col = color_per_label[l]
    axes[1].hist(df.precision, bins=bins, alpha=alpha, color=col);
    axes[1].hist(df.precision, bins=bins, label=l, histtype="step", color=col, lw=lw);
    
    axes[2].hist(df.recall, bins=bins, alpha=alpha, color=col);
    axes[2].hist(df.recall, bins=bins, label=l, histtype="step", color=col, lw=lw);
end

for ax in axes
    ax.grid(true, alpha=0.1)
    ax.set_ylabel("Num. of cells");
    ax.set_xlim(0, 1);
end

axes[1].set_xlabel("Precision")
axes[2].set_xlabel("Recall")
axes[2].legend(loc="upper left", borderpad=1, handlelength=1);
axes[1].legend(loc="upper left", borderpad=1, handlelength=1);
Plt.tight_layout();

Plt.savefig(cplotsdir("precision_recall.pdf"))
```

```julia
bins = 0.0:0.02:1.0
fig, axes = Plt.subplots(1, 2, figsize=(8, 3))
p_keys = [:frac_vals_cell_polya => "Baysor", :frac_vals_watershed_polya => "Watershed", :frac_vals_pciseq_polya => "pciSeq"]
for (k,n) in p_keys
    axes[1].hist(data[k].overlap, bins=bins, label=n, alpha=0.5, color=color_per_label[n], density=true)
end
for (k,v) in p_keys
    axes[1].hist(data[k].overlap, bins=bins, density=true, color="black", histtype="step")
end

bins = 0.0:0.04:2.0
for (k,n) in p_keys
    axes[2].hist(data[k].assigned, bins=bins, label=n, alpha=0.5, color=color_per_label[n], density=true)
end
for (k,v) in p_keys
    axes[2].hist(data[k].assigned, bins=bins, density=true, color="black", histtype="step")
end

for ax in axes
    ax.set_ylabel("Cell density")
    ax.legend(handlelength=1.0, title="Target");
    
    ymax = ax.get_ylim()[2]
    ax.vlines([0.9975], ymin=0, ymax=ymax, color="black", lw=2.0);
    ax.set_ylim(0, ymax)
end

axes[1].set_xlim(0, 1); axes[2].set_xlim(0, 2);
axes[1].set_xlabel("Overlap size / poly-A cell size");
axes[2].set_xlabel("Target cell size / poly-A cell size");
Plt.tight_layout()

Plt.savefig(cplotsdir("mol_fracs.pdf"))
```

### Poly-A brightness

```julia tags=[]
# bins = 0.0:0.016:1.0
bins = 0.0:0.05:1.0
p_cols = [:cell_watershed, :cell_polya, :cell_pciseq, :cell]
Plt.figure(figsize=(4, 4))
for s in p_cols
    t = cell_cols_aliases[s]
    Plt.hist(df_spatial.polya_brightness[df_spatial[!, s] .== 0], bins=bins, label=t, color=color_per_label[t], alpha=0.9, edgecolor="black", lw=0.2)
end

for s in p_cols
    Plt.hist(df_spatial.polya_brightness[df_spatial[!, s] .== 0], bins=bins, label="", color="black", histtype="step", lw=0.5)
end

Plt.xlim(0, 1)
Plt.xlabel("poly-A brightness"); Plt.ylabel("Num. of molecules")
Plt.legend();
Plt.tight_layout()

Plt.savefig(cplotsdir("brightness.pdf"))
```

## Watershed vs Baysor comparison on pixels


### Clustering measures

```julia tags=[]
mut_info = Dict();

samp_ids = rand(1:length(data[:watershed_labels]), 10000000);

mut_info[:watershed_pixel] = Clustering.mutualinfo(vec(data[:polya_labels])[samp_ids], vec(data[:watershed_labels])[samp_ids])
mut_info[:baysor_pixel] = Clustering.mutualinfo(vec(data[:polya_labels])[samp_ids], vec(data[:baysor_labels])[samp_ids])
mut_info[:pciseq_pixel] = Clustering.mutualinfo(vec(data[:polya_labels])[samp_ids], vec(data[:pciseq_labels])[samp_ids])

mut_info[:baysor_mol] = Clustering.mutualinfo(data[:assignment_filt][:cell], data[:assignment_filt][:cell_polya])
mut_info[:watershed_mol] = Clustering.mutualinfo(data[:assignment_filt][:cell_watershed], data[:assignment_filt][:cell_polya])
mut_info[:pciseq_mol] = Clustering.mutualinfo(data[:assignment_filt][:cell_watershed], data[:assignment_filt][:cell_pciseq])

p_df = DataFrame(
    :value => [mut_info[s] for s in [:baysor_pixel, :watershed_pixel, :pciseq_pixel, :baysor_mol, :watershed_mol, :pciseq_mol]],
    :segmentation => ["Baysor", "Watershed", "pciSeq", "Baysor", "Watershed", "pciSeq"],
    :type => ["Pixel", "Pixel", "Pixel", "Molecule", "Molecule", "Molecule"]
);
```

```julia tags=[]
plt = @df p_df StatsPlots.groupedbar(:type, :value, group=:segmentation, size=(300,400), palette=[color_per_label[s] for s in ["Baysor", "Watershed", "pciSeq"]], 
    ylabel="Mutual information", ylim=(0, 0.9), lw=0.5, xgrid=false, bg_legend=false)

Plots.savefig(cplotsdir("mutual_info.pdf"))
plt
```

### Matching-based metrics

```julia tags=[]
match_res_pix = Dict{Symbol, Any}();

# @time match_res_pix[:baysor_watershed] = BA.match_assignments(vec(data[:baysor_labels]), vec(data[:watershed_labels]));
@time match_res_pix[:baysor_polya] = BA.match_assignments(vec(data[:baysor_labels]), vec(data[:polya_labels]));
@time match_res_pix[:watershed_polya] = BA.match_assignments(vec(data[:watershed_labels]), vec(data[:polya_labels]));
@time match_res_pix[:pciseq_polya] = BA.match_assignments(vec(data[:pciseq_labels]), vec(data[:polya_labels]));
```

```julia tags=[]
ovr_counts = [B.count_array(min.(match_res_pix[s].n_overlaps[1], 3) .+ 1) for s in [:baysor_polya, :pciseq_polya, :watershed_polya]];
p_labels = ["Baysor", "pciSeq", "Watershed"]
group_var = repeat(p_labels, inner=length(ovr_counts[1]));

plt = StatsPlots.groupedbar(hcat(ovr_counts...), palette=[color_per_label[s] for s in p_labels], group=group_var,
    xlabel="Number of overlapping cells", ylabel="Num. of cells", ylim=(0, 6600), yticks=0:1500:6000, xgrid=false, lw=0.5, size=(400, 400))
Plots.hline!([maximum(data[:watershed_labels])], label="Num. poly-A cells", lw=3, color="black",
    ylims=(0, maximum(data[:watershed_labels]) + 100))
Plots.xticks!(1:4, ["0", "1", "2", "3+"])

Plots.savefig(cplotsdir("num_overlaps_pixel.pdf"))
plt
```

```julia tags=[]
frac_vals = Dict(s => Dict{Symbol, Vector}() for s in [:overlap, :assigned])

for (i,s) in zip([:baysor_labels, :watershed_labels, :pciseq_labels], [:baysor_polya, :watershed_polya, :pciseq_polya])
    overlap_mask = match_res_pix[s].n_overlaps[1] .== 1
    overlap_nums, uniq_overlap_ids = vec.(collect.(findmax(match_res_pix[s].contingency[2:end,:][overlap_mask, 2:end], dims=2)));

    frac_vals[:overlap][i] = overlap_nums ./ n_pix_per_comp[:polya_labels][getindex.(uniq_overlap_ids, 2)]
    frac_vals[:assigned][i] = n_pix_per_comp[i][overlap_mask] ./ n_pix_per_comp[:polya_labels][getindex.(uniq_overlap_ids, 2)]
end
```

```julia
fig, axes = Plt.subplots(1, 2, figsize=(8, 3))
p_keys = [:baysor_polya => "Baysor", :watershed_polya => "Watershed", :pciseq_polya => "pciSeq"]
bins = 0.0:0.04:1.0
alpha = 0.1
lw = 2
for (k, l) in p_keys
    df = BA.precision_recall_df(match_res_pix[k].contingency)
    col = color_per_label[l]
    axes[1].hist(df.precision, bins=bins, alpha=alpha, color=col);
    axes[1].hist(df.precision, bins=bins, label=l, histtype="step", color=col, lw=lw);

    axes[2].hist(df.recall, bins=bins, alpha=alpha, color=col);
    axes[2].hist(df.recall, bins=bins, label=l, histtype="step", color=col, lw=lw);
end

for ax in axes
    ax.grid(true, alpha=0.1)
    ax.set_ylabel("Num. of cells");
    ax.set_xlim(0, 1);
end

axes[1].set_xlabel("Precision")
axes[2].set_xlabel("Recall")
axes[2].legend(loc="upper left", borderpad=1, handlelength=1);
axes[1].legend(loc="upper left", borderpad=1, handlelength=1);
Plt.tight_layout();

Plt.savefig(cplotsdir("precision_recall_pixel.pdf"))
```

```julia
bins = 0.0:0.02:1.0
fig, axes = Plt.subplots(1, 2, figsize=(8, 3))
p_keys = [:baysor_labels => "Baysor", :watershed_labels => "Watershed", :pciseq_labels => "pciSeq"]
for (k,n) in p_keys
    axes[1].hist(frac_vals[:overlap][k], bins=bins, label=n, alpha=(n == "Baysor" ? 1.0 : 0.5), color=color_per_label[n], density=true)
end
for (k,v) in p_keys
    axes[1].hist(frac_vals[:overlap][k], bins=bins, density=true, color="black", histtype="step")
end

bins = 0.0:0.04:2.0
for (k,n) in p_keys
    axes[2].hist(frac_vals[:assigned][k], bins=bins, label=n, alpha=(n == "Baysor" ? 1.0 : 0.5), color=color_per_label[n], density=true)
end
for (k,v) in p_keys
    axes[2].hist(frac_vals[:assigned][k], bins=bins, density=true, color="black", histtype="step")
end

for ax in axes
    ax.set_ylabel("Cell density")
    ax.legend(handlelength=1.0, title="Target");
    
    ymax = ax.get_ylim()[2]
    ax.vlines([0.9975], ymin=0, ymax=ymax, color="black", lw=2.0);
    ax.set_ylim(0, ymax)
end

axes[1].set_xlim(0, 1); axes[2].set_xlim(0, 2);
axes[1].set_xlabel("Overlap size / poly-a cell size");
axes[2].set_xlabel("Overlap size / Target cell size");
Plt.tight_layout()

Plt.savefig(cplotsdir("mol_fracs_pixels.pdf"))
```
