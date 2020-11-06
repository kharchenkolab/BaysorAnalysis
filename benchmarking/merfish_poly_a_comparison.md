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
    display_name: Julia 1.5.2
    language: julia
    name: julia-1.5
---

```julia execution={"iopub.execute_input": "2020-11-06T12:47:06.961000+01:00", "iopub.status.busy": "2020-11-06T12:47:06.961000+01:00", "iopub.status.idle": "2020-11-06T12:47:53.141000+01:00"}
import Baysor
import Colors
import ImageMorphology
import Images
import MultivariateStats
import Plots
import Clustering
import ColorSchemes
import CSV

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

PROJECT_DIR = "/home/vpetukhov/spatial/Benchmarking/";
PLOT_DIR = "./plots/merfish/poly_a/";

Plots.theme(:wong);
Plots.default(legendfontsize=12, guidefontsize=16, tickfontsize=12, gridlinewidth=1, fg_legend=false)

cell_col_labels = ["Baysor", "poly-A", "Watershed", "Baysor Prior"]
color_per_label = Dict(Pair.(sort(cell_col_labels), Plots._all_defaults[end][:color_palette][1:length(cell_col_labels)]));
```

## Load data

```julia execution={"iopub.execute_input": "2020-11-06T12:47:53.687000+01:00", "iopub.status.busy": "2020-11-06T12:47:53.142000+01:00", "iopub.status.idle": "2020-11-06T12:48:53.891000+01:00"}
@time df_spatial, gene_names = B.load_df("/home/vpetukhov/data/spatal/merfish_moffit/merfish_coords_adj.csv");
df_spatial[!, :cell_poly_a] = denserank(df_spatial.cell);
df_spatial[!, :cell_poly_a] = ifelse.(ismissing.(df_spatial.cell_poly_a), 0, df_spatial.cell_poly_a);

@time watershed_labels = Matrix(B.load_segmentation_mask("/home/vpetukhov/data/spatal/merfish_moffit/dapi_merged_watershed.tif"));
df_spatial[!, :cell_watershed] = denserank(B.staining_value_per_transcript(df_spatial, watershed_labels)) .- 1;

df_seg = B.load_df("$PROJECT_DIR/run_results/merfish_moffit/segmentation.csv")[1];
df_spatial[!, :cell] = df_seg.cell;
df_spatial[!, :confidence] = df_seg.confidence;
gene_names_all = gene_names;

@time poly_a_arr = Float16.(Images.load("/home/vpetukhov/data/spatal/merfish_moffit/polyA_merged.tiff"));
df_spatial[!, :poly_a_brightness] = B.staining_value_per_transcript(df_spatial, poly_a_arr);

@time poly_a_labels = Matrix(B.load_segmentation_mask("/home/vpetukhov/data/spatal/merfish_moffit/moffit_poly_a_mask.tif"));
```

### Convert Baysor results to mask

```julia execution={"iopub.execute_input": "2020-11-06T12:48:53.891000+01:00", "iopub.status.busy": "2020-11-06T12:48:53.891000+01:00", "iopub.status.idle": "2020-11-06T12:53:52.891000+01:00"}
@time neighb_cm = B.neighborhood_count_matrix(df_spatial, 50);
@time color_transformation = B.gene_composition_transformation(neighb_cm, df_spatial.confidence);
@time gene_colors = B.gene_composition_colors(neighb_cm, color_transformation);

df_spatial[!, :color] = gene_colors;
```

Filter cells:

```julia execution={"iopub.execute_input": "2020-11-06T12:53:52.891000+01:00", "iopub.status.busy": "2020-11-06T12:53:52.891000+01:00", "iopub.status.idle": "2020-11-06T12:54:07.760000+01:00"}
Plots.histogram(sqrt.(B.count_array(df_spatial.cell, drop_zero=true)), bins=100, legend=false, size=(500, 300), xlabel="√(num. molecules)", ylabel="Num. cells")
Plots.vline!([sqrt(50)])
```

```julia execution={"iopub.execute_input": "2020-11-06T12:54:07.811000+01:00", "iopub.status.busy": "2020-11-06T12:54:07.811000+01:00", "iopub.status.idle": "2020-11-06T12:54:08.298000+01:00"}
filt_cells = Set(findall(B.count_array(df_spatial.cell, drop_zero=true) .< 40));
df_spatial.cell[in.(df_spatial.cell, Ref(filt_cells))] .= 0;
df_spatial.cell .= denserank(df_spatial.cell) .- 1;
```

```julia execution={"iopub.execute_input": "2020-11-06T12:54:08.329000+01:00", "iopub.status.busy": "2020-11-06T12:54:08.328000+01:00", "iopub.status.idle": "2020-11-06T12:55:31.439000+01:00"}
@time polygons = B.boundary_polygons(df_spatial, df_spatial.cell, grid_step=5.0, bandwidth=5.0, dens_threshold=1e-2, min_border_length=20);
```

```julia execution={"iopub.execute_input": "2020-11-06T12:55:31.439000+01:00", "iopub.status.busy": "2020-11-06T12:55:31.439000+01:00", "iopub.status.idle": "2020-11-06T12:55:41.240000+01:00"}
polygons_poly_a = B.extract_polygons_from_label_grid(Matrix(poly_a_labels[1:5:end, 1:5:end]'); min_border_length=50, 
    shape_method=:order, max_dev=10.0, grid_step=5.0);

polygons_watershed = B.extract_polygons_from_label_grid(Matrix(watershed_labels[1:5:end, 1:5:end]'); min_border_length=50, 
    shape_method=:order, max_dev=10.0, grid_step=5.0);
```

```julia execution={"iopub.execute_input": "2020-11-06T12:55:41.240000+01:00", "iopub.status.busy": "2020-11-06T12:55:41.240000+01:00", "iopub.status.idle": "2020-11-06T12:55:51.326000+01:00"}
Plots.plot(
    B.plot_cell_borders_polygons(@where(df_spatial, :x .< 12000, :y .< 12000, :x .> 10000, :y .> 10000), polygons, color=:color, 
        polygon_line_width=1, ms=1, ticks=false, title="Baysor"),
    B.plot_cell_borders_polygons(@where(df_spatial, :x .< 12000, :y .< 12000, :x .> 10000, :y .> 10000), polygons_poly_a, color=:color, 
        polygon_line_width=1, ms=1, ticks=false, title="Poly-A"),
    B.plot_cell_borders_polygons(@where(df_spatial, :x .< 12000, :y .< 12000, :x .> 10000, :y .> 10000), polygons_watershed, color=:color, 
        polygon_line_width=1, ms=1, ticks=false, title="Watershed"),
    size=(1200, 400), layout=(1, 3)
)
```

```julia execution={"iopub.execute_input": "2020-11-06T12:55:51.326000+01:00", "iopub.status.busy": "2020-11-06T12:55:51.326000+01:00", "iopub.status.idle": "2020-11-06T12:56:01.501000+01:00"}
@time baysor_labels = B.find_grid_point_labels_kde(B.position_data(df_spatial), denserank(df_spatial.cell) .- 1, 
    [0., 0.], Float64.(size(poly_a_labels)[[2,1]]); grid_step=5.0, bandwidth=5.0, dens_threshold=1e-2);

@time baysor_labels = ImageMorphology.closing(baysor_labels);
@time baysor_labels = B.upscale(baysor_labels', 5)[1:size(poly_a_labels, 1), 1:size(poly_a_labels, 2)];
```

### Filter segments

```julia execution={"iopub.execute_input": "2020-11-06T12:56:01.501000+01:00", "iopub.status.busy": "2020-11-06T12:56:01.501000+01:00", "iopub.status.idle": "2020-11-06T12:56:24.145000+01:00"}
labels = Dict(:baysor => baysor_labels, :poly_a => poly_a_labels, :watershed => watershed_labels);

n_pix_per_comp = Dict(n => B.count_array(vec(labels[n]), drop_zero=true) for n in keys(labels));

@showprogress for n in keys(labels)
    filt_segs = Set(findall(n_pix_per_comp[n] .< 250))
    labels[n][in.(labels[n], Ref(filt_segs))] .= 0
    vec(labels[n]) .= denserank(vec(labels[n])) .- 1
    n_pix_per_comp[n] = B.count_array(vec(labels[n]), drop_zero=true)
end
```

```julia execution={"iopub.execute_input": "2020-11-06T14:22:32.920000+01:00", "iopub.status.busy": "2020-11-06T14:22:32.920000+01:00", "iopub.status.idle": "2020-11-06T14:22:32.991000+01:00"}
bins = range(0, 250, length=100)
Plots.histogram(sqrt.(n_pix_per_comp[:baysor]), bins=bins, label="Baysor", xlabel="Num. pixels", ylabel="Num. cells")
Plots.histogram!(sqrt.(n_pix_per_comp[:poly_a]), bins=bins, label="poly-A", alpha=0.75)
Plots.histogram!(sqrt.(n_pix_per_comp[:watershed]), bins=bins, label="Watershed", alpha=0.5)
Plots.vline!([sqrt(250)], color="black", label="")
```

## Watershed vs Baysor comparison on transcripts


### Estimate matching

```julia execution={"iopub.execute_input": "2020-11-06T12:57:04.791000+01:00", "iopub.status.busy": "2020-11-06T12:57:04.791000+01:00", "iopub.status.idle": "2020-11-06T12:57:14.912000+01:00"}
qc_per_cell_dfs = B.prepare_qc_df.(Ref(df_spatial), [:cell, :cell_watershed, :cell_poly_a]; min_area=25, min_molecules_per_cell=50);
# B.plot_qc_comparison(qc_per_cell_dfs, max_quants=[0.999, 0.99, 0.99, 0.999], labels=["Baysor", "DAPI Watershed", "poly-A"])
```

```julia execution={"iopub.execute_input": "2020-11-06T12:57:14.912000+01:00", "iopub.status.busy": "2020-11-06T12:57:14.912000+01:00", "iopub.status.idle": "2020-11-06T12:57:14.936000+01:00"}
match_res = Dict{Symbol, Any}();
```

```julia execution={"iopub.execute_input": "2020-11-06T14:07:17.484000+01:00", "iopub.status.busy": "2020-11-06T14:07:17.484000+01:00", "iopub.status.idle": "2020-11-06T14:07:42.944000+01:00"}
assignments_filt = Dict(t => denserank(ifelse.(in.(df_spatial[!,s], Ref(qdf.cell_id)), df_spatial[!,s], 0)) .- 1 
    for (s, t, qdf) in zip([:cell, :cell_watershed, :cell_poly_a], [:baysor, :watershed, :poly_a], qc_per_cell_dfs));
```

```julia execution={"iopub.execute_input": "2020-11-06T14:08:10.082000+01:00", "iopub.status.busy": "2020-11-06T14:08:10.082000+01:00", "iopub.status.idle": "2020-11-06T14:08:10.334000+01:00"}
match_res[:baysor_watershed] = B.match_assignments(assignments_filt[:baysor], assignments_filt[:watershed]);
# B.plot_matching_comparison(match_res[:baysor_watershed], labels=["Baysor", "DAPI Watershed"])

match_res[:baysor_poly_a] = B.match_assignments(assignments_filt[:baysor], assignments_filt[:poly_a]);
# B.plot_matching_comparison(match_res[:baysor_poly_a], labels=["Baysor", "poly-A"])

match_res[:watershed_poly_a] = B.match_assignments(assignments_filt[:watershed], assignments_filt[:poly_a]);
# B.plot_matching_comparison(match_res[:watershed_poly_a], labels=["DAPI Watershed", "poly-A"])
```

### Plots

```julia execution={"iopub.execute_input": "2020-11-06T14:22:30.301000+01:00", "iopub.status.busy": "2020-11-06T14:22:30.301000+01:00", "iopub.status.idle": "2020-11-06T14:22:30.498000+01:00"}
n_cells_poly_a = size(qc_per_cell_dfs[3], 1)
ovr_counts = [B.count_array(min.(match_res[s].n_overlaps[1], 3) .+ 1) for s in [:baysor_poly_a, :watershed_poly_a]];
group_var = repeat(["Baysor", "Watershed"], inner=length(ovr_counts[1]));

plt = StatsPlots.groupedbar(hcat(ovr_counts...), palette=[color_per_label[s] for s in ["Baysor", "Watershed"]], group=group_var, 
    xlabel="Number of overlapping cells", ylabel="Num. of cells", ylim=(0, 6600), yticks=0:1500:6000, xgrid=false, lw=0.5, size=(400, 400))
Plots.hline!([n_cells_poly_a], label="Num. poly-A cells", lw=3, color="black", ylims=(0, n_cells_poly_a + 100))
Plots.xticks!(1:4, ["0", "1", "2", "3+"])

Plots.savefig("$PLOT_DIR/num_overlaps.pdf")
plt
```

```julia execution={"iopub.execute_input": "2020-11-06T13:41:05.498000+01:00", "iopub.status.busy": "2020-11-06T13:41:05.498000+01:00", "iopub.status.idle": "2020-11-06T13:41:05.723000+01:00"}
plt = Plots.bar(size.(qc_per_cell_dfs, 1), ylims=(0, 10000), xgrid=false, xtickfontsize=14, width=0, yticks=0:2500:10000, xlabel="Protocol",
    legend=false, ylabel="Num. of cells", color=[color_per_label[c] for c in ["Baysor", "Watershed", "poly-A"]], size=(300, 400))
# Plots.xticks!(1:3, ["Baysor", "Watershed", "poly-A"])
Plots.xticks!(Float64[])

Plots.savefig("$PLOT_DIR/num_cells.pdf")
plt
```

```julia execution={"iopub.execute_input": "2020-11-06T13:02:01.521000+01:00", "iopub.status.busy": "2020-11-06T13:02:01.521000+01:00", "iopub.status.idle": "2020-11-06T13:02:04.992000+01:00"}
frac_vals = Dict(s => Dict{Symbol, Vector}() for s in [:overlap, :assigned])

for (i, si, s) in zip(1:2, [:baysor, :watershed], [:baysor_poly_a, :watershed_poly_a])
    overlap_mask = match_res[s].n_overlaps[1] .== 1
    overlap_nums, uniq_overlap_ids = vec.(collect.(findmax(match_res[s].contingency[2:end,:][overlap_mask, 2:end], dims=2)));

    frac_vals[:overlap][si] = overlap_nums ./ qc_per_cell_dfs[3].n_transcripts[getindex.(uniq_overlap_ids, 2)]
    frac_vals[:assigned][si] = qc_per_cell_dfs[i].n_transcripts[overlap_mask] ./ qc_per_cell_dfs[3].n_transcripts[getindex.(uniq_overlap_ids, 2)]
end
```

```julia execution={"iopub.execute_input": "2020-11-06T14:22:23.129000+01:00", "iopub.status.busy": "2020-11-06T14:22:23.129000+01:00", "iopub.status.idle": "2020-11-06T14:22:23.315000+01:00"}
bins = 0.0:0.02:1.0
plt1 = Plots.histogram(frac_vals[:overlap][:baysor], bins=bins, legend=:topleft, label="Baysor", color=color_per_label["Baysor"], 
    xlabel="Overlap size / poly-a cell size", ylabel="Cell density", normalize=true, lw=0, xlim=(0, 1), ylim=(0, 3.51), legend_title="Target",
    size=(400, 300), legend_box="")
Plots.stephist!(frac_vals[:overlap][:baysor], bins=bins, label="", normalize=true, color="black")
Plots.histogram!(frac_vals[:overlap][:watershed], bins=bins, alpha=0.5, label="Watershed", normalize=true, lw=0, color=color_per_label["Watershed"])
Plots.stephist!(frac_vals[:overlap][:watershed], bins=bins, label="", normalize=true, color="black")
Plots.vline!([0.9975], label="", color="black", lw=2.5);

bins = 0.0:0.04:2.0
plt2 = Plots.histogram(frac_vals[:assigned][:baysor], bins=bins, label="Baysor", normalize=true, lw=0, color=color_per_label["Baysor"], 
    legend=false, xlabel="Target cell size / poly-A cell size", ylabel="Cell density", xlim=(0, 2), ylim=(0, 1.6), legend_title="Target",
    yticks=0:0.4:1.6, size=(400, 300))
Plots.stephist!(frac_vals[:assigned][:baysor], bins=bins, label="", normalize=true, color="black")
Plots.histogram!(frac_vals[:assigned][:watershed], bins=bins, alpha=0.5, label="Watershed", normalize=true, lw=0, color=color_per_label["Watershed"])
Plots.stephist!(frac_vals[:assigned][:watershed], bins=bins, label="", normalize=true, color="black")
Plots.vline!([1], label="", color="black", lw=2.5);

plt = Plots.plot(plt1, plt2, size=(800, 300))
Plots.savefig("$PLOT_DIR/mol_fracs.pdf")
plt
```

### Poly-A brightness

```julia execution={"iopub.execute_input": "2020-11-06T13:14:25.714000+01:00", "iopub.status.busy": "2020-11-06T13:14:25.714000+01:00", "iopub.status.idle": "2020-11-06T13:14:25.970000+01:00"}
p_dict = Dict(t => df_spatial.poly_a_brightness[df_spatial[!, s] .== 0] for (s,t) in zip([:cell, :cell_watershed, :cell_poly_a], ["Baysor", "Watershed", "poly-A"]));
```

```julia execution={"iopub.execute_input": "2020-11-06T14:22:25.554000+01:00", "iopub.status.busy": "2020-11-06T14:22:25.554000+01:00", "iopub.status.idle": "2020-11-06T14:22:26.299000+01:00"}
plt = Plots.plot(xlim=(0, 1), xlabel="poly-A brightness", ylabel="Num. of molecules", size=(400, 400))
bins = 0.0:0.016:1.0
for (s,t) in zip([:cell_watershed, :cell_poly_a, :cell], ["Watershed", "poly-A", "Baysor"])
    Plots.barhist!(df_spatial.poly_a_brightness[df_spatial[!, s] .== 0], bins=bins, label=t, color=color_per_label[t], lw=0, alpha=0.9)
end

for s in [:cell, :cell_watershed, :cell_poly_a]
    Plots.stephist!(df_spatial.poly_a_brightness[df_spatial[!, s] .== 0], bins=bins, label="", color="black", lw=0.5)
end

Plots.ylims!(0, Plots.ylims()[2])

Plots.savefig("$PLOT_DIR/brightness.pdf")

plt
```

## Watershed vs Baysor comparison on pixels


### Clustering measures

```julia execution={"iopub.execute_input": "2020-11-06T14:03:49.844000+01:00", "iopub.status.busy": "2020-11-06T14:03:49.844000+01:00", "iopub.status.idle": "2020-11-06T14:03:49.867000+01:00"}
mut_info = Dict();
```

```julia execution={"iopub.execute_input": "2020-11-06T14:03:58.327000+01:00", "iopub.status.busy": "2020-11-06T14:03:58.326000+01:00", "iopub.status.idle": "2020-11-06T14:03:58.454000+01:00"}
samp_ids = rand(1:length(watershed_labels), 10000000);
```

Poly-A vs Watershed:

```julia execution={"iopub.execute_input": "2020-11-05T14:55:24.791000+01:00", "iopub.status.busy": "2020-11-05T14:55:24.791000+01:00", "iopub.status.idle": "2020-11-05T14:55:26.286000+01:00"}
@time Clustering.randindex(vec(poly_a_labels)[samp_ids], vec(watershed_labels)[samp_ids])
```

```julia execution={"iopub.execute_input": "2020-11-06T14:04:12.641000+01:00", "iopub.status.busy": "2020-11-06T14:04:12.641000+01:00", "iopub.status.idle": "2020-11-06T14:04:13.523000+01:00"}
@time mut_info[:watershed_pixel] = Clustering.mutualinfo(vec(poly_a_labels)[samp_ids], vec(watershed_labels)[samp_ids])
```

Poly-A vs Baysor:

```julia execution={"iopub.execute_input": "2020-11-05T14:55:33.202000+01:00", "iopub.status.busy": "2020-11-05T14:55:33.202000+01:00", "iopub.status.idle": "2020-11-05T14:55:33.996000+01:00"}
@time Clustering.randindex(vec(poly_a_labels)[samp_ids], vec(baysor_labels)[samp_ids])
```

```julia execution={"iopub.execute_input": "2020-11-06T14:04:20.569000+01:00", "iopub.status.busy": "2020-11-06T14:04:20.569000+01:00", "iopub.status.idle": "2020-11-06T14:04:21.556000+01:00"}
@time mut_info[:baysor_pixel] = Clustering.mutualinfo(vec(poly_a_labels)[samp_ids], vec(baysor_labels)[samp_ids])
```

```julia execution={"iopub.execute_input": "2020-11-05T14:56:44.078000+01:00", "iopub.status.busy": "2020-11-05T14:56:44.078000+01:00", "iopub.status.idle": "2020-11-05T14:56:45.836000+01:00"}
Plots.plot([Plots.heatmap(labs[10000:5:12000, 10000:5:12000] .> 0, size=(400, 400)) for labs in [baysor_labels, poly_a_labels, watershed_labels]]..., 
    layout=(1, 3), size=(1200, 400), cbar=false)
```

```julia execution={"iopub.execute_input": "2020-11-06T14:09:01.306000+01:00", "iopub.status.busy": "2020-11-06T14:09:01.306000+01:00", "iopub.status.idle": "2020-11-06T14:09:01.599000+01:00"}
mut_info[:baysor_mol] = Clustering.mutualinfo(assignments_filt[:baysor], assignments_filt[:poly_a])
```

```julia execution={"iopub.execute_input": "2020-11-06T14:09:10.008000+01:00", "iopub.status.busy": "2020-11-06T14:09:10.008000+01:00", "iopub.status.idle": "2020-11-06T14:09:10.195000+01:00"}
mut_info[:watershed_mol] = Clustering.mutualinfo(assignments_filt[:watershed], assignments_filt[:poly_a])
```

```julia execution={"iopub.execute_input": "2020-11-06T14:12:22.333000+01:00", "iopub.status.busy": "2020-11-06T14:12:22.333000+01:00", "iopub.status.idle": "2020-11-06T14:12:22.409000+01:00"}
p_df = DataFrame(
    :value => [mut_info[s] for s in [:baysor_pixel, :watershed_pixel, :baysor_mol, :watershed_mol]], 
    :segmentation => ["Baysor", "Watershed", "Baysor", "Watershed"],
    :type => ["Pixel", "Pixel", "Molecule", "Molecule"]
);
```

```julia execution={"iopub.execute_input": "2020-11-06T14:21:45.249000+01:00", "iopub.status.busy": "2020-11-06T14:21:45.249000+01:00", "iopub.status.idle": "2020-11-06T14:21:45.457000+01:00"}
plt = @df p_df StatsPlots.groupedbar(:type, :value, group=:segmentation, size=(300,400), palette=[color_per_label[s] for s in ["Baysor", "Watershed"]], 
    ylabel="Mutual information", ylim=(0, 0.8), lw=0.5, xgrid=false, bg_legend=false)

Plots.savefig("$PLOT_DIR/pixel/mutual_info.pdf")
plt
```

### Matching-based metrics

```julia execution={"iopub.execute_input": "2020-11-06T13:57:19.680000+01:00", "iopub.status.busy": "2020-11-06T13:57:19.679000+01:00", "iopub.status.idle": "2020-11-06T13:57:47.900000+01:00"}
match_res_pix = Dict{Symbol, Any}();

@time match_res_pix[:baysor_watershed] = B.match_assignments(vec(baysor_labels), vec(watershed_labels));
@time match_res_pix[:baysor_poly_a] = B.match_assignments(vec(baysor_labels), vec(poly_a_labels));
@time match_res_pix[:watershed_poly_a] = B.match_assignments(vec(watershed_labels), vec(poly_a_labels));
```

```julia execution={"iopub.execute_input": "2020-11-06T14:22:40.410000+01:00", "iopub.status.busy": "2020-11-06T14:22:40.409000+01:00", "iopub.status.idle": "2020-11-06T14:22:40.907000+01:00"}
ovr_counts = [B.count_array(min.(match_res_pix[s].n_overlaps[1], 3) .+ 1) for s in [:baysor_poly_a, :watershed_poly_a]];
group_var = repeat(["Baysor", "Watershed"], inner=length(ovr_counts[1]));

plt = StatsPlots.groupedbar(hcat(ovr_counts...), palette=[color_per_label[s] for s in ["Baysor", "Watershed"]], group=group_var,
    xlabel="Number of overlapping cells", ylabel="Num. of cells", ylim=(0, 6600), yticks=0:1500:6000, xgrid=false, lw=0.5, size=(400, 400))
Plots.hline!([maximum(watershed_labels)], label="Num. poly-A cells", lw=3, color="black",
    ylims=(0, maximum(watershed_labels) + 100))
Plots.xticks!(1:4, ["0", "1", "2", "3+"])

Plots.savefig("$PLOT_DIR/pixel/num_overlaps.pdf")
plt
```

```julia execution={"iopub.execute_input": "2020-11-05T15:40:52.129000+01:00", "iopub.status.busy": "2020-11-05T15:40:52.129000+01:00", "iopub.status.idle": "2020-11-05T15:40:55.230000+01:00"}
frac_vals = Dict(s => Dict{Symbol, Vector}() for s in [:overlap, :assigned])

for (i,s) in zip([:baysor, :watershed], [:baysor_poly_a, :watershed_poly_a])
    overlap_mask = match_res_pix[s].n_overlaps[1] .== 1
    overlap_nums, uniq_overlap_ids = vec.(collect.(findmax(match_res_pix[s].contingency[2:end,:][overlap_mask, 2:end], dims=2)));

    frac_vals[:overlap][i] = overlap_nums ./ n_pix_per_comp[:poly_a][getindex.(uniq_overlap_ids, 2)]
    frac_vals[:assigned][i] = n_pix_per_comp[i][overlap_mask] ./ n_pix_per_comp[:poly_a][getindex.(uniq_overlap_ids, 2)]
end
```

```julia execution={"iopub.execute_input": "2020-11-06T14:22:42.752000+01:00", "iopub.status.busy": "2020-11-06T14:22:42.752000+01:00", "iopub.status.idle": "2020-11-06T14:22:42.919000+01:00"}
bins = 0.0:0.02:1.0
plt1 = Plots.histogram(frac_vals[:overlap][:baysor], bins=bins, legend=false, label="Baysor", color=color_per_label["Baysor"], 
    xlabel="Overlap size / poly-a cell size", ylabel="Cell density", normalize=true, lw=0, xlim=(0, 1), ylim=(0, 2.2), legend_title="Target",
    size=(400, 300), legend_box="")
Plots.stephist!(frac_vals[:overlap][:baysor], bins=bins, label="", normalize=true, color="black")
Plots.histogram!(frac_vals[:overlap][:watershed], bins=bins, alpha=0.5, label="Watershed", normalize=true, lw=0, color=color_per_label["Watershed"])
Plots.stephist!(frac_vals[:overlap][:watershed], bins=bins, label="", normalize=true, color="black")
Plots.vline!([0.9975], label="", color="black", lw=2.5);

bins = 0.0:0.04:2.0
plt2 = Plots.histogram(frac_vals[:assigned][:baysor], bins=bins, label="Baysor", normalize=true, lw=0, color=color_per_label["Baysor"], 
    legend=:topright, xlabel="Target cell size / poly-A cell size", ylabel="Cell density", xlim=(0, 2), ylim=(0, 1.35), legend_title="Target",
    yticks=0:0.4:1.2, size=(400, 300))
Plots.stephist!(frac_vals[:assigned][:baysor], bins=bins, label="", normalize=true, color="black")
Plots.histogram!(frac_vals[:assigned][:watershed], bins=bins, alpha=0.5, label="Watershed", normalize=true, lw=0, color=color_per_label["Watershed"])
Plots.stephist!(frac_vals[:assigned][:watershed], bins=bins, label="", normalize=true, color="black")
Plots.vline!([1], label="", color="black", lw=2.5);

plt = Plots.plot(plt1, plt2, size=(800, 300))
Plots.savefig("$PLOT_DIR/pixel/mol_fracs.pdf")
plt
```
