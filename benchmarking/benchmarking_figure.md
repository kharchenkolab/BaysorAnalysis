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
# import Gadfly
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
# GF = Gadfly;

PROJECT_DIR = "/home/vpetukhov/spatial/Benchmarking/";
PLOT_DIR = "./plots/joint_fig/";
```

```julia
Plots.theme(:wong); # :mute, :vibrant, :bright
```

```julia
R"""
library(ggplot2)
library(ggrastr)
library(ggforce)
theme_set(theme_bw())
""";

plotCorrelationEffect = R"""
function(df, frac_col, ymax=1.0, ylabel=ifelse(frac_col == "MolFrac", "Fraction of mismatching molecules", "Fraction of mismatching cells"), legend_pos=c(0,1), 
         color_pal=scales::hue_pal()(4)) {
    low.max <- max(df[[frac_col]][df$Correlation < 0.5])
    ggplot() + 
        geom_rect(aes(xmin=0.01, xmax=0.5, ymin=0.001, ymax=ym), data.frame(ym=low.max), alpha=0.5, fill=alpha("white", 0.0), color="black") +
        geom_line(aes(x=Correlation, y=.data[[frac_col]], color=Dataset, linetype=Segmentation), df) +
        theme(legend.position=legend_pos, legend.justification=legend_pos, legend.background=element_rect(fill=alpha('white', 0.2)), legend.box='horizontal') +
        guides(color=guide_legend(order=1)) +
        labs(x="Maximal correlation", y=ylabel) +
        scale_color_manual(values=color_pal) +
        scale_x_continuous(limits=c(0, 1.01), expand=c(0, 0)) +
        scale_y_continuous(limits=c(0, ymax), expand=c(0, 0))
}
""";
```

<!-- #region toc-hr-collapsed=true toc-nb-collapsed=true -->
## Load data
<!-- #endregion -->

### MERFISH

```julia
# df_spatial = CSV.read("/home/vpetukhov/data/spatal/merfish_moffit/merfish_coords_perprocessed.csv") |> DataFrame;
# df_spatial[!, :cell] = ifelse.(ismissing.(df_spatial.cell), 0, denserank(df_spatial.cell));
# CSV.write("/home/vpetukhov/data/spatal/merfish_moffit/merfish_coords_preprocessed2.csv", df_spatial);
```

```julia
@time df_spatial, gene_names = B.load_df("/home/vpetukhov/data/spatal/merfish_moffit/merfish_coords_adj.csv");
df_spatial[!, :cell_dapi] = ifelse.(ismissing.(df_spatial.cell), 0, denserank(df_spatial.cell));

df_seg = B.load_df("$PROJECT_DIR/run_results/merfish_moffit/segmentation.csv")[1];
df_spatial[!, :cell] = df_seg.cell;
df_spatial[!, :confidence] = df_seg.confidence;

df_seg_prior = B.load_df("$PROJECT_DIR/run_results/merfish_moffit_prior/segmentation.csv")[1];
df_spatial[!, :cell_prior] = df_seg_prior.cell;

@time watershed_labels = B.load_segmentation_mask("/home/vpetukhov/data/spatal/merfish_moffit/dapi_merged_watershed.tif");
df_spatial[!, :cell_watershed] = denserank(B.staining_value_per_transcript(df_spatial, watershed_labels)) .- 1;

dapi_arr = Float16.(Images.load("/home/vpetukhov/data/spatal/merfish_moffit/dapi_merged.tiff"));

watershed_labels = nothing
GC.gc()

@time paper_polys = B.boundary_polygons(df_spatial, df_spatial.cell_dapi, grid_step=5.0, bandwidth=5.0);

merfish = Dict(:df => df_spatial, :gene_names => gene_names, :min_area => 25, :name => "MERFISH", :dapi_arr => dapi_arr, :paper_polys => paper_polys);
```

### osmFISH

```julia
# import HDF5
# @time dapi = HDF5.h5read("/home/vpetukhov/data/spatal/linnarsson/Nuclei_polyT.int16.sf.hdf5", "/nuclei");
# t_quant = UInt(quantile(vec(dapi), 0.995));
# dapi[dapi .> t_quant] .= t_quant;

# @time Images.save("/home/vpetukhov/data/spatal/linnarsson/nuclei.tif", dapi);
```

```julia
@time df_spatial, gene_names = B.load_df("$PROJECT_DIR/run_results/osm_fish/segmentation.csv");
# @time df_spatial, gene_names = B.load_df("$PROJECT_DIR/run_results/test_small_cell_penalty/osm_fish/segmentation.csv");

# df_seg_prior = B.load_df("$PROJECT_DIR/run_results/osm_fish/paper_prior/segmentation.csv")[1];
# df_spatial[!, :cell_prior] = df_seg_prior.cell;

@time dapi_seg_labels = B.load_segmentation_mask("/home/vpetukhov/data/spatal/linnarsson/paper_segmentation.tiff");
df_spatial[!, :cell_dapi] = denserank(B.staining_value_per_transcript(df_spatial, dapi_seg_labels)) .- 1;

@time watershed_labels = B.load_segmentation_mask("/home/vpetukhov/data/spatal/linnarsson/nuclei_watershed.tif");
df_spatial[!, :cell_watershed] = denserank(B.staining_value_per_transcript(df_spatial, watershed_labels)) .- 1;

dapi_arr = Float16.(Images.load("/home/vpetukhov/data/spatal/linnarsson/nuclei.tif"));

df_spatial.cell .= denserank(df_spatial.cell) .- 1;

@time paper_polys = B.extract_polygons_from_label_grid(Matrix(dapi_seg_labels[1:3:end, 1:3:end]'), grid_step=3.0)

watershed_labels = dapi_seg_labels = nothing;
GC.gc()

osmfish = Dict(:df => df_spatial, :gene_names => gene_names, :min_area => 70, :name => "osmFISH", :dapi_arr => dapi_arr, :paper_polys => paper_polys);
```

### STARMAP 1020

```julia
@time df_spatial, gene_names = B.load_df("$PROJECT_DIR/run_results/star_map/vis_1020_cl0/segmentation.csv");
# @time df_spatial, gene_names = B.load_df("$PROJECT_DIR/run_results/test_small_cell_penalty/star_map_visp1020/segmentation.csv");

@time dapi_seg_labels = B.load_segmentation_mask("/home/vpetukhov/data/spatal/star_map/visual_1020_20180505_BY3_1kgenes/segmentation.tiff");
df_spatial[!, :cell_watershed] .= 0;
df_spatial[!, :cell_dapi] = denserank(B.staining_value_per_transcript(df_spatial, dapi_seg_labels)) .- 1;
df_spatial.cell .= denserank(df_spatial.cell) .- 1;

df_seg_prior = B.load_df("$PROJECT_DIR/run_results/star_map/vis_1020_prior_cl0/segmentation.csv")[1];
df_spatial[!, :cell_prior] = df_seg_prior.cell;

@time paper_polys = B.extract_polygons_from_label_grid(Matrix(dapi_seg_labels[1:3:end, 1:3:end]'), grid_step=3.0)

watershed_labels = dapi_seg_labels = nothing;
GC.gc()

starmap1020 = Dict(:df => df_spatial, :gene_names => gene_names, :min_area => 25, :name => "STARmap", :paper_polys => paper_polys);
```

```julia
# df_seg_prior = B.load_df("$PROJECT_DIR/run_results/osm_fish/paper_prior/segmentation.csv")[1];
# df_spatial[!, :cell_prior] = df_seg_prior.cell;

# starmap1020[:df][!, :cell_prior] = df_seg_prior.cell;
# datasets.starmap1020[:df][!, :cell_prior] = df_seg_prior.cell;
```

### Allen smFISH

```julia
@time df_spatial, gene_names = B.load_df("$PROJECT_DIR/run_results/spacejam2/allen_sm_fish/no_dapi/segmentation.csv");
# @time df_spatial, gene_names = B.load_df("$PROJECT_DIR/run_results/test_small_cell_penalty/allen_smfish/segmentation.csv");

df_spatial.x .-= minimum(df_spatial.x);
df_spatial.y .-= minimum(df_spatial.y);

df_spatial.x .*= 10.0;
df_spatial.y .*= 10.0;

df_seg_prior = B.load_df("$PROJECT_DIR/run_results/spacejam2/allen_sm_fish/mask_prior/segmentation.csv")[1];
df_spatial[!, :cell_prior] = df_seg_prior.cell;

@time dapi_seg_labels = B.load_segmentation_mask("/home/vpetukhov/data/spatal/SpaceJam2Full/allen_sm_fish/segmentation_labels_from_json_transposed.tiff");
df_spatial[!, :cell_dapi] = denserank(B.staining_value_per_transcript(df_spatial, dapi_seg_labels)) .- 1;

@time watershed_labels = B.load_segmentation_mask("/home/vpetukhov/data/spatal/SpaceJam2Full/allen_sm_fish/dapi_merged_watershed.tif");
df_spatial[!, :cell_watershed] = denserank(B.staining_value_per_transcript(df_spatial, watershed_labels)) .- 1;

dapi_arr = Float16.(Images.load("/home/vpetukhov/data/spatal/SpaceJam2Full/allen_sm_fish/dapi_merged.tiff"));

@time paper_polys = B.extract_polygons_from_label_grid(Matrix(dapi_seg_labels[1:3:end, 1:3:end]'), grid_step=3.0)

df_spatial.cell .= denserank(df_spatial.cell) .- 1;
watershed_labels = dapi_seg_labels = nothing;
GC.gc()

allen_smfish = Dict(:df => df_spatial, :gene_names => gene_names, :min_area => 25, :name => "Allen smFISH", :dapi_arr => dapi_arr, :paper_polys => paper_polys);
```

### ISS

```julia
@time df_spatial, gene_names = B.load_df("$PROJECT_DIR/run_results/iss_hippo/ca1_no_prior/segmentation.csv");

df_spatial[!, :cell_dapi] = df_spatial.parent_id;

df_seg_prior = B.load_df("$PROJECT_DIR/run_results/iss_hippo/ca1_paper_prior/segmentation.csv")[1];
df_spatial[!, :cell_prior] = df_seg_prior.cell;

@time watershed_labels = Matrix(B.load_segmentation_mask("/home/vpetukhov/data/spatal/iss/hippocampus/CA1/Viktor/CA1DapiBoundaries_4-3_right_watershed.tif"));
df_spatial[!, :cell_watershed] = denserank(B.staining_value_per_transcript(df_spatial, watershed_labels)) .- 1;

dapi_arr = Float16.(Images.load("/home/vpetukhov/data/spatal/iss/hippocampus/CA1/Viktor/CA1DapiBoundaries_4-3_right.tif"));

df_spatial.cell .= denserank(df_spatial.cell) .- 1;
@time paper_polys = B.boundary_polygons(df_spatial, df_spatial.cell_dapi, grid_step=1.0, bandwidth=3.0);

watershed_labels = nothing;
GC.gc()

iss = Dict(:df => df_spatial, :gene_names => gene_names, :min_area => 1, :min_mols_per_cell => 3, :name => "ISS", :dapi_arr => dapi_arr, :paper_polys => paper_polys);
```

```julia
# df_seg_prior = B.load_df("$PROJECT_DIR/run_results/iss_hippo/ca1_paper_prior/segmentation.csv")[1];
# df_spatial[!, :cell_prior] = df_seg_prior.cell;

# iss[:df][!, :cell_prior] = df_seg_prior.cell;
# datasets.iss[:df][!, :cell_prior] = df_seg_prior.cell;
```

```julia
# @time neighb_cm = Baysor.neighborhood_count_matrix(df_spatial, 50);
# @time color_transformation = Baysor.gene_composition_transformation(neighb_cm, df_spatial.confidence);
# @time gene_colors = Baysor.gene_composition_colors(neighb_cm, color_transformation; color_range=400);

# df_spatial[!, :color] = gene_colors;
```

```julia
# # df_spatial[!, :color] = B.distinguishable_colors(df_spatial.cell_watershed)[:colors];
# B.plot_comparison_for_cell(df_spatial, (1000, 2000), (2500, 3500), watershed_labels, dapi_arr, cell_col=:cell_watershed, ms=5.0, alpha=1.0, 
#     grid_step=2.0, bandwidth=2.0)
```

## Compare

```julia
datasets = deepcopy((osmfish=osmfish, merfish=merfish, allen_smfish=allen_smfish, starmap1020=starmap1020, iss=iss));
```

### Pre-process and detailed plots

```julia
# Clustering.randindex(allen_smfish[:df].cell, allen_smfish[:df].cell_dapi)[1]
```

```julia
cell_col_names = [:cell, :cell_dapi, :cell_watershed, :cell_prior];
cell_col_labels = ["Baysor", "Paper", "Watershed", "Baysor with Prior"]
for k in keys(datasets)
    d = datasets[k]
    d[:qc_per_cell_dfs] = B.prepare_qc_df.(Ref(d[:df]), cell_col_names; min_area=d[:min_area], 
        min_molecules_per_cell=get(d, :min_mols_per_cell, 50), max_elongation=15);
    
    for (i, cs) in enumerate(cell_col_names)
        d[:df][!, cs][.!in.(d[:df][!, cs], Ref(Set(d[:qc_per_cell_dfs][i].cell_id)))] .= 0
    end
    println(k)
    display(B.plot_qc_comparison(d[:qc_per_cell_dfs], max_quants=[0.995, 0.99, 0.99, 0.999, 0.999], labels=cell_col_labels))
end
```

```julia
# for k in keys(datasets)
#     d = datasets[k]
#     println(k)
#     display(B.plot_matching_comparison(d[:match_res_watershed]))
# end
```

```julia
for k in keys(datasets)
    d = datasets[k]
#     assignments_filt = [denserank(ifelse.(in.(d[:df][!,s], Ref(qdf.cell_id)), d[:df][!,s], 0)) .- 1 for (s, qdf) in zip([:cell, :cell_dapi], d[:qc_per_cell_dfs])];
    # Already filtered above
    d[:match_res] = B.match_assignments([denserank(d[:df][!, cs]) .- 1 for cs in [:cell, :cell_dapi]]...);
    d[:match_res_watershed] = B.match_assignments([denserank(d[:df][!, cs]) .- 1 for cs in [:cell, :cell_watershed]]...);
    d[:match_res_prior_dapi] = B.match_assignments([denserank(d[:df][!, cs]) .- 1 for cs in [:cell_prior, :cell_dapi]]...);
    d[:match_res_prior] = B.match_assignments([denserank(d[:df][!, cs]) .- 1 for cs in [:cell, :cell_prior]]...);

    println(k)
    display(B.plot_matching_comparison(d[:match_res]))
end
```

```julia
for k in keys(datasets)
    d = datasets[k]
    d[:stat_df] = B.build_statistics_df(d[:qc_per_cell_dfs][1:2], d[:match_res], d[:df])
    display(d[:stat_df])
end
```

```julia
@time for k in keys(datasets)
    if k == :iss
        continue
    end
    d = datasets[k]
    d[:part_cors] = [B.estimate_non_matching_part_correlation(d[:df], d[:qc_per_cell_dfs][1:2], d[:match_res], rev=r) for r in [false, true]]

    if size(d[:qc_per_cell_dfs][3], 1) > 0
        d[:part_cors_watershed] = [B.estimate_non_matching_part_correlation(d[:df], d[:qc_per_cell_dfs][[1,3]], d[:match_res_watershed], rev=r, 
                cell_cols=[:cell, :cell_watershed]) for r in [false, true]]
    end

    if size(d[:qc_per_cell_dfs][4], 1) > 0
        d[:part_cors_prior_dapi] = [B.estimate_non_matching_part_correlation(d[:df], d[:qc_per_cell_dfs][[4,2]], d[:match_res_prior_dapi], rev=r, 
                cell_cols=[:cell_prior, :cell_dapi]) for r in [false, true]]

        d[:part_cors_prior] = [B.estimate_non_matching_part_correlation(d[:df], d[:qc_per_cell_dfs][[1,4]], d[:match_res_prior], rev=r, 
                cell_cols=[:cell, :cell_prior]) for r in [false, true]]
    end

    t_bins = -0.05:0.02:1.0
    plt = Plots.histogram(d[:part_cors][1][1], bins=t_bins, widen=false, label="Baysor", legend=:topleft, 
        xlabel="Correlation", ylabel="Num. of cells", title=k);
    Plots.histogram!(d[:part_cors][2][1], bins=t_bins, label="DAPI", alpha=0.6)
    display(plt)
end
```

### Assignment stats

```julia
p_df = vcat([vcat([DataFrame(:Protocol => d[:name], :Frac => mean(d[:df][!, s] .!= 0), :Type => n) for d in datasets]...) for (s,n) in zip(cell_col_names, cell_col_labels)]...);

plt = @df p_df groupedbar(String.(:Protocol), :Frac, group=:Type, ylabel = "Fraction of molecules\nassigned to cells", xgrid=false,
    lw = 0, widen=false, ylims=(0, 1.0), size=(420, 300), bg_legend=Colors.RGBA(1.0, 1.0, 1.0, 0.5), legend=:topright, xtickfontsize=8)

Plots.savefig(plt, "$PLOT_DIR/assignment_frac.pdf")
plt
```

```julia
p_df = vcat([vcat([DataFrame(:Protocol => d[:name], :NCells => size(d[:qc_per_cell_dfs][i], 1), :Type => s) for d in datasets]...) for (i,s) in enumerate(cell_col_labels)]...);

plt = @df p_df groupedbar(String.(:Protocol), :NCells, group=:Type, ylabel = "Number of cells", xgrid=false,
    lw = 0, widen=false, ylims=(0, 11000), size=(420, 300), bg_legend=Colors.RGBA(1.0, 1.0, 1.0, 0.5), legend=:topleft, xtickfontsize=8)

Plots.savefig(plt, "$PLOT_DIR/num_cells.pdf")
plt
```

```julia
import Clustering: mutualinfo
p_df = vcat([DataFrame(
    :Rand => [mutualinfo(d[:df].cell, d[:df].cell_dapi)[1], mutualinfo(d[:df].cell, d[:df].cell_watershed)[1], mutualinfo(d[:df].cell_dapi, d[:df].cell_watershed)[1],
                mutualinfo(d[:df].cell_prior, d[:df].cell)[1], mutualinfo(d[:df].cell_prior, d[:df].cell_dapi)[1]], 
    :Type => ["Baysor vs Paper", "Watershed vs Baysor", "Watershed vs Paper", "Baysor with Prior vs Baysor", "Baysor with Prior vs Paper"],
    :Protocol => d[:name]
) for d in datasets]...)

plt = @df p_df groupedbar(String.(:Protocol), :Rand, group=:Type, ylabel = "Mutual Info", xgrid=false,
    lw = 0, widen=false, ylims=(0, 1.0), size=(570, 300), bg_legend=Colors.RGBA(1.0, 1.0, 1.0, 0.5), legend=:outerright, xtickfontsize=8)

Plots.savefig(plt, "./plots/mutual_info.pdf")
plt
```

```julia
# p_df = vcat([DataFrame(
#     :Rand => [Clustering.varinfo(d[:df].cell, d[:df].cell_dapi)[1], Clustering.varinfo(d[:df].cell, d[:df].cell_watershed)[1], Clustering.varinfo(d[:df].cell_dapi, d[:df].cell_watershed)[1],
#                 Clustering.varinfo(d[:df].cell_prior, d[:df].cell)[1], Clustering.varinfo(d[:df].cell_prior, d[:df].cell_dapi)[1]], 
#     :Type => ["Baysor vs Paper", "Baysor vs Watershed", "Watershed vs Paper", "Baysor Prior vs Baysor", "Baysor Prior vs Paper"],
#     :Protocol => d[:name]
# ) for d in datasets]...)

# plt = @df p_df groupedbar(String.(:Protocol), :Rand, group=:Type, ylabel = "Variation of Information", 
#     lw = 0, widen=false, ylims=(0, 6.0), size=(370, 300), bg_legend=Colors.RGBA(1.0, 1.0, 1.0, 0.5), legend=:topleft, xtickfontsize=8)
```

### Correlation plots

```julia
color_per_col_label = Dict(Pair.(sort(cell_col_labels), Plots._all_defaults[end][:color_palette][1:length(cell_col_labels)]));
```

```julia
pdfs1 = [DataFrame(:Correlation => d[:part_cors][1][1], :Type => d[:name]) for d in datasets if :part_cors in keys(d)];
pdfs2 = [DataFrame(:Correlation => d[:part_cors][2][1], :Type => d[:name]) for d in datasets if :part_cors in keys(d)];
t_heights = (0.9 .* hcat(size.(pdfs1, 1), size.(pdfs2, 1)) ./ max.(size.(pdfs1, 1), size.(pdfs2, 1)))[sortperm([d.Type[1] for d in pdfs2]),:]

plt = @df vcat(pdfs1...) violin(:Type, :Correlation, side=:right, label="Baysor", ylabel="Correlation between parts", 
    bg_legend=Colors.RGBA(1.0, 1.0, 1.0, 0.9), bar_width=t_heights[:,1], size=(500, 240), legend=:bottomright, xgrid=false, color=color_per_col_label["Baysor"])
@df vcat(pdfs2...) violin!(:Type, :Correlation, side=:left, label="Paper", bar_width=t_heights[:,2], color=color_per_col_label["Paper"])

Plots.savefig(plt, "$PLOT_DIR/correlations_paper.pdf")
plt
```

```julia
pdfs1 = [DataFrame(:Correlation => d[:part_cors_watershed][1][1], :Type => d[:name]) for d in datasets if :part_cors_watershed in keys(d)];
pdfs2 = [DataFrame(:Correlation => d[:part_cors_watershed][2][1], :Type => d[:name]) for d in datasets if :part_cors_watershed in keys(d)];
t_heights = (0.9 .* hcat(size.(pdfs1, 1), size.(pdfs2, 1)) ./ max.(size.(pdfs1, 1), size.(pdfs2, 1)))[sortperm([d.Type[1] for d in pdfs2]),:]

plt = @df vcat(pdfs1...) violin(:Type, :Correlation, side=:right, label="Baysor", ylabel="Correlation between parts", 
    bg_legend=Colors.RGBA(1.0, 1.0, 1.0, 0.9), bar_width=t_heights[:,1], size=(500, 240), legend=:bottomright, xgrid=false, color=color_per_col_label["Baysor"])
@df vcat(pdfs2...) violin!(:Type, :Correlation, side=:left, label="Watershed", bar_width=t_heights[:,2], color=color_per_col_label["Watershed"])

Plots.savefig(plt, "$PLOT_DIR/correlations_watershed.pdf")
plt
```

```julia
p_df = B.correlation_effect_size_df(datasets, :part_cors, cell_col_labels, [1,2]);
plts1 = plotCorrelationEffect.(Ref(p_df), ["CellFrac", "MolFrac"], [0.65, 0.32]);
p_df = B.correlation_effect_size_df(datasets, :part_cors_watershed, cell_col_labels, [1,3]);
plts2 = plotCorrelationEffect.(Ref(p_df), ["CellFrac", "MolFrac"], [0.65, 0.32]);

plt = R"cowplot::plot_grid"(plts1[1], plts2[1], plts1[2], plts2[2], ncol=2);

RCall.ijulia_setdevice(MIME("image/svg+xml"), width=8, height=5)
R"ggsave"("./plots/corr_effect.pdf", plt, width=8, height=5);

plt
```

<!-- #region toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true -->
#### Debug Prior Cors
<!-- #endregion -->

```julia
# pdfs1 = [DataFrame(:Correlation => d[:part_cors_prior][1][1], :Type => d[:name]) for d in datasets if :part_cors_prior in keys(d)];
# pdfs2 = [DataFrame(:Correlation => d[:part_cors_prior][2][1], :Type => d[:name]) for d in datasets if :part_cors_prior in keys(d)];
# t_heights = (0.9 .* hcat(size.(pdfs1, 1), size.(pdfs2, 1)) ./ max.(size.(pdfs1, 1), size.(pdfs2, 1)))[sortperm([d.Type[1] for d in pdfs2]),:]

# plt = @df vcat(pdfs1...) violin(:Type, :Correlation, side=:right, marker=(0.2, :red, stroke(0)), label="Baysor", ylabel="Correlation between parts", 
#     bg_legend=Colors.RGBA(1.0, 1.0, 1.0, 0.9), bar_width=t_heights[:,1], size=(500, 300), legend=:bottomright, xgrid=false)
# @df vcat(pdfs2...) violin!(:Type, :Correlation, side=:left, marker=(0.2, :blue, stroke(0)), label="Baysor Prior", bar_width=t_heights[:,2])

# plt
```

```julia
# pdfs1 = [DataFrame(
#     :Correlation => 
#         d[:part_cors_prior][1][1][
#             (d[:part_cors_prior][1][3] .< 0.9) .& 
#             (((1 .- d[:part_cors_prior][1][3]) .* d[:qc_per_cell_dfs][1].n_transcripts[d[:part_cors_prior][1][4]]) .> 20)
#         ], 
#     :Type => d[:name]) for d in datasets if :part_cors_prior in keys(d)];
# pdfs2 = [DataFrame(
#     :Correlation => 
#         d[:part_cors_prior][2][1][
#             (d[:part_cors_prior][2][3] .< 0.9) .& 
#             (((1 .- d[:part_cors_prior][2][3]) .* d[:qc_per_cell_dfs][4].n_transcripts[d[:part_cors_prior][2][4]]) .> 20)
#         ], 
#     :Type => d[:name]) for d in datasets if :part_cors_prior in keys(d)];
# t_heights = (0.9 .* hcat(size.(pdfs1, 1), size.(pdfs2, 1)) ./ max.(size.(pdfs1, 1), size.(pdfs2, 1)))[sortperm([d.Type[1] for d in pdfs2]),:]

# plt = @df vcat(pdfs1...) violin(:Type, :Correlation, side=:right, marker=(0.2, :red, stroke(0)), label="Baysor", ylabel="Correlation between parts", 
#     bg_legend=Colors.RGBA(1.0, 1.0, 1.0, 0.9), bar_width=t_heights[:,1], size=(500, 300), legend=:bottomright, xgrid=false)
# @df vcat(pdfs2...) violin!(:Type, :Correlation, side=:left, marker=(0.2, :blue, stroke(0)), label="Baysor Prior", bar_width=t_heights[:,2])

# plt
```

```julia
t_d = datasets.merfish
@time prior_polys = B.boundary_polygons(t_d[:df], t_d[:df].cell_prior, grid_step=5.0, bandwidth=5.0);
@time polys = B.boundary_polygons(t_d[:df], t_d[:df].cell, grid_step=5.0, bandwidth=5.0);
```

```julia
B.plot_comparison_for_cell(t_d[:df], (5870, 6290), (3545, 3906), nothing, t_d[:dapi_arr]; paper_polys=prior_polys, cell_col=:cell, 
            ms=4.0, bandwidth=5.0, ticks=true)
```

```julia
B.plot_comparison_for_cell(t_d[:df], (5870, 6290), (3545, 3906), nothing, t_d[:dapi_arr]; paper_polys=prior_polys, cell_col=:cell_prior, 
            ms=4.0, bandwidth=5.0, grid_step=5.0, ticks=true)
```

```julia
t_cd = t_d[:part_cors_prior][1]
t_qc = t_d[:qc_per_cell_dfs][1]

for ci in t_cd[2][sortperm(t_cd[1][(t_qc.n_transcripts[t_cd[4]] .* (1 .- t_cd[3])) .> 100])][1:5]
    display(B.plot_comparison_for_cell(t_d[:df], ci, nothing, t_d[:dapi_arr]; paper_polys=prior_polys, cell1_col=:cell, cell2_col=:cell_prior, 
            ms=4.0, bandwidth=5.0, ticks=true))
end
```

```julia
t_cd = t_d[:part_cors_prior][2]
t_qc = t_d[:qc_per_cell_dfs][4]

for ci in t_cd[2][sortperm(t_cd[1][(t_qc.n_transcripts[t_cd[4]] .* (1 .- t_cd[3])) .> 100])][1:5]
    display(B.plot_comparison_for_cell(t_d[:df], ci, nothing, t_d[:dapi_arr]; paper_polys=polys, cell1_col=:cell_prior, cell2_col=:cell, 
            ms=4.0, bandwidth=5.0, ticks=true))
end
```

<!-- #region toc-hr-collapsed=true toc-nb-collapsed=true -->
### Correlation examples
<!-- #endregion -->

```julia
# for d in datasets
for n in [:osmfish, :allen_smfish, :merfish]
    d = datasets[n]
    @time neighb_cm = B.neighborhood_count_matrix(d[:df], 50);
    @time color_transformation = B.gene_composition_transformation(neighb_cm, d[:df].confidence);
    @time gene_colors = B.gene_composition_colors(neighb_cm, color_transformation; color_range=400);

    d[:df][!, :gene_color] = gene_colors;
    d[:df][!, :color] = deepcopy(gene_colors);
end
```

```julia
# t_df = @where(datasets.osmfish[:df], :x .< 19070, :x .> 19045, :y .< 19585, :y .> 19550);
# B.plot_cell_borders_polygons(t_df, annotation=datasets.osmfish[:gene_names][t_df.gene], ms=8, legend=true)
```

```julia
t_d = datasets.osmfish
xls, yls = (18900, 19160), (19480, 19710)
c_polys = t_d[:paper_polys][[all((p[:,2] .< yls[2]) .& (p[:,2] .> yls[1]) .& (p[:,1] .< xls[2]) .& (p[:,1] .> xls[1])) for p in t_d[:paper_polys]]];
plt = B.plot_comparison_for_cell(t_d[:df], xls, yls, nothing, t_d[:dapi_arr]; paper_polys=c_polys, paper_poly_color="#e36200", paper_line_mult=1.5, plot_raw_dapi=false,
    cell1_col=:cell_dapi, cell2_col=:cell, ms=4.0, bandwidth=10.0, grid_step=3.0, ticks=false, alpha=0.5, dapi_alpha=0.6, polygon_line_width=3.0, polygon_alpha=0.75, 
    size_mult=1.0, ylabel="osmFISH", labelfontsize=12)

Plots.savefig(plt, "$PLOT_DIR/examples/osm_fish.png")
plt
```

```julia
t_d = datasets.allen_smfish
xls, yls = (15125, 15320), (11890, 12110)
c_polys = t_d[:paper_polys][[all((p[:,2] .< yls[2]) .& (p[:,2] .> yls[1]) .& (p[:,1] .< xls[2]) .& (p[:,1] .> xls[1])) for p in t_d[:paper_polys]]];
plt = B.plot_comparison_for_cell(t_d[:df], xls, yls, nothing, t_d[:dapi_arr]; paper_polys=c_polys, paper_poly_color="#e36200", paper_line_mult=2.0, plot_raw_dapi=false, ticks=false,
    cell1_col=:cell_dapi, cell2_col=:cell, ms=4.0, bandwidth=5.0, grid_step=3.0, polygon_line_width=3.0, alpha=0.5, dapi_alpha=0.75, polygon_alpha=0.75, size_mult=1.5,
    ylabel="Allen smFISH", labelfontsize=14)

Plots.savefig(plt, "$PLOT_DIR/examples/allen_smfish.png")
plt
```

```julia
t_d = datasets.merfish
xls, yls = (11620, 11820), (10265, 10445)
c_polys = t_d[:paper_polys][[all((p[:,2] .< yls[2]) .& (p[:,2] .> yls[1]) .& (p[:,1] .< xls[2]) .& (p[:,1] .> xls[1])) for p in t_d[:paper_polys]]];
plt = B.plot_comparison_for_cell(t_d[:df], xls, yls, nothing, t_d[:dapi_arr]; paper_polys=c_polys, paper_poly_color="#e36200", paper_line_mult=2.0, plot_raw_dapi=false, ticks=false,
    cell1_col=:cell_dapi, cell2_col=:cell, ms=4.0, bandwidth=5.0, grid_step=3.0, polygon_line_width=3.0, alpha=0.5, dapi_alpha=0.75, polygon_alpha=0.75, size_mult=1.5,
    ylabel="MERFISH", labelfontsize=14)

Plots.savefig(plt, "$PLOT_DIR/examples/merfish.png")
plt
```

```julia
# t_cd = datasets.osmfish[:part_cors][2]
# for ci in t_cd[2][sortperm(t_cd[1])][1:10]
#     display(B.plot_comparison_for_cell(datasets.osmfish[:df], ci, nothing, datasets.osmfish[:dapi_arr]; paper_polys=datasets.osmfish[:paper_polys], cell1_col=:cell_dapi, cell2_col=:cell, 
#             ms=4.0, bandwidth=5.0, ticks=true))
# end
```

```julia
# t_cd = datasets.allen_smfish[:part_cors][2]
# for ci in t_cd[2][sortperm(t_cd[1])][1:5]
#     display(B.plot_comparison_for_cell(datasets.allen_smfish[:df], ci, nothing, datasets.allen_smfish[:dapi_arr]; paper_polys=datasets.allen_smfish[:paper_polys], cell1_col=:cell_dapi, cell2_col=:cell, 
#             ms=4.0, bandwidth=5.0, ticks=true))
# end
```

```julia
# t_d = datasets.merfish
# t_cd = t_d[:part_cors][2]
# for ci in t_cd[2][sortperm(t_cd[1])][1:5]
#     display(B.plot_comparison_for_cell(t_d[:df], ci, nothing, t_d[:dapi_arr]; paper_polys=t_d[:paper_polys], cell1_col=:cell_dapi, cell2_col=:cell, 
#             ms=4.0, bandwidth=5.0, ticks=true))
# end
```

```julia
# t_d = datasets.starmap1020
# t_cd = t_d[:part_cors][2]
# for ci in t_cd[2][sortperm(t_cd[1])][1:5]
#     display(B.plot_comparison_for_cell(t_d[:df], ci, nothing, nothing; paper_polys=t_d[:paper_polys], cell1_col=:cell_dapi, cell2_col=:cell, 
#             ms=4.0, bandwidth=5.0, ticks=true))
# end
```

### Segmentation Prior

<!-- #region toc-hr-collapsed=true toc-nb-collapsed=true toc-hr-collapsed=true toc-nb-collapsed=true -->
#### Examples
<!-- #endregion -->

```julia
t_d = datasets.iss
B.plot_comparison_for_cell(t_d[:df], (5625, 5950), (2000, 2250), nothing, t_d[:dapi_arr]; ms=4.0, alpha=0.75,
    bandwidth=3.0, grid_step=1.0, ticks=true, size_mult=2.0, polygon_line_width=3.0, polygon_alpha=0.5)
```

```julia
t_d = datasets.iss
B.plot_comparison_for_cell(t_d[:df], (5625, 5950), (2000, 2250), nothing, t_d[:dapi_arr]; paper_polys=t_d[:paper_polys], ms=4.0, alpha=0.75,
    bandwidth=3.0, grid_step=1.0, ticks=true, size_mult=2.0, polygon_line_width=2.0, polygon_alpha=0.5)
```

#### Correlations

```julia
pdfs1 = [DataFrame(:Correlation => d[:part_cors_prior_dapi][1][1], :Type => d[:name]) for d in datasets if :part_cors_prior in keys(d)];
pdfs2 = [DataFrame(:Correlation => d[:part_cors_prior_dapi][2][1], :Type => d[:name]) for d in datasets if :part_cors_prior in keys(d)];
t_heights = (0.9 .* hcat(size.(pdfs1, 1), size.(pdfs2, 1)) ./ max.(size.(pdfs1, 1), size.(pdfs2, 1)))[sortperm([d.Type[1] for d in pdfs2]),:]

plt = @df vcat(pdfs1...) violin(:Type, :Correlation, side=:right, label="Baysor with Prior", ylabel="Correlation between parts", 
    bg_legend=Colors.RGBA(1.0, 1.0, 1.0, 0.9), bar_width=t_heights[:,1], size=(500, 240), legend=:bottomright, xgrid=false, color=color_per_col_label["Baysor with Prior"])
@df vcat(pdfs2...) violin!(:Type, :Correlation, side=:left, label="Paper", bar_width=t_heights[:,2], color=color_per_col_label["Paper"])

Plots.savefig(plt, "$PLOT_DIR/correlations_prior_paper.pdf")
plt
```

```julia
pdfs1 = [DataFrame(:Correlation => d[:part_cors_prior][1][1], :Type => d[:name]) for d in datasets if :part_cors_prior in keys(d)];
pdfs2 = [DataFrame(:Correlation => d[:part_cors_prior][2][1], :Type => d[:name]) for d in datasets if :part_cors_prior in keys(d)];
t_heights = (0.9 .* hcat(size.(pdfs1, 1), size.(pdfs2, 1)) ./ max.(size.(pdfs1, 1), size.(pdfs2, 1)))[sortperm([d.Type[1] for d in pdfs2]),:]

plt = @df vcat(pdfs1...) violin(:Type, :Correlation, side=:right, label="Baysor", ylabel="Correlation between parts", 
    bg_legend=Colors.RGBA(1.0, 1.0, 1.0, 0.9), bar_width=t_heights[:,1], size=(500, 240), legend=:bottomright, xgrid=false, color=color_per_col_label["Baysor"])
@df vcat(pdfs2...) violin!(:Type, :Correlation, side=:left, label="Baysor with Prior", bar_width=t_heights[:,2], 
    color=color_per_col_label["Baysor with Prior"])

Plots.savefig(plt, "$PLOT_DIR/correlations_prior.pdf")
plt
```

```julia
p_df = B.correlation_effect_size_df(datasets, :part_cors_prior_dapi, cell_col_labels, [4,2]);
plts1 = plotCorrelationEffect.(Ref(p_df), ["CellFrac", "MolFrac"], [0.65, 0.32]);

p_df = B.correlation_effect_size_df(datasets, :part_cors_prior, cell_col_labels, [1,4]);
plts2 = plotCorrelationEffect.(Ref(p_df), ["CellFrac", "MolFrac"], [0.65, 0.32]);

plt = R"cowplot::plot_grid"(plts1[1], plts2[1], plts1[2], plts2[2], ncol=2);

RCall.ijulia_setdevice(MIME("image/svg+xml"), width=8, height=5)
R"ggsave"("./plots/corr_effect_prior.pdf", plt, width=8, height=5);

plt
```

<!-- #region toc-hr-collapsed=true toc-nb-collapsed=true -->
#### Examples
<!-- #endregion -->

```julia
t_d = datasets.merfish
# @time t_d[:prior_polys] = B.boundary_polygons(t_d[:df], t_d[:df].cell_prior, grid_step=5.0, bandwidth=5.0)
@time t_d[:polys] = B.boundary_polygons(t_d[:df], t_d[:df].cell, grid_step=5.0, bandwidth=5.0);
```

```julia
t_cd = t_d[:part_cors_prior][2];
for ci in t_cd[2][sortperm(t_cd[1])[t_cd[3] .> 0.1]][1:5]
    display(B.plot_comparison_for_cell(t_d[:df], ci, nothing, nothing; paper_polys=t_d[:polys], cell1_col=:cell_prior, cell2_col=:cell, 
            ms=4.0, bandwidth=5.0, ticks=true, size_mult=1.5))
end
```

```julia
t_cd = t_d[:part_cors_prior][1];
for ci in t_cd[2][sortperm(t_cd[1])[t_cd[3] .> 0.1]][1:5]
    display(B.plot_comparison_for_cell(t_d[:df], ci, nothing, nothing; paper_polys=t_d[:prior_polys], cell1_col=:cell, cell2_col=:cell_prior, 
            ms=4.0, bandwidth=5.0, ticks=true, size_mult=1.5))
end
```

```julia
for ci in t_cd[2][sortperm(t_cd[1])[t_cd[3] .> 0.1]][1:5]
    display(B.plot_comparison_for_cell(t_d[:df], ci, nothing, nothing; paper_polys=t_d[:prior_polys], cell1_col=:cell_prior, cell2_col=:cell, 
            ms=4.0, bandwidth=5.0, ticks=true))
end
```

### Supp. plots


#### Cell stats

```julia
p_df = vcat(vcat([[@where(DataFrame(:Val => d[:qc_per_cell_dfs][i].n_transcripts, :Dataset => d[:name], :Segmentation => n), :Val .< quantile(:Val, 0.995)) 
        for d in datasets if size(d[:qc_per_cell_dfs][i], 1) > 0] for (i,n) in enumerate(cell_col_labels)]...)...);

R"""
gg <- ggplot($p_df, aes(x=Dataset, y=Val, fill=Segmentation)) + 
    geom_boxplot_jitter(outlier.jitter.width=0.05, outlier.size=0.5, outlier.alpha=0.25) +
    scale_y_continuous(minor_breaks=c(seq(1, 10, length.out=5), seq(10, 100, length.out=5), seq(100, 1000, length.out=5)), trans="log10", name='Num. of transcripts') +
    xlab("") +
    theme(legend.position=c(1, 0.05), legend.justification=c(1, 0), legend.background=element_rect(fill=alpha('white', 0.9)), 
        axis.text=element_text(size=9), axis.title.y=element_text(size=12), legend.text=element_text(size=10), legend.title=element_text(size=12))

ggsave('./plots/num_transcripts.pdf', width=4, height=4)
gg
"""
```

```julia
p_df = vcat(vcat([[DataFrame(:Val => d[:qc_per_cell_dfs][i].sqr_area, :Dataset => d[:name], :Segmentation => n) 
        for d in datasets if size(d[:qc_per_cell_dfs][i], 1) > 0] for (i,n) in enumerate(cell_col_labels)]...)...);

R"""
gg <- ggplot($p_df, aes(x=Dataset, y=Val, fill=Segmentation)) + 
    geom_boxplot_jitter(outlier.jitter.width=0.05, outlier.size=0.5, outlier.alpha=0.25) +
    scale_y_continuous(limits=c(0, 400), expand=c(0, 0), name='sqrt(Area)') +
    xlab("") +
    theme(legend.position=c(0.01, 0.99), legend.justification=c(0, 1), legend.background=element_rect(fill=alpha('white', 0.9)),
        axis.text=element_text(size=9), axis.title.y=element_text(size=12), legend.text=element_text(size=10), legend.title=element_text(size=12))

ggsave('./plots/area.pdf', width=4, height=4)
gg
"""
```

#### Similarity vs Prior Confidence 

```julia
subdirs = filter(x -> startswith(x, "ca1"), readdir("$PROJECT_DIR/run_results/iss_hippo/"));

dfs = [B.load_df("$PROJECT_DIR/run_results/iss_hippo/$(sd)/segmentation.csv")[1] for sd in subdirs];

dir_aliases = Dict(
  "ca1_paper_prior"     => "0.5",
  "ca1_paper_prior_025" => "0.25",
  "ca1_paper_prior_09"  => "0.9",
  "ca1_paper_prior_075" => "0.75",
  "ca1_paper_prior_1"   => "1.0",
  "ca1_no_prior"        => "0.0"
);
```

```julia
p_df = @orderby(vcat([DataFrame(:Val => Clustering.mutualinfo(d.cell, d.parent_id), :Type => dir_aliases[n]) for (n,d) in zip(subdirs, dfs)]...), :Type);
plt = Plots.plot(parse.(Float16, p_df.Type), p_df.Val, ylabel = "Mutual Information", legend=false, xlabel="Prior confidence", marker=true, lw=3, size=(370, 300))

Plots.savefig(plt, "./plots/mutual_info_per_confidence.pdf")
plt
```

```julia
# @orderby(vcat([DataFrame(:Val => Clustering.varinfo(d.cell, d.parent_id), :Type => dir_aliases[n]) for (n,d) in zip(subdirs, dfs)]...), :Type)
```

```julia
# t_cnts = Dict(n => counts(d.cell, d.parent_id) for (n,d) in zip(subdirs, dfs));
t_matches = [B.match_assignments(d.cell, d.parent_id) for d in dfs];
```

```julia
min_size = 1
p_df = vcat([DataFrame(:Frac => tm.max_overlaps[1][vec(sum(tm.contingency, dims=2))[2:end] .>= min_size], :Confidence => dir_aliases[n], :Type => "Baysor cells") for (n,tm) in zip(subdirs, t_matches)]...);
p_df = vcat([DataFrame(:Frac => tm.max_overlaps[2][vec(sum(tm.contingency, dims=1))[2:end] .>= min_size], :Confidence => dir_aliases[n], :Type => "Paper cells") for (n,tm) in zip(subdirs, t_matches)]..., p_df);
p_df = @orderby(p_df, :Confidence);

R"""
gg <- ggplot($p_df, aes(x=Confidence, y=Frac)) + 
    geom_boxplot_jitter(aes(fill=Type), outlier.jitter.width=0.05, outlier.size=0.2, outlier.alpha=0.1) +
    theme(legend.position=c(1, 0.05), legend.justification=c(1, 0), legend.background=element_rect(fill=alpha('white', 0.2))) +
    scale_y_continuous(limits=c(0.0, 1.01), expand=c(0, 0), name="Overlap fraction with the best matching cell")

ggsave('./plots/prior_confidence_iss.pdf', width=6, height=4)
gg
"""
```
