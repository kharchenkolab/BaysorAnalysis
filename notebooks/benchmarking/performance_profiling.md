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

```julia tags=[]
using DrWatson
quickactivate(@__DIR__)

import Baysor as B
import BaysorAnalysis as BA
import Colors
import CSV
import MultivariateStats
import Plots

using BenchmarkTools
using ProgressMeter
using OrderedCollections
using DataFrames
using DataFramesMeta
using NearestNeighbors
using Statistics
using StatsBase
```

```julia tags=[]
BenchmarkTools.DEFAULT_PARAMETERS.gcsample = true;
BenchmarkTools.DEFAULT_PARAMETERS.overhead = BenchmarkTools.estimate_overhead();
BenchmarkTools.DEFAULT_PARAMETERS.samples = 5;
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 30 * 60;
```

<!-- #region heading_collapsed="true" tags=[] -->
## Load data
<!-- #endregion -->

```julia
# @time df_spatial, gene_names = Baysor.load_df("../run_results/spacejam2/allen_sm_fish/no_dapi/segmentation.csv");
# df_spatial[!, :x] = round.(Int, 10 .* (df_spatial.x .- minimum(df_spatial.x)));
# df_spatial[!, :y] = round.(Int, 10 .* (df_spatial.y .- minimum(df_spatial.y)));
# length(gene_names)
```

```julia
@time df_spatial, gene_names = Baysor.load_df("../run_results/merfish_moffit/segmentation.csv");
length(gene_names)
```

<!-- #region heading_collapsed="true" tags=[] -->
## Molecule clustering
<!-- #endregion -->

### Baysor

```julia
bench_df = @where(df_spatial, :x .< -3300, :y .< -3300) |> deepcopy;
gn_bench = gene_names;

# confidence_nn_id = Baysor.default_param_value(:confidence_nn_id, 10);
confidence_nn_id = Baysor.default_param_value(:confidence_nn_id, 50);
@show confidence_nn_id
size(bench_df, 1)
```

```julia
bench_clust = BenchmarkGroup();
```

```julia
B.append_confidence!(bench_df, nn_id=confidence_nn_id);
bench_clust["confidence"] = @benchmarkable B.append_confidence!($bench_df, nn_id=$confidence_nn_id);
```

```julia
adjacent_points, adjacent_weights = B.build_molecule_graph(bench_df, filter=false);
bench_clust["mol_graph"] = @benchmarkable B.build_molecule_graph($bench_df, filter=false);
```

```julia
for cl in [2, 4, 6, 8, 10]
    bench_clust["clust_$cl"] = @benchmarkable B.cluster_molecules_on_mrf($bench_df.gene, $adjacent_points, $adjacent_weights, $bench_df.confidence; 
        n_clusters=$cl, max_iters=5000, n_iters_without_update=100, verbose=false);
end
```

```julia
bench_clust_res = run(bench_clust)
```

```julia
bench_res_df = vcat([DataFrame("Key" => k, "Mean time, sec" => mean(v.times) ./ 1e9, "Std time, sec" => std(v.times) ./ 1e9, 
            "Num. samples" => length(v.times)) for (k,v) in bench_clust_res]...)
```

### Leiden

```julia
using RCall
```

```julia
nm_bench = B.neighborhood_count_matrix(bench_df, 50, normalize=false);
size(nm_bench)
```

```julia
R"""
library(pagoda2)
library(conos)
library(microbenchmark)

cm <- as($nm_bench, "dgCMatrix")
rownames(cm) <- $gn_bench
colnames(cm) <- paste0("c", 1:ncol(cm))

getClusters <- function(cm, verbose=FALSE) {
    p2 <- Pagoda2$new(cm, trim=5, n.cores=1, verbose=FALSE, log.scale=verbose)
    p2$calculatePcaReduction(nPcs=50, odgenes=rownames(cm), maxit=1000, verbose=verbose, var.scale=FALSE)
    p2$makeKnnGraph(k=30, type="PCA", center=T, distance="cosine", weight.type="none", verbose=verbose)
    p2$getKnnClusters(method=conos::leiden.community, type="PCA", name="leiden", resolution=1.0)
    
    return(p2$clusters$PCA$leiden)
}

b <- microbenchmark(
    "clustering" = {getClusters(cm)},
    times=5,
    control=list(warmup=1)
)
"""
```

### Aggregate

```julia
leiden_times = rcopy(R"b").time;
```

```julia
bench_res_df
```

```julia
df1 = hcat(DataFrame("Method" => "MRF", "Num. clusters" => 2:2:10), bench_res_df[[3, 1, 5, 4, 2],2:end]);
df2 = vcat(df1, DataFrame("Method" => "Leiden", "Num. clusters" => "Any", "Mean time, sec" => mean(leiden_times) / 1e9, 
        "Std time, sec" => std(leiden_times) / 1e9, "Num. samples" => 5));

df2[:, 3:4] .= round.(df2[:, 3:4], digits=2);
df2
```

```julia
CSV.write("plots/clustering_profiling.csv", df2)
```

<!-- #region heading_collapsed="true" tags=[] -->
## Color embedding
<!-- #endregion -->

```julia
@time neighb_cm = B.neighborhood_count_matrix(df_spatial, 40);
@time color_transformation = B.gene_composition_transformation(neighb_cm, df_spatial.confidence; sample_size=20000, spread=2.0, min_dist=0.1);
@time color_emb = B.transform(color_transformation, neighb_cm);
```

```julia
bench_emb = BenchmarkGroup();
```

```julia
bench_emb["neighborhood_count_matrix_40"] = @benchmarkable B.neighborhood_count_matrix($df_spatial, 40)
```

```julia
bench_emb["gene_composition_transformation_20k"] = @benchmarkable B.gene_composition_transformation(neighb_cm, df_spatial.confidence; 
    sample_size=20000, spread=2.0, min_dist=0.1)
```

```julia
bench_emb["transform"] = @benchmarkable B.transform(color_transformation, neighb_cm)
```

```julia
bench_emb_res = run(bench_emb)
```

```julia
bench_df = vcat([DataFrame("Key" => k, "Mean time, sec" => mean(v.times) ./ 1e9, "Std time, sec" => std(v.times) ./ 1e9, 
            "Num. samples" => length(v.times)) for (k,v) in bench_emb_res]...)
```

<!-- #region heading_collapsed="true" tags=[] -->
## Segmentation
<!-- #endregion -->

```julia
bench_segmentation = BenchmarkGroup();
```

```julia
@time df_spatial, gene_names = B.load_df("../run_results/iss_hippo/ca1_no_prior/segmentation.csv");
df_spatial[!, :cell_dapi] = df_spatial.parent_id;
dapi_arr = Float16.(Images.load("/home/vpetukhov/data/spatal/iss/hippocampus/CA1/Viktor/CA1DapiBoundaries_4-3_right.tif"));
iss = Dict(:df => df_spatial, :gene_names => gene_names, :name => "ISS", :dapi_arr => dapi_arr);
```

```julia
B.append_confidence!(df_spatial, (args["prior_segmentation"]===nothing ? nothing : df_spatial.prior_segmentation), nn_id=confidence_nn_id, prior_confidence=args["prior-segmentation-confidence"])
adjacent_points, adjacent_weights = build_molecule_graph(df_spatial, filter=false)[1:2];

mol_clusts = cluster_molecules_on_mrf(df_spatial.gene, adjacent_points, adjacent_weights, df_spatial.confidence;
            n_clusters=args["n-clusters"], weights_pre_adjusted=true)

df_spatial[!, :cluster] = mol_clusts.assignment;

bm_data_arr = initial_distribution_arr(df_spatial; n_frames=args["n-frames"], scale=args["scale"], scale_std=args["scale-std"],
        n_cells_init=args["num-cells-init"], prior_seg_confidence=args["prior-segmentation-confidence"],
        min_molecules_per_cell=args["min-molecules-per-cell"], confidence_nn_id=0);

bm_data = run_bmm_parallel!(bm_data_arr, args["iters"], new_component_frac=args["new-component-fraction"], new_component_weight=args["new-component-weight"],
                            min_molecules_per_cell=args["min-molecules-per-cell"], assignment_history_depth=history_depth);
```

```julia
cur_df = deepcopy(iss[:df]);
bm_data = B.initial_distribution_arr(cur_df; n_frames=1, scale=14, scale_std="25%", min_molecules_per_cell=3)[1];
@time B.bmm!(bm_data, n_iters=350, new_component_frac=0.3, min_molecules_per_cell=3, assignment_history_depth=30, log_step=100);
cur_df[!, :cell] = B.estimate_assignment_by_history(bm_data)[1];

B.plot_comparison_for_cell(cur_df, B.val_range(cur_df.x), B.val_range(cur_df.y), nothing, iss[:dapi_arr];
    ms=2.0, bandwidth=5.0, size_mult=0.25, plot_raw_dapi=false)
```

<!-- #region heading_collapsed="true" tags=[] -->
## Full run
<!-- #endregion -->

<!-- #region heading_collapsed="true" tags=[] -->
### Run
<!-- #endregion -->

```julia
using ProgressMeter
```

```julia
dataset_paths = "/home/vpetukhov/spatial/Benchmarking/run_results/" .* 
    ["iss_hippo/ca1_no_prior", "merfish_moffit", "osm_fish", "star_map/vis_1020_cl0", "spacejam2/allen_sm_fish/no_dapi"];
param_dumps = dataset_paths .* "/segmentation_params.dump";

dataset_names = ["iss", "merfish", "osm_fish", "starmap_1020", "allen_smfish"];

param_strings = [open(p) do f readlines(f)[1][16:end-1] end for p in param_dumps];
```

```julia
baysor_path = "/home/vpetukhov/local/bin/baysor";
for i in 2:length(param_strings)
# for i in 2:2
    dataset = dataset_names[i]
    params = split(param_strings[i], ' ')

    out_path = expanduser("/home/vpetukhov/spatial/Benchmarking/run_results/profiling/$dataset/")
    mkpath(out_path)
    cmd = `/usr/bin/time -f '%e %U %P %M %t %K' -o ./profiling_output/$dataset.prof -a $baysor_path run --debug -o $out_path $params`;
#     cmd = `/usr/bin/time -f '%e %U %P %M %t %K' -o ./profiling_output/$dataset.prof -a $baysor_path run --debug --n-clusters=0 -o $out_path $params`;
    @show cmd
    
    println(dataset)
    @showprogress for ri in 1:5
        run(pipeline(cmd, stdout="./profiling_output/$dataset.log", stderr="./profiling_output/$dataset.err", append=true))
        run(pipeline(`echo -e \\n\\n\\n ----- RUN $ri ----- \\n\\n\\n`, stdout="./profiling_output/$dataset.log", append=true))
    end
end
```

<!-- #region heading_collapsed="true" tags=[] -->
### Summarize
<!-- #endregion -->

```julia
using DataFrames
using Statistics

printed_names = ["ISS", "MERFISH", "osmFISH", "STARmap 1020", "Allen smFISH"];

seg_results = dataset_paths .* "/segmentation.csv";
dataset_parameters = hcat([[size(df, 1), length(unique(df.gene))] for df in DataFrame!.(CSV.File.(seg_results))]...);
```

```julia
bench_vals = [hcat(split.(readlines("./profiling_output/$ds.prof"), ' ')...) for ds in dataset_names];
mem_vals = hcat([parse.(Float64, x[4,:]) / 1e6 for x in bench_vals]...);
cpu_vals = hcat([parse.(Float64, x[1,:]) / 60 for x in bench_vals]...);

bench_mat = round.(vcat(mean(cpu_vals, dims=1), std(cpu_vals, dims=1), mean(mem_vals, dims=1), std(mem_vals, dims=1))', digits=2);
bench_strs = [["$(r[i[1]]) ± $(r[i[2]])" for r in eachrow(bench_mat)] for i in ((1, 2), (3, 4))];
bench_df = DataFrame("Dataset" => printed_names, "Num. molecules" => dataset_parameters[1,:], "Num. genes" => dataset_parameters[2,:],
    "CPU time, min" => bench_strs[1], "Max RSS, GB" => bench_strs[2], "Num. samples" => 5)
```

```julia
CSV.write("./plots/segmentation_profiling.csv", bench_df)
```

## Parameter table

```julia tags=[]
import Pkg: TOML
using DataFrames
import CSV
```

```julia
data_paths = ["Allen smFISH" => "allen_smfish", "ISS" => "iss_hippo", "osmFISH" => "osmfish", "STARmap 1020" => "starmap_vis1020", "MERFISH Hypothalamus" => "merfish_moffit", "MERFISH Gut" => "merfish_membrane"];
prior_subfolders = ["No" => "baysor", "Paper" => "baysor_prior", "DAPI" => "baysor_dapi_prior", "Membrane" => "baysor_membrane_prior"];

p_keys = ["gene-composition-neigborhood", "scale", "prior-segmentation-confidence", "min-molecules-per-gene", "min-molecules-per-cell", "n-clusters", 
    "iters", "force-2d", "x-column", "y-column", "z-column", "gene-column", "prior_segmentation", "nuclei-genes", "cyto-genes"];
```

```julia
path_df = DataFrame([Dict(:Dataset => d, :Prior => pr, :Path => datadir("exp_pro", md, sd, "segmentation_params.dump")) for (d, md) in data_paths for (pr, sd) in prior_subfolders]);
path_df = path_df[isfile.(path_df.Path),:];

param_dicts = [OrderedDict(k => get(d, k, "NA") for k in p_keys) for d in TOML.parsefile.(path_df.Path)];
param_df = hcat(path_df[:,[:Dataset, :Prior]], vcat(DataFrame.(param_dicts)...))
CSV.write(plotsdir("parameters.csv"), param_df)
```
