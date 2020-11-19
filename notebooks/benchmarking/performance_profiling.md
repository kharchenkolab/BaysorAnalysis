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

```julia execution={"iopub.execute_input": "2020-09-30T18:33:32.842000+02:00", "iopub.status.busy": "2020-09-30T18:33:32.842000+02:00", "iopub.status.idle": "2020-09-30T18:33:32.892000+02:00"}
import Baysor
import Colors
import CSV
import MultivariateStats
import Plots

using BenchmarkTools
using ProgressMeter
using DataFrames
using DataFramesMeta
using NearestNeighbors
using Statistics
using StatsBase

B = Baysor;
```

```julia execution={"iopub.execute_input": "2020-09-30T16:01:33.073000+02:00", "iopub.status.busy": "2020-09-30T16:01:33.073000+02:00", "iopub.status.idle": "2020-09-30T16:01:33.179000+02:00"}
BenchmarkTools.DEFAULT_PARAMETERS.gcsample = true;
BenchmarkTools.DEFAULT_PARAMETERS.overhead = BenchmarkTools.estimate_overhead();
BenchmarkTools.DEFAULT_PARAMETERS.samples = 5;
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 30 * 60;
```

## Load data

```julia execution={"iopub.execute_input": "2020-09-25T19:13:42.517000+02:00", "iopub.status.busy": "2020-09-25T19:13:42.517000+02:00", "iopub.status.idle": "2020-09-25T19:13:51.412000+02:00"}
# @time df_spatial, gene_names = Baysor.load_df("../run_results/spacejam2/allen_sm_fish/no_dapi/segmentation.csv");
# df_spatial[!, :x] = round.(Int, 10 .* (df_spatial.x .- minimum(df_spatial.x)));
# df_spatial[!, :y] = round.(Int, 10 .* (df_spatial.y .- minimum(df_spatial.y)));
# length(gene_names)
```

```julia execution={"iopub.execute_input": "2020-09-30T16:01:11.758000+02:00", "iopub.status.busy": "2020-09-30T16:01:11.758000+02:00", "iopub.status.idle": "2020-09-30T16:01:22.250000+02:00"}
@time df_spatial, gene_names = Baysor.load_df("../run_results/merfish_moffit/segmentation.csv");
length(gene_names)
```

## Molecule clustering


### Baysor

```julia execution={"iopub.execute_input": "2020-09-30T16:03:23.529000+02:00", "iopub.status.busy": "2020-09-30T16:03:23.529000+02:00", "iopub.status.idle": "2020-09-30T16:03:23.667000+02:00"}
bench_df = @where(df_spatial, :x .< -3300, :y .< -3300) |> deepcopy;
gn_bench = gene_names;

# confidence_nn_id = Baysor.default_param_value(:confidence_nn_id, 10);
confidence_nn_id = Baysor.default_param_value(:confidence_nn_id, 50);
@show confidence_nn_id
size(bench_df, 1)
```

```julia execution={"iopub.execute_input": "2020-09-30T16:03:26.244000+02:00", "iopub.status.busy": "2020-09-30T16:03:26.244000+02:00", "iopub.status.idle": "2020-09-30T16:03:26.270000+02:00"}
bench_clust = BenchmarkGroup();
```

```julia execution={"iopub.execute_input": "2020-09-30T16:03:28.730000+02:00", "iopub.status.busy": "2020-09-30T16:03:28.730000+02:00", "iopub.status.idle": "2020-09-30T16:03:42.956000+02:00"}
B.append_confidence!(bench_df, nn_id=confidence_nn_id);
bench_clust["confidence"] = @benchmarkable B.append_confidence!($bench_df, nn_id=$confidence_nn_id);
```

```julia execution={"iopub.execute_input": "2020-09-30T16:03:44.226000+02:00", "iopub.status.busy": "2020-09-30T16:03:44.225000+02:00", "iopub.status.idle": "2020-09-30T16:03:47.559000+02:00"}
adjacent_points, adjacent_weights = B.build_molecule_graph(bench_df, filter=false);
bench_clust["mol_graph"] = @benchmarkable B.build_molecule_graph($bench_df, filter=false);
```

```julia execution={"iopub.execute_input": "2020-09-30T16:03:48.816000+02:00", "iopub.status.busy": "2020-09-30T16:03:48.816000+02:00", "iopub.status.idle": "2020-09-30T16:03:49.105000+02:00"}
for cl in [2, 4, 6, 8, 10]
    bench_clust["clust_$cl"] = @benchmarkable B.cluster_molecules_on_mrf($bench_df.gene, $adjacent_points, $adjacent_weights, $bench_df.confidence; 
        n_clusters=$cl, max_iters=5000, n_iters_without_update=100, verbose=false);
end
```

```julia execution={"iopub.execute_input": "2020-09-30T16:03:55.130000+02:00", "iopub.status.busy": "2020-09-30T16:03:55.129000+02:00", "iopub.status.idle": "2020-09-30T16:45:10.393000+02:00"}
bench_clust_res = run(bench_clust)
```

```julia execution={"iopub.execute_input": "2020-09-30T16:45:10.393000+02:00", "iopub.status.busy": "2020-09-30T16:45:10.393000+02:00", "iopub.status.idle": "2020-09-30T16:45:13.603000+02:00"}
bench_res_df = vcat([DataFrame("Key" => k, "Mean time, sec" => mean(v.times) ./ 1e9, "Std time, sec" => std(v.times) ./ 1e9, 
            "Num. samples" => length(v.times)) for (k,v) in bench_clust_res]...)
```

### Leiden

```julia execution={"iopub.execute_input": "2020-09-30T16:48:24.599000+02:00", "iopub.status.busy": "2020-09-30T16:48:24.599000+02:00", "iopub.status.idle": "2020-09-30T16:48:28.032000+02:00"}
using RCall
```

```julia execution={"iopub.execute_input": "2020-09-30T16:48:36.046000+02:00", "iopub.status.busy": "2020-09-30T16:48:35.550000+02:00", "iopub.status.idle": "2020-09-30T16:48:38.775000+02:00"}
nm_bench = B.neighborhood_count_matrix(bench_df, 50, normalize=false);
size(nm_bench)
```

```julia execution={"iopub.execute_input": "2020-09-30T16:48:44.712000+02:00", "iopub.status.busy": "2020-09-30T16:48:44.711000+02:00", "iopub.status.idle": "2020-09-30T17:12:26.341000+02:00"}
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

```julia execution={"iopub.execute_input": "2020-09-30T17:12:26.341000+02:00", "iopub.status.busy": "2020-09-30T17:12:26.341000+02:00", "iopub.status.idle": "2020-09-30T17:12:26.870000+02:00"}
leiden_times = rcopy(R"b").time;
```

```julia execution={"iopub.execute_input": "2020-09-30T17:12:26.870000+02:00", "iopub.status.busy": "2020-09-30T17:12:26.870000+02:00", "iopub.status.idle": "2020-09-30T17:12:27.046000+02:00"}
bench_res_df
```

```julia execution={"iopub.execute_input": "2020-09-30T17:12:27.046000+02:00", "iopub.status.busy": "2020-09-30T17:12:27.046000+02:00", "iopub.status.idle": "2020-09-30T17:12:28.091000+02:00"}
df1 = hcat(DataFrame("Method" => "MRF", "Num. clusters" => 2:2:10), bench_res_df[[3, 1, 5, 4, 2],2:end]);
df2 = vcat(df1, DataFrame("Method" => "Leiden", "Num. clusters" => "Any", "Mean time, sec" => mean(leiden_times) / 1e9, 
        "Std time, sec" => std(leiden_times) / 1e9, "Num. samples" => 5));

df2[:, 3:4] .= round.(df2[:, 3:4], digits=2);
df2
```

```julia execution={"iopub.execute_input": "2020-09-30T18:33:38.579000+02:00", "iopub.status.busy": "2020-09-30T18:33:38.578000+02:00", "iopub.status.idle": "2020-09-30T18:33:39.563000+02:00"}
CSV.write("plots/clustering_profiling.csv", df2)
```

## Color embedding

```julia execution={"iopub.execute_input": "2020-09-25T19:15:31.459000+02:00", "iopub.status.busy": "2020-09-25T19:15:31.459000+02:00", "iopub.status.idle": "2020-09-25T19:17:12.818000+02:00"}
@time neighb_cm = B.neighborhood_count_matrix(df_spatial, 40);
@time color_transformation = B.gene_composition_transformation(neighb_cm, df_spatial.confidence; sample_size=20000, spread=2.0, min_dist=0.1);
@time color_emb = B.transform(color_transformation, neighb_cm);
```

```julia execution={"iopub.execute_input": "2020-09-25T19:17:12.818000+02:00", "iopub.status.busy": "2020-09-25T19:17:12.818000+02:00", "iopub.status.idle": "2020-09-25T19:17:12.845000+02:00"}
bench_emb = BenchmarkGroup();
```

```julia execution={"iopub.execute_input": "2020-09-25T19:17:12.845000+02:00", "iopub.status.busy": "2020-09-25T19:17:12.845000+02:00", "iopub.status.idle": "2020-09-25T19:17:13.396000+02:00"}
bench_emb["neighborhood_count_matrix_40"] = @benchmarkable B.neighborhood_count_matrix($df_spatial, 40)
```

```julia execution={"iopub.execute_input": "2020-09-25T19:17:13.396000+02:00", "iopub.status.busy": "2020-09-25T19:17:13.396000+02:00", "iopub.status.idle": "2020-09-25T19:17:13.734000+02:00"}
bench_emb["gene_composition_transformation_20k"] = @benchmarkable B.gene_composition_transformation(neighb_cm, df_spatial.confidence; 
    sample_size=20000, spread=2.0, min_dist=0.1)
```

```julia execution={"iopub.execute_input": "2020-09-25T19:17:13.735000+02:00", "iopub.status.busy": "2020-09-25T19:17:13.735000+02:00", "iopub.status.idle": "2020-09-25T19:17:13.967000+02:00"}
bench_emb["transform"] = @benchmarkable B.transform(color_transformation, neighb_cm)
```

```julia execution={"iopub.execute_input": "2020-09-25T19:17:13.967000+02:00", "iopub.status.busy": "2020-09-25T19:17:13.967000+02:00", "iopub.status.idle": "2020-09-25T19:25:23.992000+02:00"}
bench_emb_res = run(bench_emb)
```

```julia execution={"iopub.execute_input": "2020-09-25T19:41:41.973000+02:00", "iopub.status.busy": "2020-09-25T19:41:41.973000+02:00", "iopub.status.idle": "2020-09-25T19:41:42.607000+02:00"}
bench_df = vcat([DataFrame("Key" => k, "Mean time, sec" => mean(v.times) ./ 1e9, "Std time, sec" => std(v.times) ./ 1e9, 
            "Num. samples" => length(v.times)) for (k,v) in bench_emb_res]...)
```

## Segmentation

```julia execution={"iopub.execute_input": "2020-09-10T13:57:52.073000+02:00", "iopub.status.busy": "2020-09-10T13:57:52.073000+02:00", "iopub.status.idle": "2020-09-10T13:57:52.106000+02:00"}
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

## Full run


### Run

```julia execution={"iopub.execute_input": "2020-09-28T10:27:13.345000+02:00", "iopub.status.busy": "2020-09-28T10:27:13.345000+02:00", "iopub.status.idle": "2020-09-28T10:27:13.369000+02:00"}
using ProgressMeter
```

```julia execution={"iopub.execute_input": "2020-09-28T10:28:06.250000+02:00", "iopub.status.busy": "2020-09-28T10:28:06.249000+02:00", "iopub.status.idle": "2020-09-28T10:28:06.390000+02:00"}
dataset_paths = "/home/vpetukhov/spatial/Benchmarking/run_results/" .* 
    ["iss_hippo/ca1_no_prior", "merfish_moffit", "osm_fish", "star_map/vis_1020_cl0", "spacejam2/allen_sm_fish/no_dapi"];
param_dumps = dataset_paths .* "/segmentation_params.dump";

dataset_names = ["iss", "merfish", "osm_fish", "starmap_1020", "allen_smfish"];

param_strings = [open(p) do f readlines(f)[1][16:end-1] end for p in param_dumps];
```

```julia execution={"iopub.execute_input": "2020-09-21T22:18:50.630000+02:00", "iopub.status.busy": "2020-09-21T22:18:50.630000+02:00", "iopub.status.idle": "2020-09-22T08:31:29.601000+02:00"}
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

### Summarize

```julia execution={"iopub.execute_input": "2020-09-28T10:29:43.277000+02:00", "iopub.status.busy": "2020-09-28T10:29:43.277000+02:00", "iopub.status.idle": "2020-09-28T10:29:46.463000+02:00"}
using DataFrames
using Statistics

printed_names = ["ISS", "MERFISH", "osmFISH", "STARmap 1020", "Allen smFISH"];

seg_results = dataset_paths .* "/segmentation.csv";
dataset_parameters = hcat([[size(df, 1), length(unique(df.gene))] for df in DataFrame!.(CSV.File.(seg_results))]...);
```

```julia execution={"iopub.execute_input": "2020-09-28T10:30:19.418000+02:00", "iopub.status.busy": "2020-09-28T10:30:19.418000+02:00", "iopub.status.idle": "2020-09-28T10:30:19.696000+02:00"}
bench_vals = [hcat(split.(readlines("./profiling_output/$ds.prof"), ' ')...) for ds in dataset_names];
mem_vals = hcat([parse.(Float64, x[4,:]) / 1e6 for x in bench_vals]...);
cpu_vals = hcat([parse.(Float64, x[1,:]) / 60 for x in bench_vals]...);

bench_mat = round.(vcat(mean(cpu_vals, dims=1), std(cpu_vals, dims=1), mean(mem_vals, dims=1), std(mem_vals, dims=1))', digits=2);
bench_strs = [["$(r[i[1]]) ± $(r[i[2]])" for r in eachrow(bench_mat)] for i in ((1, 2), (3, 4))];
bench_df = DataFrame("Dataset" => printed_names, "Num. molecules" => dataset_parameters[1,:], "Num. genes" => dataset_parameters[2,:],
    "CPU time, min" => bench_strs[1], "Max RSS, GB" => bench_strs[2], "Num. samples" => 5)
```

```julia execution={"iopub.execute_input": "2020-09-28T10:30:40.543000+02:00", "iopub.status.busy": "2020-09-28T10:30:40.543000+02:00", "iopub.status.idle": "2020-09-28T10:30:40.623000+02:00"}
CSV.write("./plots/segmentation_profiling.csv", bench_df)
```

## Parameter table

```julia execution={"iopub.execute_input": "2020-09-25T18:48:41.034000+02:00", "iopub.status.busy": "2020-09-25T18:48:41.034000+02:00", "iopub.status.idle": "2020-09-25T18:48:42.505000+02:00"}
import Pkg: TOML
using DataFrames
import CSV
```

```julia execution={"iopub.execute_input": "2020-09-25T18:48:42.645000+02:00", "iopub.status.busy": "2020-09-25T18:48:42.505000+02:00", "iopub.status.idle": "2020-09-25T18:48:43.059000+02:00"}
param_dumps = "/home/vpetukhov/spatial/Benchmarking/run_results/" .* [
    "merfish_moffit", "merfish_moffit_prior", 
    "osm_fish", "osm_fish/paper_prior",
    "star_map/vis_1020_cl0", "star_map/vis_1020_prior_cl0",
    "spacejam2/allen_sm_fish/no_dapi", "spacejam2/allen_sm_fish/mask_prior",
    "iss_hippo/ca1_no_prior", "iss_hippo/ca1_paper_prior"
    ] .* "/segmentation_params.dump";

dataset_names = repeat(["MERFISH", "osmFISH", "STARmap 1020", "Allen smFISH", "ISS"], inner=2);

p_keys = ["gene-composition-neigborhood", "scale", "prior-segmentation-confidence", 
    "min-molecules-per-gene", "min-molecules-per-cell", "n-clusters"];
```

```julia execution={"iopub.execute_input": "2020-09-25T18:48:43.059000+02:00", "iopub.status.busy": "2020-09-25T18:48:43.059000+02:00", "iopub.status.idle": "2020-09-25T18:48:45.887000+02:00"}
param_dicts = [Dict(k => get(d, k, "NA") for k in p_keys) for d in TOML.parsefile.(param_dumps)]
param_df = hcat(DataFrame("Dataset" => dataset_names, "Prior" => repeat(["No", "Yes"], outer=5)), vcat(DataFrame.(param_dicts)...));

CSV.write("plots/parameters.csv", param_df)
```
