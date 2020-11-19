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

```julia execution={"iopub.status.busy": "2020-09-08T10:05:38.694000+02:00", "iopub.execute_input": "2020-09-08T10:05:40.514000+02:00", "iopub.status.idle": "2020-09-08T10:09:09.745000+02:00"}
import Baysor
import Colors
import Plots

using Distances
using ProgressMeter
using DataFrames
using DataFramesMeta
using NearestNeighbors
using RCall
using Statistics
using StatsBase

B = Baysor;

R"""
library(ggplot2)
library(ggrastr)
library(tidyverse)
theme_set(theme_bw())
""";
```

## Load data

```julia execution={"iopub.status.busy": "2020-09-08T10:09:09.745000+02:00", "iopub.execute_input": "2020-09-08T10:09:12.882000+02:00", "iopub.status.idle": "2020-09-08T10:09:24.387000+02:00"}
@time df_spatial, gene_names = Baysor.load_df("../run_results/osm_fish/segmentation.csv");
df_spatial[!, :x] = round.(Int, 10 .* (df_spatial.x .- minimum(df_spatial.x)));
df_spatial[!, :y] = round.(Int, 10 .* (df_spatial.y .- minimum(df_spatial.y)));
length(gene_names)
```

```julia
@time neighb_cm = Baysor.neighborhood_count_matrix(df_spatial, 40);
@time color_transformation = Baysor.gene_composition_transformation(neighb_cm, df_spatial.confidence; sample_size=20000, spread=2.0, min_dist=0.1);
@time color_emb = B.transform(color_transformation, neighb_cm);
```

```julia
@time color_emb = B.transform(color_transformation, neighb_cm);
```

```julia
@time gene_colors = B.gene_composition_colors(color_emb, lrange=(10, 60));
df_spatial[!, :color] = gene_colors;
```

## LEV correlation depending on k

```julia execution={"iopub.status.busy": "2020-09-08T10:11:49.444000+02:00", "iopub.execute_input": "2020-09-08T10:11:49.445000+02:00", "iopub.status.idle": "2020-09-08T10:12:10.705000+02:00"}
@time neighb_cm = Baysor.neighborhood_count_matrix(df_spatial, 40, normalize=false);
```

```julia execution={"iopub.status.busy": "2020-09-08T10:14:06.899000+02:00", "iopub.execute_input": "2020-09-08T10:14:06.899000+02:00", "iopub.status.idle": "2020-09-08T10:14:10.324000+02:00"}
pos_data = B.position_data(df_spatial);
src_ids = B.select_ids_uniformly(pos_data', n=2000);

@time knn_info = knn(KDTree(pos_data), pos_data[:, src_ids], 1701, true);
knn_info = [(c[1:10:end], d[1:10:end]) for (c,d) in zip(knn_info...)];

plot_df = map(1:length(knn_info)) do i
    DataFrame(:k => 10 .* ((1:length(knn_info[i][1])) .- 1), :d => knn_info[i][2], 
        :Correlation => vec(cor(neighb_cm[:, src_ids[i]], neighb_cm[:, knn_info[i][1]])))
end;

plot_df = vcat(plot_df...);

R"""
# p_df <- $plot_df %>% group_by(k) %>% summarize(m=mean(Correlation), std=sd(Correlation), d=median(d)) %>% mutate(mu=m+std, ml=m-std)
p_df <- $plot_df %>% group_by(k) %>% summarize(m=mean(Correlation), mu=quantile(Correlation, 0.75), ml=quantile(Correlation, 0.25), d=mean(d))
""";
```

```julia execution={"iopub.status.busy": "2020-09-08T10:19:22.053000+02:00", "iopub.execute_input": "2020-09-08T10:19:22.053000+02:00", "iopub.status.idle": "2020-09-08T10:19:22.294000+02:00"}
R"""
gg <- ggplot(p_df, aes(x=k)) + geom_line(aes(y=m)) +
    geom_ribbon(aes(ymin=ml, ymax=mu), alpha=0.2) +
    scale_x_continuous(expand=c(0, 0), limits=c(0, 1701)) +
    labs(x="Nearest neighbor id (m)", y="Correlation of\nlocal expression vectors") +
    theme(axis.text=element_text(size=16), axis.title=element_text(size=20))

ggsave("./plots/osm_fish/cor_vs_k.pdf", width=6, height=5)
gg
"""
```

```julia execution={"iopub.status.busy": "2020-09-08T10:14:19.852000+02:00", "iopub.execute_input": "2020-09-08T10:14:19.852000+02:00", "iopub.status.idle": "2020-09-08T10:14:20.102000+02:00"}
R"""
gg <- ggplot(p_df, aes(x=d)) + geom_line(aes(y=m)) +
    geom_ribbon(aes(ymin=ml, ymax=mu), alpha=0.2) +
    scale_x_continuous(expand=c(0, 0), limits=c(0, 4201)) +
    labs(x="Eucledian distance between the vectors", y="Correlation of\nlocal expression vectors") +
    theme(axis.text=element_text(size=16), axis.title=element_text(size=20))

ggsave("./plots/osm_fish/cor_vs_distance.pdf", width=6, height=5)
gg
"""
```
