import PyPlot as Plt
import Seaborn as Sns
using DataFrames

function plot_precision_recall(pr_dfs::T where T <: AbstractDict{String, DataFrame}; col::Symbol, color_per_label::Dict{String, String}, xlabel::String,
        ax=Plt.gca(), legend::Bool=true, bins=0.0:0.05:1.0, lw=2, alpha=0.1, legend_title::String="Target", legend_loc::String="best")
    for (k, df) in pr_dfs
        color = color_per_label[k]
        ax.hist(df[!,col], bins=bins, alpha=alpha, color=color);
        ax.hist(df[!,col], bins=bins, label=k, histtype="step", color=color, lw=lw);
    end

    ax.grid(true, alpha=0.1)
    ax.set_ylabel("Num. of cells");
    ax.set_xlim(0, 1);
    ax.set_xlabel(xlabel)
    if legend
        ax.legend(title=legend_title, labelspacing=0.2, handlelength=1, loc=legend_loc);
    end
    return ax
end

function plot_correlation_violins(datasets::NamedTuple, part_cor_key::Symbol, labels::Tuple{String, String}; color_per_label::Dict{String, String}, ax=Plt.gca())
    p_df = vcat([vcat([DataFrame(:Correlation => d[part_cor_key][i][1], :Type => d[:name], :Source => t) for d in datasets if (part_cor_key in keys(d)) && (d[:name] != "ISS")]...)
        for (i,t) in enumerate(labels)]...);
    p_df = @orderby(p_df, lowercase.(:Type))

    Sns.violinplot(x=p_df.Type, y=p_df.Correlation, hue=p_df.Source, split=true, inner="quart", palette=color_per_label, bw=0.2, scale="count", saturation=1, ax=ax)
    ax.legend(title="Source", loc="lower right", frameon=false, labelspacing=0.25)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=11)
    ax.set_ylim(-1.0, 1);
    ax.set_yticks(-1.0:0.5:1.0)
    ax.set_ylabel("Correlation between parts");
    return ax
end

function plot_cell_type_clustermap(expr_per_clust, gene_names::Vector{String}, gene_subs::Vector{Int}; title::String="",
        figsize::Tuple{<:Real, <:Real}=(4, 4), width_frac::Float64=0.5, kwargs...)
    Plt.figure(figsize=(figsize[1] * width_frac, figsize[2]))
    plt, cell_ord, gene_ord = clustermap(expr_per_clust, gene_names; gene_ord=gene_subs, transpose=true, cmap="OrRd", kwargs...);
    ax = plt.axes
    ax.tick_params(axis="both", which="both",length=0)
    ax.set_xticks([]); ax.set_yticklabels(ax.get_yticklabels(), fontsize=13, rotation=0, va="baseline");
    ax.set_title(title)
    return plt
end