import Clustering
import PyPlot as Plt
import Seaborn as Sns

function clustermap(mtx::T where T <: AbstractMatrix{Float64}, gene_names::Vector{String}; gene_ord::Union{<:AbstractVector{<:Integer}, Bool}=true,
        cell_ord::Union{<:AbstractVector{<:Integer},Bool}=true, diag_genes::Bool=false, transpose::Bool=false, kwargs...)
    if typeof(cell_ord) === Bool
        if cell_ord
            cell_dists = 1 .- cor(mtx);
            cell_ord = Clustering.hclust(cell_dists, linkage=:ward).order;
        else
            cell_ord = 1:size(mtx, 2)
        end
    end

    if typeof(gene_ord) === Bool
        if gene_ord
            if diag_genes
                gene_ord = sortperm(vec(mapslices(x -> findmax(x)[2], mtx[:, cell_ord], dims=2)), rev=true);
            else
                gene_dists = 1 .- cor(mtx');
                gene_ord = Clustering.hclust(gene_dists, linkage=:ward).order;
            end
        else
            gene_ord = 1:size(mtx, 1)
        end
    end

    mtx = mtx[gene_ord, cell_ord]
    tick_fun = Plt.yticks
    if !transpose
        mtx = copy(mtx')
        tick_fun = Plt.xticks
    end

    plt = Sns.heatmap(mtx; kwargs...)
    tick_fun(1:length(gene_ord), gene_names[gene_ord], rotation=(transpose ? 0 : 90))
    return plt, cell_ord, gene_ord
end

function label_axis!(ax, label::String; x=-0.2, y=1.05, fontsize=24)
    ax.text(x, y, label, transform=ax.transAxes, fontsize=fontsize, va="top", ha="right", fontfamily="sans-serif", fontname="Helvetica")
end

function plot_embedding(embedding::Matrix{<:Real}, groups::Vector; mark_groups::Bool=true, ax=Plt.gca(),
        rasterized::Bool=true, fontsize::Int=10, palette::Union{Nothing, Vector{String}}=nothing, kwargs...)
    annot_labels = sort(unique(groups))
    ids_per_clust = B.split_ids(denserank(groups))
    kwargs = Dict{Symbol, Any}(kwargs...)
    for (i,ids) in enumerate(ids_per_clust)
        if palette !== nothing
            kwargs[:color] = palette[i]
        end

        ax.scatter(embedding[1,ids], embedding[2,ids]; label=annot_labels[i], rasterized=rasterized, kwargs...)
    end

    if mark_groups
        clust_centers = hcat([median(embedding[:, ids], dims=2) for ids in ids_per_clust]...);
        bbox_style = Dict(:facecolor=>"white", :alpha=>0.75, :edgecolor=>"none", :pad=>0.0)
        for (t, x, y) in zip(annot_labels, clust_centers[1,:], clust_centers[2,:])
            ax.annotate(t, (x,y), horizontalalignment="center", fontsize=fontsize, bbox=bbox_style)
        end
    end

    return ax
end
