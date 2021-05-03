import Clustering
import Plots

function clustermap(mtx::T where T <: AbstractMatrix{Float64}, gene_names::Vector{String}; gene_ord::Union{<:AbstractVector{<:Integer}, Bool}=true,
        cell_ord::Union{<:AbstractVector{<:Integer},Bool}=true, diag_genes::Bool=false, kwargs...)
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

    plt = Plots.heatmap(mtx[gene_ord, cell_ord]; yticks=(1:length(gene_names), gene_names[gene_ord]), kwargs...)
    return plt, cell_ord, gene_ord
end