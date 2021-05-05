import CairoMakie as MK
import Clustering
import Colors
import Random
import Plots
import SparseArrays

using LinearAlgebra
using ProgressMeter
using Statistics
using StatsPlots

function assignment_summary_df(assignments::Pair{Symbol, Vector{Int}}...; min_molecules_per_cell::Int)
    assignments = Dict(assignments)
    return DataFrame(Dict(
        "name" => collect(keys(assignments)),
        "noise_frac" => [round(mean(x .== 0), sigdigits=3) for x in values(assignments)],
        "n_cells" => [sum(B.count_array(x, drop_zero=true) .>= min_molecules_per_cell) for x in values(assignments)]
    ))[:, [:name, :n_cells, :noise_frac]]
end

function plot_subset(df_spatial::DataFrame, dapi_arr::Union{Matrix, Nothing}, (xs, xe), (ys, ye); polygons::Union{Bool, Vector{Matrix{Float64}}}=true, markersize=2.0, alpha=0.2, min_molecules_per_cell::Int=1,
        grid_step::Float64=5.0, bandwidth::Float64=grid_step, min_border_length=3, cell_col::Symbol=:cell, dapi_alpha=0.9, polygon_line_width::T1 where T1 <: Real=2, dens_threshold::Float64=1e-5,
        noise::Bool=true, size_mult=1/3, plot_raw_dapi::Bool=true, color_col::Symbol=:ncv_color, annotation_col::Union{Symbol, Nothing}=nothing, build_panel::Bool=true, grid_alpha::Float64=0.5, ticks=false,
        swap_plots::Bool=false, dapi_color::Symbol=:twilight, clims=nothing, polygon_line_color="black", plot_bg_dapi::Bool=true, kwargs...)
    df_subs = @where(df_spatial, :x .>= xs, :x .<= xe, :y .>= ys, :y .<= ye);

    if (typeof(polygons) == Bool)
        if polygons
            polygons = B.boundary_polygons(df_subs, df_subs[!, cell_col], grid_step=grid_step, min_molecules_per_cell=min_molecules_per_cell,
                bandwidth=bandwidth, min_border_length=min_border_length, dens_threshold=dens_threshold)
        else
            polygons = Matrix{Float64}[]
        end
    end

    # xs, xe, ys, ye = round.(Int, [minimum(df_subs.x), maximum(df_subs.x), minimum(df_subs.y), maximum(df_subs.y)]);
    xs, xe, ys, ye = round.(Int, [xs, xe, ys, ye]);

    xticks_vals = range(0, xe-xs, length=5)
    yticks_vals = range(0, ye-ys, length=5)

    yticks = xticks = Float64[]
    if ticks
        xticks = (xticks_vals, ["$f" for f in round.(Int, range(xs, xe, length=5))])
        yticks = (yticks_vals, ["$f" for f in round.(Int, range(ys, ye, length=5))])
    end

    plot_size = ((xe-xs), ye-ys) .* size_mult
    if plot_raw_dapi
        plot_size = (2 * plot_size[1], plot_size[2])
    end

    is_noise = noise ? (df_subs[!, cell_col] .== 0) : nothing
    annotation = (annotation_col === nothing) ? nothing : df_subs[!, annotation_col]

    polygon_kwargs = Dict(:strokewidth => polygon_line_width, :strokecolor => polygon_line_color)
    if dapi_arr === nothing
        return B.plot_molecules(df_subs, polygons; color=color_col, markersize=markersize, alpha=alpha, offset=(-xs, -ys), size=plot_size,
            is_noise=is_noise, polygon_kwargs=polygon_kwargs, noise_kwargs=Dict(:markersize => markersize), annotation=annotation, kwargs...)
    end

    dapi_subs = copy(dapi_arr[ys:ye, xs:xe]')
    plot_args = Dict(:size => plot_size, :xticks => xticks, :yticks => yticks)

    fig = MK.Figure(resolution=plot_size)
    get_axis() = ticks ?
        MK.Axis(fig, xticks=xticks, yticks=yticks) :
        MK.Axis(fig, xticklabelsvisible=false, xticksvisible=false, yticklabelsvisible=false, yticksvisible=false)

    fig[1, 1] = get_axis()

    if plot_bg_dapi
        MK.image!(Colors.GrayA.(1 .- dapi_subs ./ maximum(dapi_subs), dapi_alpha), interpolate=false)
    end

    B.plot_molecules!(df_subs, polygons; color=color_col, markersize=markersize, alpha=alpha, offset=(-xs, -ys),
        polygon_kwargs=polygon_kwargs, is_noise=is_noise,
        noise_kwargs=Dict(:markersize => markersize), annotation=annotation, kwargs...)

    if !plot_raw_dapi
        return fig
    end

    if clims === nothing
        clims = B.val_range(dapi_arr)
    else
        error("clims are not currently supported")
    end

    # Possible colorschemes: :delta, :twilight, :PuBuGn_9, :PuBu_9, :dense, :tofino, :berlin, :delta
    fig[1, 2] = get_axis()
    MK.heatmap!(dapi_subs, colormap=dapi_color)
     # TODO: p .- [xs ys] after Baysor update
    MK.poly!([MK.Point2.(eachrow(p .- [xs ys])) for p in polygons]; color="transparent", polygon_kwargs...)

    if swap_plots
        fig[1, 1:2] .= fig[1, [2, 1]]
    end

    if !build_panel
        error("build_panel is not supported")
        # return plt1, plt2
    end

    return fig
end

# function cell_coord_frame(df_spatial::DataFrame, cell_id::Int; cell_col::Symbol=:cell, offset::Float64=0.1)
#     sample_df = df_spatial[df_spatial[!, cell_col] .== cell_id,:];
#     return [round.(Int, ((1. + offset) * s - offset * e, (1. + offset) * e - offset * s)) for (s,e) in [B.val_range(sample_df.x), B.val_range(sample_df.y)]];
# end

function plot_comparison_for_cell(df_spatial::DataFrame, cell_id::Int, args...; cell1_col::Symbol=:cell, cell2_col::Symbol=:cell_paper, offset::Float64=0.1, kwargs...)
    sample_df = df_spatial[df_spatial[!, cell1_col] .== cell_id,:];
    paper_ids = setdiff(unique(sample_df[!, cell2_col]), [0]);

    if !isempty(paper_ids)
        sample_df = df_spatial[(df_spatial[!, cell1_col] .== cell_id) .| in.(df_spatial[!, cell2_col], Ref(paper_ids)),:];
    end;

    xc, yc = median(sample_df.x[sample_df[!, cell1_col] .== cell_id]), median(sample_df.y[sample_df[!, cell1_col] .== cell_id])
    xls, yls = [round.(Int, ((1. + offset) * s - offset * e, (1. + offset) * e - offset * s)) for (s,e) in [B.val_range(sample_df.x), B.val_range(sample_df.y)]];
    return plot_comparison_for_cell(df_spatial, xls, yls, args...; xc=xc, yc=yc, cell_col=cell1_col, kwargs...)
end

function plot_comparison_for_cell(df_spatial::DataFrame, xls::Tuple{T, T}, yls::Tuple{T, T}, seg_arr::Union{Matrix{<:Integer}, Nothing},
        dapi_arr::Union{Matrix{<:Real}, Nothing}; paper_polys::Array{Matrix{Float64}, 1}=Matrix{Float64}[], polygon_line_width::Float64=2.0, polygon_alpha::Float64=0.5,
        size_mult::Float64=1.0, grid_alpha::Float64=0.0, markersize::Float64=2.0, title="", center_mult::Float64=3.0, noise::Bool=false, paper_poly_color="darkred", paper_line_mult::Float64=1.0,
        plot_raw_dapi::Bool=true, kwargs...) where T <: Real

    xls, yls = max.(xls, 1), max.(yls, 1)

    if dapi_arr !== nothing
        xls = min.(xls, size(dapi_arr, 2))
        yls = min.(yls, size(dapi_arr, 1))
    end

    if seg_arr !== nothing
        paper_polys = extract_polygons_from_label_grid(copy(seg_arr[yls[1]:yls[2], xls[1]:xls[2]]'))
        # paper_polys = [Plots.Shape(pg[:,1], pg[:,2]) for pg in paper_polys]
    else
        paper_polys = [pg .- [xls[1] yls[1]] for pg in paper_polys]
    end

    if dapi_arr === nothing
        plot_raw_dapi = false
    end

    fig = plot_subset(df_spatial, dapi_arr, xls, yls; size_mult=size_mult, grid_alpha=grid_alpha, markersize=markersize, noise=noise,
        polygon_line_width=polygon_line_width, polygon_alpha=polygon_alpha, plot_raw_dapi=plot_raw_dapi, kwargs...);

    if length(paper_polys) > 0
        for i in 1:length(fig.scene.children)
            MK.poly!(fig[1, i], [MK.Point2.(eachrow(p)) for p in paper_polys]; color="transparent", strokecolor=paper_poly_color, strokewidth=polygon_line_width*paper_line_mult)
        end;
    end

    return fig
end

function joint_ordering(matrices::Matrix{T}...) where T <: Real
    joint_mtx = hcat(matrices...)
    cell_dists = Symmetric(1 .- cor(joint_mtx));
    cell_ord = Clustering.hclust(cell_dists, linkage=:ward).order;

    cell_ords = [Int[] for m in matrices]
    cum_sum = 0
    for i in 1:length(matrices)
        cell_ords[i] = cell_ord[(cell_ord .<= (cum_sum + size(matrices[i], 2))) .& (cell_ord .> cum_sum)] .- cum_sum
        cum_sum += size(matrices[i], 2)
    end

    gene_dists = Symmetric(1 .- cor(joint_mtx'));
    gene_ord = Clustering.hclust(gene_dists, linkage=:ward).order;

    return cell_ords, gene_ord
end

## Comparison of segmentations

prepare_qc_df(df_spatial::DataFrame, cell_col::Symbol=:cell; kwargs...) =
    prepare_qc_df(df_spatial, Int.(df_spatial[!, cell_col]); kwargs...)


function prepare_qc_df(df_spatial::DataFrame, assignment::Vector{<:Int}; min_area::T where T<:Real, min_molecules_per_cell::T2 where T2 <: Real,
        max_elongation::T3 where T3 <: Real = 30.0, dapi_arr::Union{Matrix{<:Real}, Nothing}=nothing)
    qc_per_cell = B.get_cell_qc_df(df_spatial, assignment, dapi_arr=dapi_arr);
    qc_per_cell[!, :cell_id] = 1:size(qc_per_cell, 1)
    qc_per_cell[!, :sqr_area] = sqrt.(qc_per_cell.area)
    return @where(qc_per_cell, :n_transcripts .>= min_molecules_per_cell, :sqr_area .>= min_area, :elongation .<= max_elongation)
end

hist_bins(vals::Vector{<:Real}...; n_bins::Int=100, min_val::T where T<: Real=0.0, m_quantile::T2 where T2<:Real=0.99, max_val_mult::Float64=1.0) =
    range(min_val, maximum([quantile(v, m_quantile) / m_quantile for v in vals if length(v) > 0])*max_val_mult, length=n_bins)

function plot_qc_comparison(qc_per_cell_dfs::Union{Array{DataFrame, 1}, Tuple{DataFrame, DataFrame}}; labels::Vector{String}=["Baysor", "DAPI"],
        max_quants::Vector{<:Real}=[0.9999, 0.99, 0.99, 0.999, 0.999], n_bins::Vector{<:Real}=[75, 100, 100, 100, 100], kwargs...)
    m_names = [:n_transcripts, :density, :elongation, :sqr_area, :mean_dapi_brightness]
    plot_titles = ["Num. of transcripts", "Density", "Elongation", "sqrt(Area)", "Mean DAPI brightness"]
    plots = Plots.Plot[]
    for (cs, xlab, mq, nb) in zip(m_names, plot_titles, max_quants, n_bins)
        if any(!(cs in propertynames(qdf)) for qdf in qc_per_cell_dfs)
            continue
        end
        t_bins = hist_bins([qdf[!, cs] for qdf in qc_per_cell_dfs]..., m_quantile=mq, n_bins=nb)
        plt = Plots.plot(widen=false, xlabel=xlab, ylabel="Num. of cells", xlims=B.val_range(t_bins))
        for (qdf, lab) in zip(qc_per_cell_dfs, labels)
            if length(unique(qdf[!, cs])) < 2
                continue
            end
            Plots.histogram!(qdf[!, cs], label=lab, bins=t_bins, alpha=0.6, kwargs...)
        end
        push!(plots, plt)
    end

    return Plots.plot(plots..., size=(400 * ceil(Int, length(plots) / 2), 600))
end

function match_assignments(qc_per_cell_dfs::Union{Array{DataFrame, 1}, Tuple{DataFrame, DataFrame}}, df_spatial::DataFrame, cell_cols::Vector{Symbol}; kwargs...)
    assignments = [denserank(ifelse.(in.(df_spatial[!,s], Ref(Set(qdf.cell_id))), df_spatial[!,s], 0)) .- 1 for (s, qdf) in zip(cell_cols, qc_per_cell_dfs)];
    return match_assignments(assignments...; kwargs...)
end

function match_assignments(assignment1::Vector{<:Integer}, assignment2::Vector{<:Integer}; bin_match_frac::Float64=0.05)
    # Basic statistics
    contingency_table = sparse_counts(assignment1 .+ 1, assignment2 .+ 1);
    ci,ri,v = SparseArrays.findnz(contingency_table)
    m_norm1 = sparse(ci, ri, v ./ sum(contingency_table, dims=1)[ri])
    m_norm2 = sparse(ci, ri, v ./ sum(contingency_table, dims=2)[ci])

    noise_fracs = Vector.((m_norm2[2:end,1], m_norm1[1,2:end]));
    max_overlaps = Vector.((maximum(m_norm2[2:end,:], dims=2)[:], maximum(m_norm1[:,2:end], dims=1)[:]));

    bin_match = dropzeros!(sparse(ci, ri, m_norm1.nzval .>= bin_match_frac))[2:end, 2:end];
    n_overlaps = (sum(bin_match, dims=2)[:], sum(bin_match, dims=1)[:])

    ctn_min = sparse(ci, ri, min.(m_norm1.nzval, m_norm2.nzval))
    match_fracs = range(0.0, 1.0, length=30)
    match_nums = [sum(any(ctn_min .>= mt, dims=1)) for mt in match_fracs];

    # Derivative statistics
    match_noise = [nf .> 0.5 for nf in noise_fracs]
    multiple_overlap = [((mo .< 0.7) .& (nf .< 0.25)) for (mo, nf) in zip(max_overlaps, noise_fracs)];
    @assert all([!any(mn .& mo) for (mn, mo) in zip(match_noise, multiple_overlap)])

    return (contingency=contingency_table, noise_fracs=noise_fracs, max_overlaps=max_overlaps, n_overlaps=n_overlaps, match_fracs=match_fracs,
        match_nums=match_nums, match_noise=match_noise, multiple_overlap=multiple_overlap)
end

function plot_matching_comparison(match_results::NamedTuple; labels::Vector{String}=["Baysor", "DAPI"], n_bins::Vector{<:Real}=[50, 30, 30], kwargs...)
    m_names = [:noise_fracs, :max_overlaps]
    plot_titles = ["Fraction of molecules, not assigned in other segmentation", "Fraction of molecules, matching to a single cell"]
    plots = Plots.Plot[]
    for (n, xlab, nb) in zip(m_names, plot_titles, n_bins)
        vals = match_results[n]
        t_bins = hist_bins(vals..., m_quantile=1.0, n_bins=nb, max_val_mult=1.02)
        plt = Plots.histogram(vals[1], label=labels[1], bins=t_bins, xlims=B.val_range(t_bins), widen=false,
            xlabel=xlab, ylabel="Num. of cells", kwargs...)
        Plots.histogram!(vals[2], label=labels[2], bins=t_bins, alpha=0.6)
        push!(plots, plt)
    end

    max_bin = maximum(maximum.(match_results.n_overlaps))
    plt = Plots.histogram(match_results.n_overlaps[1], label=labels[1], bins=0:0.5:max_bin, xticks=0:max_bin, widen=false,
        xlabel="Number of overlapping cells", ylabel="Num. of cells", kwargs...)
    Plots.histogram!(match_results.n_overlaps[2], label=labels[2], bins=0:0.5:max_bin, alpha=0.6)
    push!(plots, plt)

    plt = Plots.plot(match_results.match_fracs, match_results.match_nums, xlabel="Minimal fraction of matching molecules", ylabel="Num. of matching cells",
        xlims=(0.0, 1.0), ylims=(0, maximum(match_results.match_nums) * 1.05), widen=false, legend=:none, lw=2.0)
    push!(plots, plt)

    return Plots.plot(plots..., layout=(2, 2), size=(900, 600))
end

function build_statistics_df(qc_per_cell_dfs::Union{Array{DataFrame, 1}, Tuple{DataFrame, DataFrame}}, match_res::NamedTuple, df_spatial::DataFrame; labels::Vector{Symbol}=[:Baysor, :DAPI], sigdigits::Int=3)
    stat_df = DataFrame(Dict(:Metric => ["Num. cells", "Total num. molecules", "Fraction of noise", "Fraction of matching cells", "Fraction of matching to noise", "Fraction of multiple overlap"]))
    for i in 1:length(qc_per_cell_dfs)
        stat_df[!, labels[i]] = Any[
            size(qc_per_cell_dfs[i], 1), sum(qc_per_cell_dfs[i].n_transcripts), round(1 .- sum(qc_per_cell_dfs[i].n_transcripts) / size(df_spatial, 1), sigdigits=sigdigits),
            round(mean(.!match_res.match_noise[i] .& (.!match_res.multiple_overlap[i])), sigdigits=sigdigits), round(mean(match_res.match_noise[i]), sigdigits=sigdigits),
            round(mean(match_res.multiple_overlap[i]), sigdigits=sigdigits)
        ]
    end

    return stat_df
end

estimate_embedding(df_spatial::DataFrame, qc_per_cell_dfs::Union{Array{DataFrame, 1}, Tuple{DataFrame, DataFrame}}, cell_cols::Vector{Symbol}; kwargs...) =
    estimate_embedding([convert_segmentation_to_counts(df_spatial.gene, df_spatial[!, cq])[:, qdf.cell_id] for (cq, qdf) in zip(cell_cols, qc_per_cell_dfs)]...; kwargs...)

function estimate_embedding(count_matrices::Matrix{<:Real}...; joint::Bool=true, kwargs...)
    cm_merged = hcat(count_matrices...);
    cm_merged = cm_merged ./ sum(cm_merged, dims=1)
    ids_per_mat = B.split_ids(vcat([repeat([i], inner=size(cm, 2)) for (i,cm) in enumerate(count_matrices)]...))
    umap_merged = fit(B.UmapFit, cm_merged; kwargs...);
    return [umap_merged.embedding[:, ids] for ids in ids_per_mat], umap_merged, cm_merged
end

function plot_qc_embeddings(qc_per_cell_dfs::Union{Array{DataFrame, 1}, Tuple{DataFrame, DataFrame}}, match_res::NamedTuple, embeddings::Array{Matrix{Float64}, 1}; log_colors::Bool=true, size=(1000, 1000),
        labels::Vector{String}=["Baysor", "DAPI"], legend::Symbol=:best, ms::Float64=3.0, alpha1::Float64=0.25, plot_num_transcripts::Bool=true, layout=nothing,
        plot_partial::Bool=false, palette=Plots.distinguishable_colors(4, cchoices=[80, 100]), build_figure::Bool=true, kwargs...)
    plts = Plots.Plot[]
    c_vals = [qdf.n_transcripts for qdf in qc_per_cell_dfs]
    if log_colors
        c_vals = [log.(cv) for cv in c_vals]
    end

    factor_per_cell = String[]

    if plot_partial
        factor_per_cell = [ifelse.(match_res.multiple_overlap[i], "Multiple", ifelse.(match_res.match_noise[i], "No",
            ifelse.(match_res.max_overlaps[i] .< 0.7, "Partial", "Full"))) for i in 1:length(qc_per_cell_dfs)];
    else
        factor_per_cell = [ifelse.(match_res.multiple_overlap[i], "Multiple", ifelse.(match_res.match_noise[i], "No", "Yes")) for i in 1:length(qc_per_cell_dfs)];
    end

    for i in 1:length(qc_per_cell_dfs)
        plt = Plots.plot(format=:png, size=(600, 600), legend=legend, bg_legend=Colors.RGBA(1.0, 1.0, 1.0, 0.5), title="$(labels[i]), matching", legend_title="Match")
        for (ci,ids) in enumerate(B.split_ids(denserank(factor_per_cell[i])))
            Plots.scatter!(embeddings[i][1,ids], embeddings[i][2,ids], label=factor_per_cell[i][ids[1]], ms=ms, alpha=alpha1, markerstrokewidth=0; color=palette[ci], kwargs...)
        end
        push!(plts, plt)

        plot_num_transcripts || continue

        colors = map_to_colors(c_vals[i], lims=B.val_range(vcat(c_vals...)))[:colors]
        plt = Plots.scatter(embeddings[i][1,:], embeddings[i][2,:], color=colors, ms=ms, markerstrokewidth=0, format=:png, size=(600, 600), legend=false,
            title="$(labels[i]), num. transcripts")
        push!(plts, plt)
    end;

    if !build_figure
        return plts
    end

    return Plots.plot(plts..., size=size, format=:png, layout=something(layout, length(plts)))
end

plot_expression_vec_comparison(df_spatial::DataFrame, qc_per_cell_dfs::Union{Array{DataFrame, 1}, Tuple{DataFrame, DataFrame}}, cell_cols::Vector{Symbol}, gene_names; labels=["Baysor", "DAPI"], kwargs...) =
    plot_expression_vectors([B.count_array(vcat(split(df_spatial.gene, df_spatial[!, cs] .+ 1)[2:end][qdf.cell_id]...)) for (cs, qdf) in zip(cell_cols, qc_per_cell_dfs)]...; gene_names=gene_names, labels=labels, kwargs...)

function estimate_non_matching_part_correlation(df_spatial::DataFrame, qc_per_cell_dfs::Union{Array{DataFrame, 1}, Tuple{DataFrame, DataFrame}},
        match_res::NamedTuple; cell_cols::Vector{Symbol}=[:cell, :cell_paper], rev::Bool=false, max_overlap_borders::Tuple{Float64, Float64}=(0.25, 0.75))
    data_id = rev ? 2 : 1
    if rev
        cell_cols = reverse(cell_cols)
    end
    part_overlap_mask = (match_res.max_overlaps[data_id] .> max_overlap_borders[1]) .& (match_res.max_overlaps[data_id] .< max_overlap_borders[2]);
    cur_overlaps = match_res.max_overlaps[data_id][part_overlap_mask]
    partially_matching_ids = qc_per_cell_dfs[data_id].cell_id[part_overlap_mask];
    match_cors = Float64[]

    n_genes = maximum(df_spatial.gene)
    for t_cell_id in partially_matching_ids
        t_df = df_spatial[df_spatial[!, cell_cols[1]] .== t_cell_id,:];
        t_mcp = mode(t_df[!, cell_cols[2]][t_df[!, cell_cols[2]] .!= 0]);

        match_expr = B.prob_array(t_df.gene[t_df[!, cell_cols[2]] .== t_mcp], max_value=n_genes);
        non_match_expr = B.prob_array(t_df.gene[t_df[!, cell_cols[2]] .!= t_mcp], max_value=n_genes);

        push!(match_cors, cor(non_match_expr, match_expr))
    end

    return match_cors, partially_matching_ids, cur_overlaps, part_overlap_mask
end

function correlation_effect_size_df(datasets::NamedTuple, part_cor_key::Symbol, cell_col_labels::Vector{String}, cell_col_names::Vector{Symbol})
    dfs = DataFrame[]
    for i in 1:2
        for d in datasets
            (part_cor_key in keys(d)) || continue

            corr_res = d[part_cor_key][i]
            qc_df = d[:qc_per_cell_dfs][cell_col_names[i]]
            corr = corr_res[1];
            ord = sortperm(corr);
            corr = corr[ord]

            mismatch_frac = 1 .- corr_res[3][ord];
            n_trans_per_cell = qc_df.n_transcripts[corr_res[4]][ord];

            push!(dfs, DataFrame(
                :Correlation => corr,
                :MolFrac => cumsum(mismatch_frac .* n_trans_per_cell) ./ sum(qc_df.n_transcripts),
                :CellFrac => (1:length(corr)) ./ size(qc_df, 1),
                :Dataset => d[:name],
                :Segmentation => cell_col_labels[i]
            ))
        end
    end
    return vcat(dfs...)
end

function filter_assignment(assignment::Vector{Int}, min_mols_per_cell::Int)
    passed_ids = Set(findall(B.count_array(assignment, drop_zero=true) .>= min_mols_per_cell))
    return denserank([((a in passed_ids) ? a : 0) for a in assignment]) .- 1
end

function append_matching_statistics!(d::Dict{Symbol, Any}, cell_col_names::Vector{Symbol};
        target_cell_col_names::Union{Vector{Symbol}, Nothing}=nothing)
    cur_names = sort(intersect(cell_col_names, propertynames(d[:df])))
    if length(setdiff(cell_col_names, cur_names)) > 0
        @warn "The following col names weren't found in the data: $(setdiff(cell_col_names, cur_names))"
    end

    if target_cell_col_names === nothing
        target_cell_col_names = cur_names
    end

    d[:qc_per_cell_dfs] = Dict(k => prepare_qc_df(d[:df], k; min_area=d[:min_area],
        min_molecules_per_cell=d[:min_mols_per_cell], max_elongation=15) for k in cur_names);

    for cs in cur_names # filter cells
        d[:df][!, cs][.!in.(d[:df][!, cs], Ref(Set(d[:qc_per_cell_dfs][cs].cell_id)))] .= 0
        any(d[:df][!, cs] .!= 0) || error("All cells were filtered out for '$cs'")
    end

    d[:assignment_filt] = Dict(cs => denserank(d[:df][!, cs]) .- 1 for cs in cur_names)

    # for (i1,c1) in enumerate(cur_names)
        # for c2 in cur_names[(i1 + 1):length(cur_names)]
    @showprogress for c1 in cur_names
        for c2 in target_cell_col_names
            (c1 != c2) || continue
            parts = [get(split(String(c1), '_'), 2, "cell"), get(split(String(c2), '_'), 2, "cell")]
            parts = [p for p in parts if length(p) > 0]
            sm = Symbol(join(vcat(["match_res"], parts), "_"))
            d[sm] = match_assignments(d[:assignment_filt][c1], d[:assignment_filt][c2]);

            sel_cells = [c1, c2]
            spc = Symbol(join(vcat(["part_cors"], parts), "_"))
            qc_dfs = [d[:qc_per_cell_dfs][k] for k in sel_cells]
            d[spc] = [estimate_non_matching_part_correlation(d[:df], qc_dfs, d[sm], rev=r, cell_cols=sel_cells) for r in [false, true]]

            ssd = Symbol(join(vcat(["stat_df"], parts), "_"))
            d[ssd] = build_statistics_df([d[:qc_per_cell_dfs][k] for k in sel_cells], d[sm], d[:df]; labels=sel_cells)

            sfv = Symbol(join(vcat(["frac_vals"], parts), "_"))
            overlap_mask = d[sm].n_overlaps[1] .== 1
            overlap_nums, uniq_overlap_ids = vec.(collect.(findmax(d[sm].contingency[2:end,:][overlap_mask, 2:end], dims=2)));
            d[sfv] = (
                overlap = overlap_nums ./ d[:qc_per_cell_dfs][c2].n_transcripts[getindex.(uniq_overlap_ids, 2)],
                assigned = d[:qc_per_cell_dfs][c1].n_transcripts[overlap_mask] ./ d[:qc_per_cell_dfs][c2].n_transcripts[getindex.(uniq_overlap_ids, 2)],
                overlap_ids = findall(overlap_mask)
            )
        end
    end
end

## Quick visualizations


function plot_correlation(d::Dict{Symbol, Any}, part_cor_key::Symbol)
    pdfs1 = DataFrame(:Correlation => d[part_cor_key][1][1], :Type => d[:name]);
    pdfs2 = DataFrame(:Correlation => d[part_cor_key][2][1], :Type => d[:name]);
    t_heights = (0.9 .* (size(pdfs1, 1), size(pdfs2, 1)) ./ max(size(pdfs1, 1), size(pdfs2, 1)))

    display(("Baysor", sum(pdfs1.Correlation .> 0.5), sum(pdfs1.Correlation .< 0.5)))
    display(("Paper", sum(pdfs2.Correlation .> 0.5), sum(pdfs2.Correlation .< 0.5)))

    plt = @df pdfs1 violin(:Type, :Correlation, side=:right, label="Baysor", ylabel="Correlation between parts", legendtitle="Source",
        bg_legend=Colors.RGBA(1.0, 1.0, 1.0, 0.9), bar_width=t_heights[1], size=(300, 240), legend=:bottomright, xgrid=false) # , color=color_per_col_label["Baysor"]
    @df pdfs2 violin!(:Type, :Correlation, side=:left, label="Paper", bar_width=t_heights[2]) # , color=color_per_col_label["Paper"]

    return plt
end

function plot_correlation_lines(d::Dict{Symbol, Any}, part_cor_key::Symbol; lw=3, max_cor=0.5)
    vals = sort(d[part_cor_key][1][1])
    vals = vals[vals .<= max_cor]
    plt = Plots.plot(vals, 1:length(vals), label="Baysor", size=(300, 200), # , color=color_per_col_label["Baysor"]
        xlabel="Correlation", ylabel="Num. cells", legend=:topleft, lw=lw)
    vals = sort(d[part_cor_key][2][1])
    vals = vals[vals .<= max_cor]
    Plots.plot!(vals, 1:length(vals), label="Paper", lw=lw) # , color=color_per_col_label["Paper"]
    return plt
end

function plot_matching_hists(frac_vals::NamedTuple, match_res::NamedTuple)
    bins = 0.0:0.02:1.0
    plt1 = Plots.histogram(frac_vals.overlap, bins=bins, label=false,
        xlabel="Overlap size / poly-a cell size", ylabel="Cell density", normalize=true, lw=0, xlim=(0, 1))
    Plots.stephist!(frac_vals.overlap, bins=bins, label="", normalize=true, color="black")
    Plots.vline!([0.9975], label="", color="black", lw=2.5);

    bins = 0.0:0.04:2.0
    plt2 = Plots.histogram(min.(frac_vals.assigned, 1.99), bins=bins, normalize=true, lw=0,
        legend=false, xlabel="Target cell size / poly-A cell size", ylabel="Cell density")
    Plots.stephist!(min.(frac_vals.assigned, 1.99), bins=bins, label="", normalize=true, color="black")
    Plots.vline!([1], label="", color="black", lw=2.5);

    vals = B.prob_array(match_res.n_overlaps[1] .+ 1);
    xvals = 0:(length(vals) - 1)
    plt3 = Plots.bar(xvals, vals, xlabel="Num. overlaps", ylabel="Fraction", label=false, xticks=xvals, xgrid=false)

    plt = Plots.plot(plt1, plt2, plt3, size=(800, 200), layout=(1, 3))
    return plt
end

function estimate_jaccard_per_molecule(df::DataFrame, cell_col1::Symbol, cell_col2::Symbol; min_mols_per_cell::Int)
    mol_ids_per_cell = [Set.(B.split_ids(df[!, k] .+ 1)[2:end]) for k in [cell_col1, cell_col2]]
    cms = [B.convert_segmentation_to_counts(df.gene, df[!, k]) for k in [cell_col1, cell_col2]]

    jac_dists = Dict{Tuple{Int, Int}, Float64}();
    exp_cors = Dict{Tuple{Int, Int}, Float64}();
    dists = Float64[]
    for i in 1:size(df, 1)
        c1, c2 = df[i, cell_col1], df[i, cell_col2]
        if (c1 == 0) || (c2 == 0)
            push!(dists, 0.0)
            continue
        end

        mpc1, mpc2 = mol_ids_per_cell[1][c1], mol_ids_per_cell[2][c2]
        if (length(mpc1) < min_mols_per_cell) || (length(mpc2) < min_mols_per_cell)
            push!(dists, 0.0)
            continue
        end

        if !((c1, c2) in keys(jac_dists))
            jac_dists[(c1, c2)] = length(intersect(mpc1, mpc2)) / length(union(mpc1, mpc2))

            exp_cors[(c1, c2)] = cor(cms[1][:,c1], cms[2][:,c2])
        end
        push!(dists, jac_dists[(c1, c2)])
    end

    return (dists=dists, jac_dists=jac_dists, exp_cors=exp_cors, cms=cms)
end

function precision_recall_df(contingency::SparseArrays.SparseMatrixCSC{Int64, Int64})
    overlap_sizes, overlap_ids = findmax(contingency[2:end, 2:end], dims=1);
    overlap_sizes = collect(overlap_sizes)[:];
    overlap_ids = getindex.(overlap_ids, 1)[:];

    n_mols_per_cell_src = sum(contingency[2:end,:], dims=2)[:];
    n_mols_per_cell_targ = sum(contingency[:, 2:end], dims=1)[:];

    n_overlaps = sum(contingency[2:end, 2:end] ./ n_mols_per_cell_targ' .> 0.25, dims=1)[:]
    recall = overlap_sizes ./ n_mols_per_cell_targ;
    precision = overlap_sizes ./ n_mols_per_cell_src[overlap_ids];
    return DataFrame(Dict(:precision => precision, :recall => recall, :n_overlaps => n_overlaps, :overlap_sizes => overlap_sizes))
end