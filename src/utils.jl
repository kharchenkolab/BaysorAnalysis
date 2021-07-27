import PyPlot as Plt
using SparseArrays

function sparse_counts(vec1::Vector{<:Integer}, vec2::Vector{<:Integer})
    counts = Dict{Tuple{Int, Int}, Int}()
    for tp in zip(vec1, vec2)
        counts[tp] = get(counts, tp, 0) + 1
    end
    i1 = getindex.(keys(counts), 1)
    i2 = getindex.(keys(counts), 2)
    v = collect(values(counts))
    return sparse(i1, i2, v)
end


function vec_triu(M::AbstractMatrix{T}) where T
    m, n = size(M)
    m == n || throw(error("not square"))
    l = n*(n-1) ÷ 2
    v = zeros(T, l)
    k = 0
    @inbounds for i in 1:n
        for j in 1:(i-1)
            v[k + j] = M[j, i]
        end
        k += i-1
    end
    v
end

function method_palette(subs::Union{Vector{String}, Nothing}=nothing)
    d = Dict(
        "Baysor" => "#E69F00",
        "Watershed" => "#3093ca",
        "pciSeq" => "#1e9e00",
        "poly-A" => "#F0E442",
        "IF" => "#CC79A7",
        "Baysor, DAPI prior" => "#9566FF",
        "Baysor, IF prior" => "#6F4D0B"
    )

    (subs !== nothing) || return d
    return Dict(k => d[k] for k in subs)
end

function set_pyplot_defaults!()
    Plt.rc("axes", axisbelow=true, labelsize=16, grid=true);
    Plt.rc("axes.spines", right=false, top=false, left=true, bottom=true);
    Plt.rc("font", size=12)
    Plt.rc("grid", alpha=0.25)
    Plt.rc("legend", frameon=false, borderpad=0)
end

function upscale(image::T where T <: AbstractMatrix{TR}, ratio::Int64) where TR <: Real
    res_img = zeros(TR, size(image) .* ratio)
    ct = 1
    for c in 1:size(image, 2)
        for m1 in 1:ratio
            rt = 1
            for r in 1:size(image, 1)
                for m2 in 1:ratio
                    res_img[rt, ct] = image[r, c]
                    rt += 1
                end
            end
            ct += 1
        end
    end

    return res_img
end
