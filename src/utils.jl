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