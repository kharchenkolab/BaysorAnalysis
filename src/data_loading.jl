using DrWatson
import Images
import CSV

function load_dataset(seg_df_path::String, paper_seg_labels_path::String, watershed_path::String, dapi_path::String, prior_path::String="", pciseq_path::String="";
        name::String, min_area::Float64, min_mols_per_cell=50, paper_polygons::Bool=false, dapi::Bool=false, watershed::Bool=false, pciseq=true, prior=true,
        rescale::Bool=false, store_labels::Bool=false)
    df_spatial, gene_names = B.load_df(datadir(seg_df_path))
    res = Dict(:df => df_spatial, :gene_names => gene_names, :min_area => min_area,
        :min_mols_per_cell => min_mols_per_cell, :name => name)

    if rescale
        df_spatial.x .-= minimum(df_spatial.x);
        df_spatial.y .-= minimum(df_spatial.y);

        df_spatial.x .*= 10.0; # Specifically for Allen smFISH
        df_spatial.y .*= 10.0;
    end

    df_spatial[!, :gene_name] = gene_names[df_spatial.gene]

    paper_seg_labels = nothing
    if length(paper_seg_labels_path) > 0
        paper_seg_labels = B.load_segmentation_mask(datadir(paper_seg_labels_path));
        df_spatial[!, :cell_paper] = denserank(B.staining_value_per_transcript(df_spatial, paper_seg_labels)) .- 1;
        if store_labels
            res[:paper_labels] = paper_seg_labels
        end
    end

    # df_seg_prior = B.load_df(datadir("exp_pro/allen_smfish/baysor_prior/segmentation.csv"))[1];
    # df_spatial[!, :cell_prior] = df_seg_prior.cell;

    if watershed && length(watershed_path) > 0
        watershed_labels = B.load_segmentation_mask(datadir(watershed_path))
        df_spatial[!, :cell_watershed] = denserank(B.staining_value_per_transcript(df_spatial, watershed_labels)) .- 1
        if store_labels
            res[:watershed_labels] = watershed_labels
        end
    end

    if dapi
        res[:dapi_arr] = Float16.(Images.load(datadir(dapi_path)))
    end

    if paper_polygons && (length(paper_seg_labels_path) > 0)
        res[:paper_polys] = B.extract_polygons_from_label_grid(Matrix(paper_seg_labels[1:3:end, 1:3:end]), grid_step=3.0)
    end

    if pciseq && (length(pciseq_path) > 0)
        df_pciseq = DataFrame(CSV.File(datadir(pciseq_path)), copycols=false)
        df_spatial[!, :cell_pciseq] = df_pciseq.neighbour
    end

    if prior && (length(prior_path) > 0)
        df_prior = DataFrame(CSV.File(datadir(prior_path)), copycols=false)
        df_spatial[!, :cell_prior] = df_prior.cell
    end

    return res
end

function load_merfish(seg_df_path::String="exp_pro/merfish_moffit/baysor/segmentation.csv"; min_area::Float64=25.0, paper_polygons::Bool=false, pciseq::Bool=true, kwargs...)
    res = load_dataset(
        "exp_raw/merfish_moffit/merfish_coords_adj.csv",
        "",
        "exp_raw/merfish_moffit/dapi_merged_watershed.tif",
        "exp_raw/merfish_moffit/dapi_merged.tiff",
        "exp_pro/merfish_moffit/baysor_prior/segmentation.csv";
        paper_polygons=false, min_area=min_area, name="MERFISH", kwargs...
    )

    res[:df][!, :cell_paper] = ifelse.(ismissing.(res[:df].cell), 0, denserank(res[:df].cell));

    df_seg = B.load_df(datadir(seg_df_path))[1]
    res[:df][!, :cell] = df_seg.cell;
    res[:df][!, :confidence] = df_seg.confidence;
    res[:df][!, :cluster] = df_seg.cluster;
    res[:df][!, :ncv_color] = df_seg.ncv_color;

    if paper_polygons
        res[:paper_polys] = B.boundary_polygons(res[:df], res[:df].cell_paper, grid_step=5.0, bandwidth=5.0);
    end

    if pciseq
        df_pciseq = DataFrame!(CSV.File(datadir("exp_pro/merfish_moffit/pciseq/spots.csv")))
        res[:df] = res[:df][.!in.(res[:df].gene_name, Ref(Set(["Blank-1", "Blank-2", "Blank-3", "Blank-4", "Blank-5", "Krt90"]))),:]
        res[:df][!, :cell_pciseq] = df_pciseq.neighbour
    end

    return res
end

load_osmfish(; min_area::Float64=70.0, kwargs...) =
    load_dataset(
        "exp_pro/osmfish/baysor/segmentation.csv",
        "exp_raw/osmfish/paper_segmentation.tiff",
        "exp_raw/osmfish/nuclei_watershed.tif",
        "exp_raw/osmfish/nuclei.tif",
        "exp_pro/osmfish/baysor_prior/segmentation.csv",
        "exp_pro/osmfish/pciseq/spots.csv";
        min_area=min_area, name="osmFISH", kwargs...
    )

load_starmap1020(; min_area::Float64=25.0, kwargs...) =
    load_dataset(
        "exp_pro/starmap_vis1020/baysor/segmentation.csv",
        "exp_raw/starmap_vis1020/segmentation.tiff",
        "", "",
        "exp_pro/starmap_vis1020/baysor_prior/segmentation.csv";
        dapi=false, watershed=false,
        min_area=min_area, name="STARmap 1020", kwargs...
    )

load_allen_smfish(csv_path::String="exp_pro/allen_smfish/baysor/segmentation.csv"; min_area::Float64=25.0, kwargs...) =
    load_dataset(
        csv_path,
        "exp_raw/allen_smfish/segmentation_labels_from_json_transposed.tiff",
        "exp_raw/allen_smfish/dapi_merged_watershed.tif",
        "exp_raw/allen_smfish/dapi_merged.tiff",
        "exp_pro/allen_smfish/baysor_prior/segmentation.csv",
        "exp_pro/allen_smfish/pciseq/spots.csv";
        rescale=true, min_area=min_area, name="Allen smFISH", kwargs...
    )

function load_iss(csv_path::String="exp_pro/iss_hippo/baysor/segmentation.csv"; min_area::Float64=1.0, paper_polygons::Bool=false, kwargs...)
    res = load_dataset(
        csv_path, "",
        "exp_raw/iss_hippo/CA1DapiBoundaries_4-3_right_watershed.tif",
        "exp_raw/iss_hippo/CA1DapiBoundaries_4-3_right.tif",
        "exp_pro/iss_hippo/baysor_prior/segmentation.csv",
        "exp_pro/iss_hippo/pciseq/spots.csv";
        paper_polygons=false, min_area=min_area, min_mols_per_cell=3, name="ISS", kwargs...
    )

    res[:df][!, :cell_paper] = res[:df].parent_id;
    if paper_polygons
        res[:paper_polys] = B.boundary_polygons(res[:df], res[:df].cell_paper, grid_step=1.0, bandwidth=3.0);
    end

    return res
end