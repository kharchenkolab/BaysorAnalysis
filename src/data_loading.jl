using DrWatson
import Images

function load_dataset(seg_df_path::String, paper_seg_labels_path::String, watershed_path::String, dapi_path::String; 
        name::String, min_area::Float64, min_mols_per_cell=50, paper_polygons::Bool=false, dapi::Bool=false, watershed::Bool=false,
        rescale::Bool=false)
    df_spatial, gene_names = B.load_df(datadir(seg_df_path))

    if rescale
        df_spatial.x .-= minimum(df_spatial.x);
        df_spatial.y .-= minimum(df_spatial.y);

        df_spatial.x .*= 10.0; # Specifically for Allen smFISH
        df_spatial.y .*= 10.0;
    end

    paper_seg_labels = nothing
    if length(paper_seg_labels_path) > 0
        paper_seg_labels = B.load_segmentation_mask(datadir(paper_seg_labels_path));
        df_spatial[!, :cell_paper] = denserank(B.staining_value_per_transcript(df_spatial, paper_seg_labels)) .- 1;
    end

    # df_seg_prior = B.load_df(datadir("exp_pro/allen_smfish/baysor_prior/segmentation.csv"))[1];
    # df_spatial[!, :cell_prior] = df_seg_prior.cell;

    if watershed
        watershed_labels = B.load_segmentation_mask(datadir(watershed_path))
        df_spatial[!, :cell_watershed] = denserank(B.staining_value_per_transcript(df_spatial, watershed_labels)) .- 1
    end

    res = Dict(:df => df_spatial, :gene_names => gene_names, :min_area => min_area, :min_mols_per_cell => min_mols_per_cell, :name => name)

    if dapi
        res[:dapi_arr] = Float16.(Images.load(datadir(dapi_path)))
    end

    if paper_polygons && (length(paper_seg_labels_path) > 0)
        res[:paper_polys] = B.extract_polygons_from_label_grid(Matrix(paper_seg_labels[1:3:end, 1:3:end]), grid_step=3.0)
    end

    return res
end

function load_merfish(; min_area::Float64=25.0, paper_polygons::Bool=false, dapi::Bool=false, watershed::Bool=false)
    res = load_dataset(
        "exp_raw/merfish_moffit/merfish_coords_adj.csv",
        "",
        "exp_raw/merfish_moffit/dapi_merged_watershed.tif",
        "exp_raw/merfish_moffit/dapi_merged.tiff",
        paper_polygons=false, dapi=dapi, watershed=watershed,
        min_area=min_area, name="MERFISH"
    )

    res[:df][!, :cell_paper] = ifelse.(ismissing.(res[:df].cell), 0, denserank(res[:df].cell));

    df_seg = B.load_df(datadir("exp_pro/merfish_moffit/baysor/segmentation.csv"))[1]
    res[:df][!, :cell] = df_seg.cell;
    res[:df][!, :confidence] = df_seg.confidence;
    res[:df][!, :cluster] = df_seg.cluster;

    if paper_polygons
        res[:paper_polys] = B.boundary_polygons(res[:df], res[:df].cell_paper, grid_step=5.0, bandwidth=5.0);
    end

    return res
end

load_osmfish(; min_area::Float64=70.0, paper_polygons::Bool=false, dapi::Bool=false, watershed::Bool=false) =
    load_dataset(
        "exp_pro/osmfish/baysor/segmentation.csv",
        "exp_raw/osmfish/paper_segmentation.tiff",
        "exp_raw/osmfish/nuclei_watershed.tif",
        "exp_raw/osmfish/nuclei.tif",
        paper_polygons=paper_polygons, dapi=dapi, watershed=watershed,
        min_area=min_area, name="osmFISH"
    )

load_starmap1020(; min_area::Float64=25.0, paper_polygons::Bool=false) =
    load_dataset(
        "exp_pro/starmap_vis1020/baysor/segmentation.csv", # TODO: cl0
        "exp_raw/starmap_vis1020/segmentation.tiff",
        "", "",
        paper_polygons=paper_polygons, dapi=false, watershed=false,
        min_area=min_area, name="STARmap 1020"
    )

load_allen_smfish(; min_area::Float64=25.0, paper_polygons::Bool=false, dapi::Bool=false, watershed::Bool=false) =
    load_dataset(
        "exp_pro/allen_smfish/baysor/segmentation.csv",
        "exp_raw/allen_smfish/segmentation_labels_from_json_transposed.tiff",
        "exp_raw/allen_smfish/dapi_merged_watershed.tif",
        "exp_raw/allen_smfish/dapi_merged.tiff",
        paper_polygons=paper_polygons, dapi=dapi, watershed=watershed,
        rescale=true, min_area=min_area, name="Allen smFISH"
    )

function load_iss(; min_area::Float64=1.0, paper_polygons::Bool=false, dapi::Bool=false, watershed::Bool=false)
    res = load_dataset(
        "exp_pro/iss_hippo/baysor/segmentation.csv",
        "",
        "exp_raw/iss_hippo/CA1DapiBoundaries_4-3_right_watershed.tif",
        "exp_raw/iss_hippo/CA1DapiBoundaries_4-3_right.tif",
        paper_polygons=false, dapi=dapi, watershed=watershed,
        min_area=min_area, min_mols_per_cell=3, name="ISS"
    )

    res[:df][!, :cell_paper] = res[:df].parent_id;
    if paper_polygons
        res[:paper_polys] = B.boundary_polygons(res[:df], res[:df].cell_paper, grid_step=1.0, bandwidth=3.0);
    end

    return res
end