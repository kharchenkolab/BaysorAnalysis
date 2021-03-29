# Data

## Baysor processing

To run Baysor, the following datasets must be downloaded:

- [STARmap](https://www.starmapresources.com/data/)
  - Mouse VISp 1020 (Segmentation figure, Supp. Table 2, Bnchmarking figures). visual_1020_20180505_BY3_1kgenes
  - Mouse VISp 160 (Challenges figure). vis160_20171120_BF4_light
- osmFISH
  - [Data page](http://linnarssonlab.org/osmFISH/availability/)
  - [Dataset](http://linnarssonlab.org/osmFISH/osmFISH_SScortex_mouse_all_cells.loom)
  - [Processed data](http://pklab.med.harvard.edu/viktor/baysor/osm_fish/mRNA_coords_raw_counting.csv)
- ISS CA1
  - [Page](https://doi.org/10.6084/m9.figshare.7150760.v1)
- MERFISH
  - No data available
- seq-FISH+ NIH/3T3 cells
  - [Page](https://doi.org/10.5281/zenodo.2669683)
- [Membrane data](https://content.cruk.cam.ac.uk/jmlab/SpatialMouseAtlas2020/) from Marioni and Cai labs

## Pre-processing

### Membrane data

```r
library(dplyr)
meta <- readRDS("./exp_raw/marioni_membrane/metadata.Rds")
df <- readRDS("./exp_raw/marioni_membrane/mRNA.Rds") %>% 
  inner_join(meta[c("uniqueID", "embryo", "pos")])

dfs <- df %>% split(paste(.$embryo, .$z, sep="_"))
for (n in names(dfs)) {
  readr::write_csv(dfs[[n]], paste0("./exp_raw/marioni_membrane/molecules_", n, ".csv"))
}

dir.create("./exp_raw/marioni_membrane/tiles")
dfs <- df %>% split(paste(.$embryo, .$pos, .$z, sep="_"))
for (n in names(dfs)) {
  readr::write_csv(dfs[[n]], paste0("./exp_raw/marioni_membrane/tiles/molecules_", n, ".csv"))
}

sv <- readRDS("./exp_raw/marioni_membrane/segmentation_vertices.Rds")
readr::write_csv(sv, "./exp_raw/marioni_membrane/segmentation_vertices.csv")
```


## Baysor result analysis