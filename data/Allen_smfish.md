# Data
## Preprocessing 

Downloading the relevant files
```
mkdir -p raw_data/allen_smfish
wget https://s3.amazonaws.com/starfish.data.spacetx/spacejam2/smFISH_Allen/smFISH_MCT_CZI_Panel_0_spot_table.csv -P raw_data/allen_smfish 
```

Running Baysor (assuning Baysor is the downloaded binary)
```
mkdir -p exp_pro/allen_smfish/baysor/
Baysor run -m 70 -s 100 -i 500 -p -x x -y y --gene target -o exp_pro/allen_smfish/baysor/ smFISH_MCT_CZI_Panel_0_spot_table.csv
```

This command will create the following directory structure
```
segmentation_borders.html 
segmentation_cell_stats.csv  
segmentation_config.toml  
segmentation_counts.tsv  
segmentation.csv  
segmentation_diagnostics.html  
segmentation_log.log  
segmentation_params.dump
```

Run [notebook](https://github.com/kharchenkolab/BaysorAnalysis/blob/master/notebooks/segmentation_free/segmentation_free_allen_smfish.ipynb) that make use of `exp_pro/allen_smfish/baysor/segmentation.csv`
