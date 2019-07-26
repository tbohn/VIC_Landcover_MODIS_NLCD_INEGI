## To run this script:
## - replace $PROJECT with the path to your clone of this GitHub repo
## - replace $CONFIG with the path/filename of your config file 
## - replace $LCROOT with the location of the top-level directory for land cover classifications
## - replace $AGGROOT with the location of the top-level directory for processing aggregated MODIS observations
## - $DATA_ROOT with the location of the top-level directory for storing tables, masks, domain files, and final VIC parameters
##
## Also, you must run this script from $PROJECT/tools or add $PROJECT/tools to tyour $PATH

wrap_wrap_download_join_and_agg_MODIS_over_landcover.pl.parallel $CONFIG $LCROOT MODIS mode_PFT MCD12Q1 $PROJECT/data/CONUS_MX/lc_table.CONUS_MX.MOD_IGBP.csv CONUS_MX $PROJECT/data/CONUS_MX/mask.land.CONUS_MX.0.0625_deg.asc 2000 2016 $PROJECT/data/CONUS_MX/mask.land.CONUS_MX.10_deg.asc 30 40 -180 0 240 0.0625 $AGGROOT/MOD_IGBP/mode_PFT/aggregated veg_hist 0 &> log.wrap_wrap_download_join_and_agg_MODIS_over_landcover.pl.CONUS_MX.MODIS.mode_PFT.30_40.txt 
