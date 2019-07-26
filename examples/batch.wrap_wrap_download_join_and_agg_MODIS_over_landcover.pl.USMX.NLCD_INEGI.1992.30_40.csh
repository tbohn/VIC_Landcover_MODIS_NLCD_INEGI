### NOTE: This script is designed to generate output files for the 1992
### land cover map just covering year 2000 of MODIS observations, solely
### for the purpose generating a "Cv" variable (land cover area fractions)
### that can be used with the s1992 parameter set.
## To run this script:
## - replace $PROJECT with the path to your clone of this GitHub repo
## - replace $CONFIG with the path/filename of your config file 
## - replace $LCROOT with the location of the top-level directory for land cover classifications
## - replace $AGGROOT with the location of the top-level directory for processing aggregated MODIS observations
## - $DATA_ROOT with the location of the top-level directory for storing tables, masks, domain files, and final VIC parameters
##
## Also, you must run this script from $PROJECT/tools or add $PROJECT/tools to tyour $PATH

wrap_wrap_download_join_and_agg_MODIS_over_landcover.pl.parallel $CONFIG $LCROOT NLCD_INEGI 1992 $PROJECT/data/USMX/lc_table.USMX.NLCD_INEGI.csv USMX $PROJECT/data/USMX/mask.land.USMX.0.0625_deg.asc 2000 2000 $PROJECT/data/USMX/mask.land.USMX.10_deg.asc 30 40 -180 0 240 0.0625 $AGGROOT/NLCD_INEGI/1992/aggregated veg_hist 0 &> log.wrap_wrap_download_join_and_agg_MODIS_over_landcover.pl.USMX.NLCD_INEGI.1992.30_40.txt 
