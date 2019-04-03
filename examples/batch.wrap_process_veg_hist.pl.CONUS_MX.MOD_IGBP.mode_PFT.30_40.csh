## To run this script:
## - replace $PROJECT with the path to your clone of this GitHub repo
## - replace $LCROOT with the location of the top-level directory for land cover classifications
## - replace $AGGROOT with the location of the top-level directory for processing aggregated MODIS observations
## - $DATA_ROOT with the location of the top-level directory for storing tables, masks, domain files, and final VIC parameters
##
## Also, you must run this script from $PROJECT/tools or add $PROJECT/tools to tyour $PATH

wrap_process_veg_hist.pl $AGGROOT/MOD_IGBP $AGGROOT/MOD_IGBP mode_PFT veg_hist.30 0,1,2,3,4,5,6,7,8,9 null mode_PFT 3600 20 1 $DATA_ROOT/CONUS_MX domain.CONUS_MX.10x10 params.CONUS_MX.L2015.10x10 >& log.wrap_process_veg_hist.pl.CONUS_MX.MOD_IGBP.mode_PFT.30_40.0-9.txt

