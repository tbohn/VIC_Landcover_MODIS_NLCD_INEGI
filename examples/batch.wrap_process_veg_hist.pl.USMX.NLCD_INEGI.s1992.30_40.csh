## To run this script:
## - replace $PROJECT with the path to your clone of this GitHub repo
## - replace $LCROOT with the location of the top-level directory for land cover classifications
## - replace $AGGROOT with the location of the top-level directory for processing aggregated MODIS observations
## - $DATA_ROOT with the location of the top-level directory for storing tables, masks, domain files, and final VIC parameters
##
## Also, you must run this script from $PROJECT/tools or add $PROJECT/tools to tyour $PATH

# To run this script, replace $PROJECT with the path to your clone of this GitHub repo; $AGGROOT with the location of the top-level directory for processing aggregated MODIS observations; and $DATA_ROOT with the location of the top-level directory for storing tables, masks, domain files, and final VIC parameters
wrap_process_veg_hist.pl $AGGROOT/NLCD_INEGI $AGGROOT/NLCD_INEGI s1992 veg_hist.30 0,1,2,3,4,5,6,7,8,9 1992 s1992 3600 20 1 $DATA_ROOT/USMX domain.USMX.10x10 params.USMX.L2015.10x10 >& log.wrap_process_veg_hist.pl.USMX.NLCD_INEGI.s1992.30_40.0-9.txt

