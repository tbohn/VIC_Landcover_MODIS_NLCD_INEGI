# To run this script, replace $PROJECT with the path to your clone of this GitHub repo; and $AGGROOT with the location of the top-level directory for processing aggregated MODIS observations
wrap_process_veg_hist.pl $AGGROOT/NLCD_INEGI_MODIS $AGGROOT/NLCD_INEGI_MODIS 2011 veg_hist.30 0,1,2,3,4,5,6,7 null 2011 3600 20 1 >& log.wrap_process_veg_hist.pl.CONUS_MX.MOD_IGBP.mode_PFT.30_40.0-7.txt


wrap_process_veg_hist.pl $AGGROOT/NLCD_INEGI_MODIS $AGGROOT/NLCD_INEGI_MODIS 2011 veg_hist.30 8,9 null 2011 3600 10 1 >& log.wrap_process_veg_hist.pl.CONUS_MX.MOD_IGBP.mode_PFT.30_40.8-9.txt


