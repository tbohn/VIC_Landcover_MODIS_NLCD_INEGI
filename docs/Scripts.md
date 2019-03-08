# Scripts used in processing data for the VIC_Landcover_MODIS_NLCD_INEGI project

## 1. Downloading MODIS Data

## 2. Preparing the MOD12Q1.005 Land Cover Classification
 - download_MODIS.pl
 - find_mode_MODIS_PFT.py

 - batch.wrap_wrap_download_join_and_agg_MODIS_over_landcover.pl.CONUS_MX.30_40.csh
## 3. Aggregating MODIS Land Surface Characteristics over a Land Cover Map
 - join_and_agg_MODIS_over_landcover.py
 - wrap_download_join_and_agg_MODIS_over_landcover.pl
 - wrap_wrap_download_join_and_agg_MODIS_over_landcover.pl.parallel

 - batch.wrap_process_veg_hist.pl.PR.2001.csh
 - wrap_process_veg_hist.pl
 - process_veg_hist.single_file.pl
## 4. Gap-Filling (and other processing)
 - cleanup_albedo.py
 - compute_clim_veg_hist.py
 - compute_anom_veg_hist.py
 - gapfill_veg_hist.py
 - interp_and_agg2monthly_veg_hist.py
 - compute_clim_from_monthly.py
 - recombine_clim_anom_veg_hist.py
 - replace_vegparams_with_veghist_clim.py

## 5. Combining with Existing VIC Parameter File (Retaining the Existing Soil Properties and Replacing the Land Cover and Land Surface Properties)

## Utility Functions
 - breakup_veg_hist_into_annual_files.py
 - create_domain_file_from_asc_inputs.py
 - create_metsim_state_file_from_domain_file.py
 - disagg_veghist_monthly2hourly_nc.py
 - gapfill_vic_params.py
 - gridclip.py
 - grid_mosaic.py
 - grid_multi_clip.py
 - grid_utils.pl
 - make_runcell_consistent_with_forcing_mask.py
 - set_run_cell.py
 - subsample_forcings_over_domain.py
 - subsample_vic_params_over_domain.py
 - wrap_disagg_veghist_monthly2hourly_nc.pl
 - wrap_gridclip.pl
 - wrap_grid_multi_clip.pl

