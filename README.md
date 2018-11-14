# VIC_Landcover_MODIS_NLCD_INEGI

This project contains several versions of land cover parameters for the VIC land surface model, derived from MODIS observations and land cover classifications from MODIS, NLCD, and INEGI over the US, Mexico, and Southern Canada, as well as the scripts used to generate them.

<more info about directory structure here>

Scripts:

The processing is broken up into 10x10 degree (geographic projection) tiles.  For any given 10x10 tile, multiple MODIS sinusoidal projection tiles will overlap with the 10x10 tile (MODIS tiles are sort of rhomboidal, but they have north and south boundaries aligned with 10-degree latitude increments, so it's only the Eand W boundaries that mismatch with 10x10 tiles).

Step 1: Download and aggregate MODIS files (after which point the large MODIS files can be discarded)

 - batch.wrap_wrap_download_join_and_agg_MODIS_over_landcover.pl.CONUS_MX.30_40.csh
   - Calls wrap_wrap_download_join_and_agg_MODIS_over_landcover.pl.parallel for a specified range of 10x10 tiles

 - wrap_wrap_download_join_and_agg_MODIS_over_landcover.pl.parallel
   - Loops over the specified range of 10x10 tiles that are in the specified domain and calls a separate instance of wrap_download_join_and_agg_MODIS_over_landcover.pl in the background, i.e., in parallel, for each 10x10 tile

 - wrap_download_join_and_agg_MODIS_over_landcover.pl
   - calls wget to download the relevant MODIS tiles (with the option to ingore tiles that have already been downloaded) and calls join_and_agg_MODIS_over_landcover.py

 - join_and_agg_MODIS_over_landcover.py
   - aggregates the MODIS data over the specified land cover classification.  The two options allowed are: MODIS, which has the same gridding as the MODIS LAI and therefore has an easy 1:1 mapping with them; and NLCD, which is at 30 m resolution (reprojected to geographic and broken up into 1x1 degree tiles).  Note: this script is what checks the LAI QC codes for clouds, snow, bad retrievals, etc; it also creates urban LAI values from a prescribed NDVI-LAI relationship (since the MODIS LAI product has nulls over urban pixels!); and it computes Fcanopy from NDVI.  No land cover classes are omitted.  You have to supply a table describing, for each land cover class, how it should handle the QC flags and whether it should generate LAI from the NDVI-LAI relationship.  Note also that the MODIS land cover map is something I created - for each pixel, I took the "mode" PFT over all years that the map was produced (2001-2013).  I can give more info on that if you need.

Step 2: Gap-filling and other post-processing

 - batch.wrap_process_veg_hist.pl.PR.2001.csh
   - Calls wrap_process_veg_hist.pl for a given range of 10x10 tiles and a given domain.

 - wrap_process_veg_hist.pl
   - loops over the 10x10 tiles and calls an instance of process_veg_hist.single_file.pl for each one, in parallel. User must specify the range of 10x10 tiles, the domain, a comma-separated list of processing stages to perform, and the number of processes to run in parallel.

 - process_veg_hist.single_file.pl
   - Performs the user-specified sequence of processing steps on a given 10x10 tile.  Stages are:
     0. cleanup_albedo.py - this removes statistical outliers from the timeseries of albedo from each veg tile of each grid cell.  This because I noticed that in the time steps leading up to and immediately following the snow season, suspiciously high albedo values appeared that seemed clearly due to presence of snow that wasn't sufficient to trigger the QC flags but was clearly out of character.
     1. compute_clim_veg_hist.py - computes climatological mean and std of LAI, NDVI, Fcanopy, and albedo for each 8-day interval of the year.  Also records the count of valid observations of each for each veg tile of each grid cell.
     2. gapfill_veg_hist.py - does gap-filling in 2 stages: 1. temporally, based on a simple linear interpolation, and 2. spatially, based on a gaussian smoothing kernel.  By far, most gaps are filled by stage 1 (which is computationally less expensive).  The reason temporal takes priority is that, if you look at the data, you'll see that spatial heterogeneity is quite high compared to temporal (and compared to typical gap size in these dimensions).  So, filling spatially first results in weird smooth patches, and is more expensive (as I mentioned before).
     3. compute_anom_veg_hist.py - computes anomalies of the 8-day data relative to the climatologies computed in stages 1 and 2.  Anomalies are normal deviates (I think that's the correct term): x' = (x - x_mean)/x_std.
     4. gapfill_veg_hist.py - same script as stage 2, but now gap-fill the anomalies.  Why I've separated these two things is because climatology is much better observed/robust (with fewer gaps) than anomalies.  Interpolating the raw data would result in weird deviations from climatological behavior.  Interpolating them separately allows the climatology to be retained.
     5. recombine_clim_anom_veg_hist.py - recombines the gap-filled climatology and anomalies.
     6. interp_and_agg2monthly_veg_hist.py - interpolates the 8-day records to daily (for my use later to study veg variability); aggregates to monthly (but still a multi-year timeseries; saves these for later too); and computes climatological mean monthly cycle (saves these for using in the VIC parameter file).
     7. replace_vegparams_with_veghist_clim.py - takes an existing netcdf-format VIC parameter file over the target domain, and replaces the veg parameters with new ones.  User must supply the original VIC parameter file; a new veg library file; a file listing the root zone info for each class; and the monthly climatological LAI, Fcanopy, and albedo produced in stage 6.
     8. compute_clim_from_monthly.py - computes a climatology from the monthly files saved from stage 7, but spanning an arbitrary set of years.  I've used this for controlling more precisely which year's LAI etc are used with a given land cover map (e.g., year 2001 for the 2001 NLCD map; year 2011 for the 2011 map, etc).
     9. replace_vegparams_with_veghist_clim.py - this time, used to create veg prams using the climatology from stage 8 instead of from stage 6.

Some miscellaneous scripts that are handy utilities (for me, anyway):
 - create_domain_file_from_asc_inputs.py - I use this to create domain files
 - create_metsim_state_file_from_domain_file.py - I use this to create metsim initial state files
 - gridclip.py - clips a netcdf file to the specified lat/lon boundaries

