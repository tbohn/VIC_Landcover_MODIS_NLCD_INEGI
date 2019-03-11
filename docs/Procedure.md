# Processing Steps for the VIC_Landcover_MODIS_NLCD_INEGI Project

To use these scripts, make sure that the path to the "tools" directory is in your `$PATH` environment variable.

Quick note on MODIS files: the MODIS observations are on a sinusoidal grid, broken up into "tiles" that are roughly equal-area (see map at https://lpdaac.usgs.gov/dataset_discovery/modis). The tiles have north and south boundaries that correspond to 10-degree latitude intervals, but the east and west boundaries are slanted at various angles (relative to a geographic projection) depending on how far from the Greenwich meridian they are. So pixels in these files occur in rows that occur at regular latitude intervals, but within the rows, the pixels are not regularly spaced with respect to longitude (farther apart further from the equator).

## Stage 1: Download and aggregate MODIS files

1. Create an account with NASA Earthdata. To obtain a NASA Earthdata Login account, visit https://urs.earthdata.nasa.gov/users/new/.

2. Create a top-level directory on your machine, on a drive with sufficient space to hold the data that you will need (each MODIS sinusoidal tile for a single 8-day interval, at 500 m resolution, takes up about 8 MB for MCD15A2H.006 to 80 MB for MCD43A3.006; and there are 46 8-day intervals per year). Something like "MODIS".

  We will call the path above "MODIS" `$MODISROOT`.  So, now you have a directory called `$MODISROOT/MODIS`.

  Under `$MODISROOT/MODIS`, create the following subdirectories: "LAI", "NDVI", and "albedo".  If you also anticipate using the MCD12Q1.051 land cover classification, create a subdirectory called "PFT" (for plant functional type) to store the files for it.

3. Unfortunately, the MODIS sinusoidal tile filenames include the date/time on which they were published, which we can't just predict by any formula. I am not aware of a way to do wget with a wildcard like "*". So to download a subset of the files (say, all files for tile h08v05), we need to first download the "index.html" files that contain a list of all files in a given directory, and then parse those files with a perl script (see below) to get the specific filenames that we want to download:

   Then, run "wget" to download the index.html files.

   `wget -r -P $MODISROOT/MODIS/$SUBDIR --http-user=$USERNAME --http-password=$PASSWORD https://e4ftl01.cr.usgs.gov/$TERRA_AQUA/$PRODUCT`

   where

   `$MODISROOT` = the top-level directory for MODIS data  
   `$SUBDIR` is one of "LAI", "NDVI", "albedo", or "PFT"  
   `$TERRA_AQUA` = "MOLA" for "MOD" products, and "MOTA" for "MCD" products  
   `$PRODUCT` = product code, e.g. MOD15A2H.006  
   `$USERNAME` = your NASA Earthdata username  
   `$PASSWORD` = your NASA Earthdata password  

   This "wget" command will copy the directory structure of the data pool onto your machine, creating a subdirectory of `e4ftl01.cr.usgs.gov/$TERRA_AQUA/$PRODUCT` under your current directory, containing further subdirectories corresponding to all available 8-day intervals with a format of `$YYYY.$MM.$DD`, where `$YYYY` = 4-digit year, `$MM` = 2-digit month, and `$DD` = 2-digit day. Each of these contains and index.html file listing all available .hdf files that can be downloaded. It will unfortunately also create directories for other MODIS products; when it starts downloading index.html files from the other products, you should kill the wget command via ctrl-C (if running in the foreground) or by doing `kill -9 $PID` where `$PID` = numeric process id associated with the desired wget instance.

   The MODIS hdf filenames follow the convention

   `$PROD.A$YYYY$MM$DD.h$COLv$ROW.$COLLECTION.$FILEDATE.hdf`

   where

   `$PROD` = the first part of `$PRODUCT`, e.g., MOD15A2H  
   `$YYYY` = 4-digit year of acquisition  
   `$MM` = 2-digit month of acquisition  
   `$DD` = 2-digit day of acquisition  
   `$COL` = 2-digit column of the tile (see map at https://lpdaac.usgs.gov/dataset_discovery/modis)  
   `$ROW` = 2-digit row of the tile (see map at https://lpdaac.usgs.gov/dataset_discovery/modis)  
   `$COLLECTION` = the 2nd part of `$PRODUCT`, e.g., 006  
   `$FILEDATE` = 13-digit date/time code for the file (corresponding to when the file was generated, NOT acquisition date

   Once the downloads are finished, go to the tools directory.

4. (Optional) Downloading and preparing MCD12Q1.051 land cover classification

   If you wish to use the MCD12Q1.051 product as the underlying land cover classification, you will need to download this separately from the other MODIS products.

   Once you have created the "$MODISROOT/MODIS/PFT" directory and done an initial wget to obtain the index.html files, run the following script to download the MODIS data:

   `download_MODIS.pl $MODISROOT/MODIS/PFT/e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.051 MCD12Q1 $STARTDATE $ENDDATE $HMIN $HMAX $VMIN $VMAX $USERNAME $PASSWORD`

   where

   `$MODISROOT` = the top-level directory for MODIS data
   `$STARTDATE` = earliest desired acquisition date; for MCD12Q1, I used 2001.01.01  
   `$ENDDATE` = latest desired acquisition date; for MCD12Q1, I used 2013.01.01  
   `$HMIN` = minimum desired column of the MODIS sinusoidal tile grid  
   `$HMAX` = maximum desired column of the MODIS sinusoidal tile grid  
   `$VMIN` = minimum desired row of the MODIS sinusoidal tile grid  
   `$VMAX` = maximum desired row of the MODIS sinusoidal tile grid  
   `$USERNAME` = your NASA Earthdata username  
   `$PASSWORD` = your NASA Earthdata password  

   Once you have downloaded the MCD12Q1.051 hdf files, make sure you have a directory for storing land cover classifications. We will call it `$LCROOT`. Under `$LCROOT` there can be subdirectories for "MODIS", "NLCD_INEGI", etc.

   Run the following script to find the most frequent (mode) class for each pixel across all acquisition dates:

   `wrap_find_mode_MODIS_PFT.pl $MODISROOT/MODIS/PFT/e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.051 MCD12Q1 $STARTDATE $ENDDATE $HMIN $HMAX $VMIN $VMAX $LCROOT/MODIS/mode_PFT MCD12Q1.mode`

   where

   `$MODISROOT` = the top-level directory for MODIS data
   `$STARTDATE` = earliest desired acquisition date; for MCD12Q1, I used 2001.01.01  
   `$ENDDATE` = latest desired acquisition date; for MCD12Q1, I used 2013.01.01  
   `$HMIN` = minimum desired column of the MODIS sinusoidal tile grid  
   `$HMAX` = maximum desired column of the MODIS sinusoidal tile grid  
   `$VMIN` = minimum desired row of the MODIS sinusoidal tile grid  
   `$VMAX` = maximum desired row of the MODIS sinusoidal tile grid  
   `$LCROOT` = the top-level directory for land cover classifications

5. (Optional) Obtain the NLCD_INEGI land cover classification(s)

   You can download the NLCD_INEGI land cover classifications from [Zenodo](https://zenodo.org/record/2580428).

   Make sure you have a top-level directory for land cover classifications, which we will call `$LCROOT`. Make a subdirectory under this called "NLCD_INEGI". Store the NLCD_INEGI under there.

6. Downloading of MODIS land surface properties and aggregation over the land cover classification

   Make sure you have a top-level directory for storing and processing the aggregated MODIS observations, which we will call `$AGGROOT`.

   After you have downloaded the index.html files, you can then run the following script to download the hdf files and aggregate the data over the land cover classification:

   `wrap_wrap_download_join_and_agg_MODIS_over_landcover.pl.parallel $LCROOT $LCTYPE $LCID $LCPFX $LC_TABLE $DOMAIN $LANDMASK $STARTYEAR $ENDYEAR $COARSE_MASK $LATMIN $LATMAX $LONMIN $LONMAX $PIX_PER_DEG $OUT_RES $AGGROOT/$LCTYPE/$LCID/aggregated $OUTPFX $FORCE`

   where

   `$LCROOT` = path to top-level directory below which the land cover classifications are stored  
   `$LCTYPE` = either "MODIS" (for MCD12Q1) or "NLCD_INEGI" (for NLCD_INEGI)  
   `$LCID` = either "mode_PFT" (for MCD12Q1) or the year (2001, 2011, s1992, s2001, s2011) (for NLCD_INEGI)  
   `$LCPFX` = prefix of land cover classification files, either "MCD12Q1" or "nlcd_inegi"  
   `$LC_TABLE` = table listing land cover classes and their processing options; these can be found in the `data/$DOMAIN` directories  
   `$DOMAIN` = domain name, either "CONUS_MX" or "USMX" (NLCD_INEGI only covers the USMX domain)  
   `$LANDMASK` = domain mask at the output resolution; land mask files can be found in the `data/$DOMAIN` directories  
   `$STARTYEAR` = first year of MODIS land surface observations to process (I used 2000)  
   `$ENDYEAR` = last year of MODIS land surface observations to process (I used 2016)  
   `$COARSE_MASK` = either the filename of a mask of 10x10 degree tiles (geographic, not sinusoidal) covering the domain, or "null" to just process the whole domain as a single file; coarse masks can be found in the `data/$DOMAIN` directories  
   `$LATMIN` = southern boundary of region to process (if supplying a 10x10 tile mask, this should correspond to the southern edge of the southernmost tile that you want to process)  
   `$LATMAX` = northern boundary of region to process (if supplying a 10x10 tile mask, this should correspond to the northern edge of the northernmost tile that you want to process)  
   `$LONMIN` = western boundary of region to process (if supplying a 10x10 tile mask, this should correspond to the western edge of the westernmost tile that you want to process)  
   `$LONMAX` = eastern boundary of region to process (if supplying a 10x10 tile mask, this should correspond to the eastern edge of the easternmost tile that you want to process)  
   `$PIX_PER_DEG` = number of MODIS pixels per degree (for all products discussed here, it is 240)  
   `$OUT_RES` = output grid resolution in degrees (I used 0.0625)  
   `$AGGROOT` = path to top-level directory of the tree of output directories  
   `$LCTYPE` = either "MODIS" (for MCD12Q1) or "NLCD_INEGI" (for NLCD_INEGI)  
   `$OUTPFX` = prefix for output NetCDF files; I used "veg_hist"  
   `$FORCE` = either 0 (don't overwrite existing output files) or 1 (overwrite)  

   This script has 2 stages; 1. download the MODIS data (and figure out which MODIS tiles correspond to the region of interest); 2. aggregate over the land cover classificaton. Once stage 1 has completed successfully, running this script again will not re-run stage 1 unless `$FORCE` is set to 1.

   Variables that start with "LC" refer to the land cover classification. `$STARTYEAR` etc through `$PIX_PER_DEG` refer to the MODIS land surface observation time series. `$COARSE_MASK` is the filename of a mask over the domain at 10-degree resolution, for the purpose of breaking up large domains into smaller pieces and processing those pieces in parallel. if `$COARSE_MASK` is "null", the domain will not be divided into 10x10 degree tiles, but instead will be processed as a single region in one processing stream. `$LATMIN` through `$LONMAX` describe the geographic bounds of the region to be processed (if supplying non-null `$COARSE_MASK`, then these bounds must coincide with boundaries of 10x10 degree tiles). Output files (1 file per 10x10 degree tile if `$COARSE_MASK` is non-null) will be written to `$AGGROOT/$LCTYPE/$LCID/aggregated/`.

   If `$COARSE_MASK` is not "null", each 10x10 tile will be processed as a separate background process, so the lat/lon bounds should be chosen carefully so as not to create too many simultaneous processes. Also, each process uses a large amount of memory, which is another motivation for checking the memory usage on a single tile before submitting more than one at a time.

   This script makes at least one call to `wrap_download_join_and_agg_MODIS_over_landcover.pl`. If `$COARSE_MASK` is not "null", then for each 10x10 tile, a separate instance of `wrap_download_join_and_agg_MODIS_over_landcover.pl` will be called and run in the background, i.e., in parallel with any other tiles. If `$COARSE_MASK` is "null", then a single instance of `wrap_download_join_and_agg_MODIS_over_landcover.pl` will be called to process the entire region defined by the lat/lon bounds.

   - `wrap_download_join_and_agg_MODIS_over_landcover.pl` calls wget to download the relevant MODIS tiles (with the option to ingore tiles that have already been downloaded) and calls `join_and_agg_MODIS_over_landcover.py`  
   - `join_and_agg_MODIS_over_landcover.py` aggregates the MODIS data over the specified land cover classification.  The two options allowed are: MODIS, which has the same gridding as the MODIS LAI and therefore has an easy 1:1 mapping with them; and NLCD, which is at 30 m resolution (reprojected to geographic and broken up into 1x1 degree tiles).  This script is what checks the LAI QC codes for clouds, snow, bad retrievals, etc; it also creates urban LAI values from a prescribed NDVI-LAI relationship (since the MODIS LAI product has nulls over urban pixels); and it computes Fcanopy from NDVI.  

   Examples of running `wrap_wrap_download_join_and_agg_MODIS_over_landcover.pl.parallel` for one row of 10x10 tiles (the row with latitudes spanning 30-40 deg N) for the CONUS_MX and USMX domains, respectively, can be found in the batch files under the "examples" directory:

   - batch.wrap_wrap_download_join_and_agg_MODIS_over_landcover.pl.CONUS_MX.MODIS.mode_PFT.30_40.csh  
   - batch.wrap_wrap_download_join_and_agg_MODIS_over_landcover.pl.USMX.NLCD_INEGI.2011.30_40.csh  

   To run these batch files, replace all instances of `$PROJECT` in the files with the path to your copy of this GitHub project. Also make sure they are executable by running `chmod +x $FILENAME`, where `$FILENAME` is the path/name of the batch file that you are making executable.

   After the aggregation is complete, the large MODIS input files can be discarded.

## Stage 2: Gap-filling and other post-processing

   Before this stage can begin, we need a NetCDF VIC-5 compliant domain file. For the CONUS_MX and USMX domains, domain files are available for download on [Zenodo](https://zenodo.org/record/2564019). If you have a different domain, you can either clip out your domain from within these domain files (if your domain is inside those domains) or create a domain file from elemental inputs (mask, and DEM) with a utility script (see the "Utility Scripts" section below).

   If you wish to divide the domain into 10x10 degree boxes (for parallel processing), you will also need to divide the domain file into similar boxes.  See the "Utility Scripts" section below). If you do this, you need to name the directory containing the 10x10 tiles of the domain file as xxxx and the individual 10x10 files should have a prefix of xxx (the lat/lon range covered by the tile will be appended to the prefix).

   For domains that have been divided into 10x10 tiles, the script that manages all processing is `wrap_process_veg_hist.pl`. This script loops over the 10x10 tiles and calls an instance of `process_veg_hist.single_file.pl` for each one, in parallel. User must specify the range of 10x10 tiles, the domain, a comma-separated list of processing stages to perform, and the number of processes to run in parallel. Usage:

   `wrap_process_veg_hist.pl $AGGROOT $AGGROOT2 $LCID $PREFIX $STAGELIST $LCID_CV $LCID_OUT $TIME_OFFSET $NPARALLEL $CLEAN $DATA_ROOT $DOMAIN_PFX $PARAM_PFX $LC_SCHEME`

   where

   `$AGGROOT` = path to top-level directory of the tree of output directories  
   `$AGGROOT2` = can be equal to `$AGGROOT`; if you are storing all subsequent processing outputs in a different location, this gives you the option to do it  
   `$LCID` = either "mode_PFT" (for MCD12Q1) or the year (2001, 2011, s1992, s2001, s2011) (for NLCD_INEGI)  
   `$PREFIX` = should be same `$OUTPFX` from Stage 1, but if you are dividing into 10x10 tiles, and you want to specify a subset of files to process, you should add any information that would distinguish these files (e.g., for 30-40 latitude, specify veg_hist.30)  
   `$STAGELIST` = comma-separated list of processing stages to run, e.g., "1,2,3,4,5,6,7,8,9". Stages:  
   - 0:  
   - 1:  
   - 2:  

   `$LCID_CV` = xxx  
   `$LCID_OUT` = xxx  
   `$TIME_OFFSET` = xxx  
   `$NPARALLEL` = xxx  
   `$CLEAN` = xxx  
   `$DATA_ROOT` = xxx  
   `$DOMAIN_PFX` = xxx  
   `$PARAM_PFX` = xxx  
   `$LC_SCHEME` = xxx  
$rootdir = shift;
$archivedir = shift;
$lcid = shift;
$prefix = shift;
$stagelist = shift;
$lcid_cv = shift;
$lcid_out = shift;
$time_offset = shift; # seconds
$nParallel = shift; # number of tiles to process in parallel
$clean = shift; # 1 = delete files from previous steps as we go
$data_root = shift;
$domain_pfx = shift;
$param_pfx = shift;
$lcscheme = shift;

   If you have a small domain and want to process it as a single stream (each stage of processing deals with a single file containing the entire domain), run:

   `process_veg_hist.single_file.pl xxx`

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

   Examples of running `xxx` for one row of 10x10 tiles (the row with latitudes spanning 30-40 deg N) for the CONUS_MX and USMX domains, respectively, can be found in the batch files under the "examples" directory:

   - batch.xxx.CONUS_MX.MODIS.mode_PFT.30_40.csh  
   - batch.xxx.USMX.NLCD_INEGI.2011.30_40.csh  

   To run these batch files, replace all instances of `$PROJECT` in the files with the path to your copy of this GitHub project. Also make sure they are executable by running `chmod +x $FILENAME`, where `$FILENAME` is the path/name of the batch file that you are making executable.

## Utility Scripts

!!! Fill this out more
Some useful utility scripts:
 - create_domain_file_from_asc_inputs.py - I use this to create domain files
 - create_metsim_state_file_from_domain_file.py - I use this to create metsim initial state files
 - gridclip.py - clips a netcdf file to the specified lat/lon boundaries
 - the subsample scripts - both for subsampling and for clipping to a domain file (with or without subsampling)
