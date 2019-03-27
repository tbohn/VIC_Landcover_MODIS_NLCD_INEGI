# Processing Steps for the VIC_Landcover_MODIS_NLCD_INEGI Project

To use these scripts, make sure that the path to the "tools" directory is in your `$PATH` environment variable. These scripts require Python 3.6.4 and Perl 5 to be installed on your system.

Quick note on MODIS files: the MODIS observations are on a sinusoidal grid, broken up into "tiles" that are roughly equal-area (see map at https://lpdaac.usgs.gov/dataset_discovery/modis). The tiles have north and south boundaries that correspond to 10-degree latitude intervals, but the east and west boundaries are slanted at various angles (relative to a geographic projection) depending on how far from the Greenwich meridian they are. So pixels in these files occur in rows that occur at regular latitude intervals, but within the rows, the pixels are not regularly spaced with respect to longitude (farther apart further from the equator).

## Stages 1-2: Download and aggregate MODIS files

There are some initial steps to perform for stage 1. Then, the remainder of stage 1 and all of stage 2 are handled by a single script.

Steps:

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
   `$LCTYPE` = either "MOD_IGBP" (for MCD12Q1) or "NLCD_INEGI" (for NLCD_INEGI)  
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
   `$LCTYPE` = either "MOD_IGBP" (for MCD12Q1) or "NLCD_INEGI" (for NLCD_INEGI)  
   `$OUTPFX` = prefix for output NetCDF files; I used "veg_hist"  
   `$FORCE` = either 0 (don't overwrite existing output files) or 1 (overwrite)  

   This script has 2 stages; 1. download the MODIS data (and figure out which MODIS tiles correspond to the region of interest); 2. aggregate over the land cover classificaton. Once stage 1 has completed successfully, running this script again will not re-run stage 1 unless `$FORCE` is set to 1.

   Variables that start with "LC" refer to the land cover classification. `$STARTYEAR` etc through `$PIX_PER_DEG` refer to the MODIS land surface observation time series. `$COARSE_MASK` is the filename of a mask over the domain at 10-degree resolution, for the purpose of breaking up large domains into smaller pieces and processing those pieces in parallel. if `$COARSE_MASK` is "null", the domain will not be divided into 10x10 degree tiles, but instead will be processed as a single region in one processing stream. `$LATMIN` through `$LONMAX` describe the geographic bounds of the region to be processed (if supplying non-null `$COARSE_MASK`, then these bounds must coincide with boundaries of 10x10 degree tiles). Output files (1 file per 10x10 degree tile if `$COARSE_MASK` is non-null) will be written to `$AGGROOT/$LCTYPE/$LCID/aggregated/`.

   If `$COARSE_MASK` is not "null", each 10x10 tile will be processed as a separate background process, so the lat/lon bounds should be chosen carefully so as not to create too many simultaneous processes. Also, each process uses a large amount of memory, which is another motivation for checking the memory usage on a single tile before submitting more than one at a time.

   This script makes at least one call to `wrap_download_join_and_agg_MODIS_over_landcover.pl`. If `$COARSE_MASK` is not "null", then for each 10x10 tile, a separate instance of `wrap_download_join_and_agg_MODIS_over_landcover.pl` will be called and run in the background, i.e., in parallel with any other tiles. If `$COARSE_MASK` is "null", then a single instance of `wrap_download_join_and_agg_MODIS_over_landcover.pl` will be called to process the entire region defined by the lat/lon bounds.

   - `wrap_download_join_and_agg_MODIS_over_landcover.pl` calls wget to download the relevant MODIS tiles (with the option to ingore tiles that have already been downloaded) and calls `join_and_agg_MODIS_over_landcover.py`  
   - `join_and_agg_MODIS_over_landcover.py` aggregates the MODIS data over the specified land cover classification.  The two options allowed are: MODIS, which has the same gridding as the MODIS LAI and therefore has an easy 1:1 mapping with them; and NLCD, which is at 30 m resolution (reprojected to geographic and broken up into 1x1 degree tiles).  This script is what checks the LAI QC codes for clouds, snow, bad retrievals, etc; it also creates urban LAI values from a prescribed NDVI-LAI relationship (since the MODIS LAI product has nulls over urban pixels); and it computes Fcanopy from NDVI.  

   Examples of running `wrap_wrap_download_join_and_agg_MODIS_over_landcover.pl.parallel` for one row of 10x10 tiles (the row with latitudes spanning 30-40 deg N) for the CONUS_MX and USMX domains, respectively, can be found in the batch files under the "examples" directory:

   - `batch.wrap_wrap_download_join_and_agg_MODIS_over_landcover.pl.CONUS_MX.MODIS.mode_PFT.30_40.csh`  
   - `batch.wrap_wrap_download_join_and_agg_MODIS_over_landcover.pl.USMX.NLCD_INEGI.2011.30_40.csh`  
   - `batch.wrap_wrap_download_join_and_agg_MODIS_over_landcover.pl.USMX.NLCD_INEGI.s1992.30_40.csh`  
   - `batch.wrap_wrap_download_join_and_agg_MODIS_over_landcover.pl.USMX.NLCD_INEGI.1992.30_40.csh` - this is only for generating land cover fractions for 1992 to be used in stages 3-4 for the s1992 parameter set; MODIS observations from 2000-2016 aggregated over the 1992 map don't make much sense.  

   To run these batch files, replace all instances of `$PROJECT` in the files with the path to your copy of this GitHub project. Also make sure they are executable by running `chmod +x $FILENAME`, where `$FILENAME` is the path/name of the batch file that you are making executable.

   After the aggregation is complete, the large MODIS input files can be discarded.

## Stages 3-4: Gap-filling and other post-processing

A single script handles stages 3 and 4. Stage 3 is gap-filling. Stage 4 is the generation of VIC parameter files by combining the new land cover and time-varying land surface characteristics with the other parameters from the "source" parameter files.

Before stage 4 can begin, we need 2 other files:
 - NetCDF VIC-5 compliant domain file for the current domain
 - NetCDF VIC-5 compliant "source" parameter file (which will supply all parameters except land cover fractions, LAI, Fcanopy, and albedo) for the current domain

For the CONUS_MX and USMX domains, domain files are available for download on Zenodo [(PITRI Precipitation Disaggregation Parameters)](https://zenodo.org/record/2564019). For these same domains, "source" parameter files based on Livneh et al (2015) ("L2015" henceforth) are availabled for download on Zenodo [(MOD-LSP VIC Parameters)](https://zenodo.org/record/2612560).

If you have a different domain, you can either clip out your domain from within these domain and parameter files (if your domain is inside those domains) or create a domain file from elemental inputs (mask, and DEM) with a utility script (see the "Utility Scripts" section below).  To create a NetCDF VIC-5 compliant "source" parameter file for a domain outside L2015, you will need ascii-format VIC soil, snow, and vegetation parameters (and vegetation library file) for the domain, and you can use the tool [tonic](https://github.com/UW-Hydro/tonic).

### Option 1: if you have divided your domain into 10x10 degree boxes

#### 1.a. Clipping

If you wish to divide the domain into 10x10 degree boxes (for parallel processing), you will also need to divide (a) the domain file and (b) the "source" VIC parameter file into similar boxes. For this, use the following script:

   `grid_multi_clip.py -i $INFILE -m $COARSE_MASK -p $PREFIX -o $OUTDIR`

where

   `$INFILE` = NetCDF VIC-5 compliant input file to be clipped into boxes (either a domain file or a VIC parameter file), e.g., `$DATA_ROOT/$DOMAIN/domain.$DOMAIN.nc`, where `$DATA_ROOT` = path to top-level directory for storing tables, masks, domain files, and vic parameters; and `$DOMAIN` is either "CONUS_MX" or "USMX".  
   `$COARSE_MASK` = same 10x10 degree mask of domain as used in Stage 1.  
   `$PREFIX` = output file prefix, e.g., `domain.$DOMAIN.10x10`.  
   `$OUTDIR` = output directory, e.g., `$DATA_ROOT/domain.$DOMAIN.10x10`.  

For each output file, the lat/lon range covered by the tile will be appended to the prefix, followed by `.nc`.

#### 1.b. Processing

For domains that have been divided into 10x10 tiles, the script that manages all processing is `wrap_process_veg_hist.pl`. This script loops over the 10x10 tiles and calls an instance of `process_veg_hist.single_file.pl` for each one, in parallel. User must specify the range of 10x10 tiles, the domain, a comma-separated list of processing stages to perform, and the number of processes to run in parallel. Usage:

   `wrap_process_veg_hist.pl $AGGROOT/$LCTYPE $AGGROOT2/$LCTYPE $LCID $PREFIX $STEPLIST $LCID_CV $LCID_OUT $TIME_OFFSET $NPARALLEL $CLEAN $DATA_ROOT/$DOMAIN $DOMAINFILE_PFX $PARAM_PFX $LCTYPE`

where

   `$AGGROOT` = path to top-level directory of the tree of output directories  
   `$LCTYPE` = same `$LCTYPE` used in Stage 1 (downloading and aggregating).  
   `$AGGROOT2` = can be equal to `$AGGROOT`; if you are storing all subsequent processing outputs in a different location, this gives you the option to do it  
   `$LCID` = either "mode_PFT" (for MCD12Q1) or the year (2001, 2011, s1992, s2001, s2011) (for NLCD_INEGI)  
   `$PREFIX` = should be same `$OUTPFX` from Stage 1, but if you are dividing into 10x10 tiles, and you want to specify a subset of files to process, you should add any information that would distinguish these files (e.g., for 30-40 latitude, specify veg_hist.30)  
   `$STEPLIST` = comma-separated list of processing steps to run, e.g., "0,1,2,3,4,5,6,7,8,9". Processing steps:  
   - 0: Remove remaining albedo observations that appear impacted by snow (calls `cleanup_albedo.py`).  
   - 101: (optional) Replace land cover area fractions with those of another land cover classification. Requires specifying `$LCID_CV` which should be the `$LCID` of some other parameter set that has completed step 0 (calls `replace_cv.py`). This is primarily used for creating parameters for the NLCD_INEGI "s" (stable land cover) parameter sets. See "Examples" section below.  
   - 1: Compute the climatological (temporal) mean and standard deviation of LAI, Fcanopy, and albedo of each of the 46 8-day observations in the annual cycle, for each land cover class in each grid cell (calls `compute_clim_veg_hist.py`).  
   - 2: Fill gaps in the climatological mean and standard deviation (calls `gapfill_veg_hist.py`).  
   - 3: Compute standardized 8-day anomalies with respect to the climatological (temporal) mean and standard deviation of LAI, Fcanopy, and albedo, for each land cover class in each grid cell (calls `compute_anom_veg_hist.py`).  
   - 4: Fill gaps in the anomalies (calls `gapfill_veg_hist.py`).  
   - 5: Recombine the gap-filled climatological mean and standard deviation with the gap-filled anomalies (calls `recombine_clim_anom_veg_hist.py`).  
   - 6: Interpolate the recombined gap-filled 8-day timeseries to daily; then aggregates to monthly and computes the climatological mean annual cycle of monthly values (calls `interp_and_agg2monthly_veg_hist.py`).  
   - 7: Takes an existing VIC parameter set and replaces the land cover fractions, LAI, Fcanopy, and albedo with the land cover fractions and climatological mean monthly values of `$LCID` (calls `replace_vegparams_with_veghist_clim.py`).  
   - 8: Selects a specific year of the monthly timeseries and uses that as the monthly annual cycle (calls `compute_clim_from_monthly.py`).  
   - 9: Takes an existing VIC parameter set and replaces the land cover fractions, LAI, Fcanopy, and albedo with the land cover fractions and selected monthly values of step 8 (calls `replace_vegparams_with_veghist_clim.py`).  

   `$LCID_CV` = if step 101 selected, should be `$LCID` of the parameter set whose area fractions you want to use; otherwise should = "null".  
   `$LCID_OUT` = if step 101 selected, should be the `$LCID` that you want to call the output; otherwise should = `$LCID`.  
   `$TIME_OFFSET` = number of seconds to wait between submitting processing jobs in parallel.  
   `$NPARALLEL` = maximum number of processing jobs to run at the same time.  
   `$CLEAN` = either 1 (= delete files from previous step as we start new step) or 0 (= don't delete any files).  
   `$DATA_ROOT` = top-level directory under which the domain and parameter files for use with VIC are stored (expects `$DATA_ROOT/$DOMAIN/$DOMAINFILE_PFX/$DOMAINFILE_PFX.$LOCSTR.nc`, where `$LOCSTR` = `$LAT0_$LAT1n.$LON0_$LON1e`, where `$LAT0`...`$LON1` are the south, north, west, and east boundaries of the 10x10 tile).  
   `$DOMAINFILE_PFX` = both the name of the directory containing 10x10 domain files, and also the file prefix of those 10x10 files.  
   `$PARAM_PFX` = this refers to the "source" vic parameter dataset of which the land cover fractions and LAI, Fcanopy, and albedo annual cycles will be replaced with the ones being processed; `$PARAM_PFX` is both the name of the directory containing 10x10 files, and also the file prefix of those 10x10 files.  

This script calls `process_veg_hist.single_file.pl` for each 10x10 tile.

#### 1.c. Mosaicking

After running `wrap_process_veg_hist.pl` on all 10x10 tiles, the end result can be mosaicked back together into a single file covering the entire domain via:

   `grid_mosaic.py -d $DOMAIN_FILE -i $AGGROOT2/$LCTYPE/$LCID/$SUBDIR -p $PREFIX -o $OUTFILE`

where

   `$SUBDIR` = subdirectory containing the files you want to mosaic together. For the output of stage 6 (17-year monthly timeseries of LAI etc.), `$SUBDIR` = "monthly". For the output of stage 7 (vic parameters with monthly annual cycle derived from climatological mean between start and end years), `$SUBDIR` = "vic_params.allyears". For the output of stage 9 (vic parameters with monthly annual cycle derived from a single year `$YEAR`), `$SUBDIR` = `vic_params.$YEAR_$YEAR`.  
   `$STARTYEAR` and `$ENDYEAR` = first and last years of period used for computing monthly annual cycle; for parameter sets using climatological mean over the period 2000-2016, these are 2000 and 2016; for parameter sets using the monthly values from a single year, `$STARTYEAR` and `$ENDYEAR` are both set to that year.  
   `$OUTFILE` = output path/filename; filename should be basically `$PREFIX.nc`.  

#### Examples

Examples of running `wrap_process_veg_hist.pl` for one row of 10x10 tiles (the row with latitudes spanning 30-40 deg N) for the CONUS_MX and USMX domains, respectively, can be found in the batch files under the "examples" directory:

 - `batch.wrap_process_veg_hist.pl.CONUS_MX.MODIS.mode_PFT.30_40.csh`  
 - `batch.wrap_process_veg_hist.pl.USMX.NLCD_INEGI.2011.30_40.csh`  
 - `batch.wrap_process_veg_hist.pl.USMX.NLCD_INEGI.s1992.30_40.csh` - note: for s1992, you need to have run batch scripts for Stages 1-2 for both s1992 and 1992; the land cover fractions from 1992 will be used in s1992.  

To run these batch files, replace all instances of `$PROJECT` in the files with the path to your copy of this GitHub project. Also make sure they are executable by running `chmod +x $FILENAME`, where `$FILENAME` is the path/name of the batch file that you are making executable.

### Option 2: if you are processing your domain as a single stream and file

#### Processing

If you have a small domain and want to process it as a single stream (each stage of processing deals with a single file containing the entire domain), run:

   `process_veg_hist.single_file.pl $AGGROOT $AGGROOT2 $LCID $AGGFILE $STAGELIST $LCID_CV $LCID_OUT $CLEAN $DATA_ROOT $DOMAIN_PFX $PARAM_PFX $LCTYPE`

where

   `$AGGFILE` = the path/name of the output file of Stage 1 (i.e., aggregated MODIS observations over land cover fractions).  
   all other `$*` variables = same definitions as under Option 1 (the 10x10 tile option).  

#### Examples

Examples of running `process_veg_hist.single_file.pl` for a single 10x10 tile for latitudes 30-40 deg N and longitudes -130--120 deg E from the CONUS_MX and USMX domains, respectively, can be found in the batch files under the "examples" directory:

 - batch.process_veg_hist.single_file.pl.CONUS_MX.MODIS.mode_PFT.30_40n.-130_-120e.csh  
 - batch.process_veg_hist.single_file.pl.USMX.NLCD_INEGI.2011.30_40n.-130_-120e.csh  
 - batch.process_veg_hist.single_file.pl.USMX.NLCD_INEGI.s1992.30_40n.-130_-120e.csh  

To run these batch files, replace all instances of `$PROJECT` in the files with the path to your copy of this GitHub project. Also make sure they are executable by running `chmod +x $FILENAME`, where `$FILENAME` is the path/name of the batch file that you are making executable.

## Stage 5: (Optional) Prepare monthly timeseries of MODIS Observations

This stage is only necessary if you wish to drive VIC with an explicit 17-year timeseries of LAI, Fcanopy, and albedo.

1. Break up the time series of monthly-varying files (an output of processing step 6) into 1-year files:

   `breakup_veg_hist_into_annual_files.py -i $INFILE -v $VARNAMELIST -o $OUTDIR -p $PREFIX`

   where

   `$INFILE` = path/filename of the monthly (mosaicked, if you divided the domain into 10x10s) file.  
   `$VARNAMELIST` = comma-separated list of the monthly variables in the file (e.g., "LAI,Fcanopy,albedo").  
   `$OUTDIR` = output directory.  
   `$PREFIX` = prefix of output filenames.  

2. Disaggregate the monthly timeseries to daily:

   `wrap_disagg_veghist_monthly2hourly_nc.pl $INDIR $PREFIX $OUTDIR`

   where

   `$INDIR` = input directory (=output directory of previous step).  
   `$PREFIX` = prefix of input/output filenames.  
   `$OUTDIR` = output directory.  

## Utility Scripts

Some useful utility scripts:
 - To create a domain file for a new domain:

   `create_domain_file_from_asc_inputs.py -m $MASK -f $FRACMASK -d $DEM -p $PRCP_PARAMDIR -o $OUTFILE`

   where

   `$MASK` = ESRI ascii-format mask delineating the domain (1 = in domain; 0 = outside)  
   `$FRACMASK` = ESRI ascii-format grid showing the fraction of each cell's area that lies within the domain  
   `$DEM` = ESRI ascii-format grid of elevation values (Digital Elevation Model)  
   `$PRCP_PARAMDIR` = directory containing ESRI ascii-format grid files of precipitation parameters (2 parameters, "duration" and "peak_time", with 12 monthly mean values each, one file per month)  
   `$OUTFILE` = output path/filename

 - To create initial state files for MetSim:

   `create_metsim_state_file_from_domain_file.py -i $INFILE -s $STATEDATE -o $OUTFILE`

   where

   `$INFILE` = domain file  
   `$STATEDATE` = date of 1 day prior to the desired start; format `$YYYY-$MM-$SDD`  
   `$OUTFILE` = output MetSim state file  

 - To clip out a spatial subset (lat/lon box) of a netcdf file:

   `gridclip.py -i $INFILE -s $SOUTH -n $NORTH -w $WEST -e $EAST -o $OUTFILE`

   where

   `$INFILE` = input file  
   `$SOUTH`,`$NORTH`,`$WEST`,`$EAST` = geographic boundaries of subset  
   `$OUTFILE` = output file  

 - To subsample VIC parameters to a finer resolution AND/OR clip to a new irregularly-shaped domain (via the mask of a domain file):

   `subsample_vic_params_over_domain.py -i $INFILE -d $DOMAINFILE [-n] [-s] [-f] -o $OUTFILE`

   where

   `$INFILE` = input file  
   `$DOMAINFILE` = domain file to clip to (at either the same or finer resolution)  
   `-n` = optional; if specified, use default fill ("nodata") values  
   `-s` = optional; if specified, `$INFILE` contains snow bands  
   `-f` = optional; if specified, `$INFILE` contains Fcanopy  
   `$OUTFILE` = output file  

 - To subsample VIC gridded daily forcings to a finer resolution AND/OR clip to a new irregularly-shaped domain (via the mask of a domain file):

   `subsample_forcings_over_domain.py -i $INFILE -d $DOMAINFILE -v $VARNAMELIST -o $OUTFILE`

   where

   `$INFILE` = input file  
   `$DOMAINFILE` = domain file to clip to (at either the same or finer resolution)  
   `$VARNAMELIST` = (optional) comma-separated list of variables to process; default = "Prec,Tmin,Tmax,wind"  
   `$OUTFILE` = output file  

 - To set the "run_cell" variable equal to the mask of a forcing file (to ensure that VIC doesn't try to run where no forcings are available):

   `set_run_cell.py -p $PARAMFILE -r $RUNCELLFILE -v $VARNAME -o $OUTFILE`

   where

   `$PARAMFILE` = VIC parameter file  
   `$RUNCELLFILE` = file from which to get the values for setting "run_cell", e.g., a domain file, or maybe another VIC parameter file  
   `$VARNAME` = name of variable in `$RUNCELLFILE` the values of which will be assigned to "run_cell"  
   `$OUTFILE` = output file
