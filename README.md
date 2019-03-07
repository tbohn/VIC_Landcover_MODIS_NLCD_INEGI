# VIC_Landcover_MODIS_NLCD_INEGI

This project contains the scripts used to created MODIS-based land surface parameters [(MOD-LSP)](https://zenodo.org/record/2559631) for the VIC model (Liang et al., 1994) release 5.0 and later (Hamman et al., 2018). These parameters are intended for use in VIC's "image" driver.

## This project had several goals:
1. update the land cover classification underlying the land surface parameters widely used in VIC modeling over North America
2. provide a set of parameters that can be used in land cover change studies
3. supply estimates of canopy fraction for the recently developed canopy fraction feature in VIC (Bohn and Vivoni, 2016)
4. provide improved, spatially explicit estimates of time-varying surface properties (LAI, canopy fraction, albedo) over the domain
5. provide values of these properties for a wider range of classes such as urban, which had been ignored in prior parameter sets

The VIC land surface parameters consist of gridded values of: soil properties, fractional area coverage of land cover classes, and physical properties of the surface and vegetation structure for each class. Thus, to achieve the project's goals required new land cover classifications (both to update the map and to obtain multi-decadal snapshots of land cover change) and gridded observations of surface properties. These surface properties needed to be aggregated over the pixels of the land cover maps to yield separate spatial average values for each land cover class in each grid cell.

The chosen domain was the continental US, Mexico, and southern Canada (the CONUS + Mexico, or CONUS_MX domain), at 1/16 degree (6 km) spatial resolution. This is the domain used in the Livneh et al (2015) (L2015 hereafter) gridded daily meteorology dataset, and these parameters are designed to be compatible with that dataset.

## Processing Steps
The processing steps covered by the scripts in this dataset are:
1. Downloading of MODIS data
2. Aggregation of MODIS land surface properties over the desired land cover classification
3. Gap-filling
4. Replacement of the L2015 land cover fractions and land surface properties with the MOD-LSP values

In addition, there are various utility scripts (e.g., clipping to smaller region of interest, subsampling, etc).

## Inputs:
### Land Cover Classifications:
 - MOD12Q1.005 MODIS-based classification of Friedl et al (2010) at 500-m resolution, for years 2001-2013.  This cand be downloaded from the [MODIS site](https://lpdaac.usgs.gov/dataset_discovery/modis).
 - NLCD_INEGI land cover classifications (Bohn and Vivoni 2019a) for years 1992, 2001, and 2011. These are a combination of the National Land Cover Database (NLCD; Homer et al., 2015) over the United States and a modification of the INEGI Uso del Suelo y Vegetacion (INEGI, 2014) classification over Mexico that was harmonized with the NLCD legend. These are available for download at [Zenodo](https://zenodo.org/record/2580428).

### Land Surface Properties:
MODIS phenology available from [MODIS site](https://lpdaac.usgs.gov/dataset_discovery/modis)
 - MOD15A2H.006 (2000-02-18 to 2001-06-26) and MCD15A2H.006 (2001-07-04 to end of 2016) Leaf Area Index (LAI) (Myneni et al., 2002)
 - MOD43A3.006 White Sky Albedo from shortwave broadband (Schaaf et al., 2002)
 - MOD13A1.006 Normalized Different Vegetation Index (NDVI) (Huete et al., 2002)

### Miscellaneous Files:
 - NetCDF format VIC-5 image driver parameter file from the L2015 simulations (soil parameters were taken directly from this file) available from [Zenodo](https://zenodo.org/record/2564019)
 - Land cover class tables listing the classes in the underlying classification, along with ID codes and flags indicating how each class should be processed (examples of which are included in the ./examples directory)

## Outputs:
 - Several versions of land cover parameters for the VIC land surface model, available for download from [Zenodo](https://zenodo.org/record/2559631).

These parameter files were designed for use with VIC 5 (image driver). VIC 5 image driver requires a "domain" file to accompany the parameter file. This domain file is also necessary for disaggregating the daily gridded meteorological forcings to hourly for input to VIC via the disaggregating tool MetSim (https://github.com/UW-Hydro/MetSim).  We have provided a domain file compatible with the L2015 forcings and the MOD-LSP parameters, on Zenodo (https://zenodo.org/record/2564019).


## Directory Structure:

./
 - License.txt - GNU public license
 - README.md - this file

docs/
 - Overview.md - overview of whole project
 - Scripts.md - list of scripts, including description and usage
 - Procedure.md - list of processing steps

examples/ - the inputs and intermediate files for an example 10x10 degree "tile"
 - batch.process_one_tile.40_50n.130_120e.sh - linux shell script that contains the commands for processing an example 10x10 tile
 - inputs/ - input files for the example
 - outputs/ - intermediate and final output files for the example

tools/ - processing scripts are stored here; for details, see the files under docs/
 - scripts for aggregating MODIS data over land cover maps and generation of VIC parameter files
