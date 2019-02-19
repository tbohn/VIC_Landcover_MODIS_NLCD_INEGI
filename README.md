# VIC_Landcover_MODIS_NLCD_INEGI

This project contains several versions of land cover parameters for the VIC land surface model, derived from MODIS observations and land cover classifications from MODIS, NLCD, and INEGI over the US, Mexico, and Southern Canada, as well as the scripts used to generate them.

## Directory Structure:

./
 - License.txt - GNU public license
 - README.md - this file

docs/
 - Overview.md - overview of whole project
 - INEGI/ - processing of INEGI Usos de Suelo y Vegetacion maps

 ...- "Processing of INEGI USOSV dataset.docx" - describes the processing of INEGI maps, including the names of the scripts involved in each processing step

 - NLCD_INEGI/ - further processing and merging of NLCD and INEGI maps

 ...- Scripts.md - list of scripts, including description and usage
 ...- Procedure.md - list of processing steps

 - MODIS/ - aggregation of MODIS data over land cover maps and generation of VIC parameter files

 ...- Scripts.md - list of scripts, including description and usage
 ...- Procedure.md - list of processing steps

examples/ - the inputs and intermediate files for an example 10x10 degree "tile"
 - batch.process_one_tile.40_50n.130_120e.sh - linux shell script that contains the commands for processing an example 10x10 tile
 - inputs/ - input files for the example
 - outputs/ - intermediate and final output files for the example

tools/ - processing scripts are stored here; for details, see the files under docs/
 - INEGI/ - ArcPy scripts (used in conjunction with ESRI ArcGIS) for processing the INEGI Usos de Suelo y Vegetacion maps
 - NLCD_INEGI/ - scripts for further processing and merging NLCD and INEGI maps
 - MODIS/ - scripts for aggregating MODIS data over land cover maps and generation of VIC parameter files
