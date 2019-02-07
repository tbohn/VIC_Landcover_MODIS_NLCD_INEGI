# VIC_Landcover_MODIS_NLCD_INEGI

This project contains several versions of land cover parameters for the VIC land surface model, derived from MODIS observations and land cover classifications from MODIS, NLCD, and INEGI over the US, Mexico, and Southern Canada, as well as the scripts used to generate them.

## Directory Structure:

./
 - License.txt - GNU public license
 - README.md - this file

docs/
 - Scripts.md - list of scripts, including description and usage
 - Procedure.md - list of processing steps

examples/ - the inputs and intermediate files for an example 10x10 degree "tile"
 - batch.process_one_tile.40_50n.130_120e.sh - linux shell script that contains the commands for processing an example 10x10 tile
 - inputs/ - input files for the example
 - outputs/ - intermediate and final output files for the example

tools/
 - processing scripts are stored here; for list and descriptions, see the files under docs/
