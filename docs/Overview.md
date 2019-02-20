# Overview of processing

Creating VIC land surface parameters requires (1) enumerating the area fractions of land cover classes within each grid cell and (2) assigning physical properties to those land cover classes.

To accomplish this, we need to draw on (1) land cover maps at finer resolution than the model grid; (2) spatial datasets of physical properties, also at finer resolution than the model grid; and (3) literature-based values of physical properties on a per-class basis (uniform in space).

This project focused on improving (1) land cover maps and (2) spatial datasets of physical properties. It draws on several land cover maps over North America and MODIS observations of surface properties. For land cover maps, we had two main sources: (a) the MOD12Q1.051 MODIS-based classification of Friedl et al.  (2010), at 500-m resolution, covering the period 2000-2013 with global coverage; and (b) LANDSAT- based 30-m classifications of the National Land Cover Database (NLCD; Homer et al., 2015) over the US and the Usos del Suelo y Vegetacion dataset of Mexico's Instituto Nacional de Estadistica y Geografia (INEGI, 2015), with "snapshots" from years 1992/3, 2001/2, and 2011.

Accordingly, processing fell into two main stages: (1) processing of land cover maps and (2) processing of MODIS products, relative to the land cover maps. Spatially uniform properties of each class were copied from the Livneh et al (2015) VIC parameter set.

Processing stages:
 1. Processing of INEGI Usos del Suelo y Vegetacion maps to (a) reconcile their classification schemes to the scheme of NLCD 2011; (b) reproject to geographic and rasterize them at approximately 30m resolution. **This assumes that you have the original INEGI shapefiles on your system.** The INEGI shapefiles can be found at xxx.
 2. Further processing and merging of the NLCD and INEGI datasets into a single geographic projection raster at approximately 30m resolution. **This assumes that you have the original NLCD maps and the rasterized INEGI maps on your system.** The NLCD maps can be found at xxx. The rasterized INEGI maps can be found at xxx.
 3. Aggregation of MODIS 500-m LAI, NDVI, and albedo products over either (a) the MOD12Q1.051 land cover map or the processed NLCD_INEGI maps from various years. **This assumes that you have (a) all necessary MODIS products and (b) either the MOD12Q1.051 land cover product or the processed NLCD and INEGI maps on your system.** MODIS products can be found at xxx. NLCD and INEGI maps must be the result of processing stage 2.

There are separate directories under docs/ corresponding to these 3 stages (INEGI, NLCD_INEGI, and MODIS, respectively), and corresponding directories under tools/ pertaining to these same stages. Detailed documentation of each stage can be found in the docs/*/ subdirectories.

