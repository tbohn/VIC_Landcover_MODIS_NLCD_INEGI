import arcpy
import os

# If output files already exist, go ahead and overwrite them
# (otherwise we'd get an error)
arcpy.env.overwriteOutput = True

# Get user parameters
year = arcpy.GetParameterAsText(0)

# Root directory where input and output files are located
rootDir = "Z:/ArcGIS/NLCD/" + year

# Areas to process
areaList = ["02","04","05","07"]

# Loop over areas
for area in areaList:

    rasterIn = rootDir + "/clip/area" + area + ".img"
    rasterOut = rootDir + "/geo/area" + area + ".img"
    arcpy.ProjectRaster_management(rasterIn,rasterOut,'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')

