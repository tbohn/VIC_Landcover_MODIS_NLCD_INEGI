import arcpy
import os

# If output files already exist, go ahead and overwrite them
# (otherwise we'd get an error)
arcpy.env.overwriteOutput = True

# This script requires user to input one parameter: the year
year = arcpy.GetParameterAsText(0)

# Root directory where input and output files are located
rootDir = "Z:/ArcGIS/NLCD/" + year

# Areas to process
areaList = ["02","04","05","07"]

# Loop over areas
for area in areaList:

    rasterIn = rootDir + "/geo/area" + area + ".img"
    asciiOut = rootDir + "/asc.geo/area" + area + ".asc"
    arcpy.RasterToASCII_conversion(rasterIn,asciiOut)
