import arcpy
import os

# If output files already exist, go ahead and overwrite them
# (otherwise we'd get an error)
arcpy.env.overwriteOutput = True

# This script requires user to input one parameter: the year
year = arcpy.GetParameterAsText(0)

# Root directory where input and output files are located
rootDir = "Z:/ArcGIS"
rootDirNLCD = "Z:/ArcGIS/NLCD"
prefix = "nlcd_" + year + "_landcover_2011_edition_2014_10_10"

# Areas to clip
areaList = ["02","04","05","07"]

# Loop over areas
for area in areaList:

    rasterIn = rootDirNLCD + "/" + prefix + "/" + prefix + ".img"
    clipPoly = rootDir + "/bound" + area + ".shp"
    rasterOut = rootDirNLCD + "/" + year + "/clip/area" + area + ".img"
    arcpy.Clip_management(rasterIn,"#",rasterOut,clipPoly,"0","ClippingGeometry")

