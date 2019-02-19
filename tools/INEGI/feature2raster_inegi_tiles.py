import arcpy
import os

# If output files already exist, go ahead and overwrite them
# (otherwise we'd get an error)
arcpy.env.overwriteOutput = True

series = arcpy.GetParameterAsText(0)
suffix = arcpy.GetParameterAsText(1)

rootDir = "Z:/ArcGIS/INEGI/SERIE_" + series
inDir = "geo"
if (series == "V"):
    inDir = "geo.clip"
outDir = "raster"
ascDir = "asc"

metatiles = ['ghi11','hi12','fg12','h13','g13','ef13','gh14','f14','de14','def15','ef16']

for k in metatiles:
    featureclass = rootDir + "/" + inDir + "/" + k + suffix + ".shp"
    rasterOut = rootDir + "/" + outDir + "/" + k + suffix + ".tif"
    asciiOut = rootDir + "/" + ascDir + "/" + k + suffix + ".asc"
    arcpy.FeatureToRaster_conversion(featureclass,"ClassNum",rasterOut,0.0004166666666667)
    arcpy.RasterToASCII_conversion(rasterOut, asciiOut)
    

