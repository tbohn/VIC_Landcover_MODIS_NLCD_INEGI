import arcpy
import os

series = arcpy.GetParameterAsText(0)
suffix = arcpy.GetParameterAsText(1)

rootDir = "Z:/ArcGIS/INEGI/SERIE_" + series
rootDir3 = "Z:/ArcGIS/INEGI/SERIE_III"

metatiles = ['ghi11','hi12','fg12','h13','g13','ef13','gh14','f14','de14','def15','ef16']

for k in metatiles:
    featureclass = rootDir + "/merge/" + k + suffix + ".shp"
    prjfile = rootDir3 + "/merge/" + k  + "v.prj"
    arcpy.DefineProjection_management(featureclass,prjfile)
    

