import arcpy
import os

# If output files already exist, go ahead and overwrite them
# (otherwise we'd get an error)
arcpy.env.overwriteOutput = True

series = arcpy.GetParameterAsText(0)
suffix = arcpy.GetParameterAsText(1)

rootDir = "Z:/ArcGIS/INEGI/SERIE_" + series

metatiles = ['ghi11','hi12','fg12','h13','g13','ef13','gh14','f14','de14','def15','ef16']

for k in metatiles:
    featureclass = rootDir + "/cve_union/" + k + suffix + ".shp"
    featureclassOut = rootDir + "/geo/" + k + suffix + ".shp"
    arcpy.Project_management(featureclass,featureclassOut,'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')
    

