import arcpy
import os

# If output files already exist, go ahead and overwrite them
# (otherwise we'd get an error)
arcpy.env.overwriteOutput = True

# This script requires user to input one parameter: the series number
# (which needs to be a roman numeral, i.e. "I", "II", "III", or "IV")
series = arcpy.GetParameterAsText(0)

# Root directory where input and output files are located
# (directory name contains the series number)
rootDir = "Z:/ArcGIS/INEGI/SERIE_" + series

# List of the metatiles we're going to process
metatiles = ['ghi11','hi12','fg12','h13','g13','ef13','gh14','f14','de14','def15','ef16']

# Loop over metatiles
for k in metatiles:

    # Define the input and output files
    featureclassV = rootDir + "/merge.tmp/" + k + ".shp"
    featureclassOut = rootDir + "/merge/" + k + "v"

    # Copy input to output
    arcpy.CopyFeatures_management(featureclassV,featureclassOut)

    # Extra "g" file for series 3 and 4
    if (series == "III" or series == "IV"):
        featureclassG = rootDir + "/merge.tmp/" + k + "g.shp"
        featureclassOut = rootDir + "/merge/" + k + "g"
        arcpy.CopyFeatures_management(featureclassG,featureclassOut)
