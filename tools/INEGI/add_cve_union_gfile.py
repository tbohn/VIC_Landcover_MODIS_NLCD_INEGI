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

# Max allowable string lengths for the fields we're adding to the shapefile
fieldLength = 7
fieldLength2 = 10

# List of the metatiles we're going to process
metatiles = ["ghi11","hi12","fg12","h13","g13","ef13","gh14","f14","de14","def15","ef16"]

# Loop over metatiles
for k in metatiles:

    # Define the input "v" file and output file
    featureclassG = rootDir + "/merge/" + k + "g.shp"
    featureclassOut = rootDir + "/cve_union/" + k + "g.shp"

    # Copy the "g" file to the output file
    arcpy.CopyFeatures_management(featureclassG,featureclassOut)

    # Add fields to the output file
    arcpy.AddField_management(featureclassOut,"ClassName","TEXT","","",fieldLength2)
    arcpy.AddField_management(featureclassOut,"ClassNum","LONG")

    # Read all of the polygons of the output file into memory as a list called "rows"
    # Here we use "UpdateCursor" function so that if we modify "rows", it modifies
    # the output file
    rows = arcpy.UpdateCursor(featureclassOut)

    # Loop over all the polygons. Note that row is a copy of the
    # current polygon in rows.  Modifying row won't have any effect
    # on rows, without an updateRow command.
    for row in rows:
        
        # Store the value of this row's CLAVE field in the variable "code"
        if (series == "III"):
            fieldname = "CLAVE"
        elif (series == "IV"):
            fieldname = "CVE_G"
        code = row.getValue(fieldname)

        # Compute values of ClassName and ClassNum
        ClassName = "Missing"
        ClassNum = -1
        if (code == "P/E"):
            ClassName = "Missing"
            ClassNum = -1
        elif (code == "H2O"):
            ClassName = "Water"
            ClassNum = 0
        elif (code == "AH" or code == "ZU"):
            ClassName = "Urban"
            ClassNum = 2
        elif (code == "ADV" or code == "DV"):
            ClassName = "Bare"
            ClassNum = 3
        elif (code[0:1] == "P"):
            if (code == "PH" or code == "PN"):
                ClassName = "Grassland"
                ClassNum = 6
            else:
                ClassName = "Pasture"
                ClassNum = 7
        elif (code == "[R]" or code[0:1] == "H" or code[0:1] == "R" or code[0:1] == "T" or code == "IAPF"):
            ClassName = "Cropland"
            ClassNum = 8
        elif (code == "ACUI"):
            ClassName = "Wetland"
            ClassNum = 9

        # Store the ClassName and ClassNum
        row.ClassName = ClassName
        row.ClassNum = ClassNum

        # Transfer the modified row to the rows list (which is what gets written)
        rows.updateRow(row)
