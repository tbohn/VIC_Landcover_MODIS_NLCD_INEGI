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
if (series != "V"):
    metatiles = ["ghi11","hi12","fg12","h13","g13","ef13","gh14","f14","de14","def15","ef16"]
else:
    metatiles = ["usv250s5_union"]

# Loop over metatiles
for k in metatiles:

    # Define the input "v" file and output file
    if (series != "V"):
        featureclassV = rootDir + "/merge/" + k + "v.shp"
        featureclassOut = rootDir + "/cve_union/" + k + ".shp"
        inField = "CLAVEFOT"
    else:
        featureclassV = rootDir + "/cve_union.old/" + k + ".shp"
        featureclassOut = rootDir + "/cve_union/" + k + ".shp"
        inField = "CVE_UNION"

    # Copy the "v" file to the output file
    arcpy.CopyFeatures_management(featureclassV,featureclassOut)

    # Add fields to the output file
    if (series != "V"):
        arcpy.AddField_management(featureclassOut,"CVE_UNION","TEXT","","",fieldLength)
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
        
        # Store the value of this row's CLAVEFOT field in the variable "code"
        code = row.getValue(inField)

        # Strip out the erosion flag
        if (code[0:2] == "E-"):
            code = code[2:]

        # Move the "VS*" codes to the front and remove the "M*" codes
        if (code != "P/E"):
            words = code.split("/")
            if (len(words) > 1 and (words[1] == "VSa" or words[1] == "VSA" or words[1] == "VSh")):
               code = words[1] + "/" + words[0]
            elif (words[0] == "VSa" or words[0] == "VSA" or words[0] == "VSh"):
               code = words[0] + "/" + words[1]
            else:
               code = words[0]

        # Convert Re code to R
        if (code[0:2] == "Re"):
            code = "R" + code[2:]

        # Convert the 4-character agric. codes to 3-characters
        if (code[0:1] == "R" or (code[0:1] == "H" and code[1:1] != "2") or code[0:1] == "T"):
            code = code[0:3]
            if (code[1:2] == "PA"):
               code = code[0:1] + "AP"
            elif (code[1:2] == "SA"):
               code = code[0:1] + "AS"
            elif (code[1:2] == "PS"):
                code = code[0:1] + "SP"

        if (series != "V"):
            # Store the new codes in the "CVE_UNION" field
            row.CVE_UNION = code
            #row.setValue("CVE_UNION",code) # this also works

        # Compute values of ClassName and ClassNum
        ClassName = "Missing"
        ClassNum = -1
        if (code[0:3] == "VSa" or code[0:3] == "VSA" or code[0:3] == "VSh"):
            words = code.split("/")
            code = words[1]
        if (code == "P/E"):
            ClassName = "Missing"
            ClassNum = -1
        elif (code == "H2O"):
            ClassName = "Water"
            ClassNum = 0
        elif (code == "AH" or code == "ZU"):
            ClassName = "Urban"
            ClassNum = 2
        elif (code == "ADV" or code == "DV" or code == "VD"):
            ClassName = "Bare"
            ClassNum = 3
        elif (code[0:1] == "B" or code[0:1] == "S" or code == "MST"):
            ClassName = "Forest"
            ClassNum = 4
        elif (code[0:1] == "M"):
            ClassName = "Shrubland"
            ClassNum = 5
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
        elif (code == "ACUI" or code == "VA" or code == "VG"
                 or code == "VH" or code == "VHH" or code == "VM"
                 or code == "VT"):
            ClassName = "Wetland"
            ClassNum = 9
        elif (code == "VU" or code == "VY"):
            ClassName = "Shrubland"
            ClassNum = 5
        elif (code[0:2] == "VP" or code [0:2] == "VS"):
            ClassName = "Forest"
            ClassNum = 4
        elif (code == "VW"):
            ClassName = "Grassland"
            ClassNum = 6

        # Store the ClassName and ClassNum
        row.ClassName = ClassName
        row.ClassNum = ClassNum

        # Transfer the modified row to the rows list (which is what gets written)
        rows.updateRow(row)
