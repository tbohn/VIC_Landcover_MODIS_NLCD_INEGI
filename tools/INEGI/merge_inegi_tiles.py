import arcpy
import os

series = arcpy.GetParameterAsText(0)

rootDir = "Z:/ArcGIS/INEGI/SERIE_" + series

if (series == "I"):
    seriesNum = "1"
else if (series == "II"):
    seriesNum = "2"
else if (series == "III"):
    seriesNum = "3"
else if (series == "IV"):
    seriesNum = "4"

if (series == "I" or series == "II"):
    datum = "n27"
    suffixes = ['v']
    rootDirPRJ = "Z:/ArcGIS/INEGI/SERIE_I"
else:
    datum = "i92"
    suffixes = ['g','v']
    rootDirPRJ = "Z:/ArcGIS/INEGI/SERIE_III"

dirpfx = rootDir + "/conjuntoUsoSueloVegetacion_RF_250k_" + datum + "_SHP_" + series + "_"

metatiles = {'ghi11':['i1111','i1112','h1102','h1103','h1105','h1106','h1109','h1112','g1103','g1106'],
             'hi12':['i1210','h1201','h1202','h1203','h1204','h1205','h1206','h1207','h1208','h1209','h1210','h1211','h1212'],
             'fg12':['g1201','g1202','g1203','g1204','g1205','g1206','g1207','g1208','g1209','g1210','g1211','g1212','f1202','f1203','f1205','f1206'],
             'h13':['h1301','h1302','h1304','h1305','h1307','h1308','h1309','h1310','h1311','h1312'],
             'g13':['g1301','g1302','g1303','g1304','g1305','g1306','g1307','g1308','g1309','g1310','g1311','g1312'],
             'ef13':['f1301','f1302','f1303','f1304','f1305','f1306','f1307','f1308','f1309','f1311','f1312','e1302','e1303','e1305','e1306','e1309'],
             'gh14':['h1407','h1410','h1411','g1401','g1402','g1404','g1405','g1406','g1407','g1408','g1409','g1410','g1411','g1412'],
             'f14':['f1401','f1402','f1403','f1404','f1405','f1406','f1407','f1408','f1409','f1410','f1411','f1412'],
             'de14':['e1401','e1402','e1403','e1404','e1405','e1406','e1407','e1408','e1409','e1410','e1411','e1412','d1403'],
             'def15':['f1509','f1512','e1501','e1503','e1504','e1505','e1506','e1507','e1508','e1509','e1510','e1511','e1512','d1501','d1502','d1503','d1505'],
             'ef16':['f1607','f1608','f1610','f1611','e1601','e1602','e1604','e1605','e1607']}

for suffix in suffixes:
    for k, v in metatiles.iteritems():
        featureclassOut = rootDir + "/merge/" + k + suffix + ".shp"
        featureclassList = []
        for tile in metatiles[k]:
            filename = dirpfx + tile.upper() + "/" + tile + "u" + seriesNum + suffix + "/shape/" + tile + "u" + seriesNum + suffix + ".shp"
            featureclassList.append(filename)
        arcpy.Merge_management(featureclassList,featureclassOut)
        prjfile = rootDirPRJ + "/" + k + "v.prj"
        arcpy.DefineProjection_management(featureclassOut,prjfile)
    

