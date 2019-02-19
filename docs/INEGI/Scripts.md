# Scripts for processing original INEGI Usos del Suelo maps

Original INEGI Usos del Suelo maps are in ESRI polygon shapefile format, and
the polygons are classified according to slightly different schemes in each
year, each with many classes. Processing scripts first reconcile the
classifications with the NLCD 2011 classification, then reproject to a
geographic projection, and then convert the polygons to raster maps at
approximately the same resolution as the NLCD maps.

add_cve_union_gfile.py
add_cve_union.py
add_cve_union_tmp.py
clip_areas_nlcd.py
copy_shapefiles.py
defproj_inegi_tiles.py
feature2raster_inegi_tiles.py
merge_inegi_tiles.py
raster2ascii_nlcd.py
reproj_inegi_tiles_geo.py
reproj_nlcd_tiles_geo.py
