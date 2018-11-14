#!/usr/bin/env python

# NOTE this script assumes land cover data are in geographic projection, broken up into 1x1 degree tiles

import os.path
import numpy as np
from pyhdf.SD import SD, SDC
from netCDF4 import Dataset
import sys, getopt
import re
import math
import time

def compute_day_of_year (year, month, day):

  # month and day are numeric but can be strings (0-padded numbers)
  # month should range from 1 (jan) to 12 (dec)

  month_days = [31,28,31,30,31,30,31,31,30,31,30,31]
  if year % 4 == 0:
    month_days[2] = 29

  yday = int(day)
  mmax = int(month) - 1
  for m in range(mmax):
    yday += month_days[m]

  return yday


def modis_sinusoidal_to_latlon (h,v,npix_per_deg_y,row,col):

  PI = 3.14159265358979328462
  minlat = 80 - 10*v
  maxlat = 90 - 10*v
  cellsize_y = 1/npix_per_deg_y
  latcenter = maxlat - (row/npix_per_deg_y + 0.5*cellsize_y)
  if (v < 9):
    lat_for_x = latcenter - 0.5*cellsize_y
  else:
    lat_for_x = latcenter + 0.5*cellsize_y
  npix_per_deg_x = npix_per_deg_y*math.cos(PI/180*lat_for_x)
  # number of columns between 0 longitude and 180 longitude for this row
  ncols_map_half_width = 180*npix_per_deg_x
  cellsize_x = 1/(npix_per_deg_x)

  if (h < 18):
    minlon = 10*(h-18)*npix_per_deg_y*cellsize_x
    loncenter = minlon + (col+0.5)*cellsize_x
  else:
    minlon = 10*(h-18)*npix_per_deg_y*cellsize_x
    loncenter = minlon + (col+0.5)*cellsize_x

  return [latcenter,loncenter,cellsize_y,cellsize_x]


def latlon_to_modis_sinusoidal (lat, lon, npix_per_deg_y):

  PI = 3.14159265358979328462
  v = int((90-lat)/10)
  row = int((90-(v*10)-lat)*npix_per_deg_y)
  npix_per_deg_x = npix_per_deg_y*math.cos(PI/180*lat)
  # number of columns between 0 longitude and 180 longitude for this row
  ncols_map_half_width = int(180*npix_per_deg_x)
  if (lon >= 0):
    col = int(ncols_map_half_width*lon/180 + 0.5)
  else:
    col = -1 - int(ncols_map_half_width*math.fabs(lon)/180 + 0.5)
  col = col + 180*npix_per_deg_y
  h = int(col/(npix_per_deg_y*10))
  col = col-h*10*npix_per_deg_y

  return [h,v,row,col]


def compute_LAI (ndvi):

  minNDVI = 0.1
  critNDVI = 0.6
  # A = 2 / (critNDVI-minNDVI)**2, i.e., when NDVI=critNDVI, LAI=2
  A = 8

  ndvi = max(ndvi, minNDVI)
  LAI = A * (ndvi - minNDVI) * (ndvi - minNDVI)

  return LAI


def compute_fcan (ndvi):

  minNDVI = 0.1
  maxNDVI = 0.8
  minFcan = 0.01

  ndvi = max(ndvi, minNDVI)
  ndvi = min(ndvi, maxNDVI)

  fcan = (ndvi-minNDVI)/(maxNDVI-minNDVI)
  fcan = fcan*fcan

  fcan = max(fcan, minFcan)

  return fcan


def limit_fcan_for_reasonable_lai (fcan, lai):

  # specific_lai = lai/fcan; specific_lai must be less than 4 + 16 * fcan; i.e. specific_lai limited to 4 for extremely sparse veg and 20 for completely dense veg
  specific_lai = 4 + 16 * fcan

  if (lai > specific_lai * fcan):
    # solve for fcan
    fcan = (math.sqrt(1 + 4*lai) - 1) / 8

  return fcan


def main():
  infile = ''
  maskfile = ''
  outfile = ''
  try:
    opts, args = getopt.getopt(sys.argv[1:],"hc:i:l:t:m:d:f:b:s:n:w:e:r:o:",["lctype=","lcid=","lcmaplist_file=","lctable=","mask=","domain=","filelist_file=","hv_boundlist=","south=","north=","west=","east=","resolution=","outpathpfx="])
  except getopt.GetoptError:
    print(sys.argv[0], ' -c <lctype> -i <lcid> -l <lcmaplist_file> -t <lctable> -m <landmask> -d <domain> -f <filelist_file> -b <hv_boundlist> -s <south> -n <north> -w <west> -e <east> -r <resolution> -o <outpathpfx>')
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print(sys.argv[0], ' -c <lctype> -i <lcid> -l <lcmaplist_file> -t <lctable> -m <landmask> -d <domain> -f <filelist_file> -b <hv_boundlist> -s <south> -n <north> -w <west> -e <east> -r <resolution> -o <outpathpfx>')
      sys.exit()
    elif opt in ("-c", "--lctype"):
      lctype = arg
    elif opt in ("-i", "--lcid"):
      lcid = arg
    elif opt in ("-l", "--lcmaplist_file"):
      lc_maplist_file = arg
    elif opt in ("-t", "--lctable"):
      lc_table = arg
    elif opt in ("-m", "--mask"):
      landmask = arg
    elif opt in ("-d", "--domain"):
      domain = arg
    elif opt in ("-f", "--filelist_file"):
      filelist_file = arg
    elif opt in ("-b", "--hv_boundlist"):
      hv_boundlist = arg
    elif opt in ("-s", "--south"):
      if (len(arg.split('.')) > 1):
        minlat = float(arg)
      else:
        minlat = int(arg)
    elif opt in ("-n", "--north"):
      if (len(arg.split('.')) > 1):
        maxlat = float(arg)
      else:
        maxlat = int(arg)
    elif opt in ("-w", "--west"):
      if (len(arg.split('.')) > 1):
        minlon = float(arg)
      else:
        minlon = int(arg)
    elif opt in ("-e", "--east"):
      if (len(arg.split('.')) > 1):
        maxlon = float(arg)
      else:
        maxlon = int(arg)
    elif opt in ("-r", "--resolution"):
      resolution = float(arg)
    elif opt in ("-o", "--outpathpfx"):
      outpathpfx = arg


  # Parse hv_boundlist
  h_list = []
  v_list = []
  tmpList = hv_boundlist.lstrip().rstrip().split(',')
  hmin = int(tmpList[0])
  hmax = int(tmpList[1])
  vmin = int(tmpList[2])
  vmax = int(tmpList[3])
  nh = hmax-hmin+1
  nv = vmax-vmin+1
  nTiles = nh*nv
  for i in range(nh):
    h = hmin+i
    for j in range(nv):
      v = vmin+j
      h_list.append(h)
      v_list.append(v)

  # Define some constants
  LAI_varname = 'Lai_500m'
  LAI_qc_varname = 'FparLai_QC'
  LAI_qc2_varname = 'FparExtra_QC'
  NDVI_varname = '500m 16 days NDVI'
  alb_varname = 'Albedo_WSA_shortwave'
  fill_value_LAI = 255.0
  fill_value_LAI_qc = 255
  fill_value_NDVI = -3000.0
  fill_value_alb = 32767.0

  # Read lc_table
  lc_classes = {}
  lc_LAI_use_LAI_qc = {}
  lc_NDVI_use_LAI_qc = {}
  lc_use_LAI_NDVI_rel = {}
  with open (lc_table) as f:
    for line in f:
      tmpList = line.lstrip().rstrip().split(',')
      if tmpList[0][0] == '#':
        continue
      lc_class = int(tmpList[0])
      lc_classes[lc_class] = tmpList[4]
      lc_LAI_use_LAI_qc[lc_class] = bool(int(tmpList[5]))
      lc_NDVI_use_LAI_qc[lc_class] = bool(int(tmpList[6]))
      lc_use_LAI_NDVI_rel[lc_class] = bool(int(tmpList[7]))
  nClasses = len(lc_classes)

  pNrows = re.compile('nrows',re.I)
  pNcols = re.compile('ncols',re.I)
  pXllc = re.compile('xllcorner',re.I)
  pYllc = re.compile('yllcorner',re.I)
  pCellsize = re.compile('cellsize',re.I)
  pNodata = re.compile('nodata(_val)?',re.I)


  # Initialize output data structures
  nrows_out = int(math.ceil((maxlat-minlat)/resolution))
  ncols_out = int(math.ceil((maxlon-minlon)/resolution))
  data_out = {}
  nodata_out = {}
  keys_out = ['mask','Cv','LAI','NDVI','albedo','fcanopy']
  for key in keys_out:
    if (key == 'mask'):
      nodata_out[key] = int(0)
      data_out[key] = np.full([nrows_out,ncols_out], nodata_out[key], dtype=np.int)
    else:
      nodata_out[key] = float(-1)
      data_out[key] = {}
      if (key == 'Cv'):
        data_out[key]['count'] = np.zeros([nrows_out,ncols_out]).astype(int)
      for lc_class in sorted(lc_classes.keys()):
        data_out[key][lc_class] = {}
        data_out[key][lc_class]['count'] = np.zeros([nrows_out,ncols_out]).astype(int)
        data_out[key][lc_class]['mean'] = np.zeros([nrows_out,ncols_out]).astype(float)
  mask_varname_out = 'LandMask'
  Cv_varname_out = 'Cv'
  LAI_varname_out = 'LAI'
  NDVI_varname_out = 'NDVI'
  fcan_varname_out = 'fcanopy'
  alb_varname_out = 'albedo'

  # Initial time for benchmarking
  clock_new = time.clock()
  time_new = time.time()
  clock_old = clock_new
  time_old = time_new

  # Read land mask
  if (landmask != 'null'):
    print('reading landmask file',landmask)
    with open (landmask) as f:

      # read file
      row_mask = -6
      for line in f:
        tmpList = line.lstrip().rstrip().split(' ')

        # read header to compute lat/lons of rows/cols
        if row_mask < 0:
          if re.match(pNrows,tmpList[0]):
            nrows_mask = int(tmpList[1])
          elif re.match(pNcols,tmpList[0]):
            ncols_mask = int(tmpList[1])
          elif re.match(pXllc,tmpList[0]):
            xllcorner_mask = float(tmpList[1])
          elif re.match(pYllc,tmpList[0]):
            yllcorner_mask = float(tmpList[1])
          elif re.match(pCellsize,tmpList[0]):
            cellsize_mask = float(tmpList[1])
          elif re.match(pNodata,tmpList[0]):
            nodata_mask = tmpList[1]

        # read data lines
        else:
          lat_mask = yllcorner_mask + (nrows_mask - 1 - row_mask + 0.5) * cellsize_mask
          row_out = nrows_out - 1 - int((maxlat - lat_mask)/resolution)

          for col_mask in range(len(tmpList)):

            lon_mask = xllcorner_mask + (col_mask + 0.5) * cellsize_mask
            col_out = int((lon_mask - minlon)/resolution)
            if (lat_mask >= minlat and lat_mask <= maxlat and lon_mask >= minlon and lon_mask <= maxlon):
              tmp_mask_data = int(tmpList[col_mask])
              if (tmp_mask_data != int(nodata_mask)):
                data_out['mask'][row_out,col_out] = tmp_mask_data

        row_mask += 1

    # Benchmark time
    clock_new = time.clock()
    time_new = time.time()
    print('seconds elapsed:','clock:',clock_new-clock_old,'wall:',time_new-time_old)
    clock_old = clock_new
    time_old = time_new

  # Read lc_maplist_file
  with open (lc_maplist_file) as f:
    for line in f:
      PFTfiles_list = line.lstrip().rstrip().split(',')
      break
  nPFTfiles = len(PFTfiles_list)


  # Read filelist_file
  date_list = []
  LAIfiles = {}
  NDVIfiles = {}
  albfiles = {}
  with open (filelist_file) as f:
    for line in f:
      tmpList = line.lstrip().rstrip().split(' ')
      date_list.append(tmpList[0])
      LAIfiles[tmpList[0]] = tmpList[1]
      NDVIfiles[tmpList[0]] = tmpList[2]
      albfiles[tmpList[0]] = tmpList[3]
  year0 = date_list[0][0:4]
  ndates = len(date_list)

  # Open output file and write header
  outfile = outpathpfx + '.' + lcid + '.nc'
  print('creating output file',outfile)
  newrootgrp = Dataset(outfile, 'w', format='NETCDF4')
  timedim = newrootgrp.createDimension('time', None)
  classdim = newrootgrp.createDimension('veg_class', nClasses)
  latdim = newrootgrp.createDimension('lat', nrows_out)
  londim = newrootgrp.createDimension('lon', ncols_out)
  timevar = newrootgrp.createVariable('time', 'i4', ('time'))
  timevar.units = 'days since ' + year0 + '-01-01 00:00:00.0'
  yearvar = newrootgrp.createVariable('year', 'i4', ('time'))
  yearvar.units = 'year'
  dayvar = newrootgrp.createVariable('day_of_year', 'i4', ('time'))
  dayvar.units = 'day of year'
  latvar = newrootgrp.createVariable('lat', 'f4', ('lat'))
  latvar.units = 'degrees_north'
  latvar.standard_name = 'latitude'
  latvar.long_name = 'latitude of grid cell center'
  lonvar = newrootgrp.createVariable('lon', 'f4', ('lon'))
  lonvar.units = 'degrees_east'
  lonvar.standard_name = 'longitude'
  lonvar.long_name = 'longitude of grid cell center'
  latsvar = newrootgrp.createVariable('lats', 'f4', ('lat','lon'))
  latsvar.units = 'degrees'
  latsvar.standard_name = 'lats'
  latsvar.long_name = 'lats'
  latsvar.description = 'Latitude of grid cell center'
  lonsvar = newrootgrp.createVariable('lons', 'f4', ('lat','lon'))
  lonsvar.units = 'degrees'
  lonsvar.standard_name = 'lons'
  lonsvar.long_name = 'lons'
  lonsvar.description = 'Longitude of grid cell center'
  classvar = newrootgrp.createVariable('veg_class', 'i2', ('veg_class'))
  classvar.units = '-'
  classvar.standard_name = 'land cover class ID code'
  classvar.long_name = 'land cover class ID code'
  classname = newrootgrp.createVariable('class_name', 'str', ('veg_class'))
  classname.units = '-'
  classname.standard_name = 'land cover class name'
  classname.long_name = 'land cover class name'

  # Create data variables
  mask_out = newrootgrp.createVariable(mask_varname_out, 'i2', ('lat','lon'), fill_value=nodata_out['mask'])
  mask_out.units = '-'
  mask_out.long_name = 'Land Mask'
  mask_out.missing_value = nodata_out['mask']
  Cv_out = newrootgrp.createVariable(Cv_varname_out, 'f4', ('veg_class','lat','lon'), fill_value=nodata_out['Cv'])
  Cv_out.units = '-'
  Cv_out.long_name = 'Area Fraction'
  Cv_out.missing_value = nodata_out['Cv']
  LAI_out = newrootgrp.createVariable(LAI_varname_out, 'f4', ('time','veg_class','lat','lon'), fill_value=nodata_out['LAI'])
  LAI_out.units = 'm2/m2'
  LAI_out.long_name = 'Leaf Area Index'
  LAI_out.missing_value = nodata_out['LAI']
  NDVI_out = newrootgrp.createVariable(NDVI_varname_out, 'f4', ('time','veg_class','lat','lon'), fill_value=nodata_out['NDVI'])
  NDVI_out.units = '-'
  NDVI_out.long_name = 'Normalized Difference Vegetation Index'
  NDVI_out.missing_value = nodata_out['NDVI']
  fcan_out = newrootgrp.createVariable(fcan_varname_out, 'f4', ('time','veg_class','lat','lon'), fill_value=nodata_out['fcanopy'])
  fcan_out.units = '-'
  fcan_out.long_name = 'Canopy fraction'
  fcan_out.missing_value = nodata_out['fcanopy']
  alb_out = newrootgrp.createVariable(alb_varname_out, 'f4', ('time','veg_class','lat','lon'), fill_value=nodata_out['albedo'])
  alb_out.units = '-'
  alb_out.long_name = 'Albedo'
  alb_out.missing_value = nodata_out['albedo']

  # Populate dimension variables except time
  for row in range(nrows_out):
    latvar[row] = minlat + (row + 0.5) * resolution
  for col in range(ncols_out):
    lonvar[col] = minlon + (col + 0.5) * resolution
  for row in range(nrows_out):
    for col in range(ncols_out):
      latsvar[row,col] = latvar[row]
      lonsvar[row,col] = lonvar[col]
  i = 0
  for lc_class in sorted(lc_classes.keys()):
    classvar[i] = lc_class
    classname[i] = lc_classes[lc_class]
    i += 1
  if (landmask != 'null'):
    mask_out[:,:] = data_out['mask'][:,:]



  # Benchmark time
  clock_new = time.clock()
  time_new = time.time()
  print('seconds elapsed:','clock:',clock_new-clock_old,'wall:',time_new-time_old)
  clock_old = clock_new
  time_old = time_new


  # Write dates to output file
  offset = 0
  year_save = -1
  for t in range(ndates):
    year = int(date_list[t][0:4])
    day = int(date_list[t][4:7])
    yearvar[t] = year
    dayvar[t] = day
    if (year != year_save and year_save >= 0):
      offset += 365
      if (year_save % 4 == 0): 
        offset += 1
    timevar[t] = day + offset - 1 # subtract 1 since it's days since Jan 1 of first year
    year_save = year
 
  # Benchmark time
  clock_new = time.clock()
  time_new = time.time()
  print('seconds elapsed:','clock:',clock_new-clock_old,'wall:',time_new-time_old)
  clock_old = clock_new
  time_old = time_new

  nrows_modis_save = 2400
  ncols_modis_save = 2400
  day_save = -7

  # Loop over dates
  for t in range(len(date_list)):

    # Re-initialize data_out
    for key in keys_out:
      if (key != 'mask'):
        if (key == 'Cv'):
          data_out[key]['count'] = np.zeros([nrows_out,ncols_out]).astype(int)
        for lc_class in sorted(lc_classes.keys()):
          data_out[key][lc_class]['count'] = np.zeros([nrows_out,ncols_out]).astype(int)
          data_out[key][lc_class]['mean'] = np.zeros([nrows_out,ncols_out]).astype(float)

    # Split the 3 file lists into arrays
    LAIfiles_list = LAIfiles[date_list[t]].split(",")
    NDVIfiles_list = NDVIfiles[date_list[t]].split(",")
    albfiles_list = albfiles[date_list[t]].split(",")
    nTiles = len(LAIfiles_list)

    # Open first MODIS file to get dimensions
    # this whole script would benefit from creating a dict for files and filelists where LAI, NDVI, alb are keys, so that we can reuse code
    k = 0
    nrows_modis = -1
    while nrows_modis < 0 and k < nTiles:
      if (LAIfiles_list[k] != 'null'):
        print("reading file to get dimensions",LAIfiles_list[k])
        hdf = SD(LAIfiles_list[k],SDC.READ)
        nrows_modis = hdf.datasets()[LAI_varname][1][0]
        ncols_modis = hdf.datasets()[LAI_varname][1][1]
        hdf.end()
        tmpList = re.split('\/',LAIfiles_list[k])
        tmpList2 = re.split('\.',tmpList.pop())
        day = tmpList2[1][5:8]
        ndays = int(day)
      elif (NDVIfiles_list[k] != 'null'):
        print("reading file to get dimensions",NDVIfiles_list[k])
        hdf = SD(NDVIfiles_list[k],SDC.READ)
        nrows_modis = hdf.datasets()[NDVI_varname][1][0]
        ncols_modis = hdf.datasets()[NDVI_varname][1][1]
        hdf.end()
        tmpList = re.split('\/',NDVIfiles_list[k])
        tmpList2 = re.split('\.',tmpList.pop())
        day = tmpList2[1][5:8]
        ndays = int(day)
      elif (albfiles_list[k] != 'null'):
        print("reading file to get dimensions",albfiles_list[k])
        hdf = SD(albfiles_list[k],SDC.READ)
        nrows_modis = hdf.datasets()[alb_varname][1][0]
        ncols_modis = hdf.datasets()[alb_varname][1][1]
        hdf.end()
        tmpList = re.split('\/',albfiles_list[k])
        tmpList2 = re.split('\.',tmpList.pop())
        day = tmpList2[1][5:8]
        ndays = int(day)
      k += 1

    if (nrows_modis < 0):
      nrows_modis = nrows_modis_save
      ncols_modis = ncols_modis_save
      if (k == 0):
        day = day_save + 8
      else:
        day = day_save
    else:
      nrows_modis_save = nrows_modis
      ncols_modis_save = ncols_modis
      day_save = day

    # Dimensions - assume same across all files
    npix_per_deg_y = int(nrows_modis/10)
    LAI_data = np.full([nTiles, nrows_modis, ncols_modis], fill_value_LAI, dtype=np.single)
    LAI_qc_data = np.full([nTiles, nrows_modis, ncols_modis], fill_value_LAI_qc, dtype=np.int)
    LAI_qc2_data = np.full([nTiles, nrows_modis, ncols_modis], fill_value_LAI_qc, dtype=np.int)
    NDVI_data = np.full([nTiles, nrows_modis, ncols_modis], fill_value_LAI, dtype=np.single)
    alb_data = np.full([nTiles, nrows_modis, ncols_modis], fill_value_LAI, dtype=np.single)

    # Open all MODIS files and store in memory
    for k in range(0,nTiles):
 
      # LAI
      if (LAIfiles_list[k] != 'null'):
        print("reading file",LAIfiles_list[k])
        hdf = SD(LAIfiles_list[k],SDC.READ)
        tmpdata = hdf.select(LAI_varname)
        scale_factor_LAI = tmpdata.attributes()['scale_factor']
        add_offset_LAI = tmpdata.attributes()['add_offset']
        valid_range_LAI = tmpdata.attributes()['valid_range']
        fill_value_LAI = tmpdata.attributes()['_FillValue']
        LAI_data[k,:,:] = tmpdata[:,:]
        for val in [252,253,254]:
          LAI_data[k, LAI_data[k] == val] = 0
        for val in [249,250,251]:
          LAI_data[k, LAI_data[k] == val] = fill_value_LAI
        LAI_data[k, LAI_data[k] < 0] = fill_value_LAI
        LAI_data[k, LAI_data[k] != fill_value_LAI] -= add_offset_LAI
        LAI_data[k, LAI_data[k] != fill_value_LAI] *= scale_factor_LAI
        LAI_data[k, LAI_data[k] < 0] = fill_value_LAI
        tmpdata = hdf.select(LAI_qc_varname)
        LAI_qc_data[k,:,:] = tmpdata[:,:]
        tmpdata = hdf.select(LAI_qc2_varname)
        LAI_qc2_data[k,:,:] = tmpdata[:,:]
        hdf.end()
 
      # NDVI
      if (NDVIfiles_list[k] != 'null'):
        print("reading file",NDVIfiles_list[k])
        hdf = SD(NDVIfiles_list[k],SDC.READ)
        tmpdata = hdf.select(NDVI_varname)
        scale_factor_NDVI = tmpdata.attributes()['scale_factor']
        add_offset_NDVI = tmpdata.attributes()['add_offset']
        valid_range_NDVI = tmpdata.attributes()['valid_range']
        fill_value_NDVI = tmpdata.attributes()['_FillValue']
        NDVI_data[k,:,:] = tmpdata[:,:]
        NDVI_data[k, NDVI_data[k] < 0] = fill_value_NDVI
        NDVI_data[k, NDVI_data[k] != fill_value_NDVI] -= add_offset_NDVI
#        NDVI_data[k, NDVI_data[k] != fill_value_NDVI] *= scale_factor_NDVI
        NDVI_data[k, NDVI_data[k] != fill_value_NDVI] /= scale_factor_NDVI
        NDVI_data[k, NDVI_data[k] < 0] = fill_value_NDVI
        hdf.end()
      else:
        NDVI_data[k,:,:] = np.full([nrows_modis,ncols_modis], fill_value_NDVI, dtype=np.float)
 
      # Albedo
      if (albfiles_list[k] != 'null'):
        print("reading file",albfiles_list[k])
        hdf = SD(albfiles_list[k],SDC.READ)
        tmpdata = hdf.select(alb_varname)
        scale_factor_alb = tmpdata.attributes()['scale_factor']
        add_offset_alb = tmpdata.attributes()['add_offset']
        valid_range_alb = tmpdata.attributes()['valid_range']
        fill_value_alb = tmpdata.attributes()['_FillValue']
        alb_data[k,:,:] = tmpdata[:,:]
        alb_data[k, alb_data[k] < 0] = fill_value_alb
        alb_data[k, alb_data[k] != fill_value_alb] -= add_offset_alb
        alb_data[k, alb_data[k] != fill_value_alb] *= scale_factor_alb
        alb_data[k, alb_data[k] < 0] = fill_value_alb
        hdf.end()
      else:
        alb_data[k,:,:] = np.full([nrows_modis,ncols_modis], fill_value_alb, dtype=np.float)

    # Compute relationship between modis row/col and out row/col
    row_out_of_row_modis = np.empty([nTiles,nrows_modis]).astype(int)
    col_out_of_rowcol_modis = np.empty([nTiles,nrows_modis,ncols_modis]).astype(int)
    for k in range(nTiles):
      for row_modis in range(nrows_modis):
        returnList = modis_sinusoidal_to_latlon(h_list[k],v_list[k],npix_per_deg_y,row_modis,0)
        lat_modis = returnList[0]
        row_out_of_row_modis[k,row_modis] = nrows_out - 1 - int((maxlat - lat_modis)/resolution)
        for col_modis in range(ncols_modis):
          returnList = modis_sinusoidal_to_latlon(h_list[k],v_list[k],npix_per_deg_y,row_modis,col_modis)
          lat_modis = returnList[0]
          lon_modis = returnList[1]
          col_out_of_rowcol_modis[k,row_modis,col_modis] = int((lon_modis - minlon)/resolution)

    # Benchmark time
    clock_new = time.clock()
    time_new = time.time()
    print('seconds elapsed:','clock:',clock_new-clock_old,'wall:',time_new-time_old)
    clock_old = clock_new
    time_old = time_new


    # Read land cover map and summarize at MODIS resolution
    if (t == 0):

      # Define and initialize land cover map structure
      lc_class_count = {}
      for k in range(nTiles):
        lc_class_count[k] = {}
        for row_modis in range(nrows_modis):
          lc_class_count[k][row_modis] = {}
          for col_modis in range(ncols_modis):
            lc_class_count[k][row_modis][col_modis] = {}

      # loop over lc files
      first = 1
      ind = []
      ind_dict = {}
      for k in range(nPFTfiles):

        if (PFTfiles_list[k] != 'null'):

          print("reading landcover file",PFTfiles_list[k])

          if lctype == 'modis':

            # HDF files with data on MODIS grid
            PFT_data = np.full([nPFTfiles, nrows_modis, ncols_modis], fill_value_PFT, dtype=np.int)
            hdf = SD(PFTfiles_list[k],SDC.READ)
            tmpdata = hdf.select(PFT_varname)
            PFT_data[k,:,:] = tmpdata[:,:]
            hdf.end()

            for row_modis in range(nrows_modis):
              for col_modis in range(ncols_modis):

                lc_class = PFT_data[k,row_modis,col_modis]
                # HACK!!!!!!!!!!!!!!!!!!!!
                if (lc_class == 254):
                  lc_class = 16
                row_out = row_out_of_row_modis[k,row_modis]
                col_out = col_out_of_rowcol_modis[k,row_modis,col_modis]

                if (row_out >= 0 and row_out < nrows_out and col_out >= 0 and col_out < ncols_out and lc_class != nodata_lc and lc_class >= 0 and (landmask == 'null' or data_out['mask'][row_out,col_out] != nodata_out['mask'])):

                  ind.append([k,row_modis,col_modis])

                  if col_out in lc_class_count[k][row_modis][col_modis]:
                    if lc_class in lc_class_count[k][row_modis][col_modis][col_out]:
                      lc_class_count[k][row_modis][col_modis][col_out][lc_class] += 1
                    else:
                      lc_class_count[k][row_modis][col_modis][col_out][lc_class] = 1
                    lc_class_count[k][row_modis][col_modis][col_out]['total'] += 1
                  else:
                    lc_class_count[k][row_modis][col_modis][col_out] = {}
                    lc_class_count[k][row_modis][col_modis][col_out][lc_class] = 1
                    lc_class_count[k][row_modis][col_modis][col_out]['total'] = 1
          elif lctype == 'asc':

            # Ascii grid files with data on lat/lon grid
            with open (PFTfiles_list[k]) as f:

              # read file
              row_tmp = -6
              for line in f:
                tmpList = line.lstrip().rstrip().split(' ')

                # read header to compute lat/lons of rows/cols
                if row_tmp < 0:
                  if re.match(pNrows,tmpList[0]):
                    nrows_tmp = int(tmpList[1])
                  elif re.match(pNcols,tmpList[0]):
                    ncols_tmp = int(tmpList[1])
                  elif re.match(pXllc,tmpList[0]):
                    xllcorner_tmp = float(tmpList[1])
                  elif re.match(pYllc,tmpList[0]):
                    yllcorner_tmp = float(tmpList[1])
                  elif re.match(pCellsize,tmpList[0]):
                    cellsize_lc = float(tmpList[1])
                  elif re.match(pNodata,tmpList[0]):
                    nodata_lc = int(tmpList[1])
                    if (first):
                      nrows_lc = int((maxlat-minlat)/cellsize_lc + 0.5)
                      ncols_lc = int((maxlon-minlon)/cellsize_lc + 0.5)

                # read data lines
                else:

                  lat_lc = yllcorner_tmp + (nrows_tmp - 1 - row_tmp + 0.5) * cellsize_lc
                  row_out = nrows_out - 1 - int((maxlat - lat_lc)/resolution)
                  if (row_out < 0):
                    row_out = 0
                  elif (row_out > nrows_out-1):
                    row_out = nrows_out-1

                  for col_tmp in range(len(tmpList)):

                    lon_lc = xllcorner_tmp + (col_tmp + 0.5) * cellsize_lc
                    col_out = int((lon_lc - minlon)/resolution)
                    if (col_out < 0):
                      col_out = 0
                    elif (col_out > ncols_out-1):
                      col_out = ncols_out-1

                    returnList = latlon_to_modis_sinusoidal(lat_lc,lon_lc,npix_per_deg_y)

                    this_h = returnList[0]
                    this_v = returnList[1]
                    row_modis = returnList[2]
                    col_modis = returnList[3]
                    if (row_modis < 0):
                      row_modis = 0
                    elif (row_modis > nrows_modis-1):
                      row_modis = nrows_modis-1
                    if (col_modis < 0):
                      col_modis = 0
                    elif (col_modis > ncols_modis-1):
                      col_modis = ncols_modis-1

                    lc_class = int(tmpList[col_tmp])

                    this_k = -1
                    for k in range(nTiles):
                      if (h_list[k] == this_h and v_list[k] == this_v):
                        this_k = k
                        break
                    if (this_k == -1):
                      continue

                    if (lc_class != nodata_lc and lc_class >= 0 and (landmask == 'null' or data_out['mask'][row_out,col_out] != nodata_out['mask'])):

                      if k not in ind_dict.keys():
                        ind_dict[k] = {}
                      if row_modis not in ind_dict[k].keys():
                        ind_dict[k][row_modis] = {}
                      if col_modis not in ind_dict[k][row_modis].keys():
                        ind_dict[k][row_modis][col_modis] = 1
                        ind.append([k,row_modis,col_modis])

                      if col_out in lc_class_count[this_k][row_modis][col_modis]:
                        if lc_class in lc_class_count[this_k][row_modis][col_modis][col_out]:
                          lc_class_count[this_k][row_modis][col_modis][col_out][lc_class] += 1
                        else:
                          lc_class_count[this_k][row_modis][col_modis][col_out][lc_class] = 1
                        lc_class_count[this_k][row_modis][col_modis][col_out]['total'] += 1
                      else:
                        lc_class_count[this_k][row_modis][col_modis][col_out] = {}
                        lc_class_count[this_k][row_modis][col_modis][col_out][lc_class] = 1
                        lc_class_count[this_k][row_modis][col_modis][col_out]['total'] = 1

                row_tmp += 1

          first = 0


    # Benchmark time
    clock_new = time.clock()
    time_new = time.time()
    print('seconds elapsed:','clock:',clock_new-clock_old,'wall:',time_new-time_old)
    clock_old = clock_new
    time_old = time_new


    # Initialize output data structures
    print('initializing output data structures')
    for key in keys_out:
      if (key != 'mask'):
        if (t == 0):
          data_out[key] = {}
          if (key == 'Cv'):
            data_out[key]['count'] = np.zeros([nrows_out,ncols_out]).astype(int)
        for lc_class in sorted(lc_classes.keys()):
          if (t == 0):
            data_out[key][lc_class] = {}
          data_out[key][lc_class]['count'] = np.zeros([nrows_out,ncols_out]).astype(int)
          data_out[key][lc_class]['mean'] = np.zeros([nrows_out,ncols_out]).astype(float)

    # Benchmark time
    clock_new = time.clock()
    time_new = time.time()
    print('seconds elapsed:','clock:',clock_new-clock_old,'wall:',time_new-time_old)
    clock_old = clock_new
    time_old = time_new

    # Loop over MODIS grid
    print('looping over MODIS grid')
    for k,row_modis,col_modis in ind:
      row_out = row_out_of_row_modis[k,row_modis]
      for col_out in sorted(lc_class_count[k][row_modis][col_modis].keys()):

        # Populate running stats but filter out nodatas
        lc_class_count_total = lc_class_count[k][row_modis][col_modis][col_out]['total']
        if (lc_class_count_total > 0 and data_out['mask'][row_out,col_out] != nodata_out['mask']):

          for lc_class_tmp in (lc_class_count[k][row_modis][col_modis][col_out].keys()):

            if (lc_class_tmp == 'total'):
              if (t == 0):
                data_out['Cv']['count'][row_out,col_out] += lc_class_count_total
            else:
              lc_class = int(lc_class_tmp)
              if (lc_class not in lc_classes.keys()):
                continue
              lc_class_count_per_class = lc_class_count[k][row_modis][col_modis][col_out][lc_class]

              if (t == 0):
                data_out['Cv'][lc_class]['count'][row_out,col_out] += lc_class_count_per_class
                if (landmask == 'null' and lc_class_count_per_class > 0):
                  data_out['mask'][row_out,col_out] = 1
              bin = '{0:08b}'.format(LAI_qc_data[k,row_modis,col_modis])
              bin2 = '{0:08b}'.format(LAI_qc2_data[k,row_modis,col_modis])
              goodqc = True
              if (bin[7] == '1' or bin[5] == '1' or bin[3:5] == '01' or bin[3:5] == '10' or bin2[1] == '1' or bin2[5] == '1' or bin2[3] == '1' or bin2[2] == '1' or bin2[1] == '1'):
                goodqc = False
              LAI_tmp = fill_value_LAI
              if (LAI_data[k,row_modis,col_modis] != fill_value_LAI
                  and (goodqc or lc_LAI_use_LAI_qc[lc_class] == False)
                  and (lc_use_LAI_NDVI_rel[lc_class] == False)):
                LAI_tmp = LAI_data[k,row_modis,col_modis]
                data_out['LAI'][lc_class]['count'][row_out,col_out] += lc_class_count_per_class
                data_out['LAI'][lc_class]['mean'][row_out,col_out] += LAI_tmp * lc_class_count_per_class
              if (NDVI_data[k,row_modis,col_modis] != fill_value_NDVI
                  and (goodqc or lc_NDVI_use_LAI_qc[lc_class] == False)):
                data_out['NDVI'][lc_class]['count'][row_out,col_out] += lc_class_count_per_class
                data_out['NDVI'][lc_class]['mean'][row_out,col_out] += NDVI_data[k,row_modis,col_modis] * lc_class_count_per_class
                if lc_use_LAI_NDVI_rel[lc_class]:
                  LAI_tmp = compute_LAI(NDVI_data[k,row_modis,col_modis])
                  data_out['LAI'][lc_class]['count'][row_out,col_out] += lc_class_count_per_class
                  data_out['LAI'][lc_class]['mean'][row_out,col_out] += LAI_tmp * lc_class_count_per_class
                if (LAI_tmp != fill_value_LAI):
                  tmp_fcan = compute_fcan(NDVI_data[k,row_modis,col_modis])
                  tmp_fcan = limit_fcan_for_reasonable_lai(tmp_fcan, LAI_tmp)
                  data_out['fcanopy'][lc_class]['count'][row_out,col_out] += lc_class_count_per_class
                  data_out['fcanopy'][lc_class]['mean'][row_out,col_out] += tmp_fcan * lc_class_count_per_class
              if (alb_data[k,row_modis,col_modis] != fill_value_alb
                  and (goodqc or lc_NDVI_use_LAI_qc[lc_class] == False)):
                data_out['albedo'][lc_class]['count'][row_out,col_out] += lc_class_count_per_class
                data_out['albedo'][lc_class]['mean'][row_out,col_out] += alb_data[k,row_modis,col_modis] * lc_class_count_per_class

    if (landmask == 'null'):
      mask_out[:,:] = data_out['mask'][:,:]

    # Benchmark time
    clock_new = time.clock()
    time_new = time.time()
    print('seconds elapsed:','clock:',clock_new-clock_old,'wall:',time_new-time_old)
    clock_old = clock_new
    time_old = time_new

    # Compute means
    print('compute means')
    for key in keys_out:
      if (key == 'mask'):
        continue
      for lc_class in sorted(lc_classes.keys()):
        for row in range(nrows_out):
          for col in range(ncols_out):
            if (key == 'Cv'):
              if (t == 0):
                if (data_out[key]['count'][row,col] > 0):
                  data_out[key][lc_class]['mean'][row,col] = data_out[key][lc_class]['count'][row,col] / data_out[key]['count'][row,col]
                else:
                  data_out[key][lc_class]['mean'][row,col] = nodata_out[key]
            else:
              if (data_out[key][lc_class]['count'][row,col] > 0):
                data_out[key][lc_class]['mean'][row,col] = data_out[key][lc_class]['mean'][row,col] / data_out[key][lc_class]['count'][row,col]
              else:
                data_out[key][lc_class]['mean'][row,col] = nodata_out[key]

    # Benchmark time
    clock_new = time.clock()
    time_new = time.time()
    print('seconds elapsed:','clock:',clock_new-clock_old,'wall:',time_new-time_old)
    clock_old = clock_new
    time_old = time_new

    # Write to output files
    print('writing to',outfile)
    i = 0
    for lc_class in sorted(lc_classes.keys()):
      if (t == 0):
        Cv_out[i,:,:] = data_out['Cv'][lc_class]['mean'][:,:]
      LAI_out[t,i,:,:] = data_out['LAI'][lc_class]['mean'][:,:]
      NDVI_out[t,i,:,:] = data_out['NDVI'][lc_class]['mean'][:,:]
      fcan_out[t,i,:,:] = data_out['fcanopy'][lc_class]['mean'][:,:]
      alb_out[t,i,:,:] = data_out['albedo'][lc_class]['mean'][:,:]
      i += 1

    # Benchmark time
    clock_new = time.clock()
    time_new = time.time()
    print('seconds elapsed:','clock:',clock_new-clock_old,'wall:',time_new-time_old)
    clock_old = clock_new
    time_old = time_new


  # close output file
  newrootgrp.close()

if __name__ == "__main__":
  main()

