#!/usr/bin/env python

import os.path
import numpy as np
from pyhdf.SD import SD, SDC
import sys, getopt
import re
import math
import time

def decide_pft(x):
  x = np.where((x == 254), 17, x)
  x = np.where((x == 255), 18, x)
  meta_classes = {
    0: 0,
    1: 5,
    2: 5,
    3: 5,
    4: 5,
    5: 5,
    6: 7,
    7: 7,
    8: 9,
    9: 9,
    10: 10,
    11: 11,
    12: 12,
    13: 13,
    14: 12,
    15: 15,
    16: 16,
    17: 17, # 254
    18: 18, # 255
  }
  nClasses = len(meta_classes)
  freq = np.full(nClasses,0,dtype='float32')
  freq_meta = np.full(nClasses,0,dtype='float32')
  N = len(x)
  tot = 0
  flag = 0
  for i in range(N):
    if (x[i] < nClasses-2):
      pft = int(x[i])
      freq[pft] += 1
      tot += 1
    else:
      flag = x[i]
  if (tot > 0):
    for i in range(nClasses):
      freq[i] /= tot
    if (np.max(freq) > 0.5):
      pft_out = np.where(freq > 0.5)[0]
    else: # no clear majority
      for i in range(nClasses):
        if (freq[i] > 0):
          freq_meta[meta_classes[i]] += freq[i]
      if (np.max(freq_meta) > 0.5): # does a lumped class have a majority?
        pft_out = np.where(freq_meta > 0.5)[0]
        flag = 1
      else: # look for mode
        flag = 2
        tmp = np.where(freq_meta == np.max(freq_meta))[0]
        pft_out = tmp[0]
        # Hack to favor crops in event of a tie with other covers
        if (freq_meta[12] == np.max(freq_meta)):
          pft_out = 12
        if (pft_out == 5 or pft_out == 7 or pft_out == 9 or pft_out == 10):
          flag = 3
          freq_trees = freq_meta[5] + freq_meta[9]
          freq_grass_shrub = freq_meta[7] + freq_meta[10]
          freq_sav_grass_shrub = freq_meta[7] + freq_meta[9] + freq_meta[10]
          freq_trees_sav_grass_shrub = freq_meta[5] + freq_meta[7] + freq_meta[9] + freq_meta[10]
          if (freq_trees >= 0.5):
            pft_out = 8
          elif (freq_sav_grass_shrub >= 0.5):
            if (freq_meta[9] >= freq_meta[7]+freq_meta[10]):
              pft_out = 9
            elif (freq_meta[7] >= freq_meta[10] ):
              pft_out = 7
            else:
              pft_out = 10
          elif (freq_trees_sav_grass_shrub > 0.5):
            pft_out = 9
          else: # take mode of all
            tmp = np.where(freq == np.max(freq))[0]
            pft_out = tmp[0]
  else:
    pft_out = nClasses-1

  if (pft_out == nClasses-1):
    pft_out = 255
  elif (pft_out == nClasses-2):
    pft_out = 254
  if (flag == nClasses-1):
    flag = 255
  elif (flag == nClasses-2):
    flag = 254

  return (pft_out,flag)



def main():
  infile = ''
  maskfile = ''
  outfile = ''
  try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:f:o:",["indir=","filelist=","outfile="])
  except getopt.GetoptError:
    print(sys.argv[0], ' -i <indir> -f <filelist> -o <outfile>')
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print(sys.argv[0], ' -i <indir> -f <filelist> -o <outfile>')
      sys.exit()
    elif opt in ("-i", "--indir"):
      indir = arg
    elif opt in ("-f", "--filelist"):
      filelist = arg
    elif opt in ("-o", "--outfile"):
      outfile = arg

  # Define some constants
  PFT_varname = 'Land_Cover_Type_1'

  # Parse filelist
  fileList = []
  fileList = filelist.lstrip().rstrip().split(',')
  nFiles = len(fileList)

  # Read files
  for k in range(nFiles):
    filename = indir + '/' + fileList[k]
    if (k == 0):
      print("reading file to get dimensions",filename)
      hdf = SD(filename,SDC.READ)
      nrows_modis = hdf.datasets()[PFT_varname][1][0]
      ncols_modis = hdf.datasets()[PFT_varname][1][1]
      hdf.end()
      data = np.full([nrows_modis,ncols_modis,nFiles],255,dtype='int16')
#      rowvals = np.full([nrows_modis,ncols_modis]),0,dtype='int32')
#      colvals = np.full([nrows_modis,ncols_modis]),0,dtype='int32')
#      for i in range(nrows_modis):
#        rowvals[i,:] = i
#      for j in range(ncols_modis):
#        colvals[:,j] = j
    print("reading file",filename)
    try:
      hdf = SD(filename,SDC.READ)
      tmpdata = hdf.select(PFT_varname)
      data[:,:,k] = tmpdata[:,:]
      hdf.end()
    except:
      data[:,:,k] = data[:,:,k-1] # repeat previous (doesn't work if first file has an error)

  flag = np.full([nrows_modis,ncols_modis],255,dtype='int16')
  pft = np.full([nrows_modis,ncols_modis],255,dtype='int16')
  for i in range(nrows_modis):
    print('processing row',i)
    for j in range(ncols_modis):
      (pft_tmp,flag_tmp) = decide_pft(data[i,j,:])
      pft[i,j] = pft_tmp
      flag[i,j] = flag_tmp
#      if (i < 10 and pft[i,j] != data[i,j,0]):
#        print(i,j,'data',data[i,j,:],'pft',pft[i,j],'flag',flag[i,j])

  hdf = SD(outfile, SDC.WRITE | SDC.CREATE)
  sds_pft = hdf.create(PFT_varname, SDC.INT16, (nrows_modis, ncols_modis))
  sds_flag = hdf.create('Data_Flag', SDC.INT16, (nrows_modis, ncols_modis))
  sds_pft.setfillvalue(255)
  sds_flag.setfillvalue(255)
  dim1 = sds_pft.dim(0)
  dim1.setname('row')
  dim2 = sds_pft.dim(1)
  dim2.setname('col')
  dim3 = sds_flag.dim(0)
  dim3.setname('row')
  dim4 = sds_flag.dim(1)
  dim4.setname('col')
  sds_pft[:,:] = pft[:,:]
  sds_flag[:,:] = flag[:,:]
  sds_pft.endaccess()
  sds_flag.endaccess()
  hdf.end()

if __name__ == "__main__":
  main()

