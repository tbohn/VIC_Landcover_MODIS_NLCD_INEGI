#!/usr/bin/env python

import os.path
import numpy as np
import xarray as xr
import sys, getopt

def main():
  infile = ''
  outfile = ''
  try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:s:n:w:e:o:",["infile=","south=","north=","west=","east=","outfile="])
  except getopt.GetoptError:
    print(sys.argv[0], ' -i <infile> -s <south> -n <north> -w <west> -e <east> -o <outfile>')
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print(sys.argv[0], ' -i <infile> -s <south> -n <north> -w <west> -e <east> -o <outfile>')
      sys.exit()
    elif opt in ("-i", "--infile"):
      infile = arg
    elif opt in ("-s", "--south"):
      south = float(arg)
    elif opt in ("-n", "--north"):
      north = float(arg)
    elif opt in ("-w", "--west"):
      west = float(arg)
    elif opt in ("-e", "--east"):
      east = float(arg)
    elif opt in ("-o", "--outfile"):
      outfile = arg


  ds = xr.open_dataset(infile)

  latvar_tmp = ds['lat']
  lonvar_tmp = ds['lon']

  nLat = len(latvar_tmp)
  nLon = len(lonvar_tmp)
  minLat = round(np.asscalar(latvar_tmp[0]),6)
  maxLat = round(np.asscalar(latvar_tmp[-1]),6)
  minLon = round(np.asscalar(lonvar_tmp[0]),6)
  maxLon = round(np.asscalar(lonvar_tmp[-1]),6)
  resolution = (maxLat-minLat)/(nLat-1)
  minLat -= 0.5*resolution
  maxLat += 0.5*resolution
  minLon -= 0.5*resolution
  maxLon += 0.5*resolution
  y0 = int((south-minLat)/resolution)
  y1 = int((north-minLat)/resolution)
  x0 = int((west-minLon)/resolution)
  x1 = int((east-minLon)/resolution)

  ds_out = ds.isel(lat=slice(y0,y1),lon=slice(x0,x1))

  # Write to output file
  ds_out.to_netcdf(outfile, engine='scipy')

  ds.close()

if __name__ == "__main__":
  main()

