#!/usr/bin/env python

import os.path
import numpy as np
import xarray as xr
import sys, getopt

def main():
  infile = ''
  outfile = ''
  try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:c:o:",["infile=","cvfile=","outfile="])
  except getopt.GetoptError:
    print(sys.argv[0], ' -i <infile> -c <cvfile> -o <outfile>')
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print(sys.argv[0], ' -i <infile> -c <cvfile> -o <outfile>')
      sys.exit()
    elif opt in ("-i", "--infile"):
      infile = arg
    elif opt in ("-c", "--cvfile"):
      cvfile = arg
    elif opt in ("-o", "--outfile"):
      outfile = arg

  ds = xr.open_dataset(infile)
  ds_cv = xr.open_dataset(cvfile)

  ds_out = ds.copy()

  ds_out['Cv'][:] = ds_cv['Cv'][:]

  # Write to output file
  ds_out.to_netcdf(outfile)

  ds.close()

if __name__ == "__main__":
  main()

