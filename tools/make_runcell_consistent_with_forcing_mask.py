#!/usr/bin/env python

import os.path
import numpy as np
import xarray as xr
import sys, getopt

def main():
    paramfile = ''
    maskfile = ''
    outfile = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hp:m:o:",["paramfile=","maskfile=","outfile="])
    except getopt.GetoptError:
        print(sys.argv[0], ' -p <paramfile> -m <maskfile> -o <outfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -p <paramfile> -m <maskfile> -o <outfile>')
            sys.exit()
        elif opt in ("-p", "--paramfile"):
            paramfile = arg
        elif opt in ("-m", "--maskfile"):
            maskfile = arg
        elif opt in ("-o", "--outfile"):
            outfile = arg

    # Open mask file
    ds_mask = xr.open_dataset(maskfile)
    mask = np.where(ds_mask['mask'] == 1, 1, 0).astype(np.int32)
    ds_mask.close()

    # Open input param file
    ds = xr.open_dataset(paramfile)

    fill_value = ds['run_cell'].encoding['_FillValue']
    mask_tmp = np.where(ds['mask'] == 1, 1, 0).astype(np.int32)
    mask_intersection = mask * mask_tmp
    mask_intersection[mask_intersection == 0] = fill_value
    ds['run_cell'][:] = mask_intersection[:]

    # Write to output file
    ds.to_netcdf(outfile)

    ds.close()

if __name__ == "__main__":
    main()

