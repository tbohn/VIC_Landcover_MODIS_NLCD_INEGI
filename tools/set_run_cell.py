#!/usr/bin/env python

import os.path
import numpy as np
import xarray as xr
import sys, getopt

def main():
    paramfile = ''
    runcellfile = ''
    varname = ''
    outfile = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hp:r:v:o:",["paramfile=","runcellfile=","varname=","outfile="])
    except getopt.GetoptError:
        print(sys.argv[0], ' -p <paramfile> -r <runcellfile> -v <varname> -o <outfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -p <paramfile> -r <runcellfile> -v <varname> -o <outfile>')
            sys.exit()
        elif opt in ("-p", "--paramfile"):
            paramfile = arg
        elif opt in ("-r", "--runcellfile"):
            runcellfile = arg
        elif opt in ("-v", "--varname"):
            varname = arg
        elif opt in ("-o", "--outfile"):
            outfile = arg

    # Open runcell file
    ds_run = xr.open_dataset(runcellfile)
    runcell = np.where(ds_run[varname] == 1, 1, 0).astype(np.int32)
    ds_run.close()

    # Open input param file
    ds = xr.open_dataset(paramfile)

    fill_value = ds['run_cell'].encoding['_FillValue']
    mask_tmp = np.where(ds['mask'] == 1, 1, 0).astype(np.int32)
    runcell[mask_tmp == 0] = fill_value
    ds['run_cell'][:] = runcell[:]

    # Write to output file
    ds.to_netcdf(outfile)

    ds.close()

if __name__ == "__main__":
    main()

