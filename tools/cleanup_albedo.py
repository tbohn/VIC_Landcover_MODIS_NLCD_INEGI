#!/usr/bin/env python

import os.path
import numpy as np
import xarray as xr
import sys, getopt
import math
import time
import pandas as pd

def main():
    infile = ''
    outfile = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:o:",["infile=","outfile="])
    except getopt.GetoptError:
        print(sys.argv[0], ' -i <infile> -o <outfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -i <infile> -o <outfile>')
            sys.exit()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-o", "--outfile"):
            outfile = arg

    # Open input file
    ds = xr.open_dataset(infile)

    latvar = ds['lat']
    lonvar = ds['lon']
    classvar = ds['veg_class']
    classname = ds['class_name']
    timevar = ds['time']
    yearvar = ds['year']
    dayvar = ds['day_of_year']
    mask = ds['LandMask']
    Cv = ds['Cv']

    varname = 'albedo'
    data = ds[varname]
    nTime = data.shape[0]
    nClass = data.shape[1]
    nLat = data.shape[2]
    nLon = data.shape[3]
#    fill_value = ds[varname].attrs['_FillValue']
    fill_value = np.asscalar(np.full([1],-1,dtype=np.single))

    # For each class, compute albedo stats and filter out positive outliers
    for c in range(nClass):
        tmp = np.empty([nTime,nLat,nLon],dtype=np.single)
        tmp2 = np.empty([nTime,nLat,nLon],dtype=np.single)
        tmp[:] = data[:,c]
        thresh = np.empty([nLat,nLon],dtype=np.single)
        thresh = np.nanpercentile(tmp,75,axis=0)
        tmp2[:] = tmp[:]
        tmp2[tmp2>thresh] = np.nan
        mean = np.nanmean(tmp2,axis=0)
        std = np.nanstd(tmp2,axis=0)
        thresh = mean+4*std
        tmp2[:] = tmp[:]
        tmp2[tmp2>thresh] = fill_value
        data[:,c] = tmp2

    ds[varname] = data

    for varname in ['LAI','fcanopy','albedo']:
        ds[varname].encoding['zlib'] = True

    # Write to output file
    ds.to_netcdf(outfile, engine='scipy')

    ds.close()

if __name__ == "__main__":
    main()

