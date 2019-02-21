#!/usr/bin/env python

import os.path
import numpy as np
import xarray as xr
import sys, getopt

def compute_ndvi (fcan):

    nMin = 0.1
    nMax = 0.8

    ndvi = nMin + np.sqrt(fcan) * (nMax - nMin)

    return ndvi


def compute_fcan (ndvi):

    nMin = 0.1
    nMax = 0.8
    nThresh = nMin + 0.75 * (nMax - nMin)
    fcanMin = 0.01

    ndvi = max(ndvi, nMin)
    ndvi = min(ndvi, nMax)

    if ndvi >= nThresh:
        # Parabolic
        fcan = (ndvi - nMin) / (nMax - nMin)
        fcan = fcan*fcan
    else:
        # Linear
        tmp = (nThresh - nMin) / (nMax - nMin)
        fcan = tmp*tmp * (ndvi - nMin) / (nThresh - nMin)

    fcan = max(fcan, fcanMin)

    return fcan


def main():
    infile = ''
    outfile = ''
    recompute_fcan = False
    try:
        opts, args = getopt.getopt(sys.argv[1:],"ha:c:v:fo:",["anomfile=","climfile=","varnamelist=","outfile="])
    except getopt.GetoptError:
        print(sys.argv[0], ' -a <anomfile> -c <climfile> -v <varnamelist> [-f] -o <outfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -a <anomfile> -c <climfile> -v <varnamelist> [-f] -o <outfile>')
            sys.exit()
        elif opt in ("-a", "--anomfile"):
            anomfile = arg
        elif opt in ("-c", "--climfile"):
            climfile = arg
        elif opt in ("-v", "--varnamelist"):
            varnamelist = arg
            varnames = varnamelist.lstrip().rstrip().split(',')
        elif opt in ("-f"):
            recompute_fcan = True
        elif opt in ("-o", "--outfile"):
            outfile = arg


    # Valid ranges
    valid_min = {
        'LAI': 0,
        'NDVI': -1,
        'fcanopy': 0,
        'albedo': 0,
    }
    valid_max = {
        'LAI': 20,
        'NDVI': 1,
        'fcanopy': 1,
        'albedo': 1,
    }

    ds = xr.open_dataset(anomfile)
    # set day_of_year to be a coordinate variable
    ds = ds.set_coords('day_of_year')

    latvar = ds['lat']
    lonvar = ds['lon']
    latsvar = ds['lats']
    lonsvar = ds['lons']
    classvar = ds['veg_class']
    classname = ds['class_name']
    timevar = ds['time']
    yearvar = ds['year']
    dayvar = ds['day_of_year']
    mask = ds['LandMask']
    Cv = ds['Cv']

    ds2 = xr.open_dataset(climfile)

    # Loop over variables
    for varname in varnames:

        tmpvar = varname + '_anom'
        print('processing',tmpvar)
        anom = ds[tmpvar]

        nTime = anom.shape[0]
        nClasses = anom.shape[1]
        nLat = anom.shape[2]
        nLon = anom.shape[3]

        print('reading mean')
        tmpvar = varname + '_mean'
        mean = ds2[tmpvar]
        print('reading std')
        tmpvar = varname + '_std'
        std = ds2[tmpvar]

        nDays = mean.shape[0]

        # Reconstruct original time series from anomaly, mean, and std
        print('combining anom and mean')
        data = np.full([nTime,nClasses,nLat,nLon],np.nan,dtype=np.single)
        for t in range(nTime):
            d = t % nDays
            data[t] = anom[t] * std[d] + mean[d]

        # Keep values within valid ranges
        data = np.where((data < valid_min[varname]), valid_min[varname], data)
        data = np.where((data > valid_max[varname]), valid_max[varname], data)

        tmpvar = varname + '_anom'
        ds[varname] = (['time','veg_class','lat','lon'], data)
        ds[varname].attrs['long_name'] = ds[tmpvar].attrs['long_name'][:-23]
        ds = ds.drop(tmpvar)
        tmpvar = varname + '_count'
        ds = ds.drop(tmpvar)
        tmpvar = varname + '_anom_gapfill_flag'
        ds = ds.drop(tmpvar)

    if recompute_fcan:
        fcanopy = np.empty(ds['fcanopy'].shape, dtype=np.single)
        NDVI = np.empty(ds['fcanopy'].shape, dtype=np.single)
        fcanopy[:] = ds['fcanopy'][:]
        vcompute_ndvi = np.vectorize(compute_ndvi)
        NDVI = vcompute_ndvi(fcanopy)
        vcompute_fcan = np.vectorize(compute_fcan)
        fcanopy = vcompute_fcan(NDVI)
        ds['fcanopy'][:] = fcanopy[:]

    # Write to output file
    ds.to_netcdf(outfile)

    ds.close()
    ds2.close()

if __name__ == "__main__":
    main()

