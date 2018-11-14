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
        opts, args = getopt.getopt(sys.argv[1:],"hi:v:o:",["infile=","varnamelist=","outfile="])
    except getopt.GetoptError:
        print(sys.argv[0], ' -i <infile> -v <varnamelist> -o <outfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -i <infile> -v <varnamelist> -o <outfile>')
            sys.exit()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-v", "--varnamelist"):
            varnamelist = arg
            varnames = varnamelist.lstrip().rstrip().split(',')
        elif opt in ("-o", "--outfile"):
            outfile = arg


    # Open input file
    ds = xr.open_dataset(infile)
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

    dayvar_out = dayvar.groupby('day_of_year').mean('time').astype(int)

    # Output dataset
    ds_out = xr.Dataset(
        {
            'LandMask': (['lat','lon'], mask),
            'Cv': (['veg_class','lat','lon'], Cv),
        },
        coords={
            'lat': (['lat'], latvar),
            'lon': (['lon'], lonvar),
            'day_of_year': (['time'], dayvar_out),
            'veg_class': (['veg_class'], classvar),
            'class_name': (['veg_class'], classname),
        }
    )
  
    for tmpvar in ['lat','lon','day_of_year','veg_class','class_name','LandMask','Cv']:
        ds_out[tmpvar].attrs = ds[tmpvar].attrs


    # Loop over variables
    for varname in varnames:

        print('processing',varname)
        data = ds[varname]

        # Compute climatological summaries
        mean = data.groupby('day_of_year').mean('time')
        std = data.groupby('day_of_year').std('time')
        count = data.groupby('day_of_year').count('time').where(~np.isnan(mean))
        nDay = mean.shape[0]
        nClasses = mean.shape[1]
        nLat = mean.shape[2]
        nLon = mean.shape[3]
        for d in range(nDay):
            for c in range(nClasses):
                if (np.sum(count[d, c, :, :]) > 0):
                    tmp = Cv[c, :, :] * np.isnan(mean[d, c, :, :]).astype(float)
                    a = np.where(tmp > 0)[0]
                    b = np.where(tmp > 0)[1]
                    for i,j in zip(a,b):
                        count[d, c, i, j] = 0

        tmpvar = varname + '_mean'
        ds_out[tmpvar] = (['day_of_year','veg_class','lat','lon'], mean)
        ds_out[tmpvar].attrs = ds[varname].attrs
        ds_out[tmpvar].attrs['long_name'] = ds[varname].attrs['long_name'] + ' Climatological Mean'
        ds_out[tmpvar].encoding['zlib'] = True
        tmpvar = varname + '_std'
        ds_out[tmpvar] = (['day_of_year','veg_class','lat','lon'], std)
        ds_out[tmpvar].attrs = ds[varname].attrs
        ds_out[tmpvar].attrs['long_name'] = ds[varname].attrs['long_name'] + ' Climatological Standard Deviation'
        ds_out[tmpvar].encoding['zlib'] = True
        tmpvar = varname + '_count'
        ds_out[tmpvar] = (['day_of_year','veg_class','lat','lon'], count)
        ds_out[tmpvar].attrs['long_name'] = ' Number of observations in climatological stats'
        ds_out[tmpvar].encoding['zlib'] = True


    # Write to output file
    ds_out.to_netcdf(outfile, engine='scipy')

    ds_out.close()
    ds.close()

if __name__ == "__main__":
    main()

