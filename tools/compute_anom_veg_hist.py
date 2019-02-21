#!/usr/bin/env python

import os.path
import numpy as np
import pandas as pd
import xarray as xr
import sys, getopt
from datetime import datetime

def main():
    infile = ''
    outfile = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:c:v:o:",["infile=","climfile=","varnamelist=","outfile="])
    except getopt.GetoptError:
        print(sys.argv[0], ' -i <infile> -c <climfile> -v <varnamelist> -o <outfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -i <infile> -c <climfile> -v <varnamelist> -o <outfile>')
            sys.exit()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-c", "--climfile"):
            climfile = arg
        elif opt in ("-v", "--varnamelist"):
            varnamelist = arg
            varnames = varnamelist.lstrip().rstrip().split(',')
        elif opt in ("-o", "--outfile"):
            outfile = arg


    # Open input file
    ds = xr.open_dataset(infile)

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

    # Open climfile
    ds2 = xr.open_dataset(climfile)

    # Output dataset
    ds_out = xr.Dataset(
        coords={
            'time': (['time'], timevar),
            'year': (['time'], yearvar),
            'day_of_year': (['time'], dayvar),
            'lat': (['lat'], latvar),
            'lon': (['lon'], lonvar),
            'lats': (['lat','lon'], latsvar),
            'lons': (['lat','lon'], lonsvar),
            'veg_class': (['veg_class'], classvar),
            'class_name': (['veg_class'], classname),
        },
        data_vars={
            'LandMask': (['lat','lon'], mask),
            'Cv': (['veg_class','lat','lon'], Cv),
        },
    )
  
    for tmpvar in ['time','year','day_of_year','lat','lon','lats','lons','veg_class','class_name','LandMask','Cv']:
        ds_out[tmpvar].attrs = ds[tmpvar].attrs

    # Loop over variables
    for varname in varnames:

        print('processing',varname)
        data = ds[varname]

        nTime = data.shape[0]
        nClasses = data.shape[1]
        nLat = data.shape[2]
        nLon = data.shape[3]

        # Read climatological mean and std from climfile
        print('reading mean')
        tmpvar = varname + '_mean'
        mean = ds2[tmpvar]
        print('reading std')
        tmpvar = varname + '_std'
        std = ds2[tmpvar]
        nDays = mean.shape[0]
        nYears = int(nTime / nDays)

        # Compute normalized anomalies (z-transform)
        print('computing anom')
        anom = np.empty([nYears,nDays,nClasses,nLat,nLon],dtype=np.single)
        data_tmp = np.empty([nTime,nClasses,nLat,nLon],dtype=np.single)
        data_tmp[:] = data[:]
        data_tmp = data_tmp.reshape([nYears,nDays,nClasses,nLat,nLon])
        for y in range(nYears):
            anom[y] = (data_tmp[y] - mean) / std
        anom = np.where(np.isinf(anom), 0, anom)
        anom = anom.reshape([nTime,nClasses,nLat,nLon])

        # Define count
        count = np.zeros([nTime,nClasses,nLat,nLon],np.int)
        count = np.where(~np.isnan(anom), 1, 0)

        tmpvar = varname + '_anom'
        ds_out[tmpvar] = (['time','veg_class','lat','lon'], anom)
        ds_out[tmpvar].attrs = ds[varname].attrs
        ds_out[tmpvar].attrs['long_name'] = ds[varname].attrs['long_name'] + ' Climatological Anomaly'
        ds_out[tmpvar].encoding['zlib'] = True
        tmpvar = varname + '_count'
        ds_out[tmpvar] = (['time','veg_class','lat','lon'], count)
        ds_out[tmpvar].attrs['long_name'] = ' Number of observations in climatological stats'
        ds_out[tmpvar].encoding['zlib'] = True

    # Write to output file
    ds_out.to_netcdf(outfile)

    ds.close()
    ds2.close()
    ds_out.close()

if __name__ == "__main__":
    main()

