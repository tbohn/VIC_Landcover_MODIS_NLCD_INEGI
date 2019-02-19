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
        opts, args = getopt.getopt(sys.argv[1:],"hi:v:s:e:o:",["infile=","varnamelist=","startyearidx=","endyearidx=","outfile="])
    except getopt.GetoptError:
        print(sys.argv[0], ' -i <infile> -v <varnamelist> -s <startyear> -e <endyear> -o <outfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -i <infile> -v <varnamelist> -s <startyear> -e <endyear> -o <outfile>')
            sys.exit()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-v", "--varnamelist"):
            varnamelist = arg
            varnames = varnamelist.lstrip().rstrip().split(',')
        elif opt in ("-s", "--startyearidx"):
            startyearidx = int(arg)
        elif opt in ("-e", "--endyearidx"):
            endyearidx = int(arg)
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
    mask = ds['LandMask']
    Cv = ds['Cv']
    nTime = len(timevar)
    nClass = len(classvar)
    nLat = len(latvar)
    nLon = len(lonvar)

    # Output dataset
    nMonth = 12
    month = np.arange(nMonth).astype(np.int32)
    nYears = int(nTime/nMonth)
    ds_out = xr.Dataset(
        {
            'LandMask': (['lat','lon'], mask),
            'Cv': (['veg_class','lat','lon'], Cv),
        },
        coords={
            'lat': (['lat'], latvar),
            'lon': (['lon'], lonvar),
            'month_of_year': (['month_of_year'], month),
            'veg_class': (['veg_class'], classvar),
            'class_name': (['veg_class'], classname),
        }
    )
  
    for tmpvar in ['lat','lon','veg_class','class_name','LandMask','Cv']:
        ds_out[tmpvar].attrs = ds[tmpvar].attrs
    ds_out['month_of_year'].attrs['units'] = 'month'
    ds_out['month_of_year'].attrs['long_name'] = 'month of climatological year'


    # Loop over variables
    for varname in varnames:

        print('processing',varname)

        # Compute climatological summaries
        data = np.empty([nTime,nClass,nLat,nLon],dtype=np.single)
        data[:] = ds[varname][:]
        data = data.reshape(nYears,nMonth,nClass,nLat,nLon)
        if endyearidx > startyearidx:
            mean = np.mean(data[startyearidx:endyearidx], 0)
        else:
            mean = data[startyearidx]

        tmpvar = varname + '_mean'
        ds_out[tmpvar] = (['month_of_year','veg_class','lat','lon'], mean)
        ds_out[tmpvar].attrs = ds[varname].attrs
#        ds_out[tmpvar].encoding = ds[tmpvar].encoding
        ds_out[tmpvar].encoding['zlib'] = True


    # Write to output file
    ds_out.to_netcdf(outfile)

if __name__ == "__main__":
    main()

