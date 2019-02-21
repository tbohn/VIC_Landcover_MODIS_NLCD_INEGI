#!/usr/bin/env python

import os.path
import numpy as np
import xarray as xr
import pandas as pd
import sys, getopt

def main():
    infile = ''
    outfile = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:s:e:o:",["infile=","startyear=","endyear=","outfile="])
    except getopt.GetoptError:
        print(sys.argv[0], ' -i <infile> -s <startyear> -e <endyear> -o <outfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -i <infile> -s <startyear> -e <endyear> -o <outfile>')
            sys.exit()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-s", "--startyear"):
            startyear = int(arg)
        elif opt in ("-e", "--endyear"):
            endyear = int(arg)
        elif opt in ("-o", "--outfile"):
            outfile = arg


    ds = xr.open_dataset(infile)

    nTime_in = len(ds['time'])
    nClass = len(ds['veg_class'])
    nLat = len(ds['lat'])
    nLon = len(ds['lon'])
    lat = np.empty([nLat])
    lon = np.empty([nLon])
    lat[:] = ds['lat'][:]
    lon[:] = ds['lon'][:]
    startdate = '{:04d}-01-01 00:00:00'.format(startyear)
    enddate = '{:04d}-12-31 23:59:00'.format(endyear)
    time_in = pd.date_range(start=startdate,end=enddate,freq='M')
    month_in = time_in.month
    time_out = pd.date_range(start=startdate,end=enddate,freq='H')
    month_out = time_out.month
    nTime_out = len(time_out)

    varnames = ['LAI','fcanopy','albedo']
    outvars = {}
    for varname in varnames: 
        outvars[varname] = np.empty([nTime_out,nClass,nLat,nLon])
        for m in range(nTime_in):
            tlist = np.where(month_out == month_in[m])
            outvars[varname][tlist] = ds[varname][m]

    ds.close()

    time = pd.date_range(start=startdate,end=enddate,freq='M')
    ds_out = xr.Dataset(
        coords = {
            'time': (['time'], time_out),
            'veg_class': (['veg_class'], ds['veg_class']),
            'lat': (['lat'], lat),
            'lon': (['lon'], lon),
        }
    )
    for varname in varnames:
        ds_out[varname] = (['time','veg_class','lat','lon'], outvars[varname])

    ds_out.to_netcdf(outfile)
    ds_out.close()

if __name__ == "__main__":
    main()

