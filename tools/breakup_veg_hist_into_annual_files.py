#!/usr/bin/env python

import os.path
import numpy as np
import xarray as xr
import pandas as pd
import sys, getopt


def main():
    infile = ''
    outdir = ''
    prefix = ''
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:v:o:p:",["infile=","varnamelist=","outdir=","prefix="])
    except getopt.GetoptError:
        print(sys.argv[0], ' -i <infile> -v <varnamelist> -o <outdir> -p <prefix>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -i <infile> -v <varnamelist> -o <outdir> -p <prefix>')
            sys.exit()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-v", "--varnamelist"):
            varnamelist = arg
            varnames = varnamelist.lstrip().rstrip().split(',')
        elif opt in ("-o", "--outdir"):
            outdir = arg
        elif opt in ("-p", "--prefix"):
            prefix = arg

    ds = xr.open_dataset(infile)

    yearvar = ds['year']
    timevar = ds['time']
    latvar = ds['lat']
    lonvar = ds['lon']
    classvar = ds['veg_class']
    classname = ds['class_name']
    mask = ds['LandMask']
    Cv = ds['Cv']

    nTime = len(timevar)
    nClass = len(classvar)
    nLat = len(latvar)
    nLon = len(lonvar)

    year0 = np.asscalar(yearvar[0])
    year1 = np.asscalar(yearvar[-1])

    nYears = year1 - year0 + 1

    startdate = str(year0) + '-01-01'
    time_out = pd.date_range(start=startdate,periods=nTime,freq='M')

    # Loop over years
    for y in range(nYears):
        year = year0 + y
        outfile = outdir + '/' + prefix + '.' + str(year) + '.nc'
        nTime = 12
        t0 = y * nTime
        t1 = t0 + nTime

        ds_out = xr.Dataset(
            {
                'LandMask': (['lat','lon'], mask),
                'Cv': (['veg_class','lat','lon'], Cv),
            },
            coords={
                'time': (['time'], time_out[t0:t1]),
                'lat': (['lat'], latvar),
                'lon': (['lon'], lonvar),
                'veg_class': (['veg_class'], classvar),
                'class_name': (['veg_class'], classname),
            },
        )

        # Loop over variables
        for varname in varnames:
            tmp = np.empty([nTime,nClass,nLat,nLon])
            tmp[:] = ds[varname][t0:t1]
            ds_out[varname] = (['time','class','lat','lon'], tmp)

        # Assign attributes
        for varname in ds_out.variables.keys():
            ds_out[varname].attrs = ds[varname].attrs
            ds_out[varname].encoding = ds[varname].encoding

        # Write to output file
        ds_out.to_netcdf(outfile)
        ds_out.close()

    ds.close()

if __name__ == "__main__":
    main()

