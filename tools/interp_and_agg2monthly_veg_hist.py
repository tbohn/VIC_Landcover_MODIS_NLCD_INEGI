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
        opts, args = getopt.getopt(sys.argv[1:],"hi:v:o:p:m:c:",["infile=","varnamelist=","outdir=","prefix=","monthlyfile=","climfile="])
    except getopt.GetoptError:
        print(sys.argv[0], ' -i <infile> -v <varnamelist> -o <outdir> -p <prefix> -m <monthlyfile> -c <climfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -i <infile> -v <varnamelist> -o <outdir> -p <prefix> -m <monthlyfile> -c <climfile>')
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
        elif opt in ("-m", "--monthlyfile"):
            monthlyfile = arg
        elif opt in ("-c", "--climfile"):
            climfile = arg


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

    data = {}
    first = 1
    for varname in varnames:
        print('reading',varname)
        data[varname] = ds[varname]
        if first:
            nTime = data[varname].shape[0]
            nClasses = data[varname].shape[1]
            nLat = data[varname].shape[2]
            nLon = data[varname].shape[3]
        first = 0
    year_dict = {}
    dayofyear_dict = {}
    for t in range(nTime):
        year_dict[np.asscalar(yearvar[t])] = 1
        dayofyear_dict[np.asscalar(dayvar[t])] = 1
    nYear = len(year_dict)
    nDayofyear = len(dayofyear_dict)
    dDay = int(np.ceil(365/nDayofyear))
    nDayspyear = {}
    for year in year_dict.keys():
        if year % 4 == 0:
            nDayspyear[year] = 366
        else:
            nDayspyear[year] = 365

    # Add record corresponding to final day of final year
    # Create temporary dataset for this; not for output
    endyear = int(yearvar[-1])
    enddatestr = str(endyear) + '-12-31'
    enddate = np.datetime64(enddatestr, 'ns')
    nTime_tmp = nTime + 1
    timevar_tmp = np.empty([nTime+1],'M8[ns]')
    timevar_tmp[:nTime] = timevar
    timevar_tmp[-1] = enddate
    yearvar_tmp = timevar_tmp.astype('M8[Y]')
    tmp = pd.Series(pd.to_datetime(timevar_tmp)).dt.dayofyear
    dayvar_tmp = tmp.values
    ds_tmp = xr.Dataset(
        coords={
            'time': (['time'], timevar_tmp),
            'year': (['time'], yearvar_tmp),
            'day_of_year': (['time'], dayvar_tmp),
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

    data_tmp = {}
    for varname in varnames:
        print('adding final record to',varname)
        data_tmp[varname] = np.empty([nTime_tmp,nClasses,nLat,nLon],dtype=np.single)
        data_tmp[varname][:nTime] = data[varname]
        data_tmp[varname][-1] = data_tmp[varname][-2]
        ds_tmp[varname] = (['time','veg_class','lat','lon'], data_tmp[varname])

    ds.close()

    #
    # Resample to daily with linear interpolation
    #

    # Daily time variables
    year = pd.to_datetime(timevar_tmp).year
    startyear = year[0]
    endyear = year[-1]
    startstr = str(startyear) + '-01-01'
    endstr = str(endyear) + '-12-31'
    timevar_daily = pd.date_range(start=startstr,end=endstr,freq='D')
    yearvar_daily = timevar_daily.year
    dayvar_daily = timevar_daily.dayofyear
    leapday = np.where(year % 4 == 0, 1, 0)

    # Monthly time variables
    timevar_monthly = pd.date_range(start=startstr,end=endstr,freq='M')
    yearvar_monthly = timevar_monthly.year
    monthvar_monthly = timevar_monthly.month
    nTime_monthly = len(timevar_monthly)

    # Allocate space for data vars
    data_daily = {}
    data_monthly = {}
    for varname in varnames:
        print('allocating daily and monthly',varname)
        data_daily[varname] = np.empty([nDayofyear,dDay,nClasses,nLat,nLon],dtype=np.single)
        data_monthly[varname] = np.empty([nTime_monthly,nClasses,nLat,nLon],dtype=np.single)

    # Loop over years and write daily outputs to separate yearly files
    offset = 0
    offset_monthly = 0
    for year in range(startyear,endyear+1):
        yearstr = str(year)
        print('computing daily data for year', year)
        idx_out = np.where(yearvar_daily.astype(int) == year)
        timevar_out = timevar_daily[idx_out]
        yearvar_out = yearvar_daily[idx_out]
        dayvar_out = dayvar_daily[idx_out]
        nTime_out = len(timevar_out)
        ds_daily = xr.Dataset(
            coords={
                'time': (['time'], timevar_out),
                'year': (['time'], yearvar_out),
                'day_of_year': (['time'], dayvar_out),
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
            ds_daily[tmpvar].attrs = ds[tmpvar].attrs

        # Resample to daily and interpolate
        for varname in varnames:
            print('interpolating', varname)
            X0 = data_tmp[varname][offset:offset+nDayofyear-1]
            X1 = data_tmp[varname][offset+1:offset+nDayofyear]
            for j in range(dDay):
                data_daily[varname][:-1,j] = X0 + (X1 - X0) / dDay * j
            # Account for partial interval at end of year
            dRemaining = nTime_out - (nDayofyear-1)*dDay
            X0 = data_tmp[varname][offset+nDayofyear-1]
            if np.asscalar(yearvar_out[0]) != endyear:
                X1 = data_tmp[varname][offset+nDayofyear]
            else:
                X1 = data_tmp[varname][-1]
            for j in range(dRemaining):
                data_daily[varname][-1,j] = X0 + (X1 - X0) / dRemaining * j

            ds_daily[varname] = (['time','veg_class','lat','lon'], data_daily[varname].reshape([nDayofyear*dDay,nClasses,nLat,nLon])[:nTime_out])
            ds_daily[varname].attrs = ds[varname].attrs
            ds_daily[varname].encoding['zlib'] = True

            print('aggregating to monthly',varname)
#            data_monthly[varname][offset_monthly:offset_monthly+12] = data_daily[varname][:nTime_out+1].reshape([nTime_out,nClasses,nLat,nLon]).resample(freq='M',dim='time')
            data_monthly[varname][offset_monthly:offset_monthly+12] = ds_daily[varname].resample(freq='M',dim='time')

        outfile = outdir + '/' + prefix + '.' + yearstr + '.nc'
        print('writing to', outfile)
        ds_daily.to_netcdf(outfile)
        ds_daily.close()

        offset += nDayofyear
        offset_monthly += 12


    # set up monthly output file
    print('setting up monthly file',monthlyfile)
    ds_monthly = xr.Dataset(
        coords={
            'time': (['time'], timevar_monthly),
            'year': (['time'], yearvar_monthly),
            'month': (['time'], monthvar_monthly),
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

    for tmpvar in ['time','year','lat','lon','lats','lons','veg_class','class_name','LandMask','Cv']:
        ds_monthly[tmpvar].attrs = ds[tmpvar].attrs
    ds_monthly['month'].attrs['long_name'] = 'month of calendar year'

    for varname in varnames:
        ds_monthly[varname] = (['time','veg_class','lat','lon'], data_monthly[varname])
        ds_monthly[varname].attrs = ds[tmpvar].attrs
        ds_monthly[varname].encoding['zlib'] = True

    ds_monthly.to_netcdf(monthlyfile)
    ds_monthly.close()

    # Compute monthly climatology
    print('computing monthly climatology')
    nTime_clim = 12
    nYears = np.ceil(nTime_monthly / nTime_clim).astype(int)
    timevar_clim = np.arange(nTime_clim).astype(int)

    mean = {}
    std = {}
    for varname in varnames:
        mean[varname] = ds_monthly[varname].groupby('month').mean('time')
        std[varname] = ds_monthly[varname].groupby('month').std('time')

    # Compute monthly normalized anomalies (z-transform)
    print('computing monthly anomalies')
    anom = {}
    for varname in varnames:
        anom[varname] = np.empty([nYears,nTime_clim,nClasses,nLat,nLon],dtype=np.single)
        data_tmp = np.empty([nTime_monthly,nClasses,nLat,nLon],dtype=np.single)
        data_tmp[:] = data_monthly[varname][:]
        data_tmp = data_tmp.reshape([nYears,nTime_clim,nClasses,nLat,nLon])
        for y in range(nYears):
            anom[varname][y] = (data_tmp[y] - mean[varname]) / std[varname]
        anom[varname] = np.where(np.isinf(anom[varname]), 0, anom[varname])
        anom[varname] = anom[varname].reshape([nTime_monthly,nClasses,nLat,nLon])

    # construct output dataset
    ds_clim = xr.Dataset(
        coords={
            'month_of_year': (['month_of_year'], timevar_clim),
            'time': (['time'], timevar_monthly),
            'year': (['time'], yearvar_monthly),
            'month': (['time'], monthvar_monthly),
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

    for tmpvar in ['time','year','lat','lon','lats','lons','veg_class','class_name','LandMask','Cv']:
        ds_clim[tmpvar].attrs = ds[tmpvar].attrs
    ds_clim['month_of_year'].attrs['units'] = 'month'
    ds_clim['month_of_year'].attrs['long_name'] = 'month of climatological year'
    ds_clim['month'].attrs['long_name'] = 'month of calendar year'

    for varname in varnames:
        tmpvar = varname + '_mean'
        ds_clim[tmpvar] = (['month_of_year','veg_class','lat','lon'], mean[varname])
        ds_clim[tmpvar].attrs['long_name'] = ds[varname].attrs['long_name'] + ' Climatological Mean'
        ds_clim[tmpvar].encoding['zlib'] = True
        tmpvar = varname + '_std'
        ds_clim[tmpvar] = (['month_of_year','veg_class','lat','lon'], std[varname])
        ds_clim[tmpvar].attrs['long_name'] = ds[varname].attrs['long_name'] + ' Climatological Standard Deviation'
        ds_clim[tmpvar].encoding['zlib'] = True
        tmpvar = varname + '_anom'
        ds_clim[tmpvar] = (['time','veg_class','lat','lon'], anom[varname])
        ds_clim[tmpvar].attrs['long_name'] = ds[varname].attrs['long_name'] + ' Climatological Anomaly'
        ds_clim[tmpvar].encoding['zlib'] = True

    ds_clim.to_netcdf(climfile)
    ds_clim.close()

    ds_tmp.close()

if __name__ == "__main__":
    main()

