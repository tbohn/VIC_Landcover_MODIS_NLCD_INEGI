#!/usr/bin/env python

import os.path
import numpy as np
import xarray as xr
import sys, getopt
import math
import time

def main():
    infile = ''
    outfile = ''
    has_nsa = False
    has_irr = False
    has_imp = False
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:nrmo:",["infile=","outfile="])
    except getopt.GetoptError:
        print(sys.argv[0], ' -i <infile> [-n] [-r] [-m] -o <outfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -i <infile> [-n] [-r] [-m] -o <outfile>')
            sys.exit()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-n", "--has_nsa"):
            has_nsa = True
        elif opt in ("-r", "--has_irr"):
            has_irr = True
        elif opt in ("-m", "--has_imp"):
            has_imp = True
        elif opt in ("-o", "--outfile"):
            outfile = arg

    nodata_default = True

    # Define nodata values
    fill_value_mask = np.int32(0)
    fill_value_int = np.int32(-9999)
    fill_value_str = ''
    fill_value_float = np.nan

    # Variables
    varnames_2d_soil = [
                        'run_cell',
                        'gridcell',
                        'infilt',
                        'Ds',
                        'Dsmax',
                        'Ws',
                        'c',
                        'elev',
                        'avg_T',
                        'dp',
                        'off_gmt',
                        'rough',
                        'snow_rough',
                        'annual_prec',
                        'fs_active',
                       ]
    if (has_nsa):
        varnames_2d_soil.append('new_snow_albedo')

    varnames_2d_veg = [
                       'Nveg',
                      ]

    varnames_3d_soil = [
                        'expt',
                        'Ksat',
                        'phi_s',
                        'init_moist',
                        'depth',
                        'bubble',
                        'quartz',
                        'bulk_density',
                        'soil_density',
                        'Wcr_FRACT',
                        'Wpwp_FRACT',
                        'resid_moist',
                       ]

    varnames_3d_snow = [
                        'AreaFract',
                        'Pfactor',
                        'elevation',
                       ]

    varnames_3d_veg = [
                       'Cv',
                       'overstory',
                       'rarc',
                       'rmin',
                       'wind_h',
                       'RGL',
                       'rad_atten',
                       'wind_atten',
                       'trunk_ratio',
                      ]

    if (has_irr):
        varnames_3d_veg.append('irr_active')
        varnames_3d_veg.append('ithresh')
        varnames_3d_veg.append('itarget')
    if (has_imp):
        varnames_3d_veg.append('fimp')
        varnames_3d_veg.append('feffimp')

    varnames_4d_root = [
                        'root_depth',
                        'root_fract',
                       ]

    varnames_4d_clim = [
                        'LAI',
                        'fcanopy',
                        'albedo',
                        'veg_rough',
                        'displacement',
                       ]
    if (has_irr):
        varnames_4d_clim.append('fcrop')
        varnames_4d_clim.append('firr')
        varnames_4d_clim.append('irr_clim')

    varnames_int = [
                    'layer',
                    'run_cell',
                    'gridcell',
                    'fs_active',
                    'Nveg',
                    'overstory',
                   ]

    varnames_str = []

    ds = xr.open_dataset(infile)

    nClass = ds['Cv'].shape[0]
    nLat = ds['Cv'].shape[1]
    nLon = ds['Cv'].shape[2]
    nMonth = ds['LAI'].shape[1]
    nRoot = ds['root_depth'].shape[1]
    nLayer = ds['layer'].shape[0]

    Cv = np.empty(ds['Cv'].shape)
    Cv[:] = ds['Cv'][:]

    # Write output file
    ds_out = xr.Dataset(
                        {
                         'mask': (['lat','lon'],ds['mask']),
                        },
                        coords={
                                'lat': (['lat'],ds['lat']),
                                'lon': (['lon'],ds['lon']),
                                'month': (['month'],ds['month']),
                                'layer': (['nlayer'],ds['layer']),
                                'root_zone': (['root_zone'],ds['root_zone']),
                                'veg_class': (['veg_class'],ds['veg_class']),
                                'veg_descr': (['veg_class'],ds['veg_descr']),
                               }
                       )

    if 'snow_band' in ds.coords.keys():
        nBand = len(ds['snow_band'])
        snow_band = np.arange(nBand).astype(int)
        ds_out['snow_band'] = (['snow_band'],snow_band)
        ds_out.set_coords('snow_band', inplace=True)
    for varname in ds_out.coords.keys():
        ds_out[varname].attrs = ds[varname].attrs
        ds_out[varname].encoding = ds[varname].encoding
    ds_out['mask'].attrs = ds['mask'].attrs
    if nodata_default:
        ds_out['mask'].encoding['_FillValue'] = fill_value_mask
    else:
        ds_out['mask'].encoding['_FillValue'] = ds['mask'].encoding['_FillValue']

    # Copy data and gapfill where appropriate
    out_dict = {}

    # Special cases
    # Copy Cv
    out_dict['Cv'] = np.empty([nClass,nLat,nLon])
    out_dict['Cv'][:] = ds['Cv'][:]
    Cv_sum = np.sum(out_dict['Cv'],axis=0)

    # Compute Nveg
    Cv_mask = np.where(out_dict['Cv'] > 0, 1, 0).astype(np.int32)
    Nveg = np.sum(Cv_mask, axis=0, dtype=np.int32).astype(np.int32)
    Nveg[Nveg == 0] = fill_value_int
    out_dict['Nveg'] = Nveg

    # Loop over all other variables
    for varname in ds.variables.keys():

        if varname in ds.coords.keys():
            continue
        print('varname',varname)
        shape_in = ds[varname].shape
        nDims = len(shape_in)
        if (varname in varnames_2d_soil or
            varname in varnames_2d_veg or
            varname in varnames_3d_soil or
            varname in varnames_3d_snow):

            if varname in varnames_2d_soil:
                print('2d var:',varname)
                shape = [nLat,nLon]
                dimstr = ['lat','lon']
            elif varname in varnames_2d_veg:
                print('2d var:',varname)
                shape = [nLat,nLon]
                dimstr = ['lat','lon']
            elif varname in varnames_3d_soil:
                print('3d_soil var:',varname)
                shape = [nLayer,nLat,nLon]
                dimstr = ['nlayer','lat','lon']
            elif varname in varnames_3d_snow:
                print('3d_snow var:',varname)
                shape = [nBand,nLat,nLon]
                dimstr = ['snow_band','lat','lon']

            out_dict[varname] = ds[varname]

        elif (varname in varnames_3d_veg or
              varname in varnames_4d_root or
              varname in varnames_4d_clim):

            if varname in varnames_3d_veg:
                print('3d_veg var:',varname)
                shape = [nClass,nLat,nLon]
                dimstr = ['veg_class','lat','lon']
            elif varname in varnames_4d_root:
                print('4d_root var:',varname)
                shape = [nClass,nRoot,nLat,nLon]
                dimstr = ['veg_class','root_zone','lat','lon']
            elif varname in varnames_4d_clim:
                print('4d_clim var:',varname)
                shape = [nClass,nMonth,nLat,nLon]
                dimstr = ['veg_class','month','lat','lon']

            if varname in varnames_int:
                dtype = np.int32
                fill_value = fill_value_int
            elif varname in varnames_str:
                dtype = np.str
                fill_value = fill_value_str
            else:
                dtype = np.single
                fill_value = fill_value_float

            # Determine where gaps exist
            # gapmask = 1 where data should exist but doesn't; = 0 otherwise
            should_exist = Cv_mask
            exist_before = np.where(~np.isnan(ds[varname]), 1, 0)
            if varname in varnames_3d_veg:
                gapmask = should_exist - exist_before
            elif varname in varnames_4d_root:
                gapmask = np.empty([nClass,nRoot,nLat,nLon])
                for c in range(nClass):
                    for r in range(nRoot):
                        gapmask[c,r] = should_exist[c] - exist_before[c,r]
            elif varname in varnames_4d_clim:
                gapmask = np.empty([nClass,nMonth,nLat,nLon])
                for c in range(nClass):
                    for t in range(nMonth):
                        gapmask[c,t] = should_exist[c] - exist_before[c,t]

            # Do gapfilling
            if (varname in varnames_3d_veg and
                varname not in ['Cv','Nveg']):
                data = np.empty([nClass,nLat,nLon])
                data[:] = ds[varname].copy()
                tmpmean = np.nanmean(data, axis=(1,2))
                for c in range(nClass):
                    data[c,gapmask[c] == 1] = tmpmean[c]
                out_dict[varname] = np.empty([nClass,nLat,nLon])
                out_dict[varname][:] = data[:]
            elif varname in varnames_4d_root:
                data = np.empty([nClass,nRoot,nLat,nLon])
                data[:] = ds[varname].copy()
                tmpmean = np.nanmean(data, axis=(2,3))
                for c in range(nClass):
                    for r in range(nRoot):
                        data[c,r,gapmask[c,r] == 1] = tmpmean[c,r]
                out_dict[varname] = np.empty([nClass,nRoot,nLat,nLon])
                out_dict[varname][:] = data[:]
            elif varname in varnames_4d_clim:
                data = np.empty([nClass,nMonth,nLat,nLon])
                data[:] = ds[varname].copy()
                tmpmean = np.nanmean(data, axis=(2,3))
                for c in range(nClass):
                    for t in range(nMonth):
                        data[c,t,gapmask[c,t] == 1] = tmpmean[c,t]
                out_dict[varname] = np.empty([nClass,nMonth,nLat,nLon])
                out_dict[varname][:] = data[:]

            # null out values where Cv==0
            if (varname in varnames_3d_veg and
                varname not in ['Cv', 'Nveg']):
                tmp = np.empty(out_dict[varname].shape)
                tmp[:] = out_dict[varname][:]
                tmp[Cv == 0] = np.nan
                out_dict[varname][:] = tmp[:]
            elif varname in varnames_4d_root:
                for c in range(nClass):
                    for r in range(nRoot):
                        out_dict[varname][c,r,Cv[c] == 0] = np.nan
            elif varname in varnames_4d_clim:
                for c in range(nClass):
                    for t in range(nMonth):
                        out_dict[varname][c,t,Cv[c] == 0] = np.nan

            if varname in varnames_int:
                tmp = np.empty(out_dict[varname].shape)
                tmp[:] = out_dict[varname][:]
                tmp[np.isnan(tmp)] = fill_value_int
                out_dict[varname] = tmp.astype(np.int32).copy()

        else:
            print('no case for',varname)
            continue

        # Populate ds_out
        ds_out[varname] = (dimstr, out_dict[varname])
        ds_out[varname].attrs = ds[varname].attrs
        if nodata_default:
            if varname in varnames_int:
                ds_out[varname].encoding['_FillValue'] = fill_value_int
            elif varname not in varnames_str:
                ds_out[varname].encoding['_FillValue'] = fill_value_float
        else:
            ds_out[varname].encoding = ds[varname].encoding

    # Write to output file
    ds_out.to_netcdf(outfile)

    ds.close()
    ds_out.close()

if __name__ == "__main__":
    main()

