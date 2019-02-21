#!/usr/bin/env python

import os.path
import numpy as np
import xarray as xr
import sys, getopt

def subsample_variable(data_in, shape, dtype, llcorner_in,
                       cellsize_in, llcorner, cellsize, mask, fill_value):

    # Replace nans with fill_value
    if ((dtype == 'float' or dtype == 'np.single' or dtype == 'double' or
         dtype == 'int' or dtype == 'np.int32' or dtype == 'np.int64') and
        np.any(np.isnan(data_in))):
        tmp = np.where(np.isnan(data_in),fill_value,data_in)
        data_in = tmp

    # Copy values to output data array
    # Using an explicit loop to account for different grid resolutions
    shape_in = data_in.shape
    data = np.full(shape, fill_value, dtype=dtype)
    for y in range(shape[-2]):
        ctr_lat = llcorner[0] + (y + 0.5) * cellsize
        y_in = int((ctr_lat - llcorner_in[0]) / cellsize_in)
        if y_in < 0 or y_in >= shape_in[-2]:
            continue
        for x in range(shape[-1]):
            ctr_lon = llcorner[1] + (x + 0.5) * cellsize
            x_in = int((ctr_lon - llcorner_in[1]) / cellsize_in)
            if x_in < 0 or x_in >= shape_in[-1]:
                continue
            if len(shape) == 2:
              data[...,y,x] = np.asscalar(data_in[...,y_in,x_in].astype(dtype))
            else:
              data[...,y,x] = data_in[...,y_in,x_in].astype(dtype)

    # Overwrite values from data_in with fill_value where mask is 0
    data[..., mask != 1] = fill_value

    return data


def main():
    infile = ''
    domainfile = ''
    outfile = ''
    snow_bands = False
    new_snow_alb_supplied = False
    fcanopy = False
    irrigation = False
    nodata_default = False
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:d:nsafro:",["infile=","domainfile=","nodata_default","snow_bands","new_snow_albedo","fcanopy","irrigation","outfile="])
    except getopt.GetoptError:
        print(sys.argv[0], ' -i <infile> -d <domainfile> [-n] [-s] [-a] [-f] [-r] -o <outfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -i <infile> -d <domainfile> [-n] [-s] [-a] [-f] [-r] -o <outfile>')
            sys.exit()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-d", "--domainfile"):
            domainfile = arg
        elif opt in ("-n", "--nodata_default"):
            nodata_default = True
        elif opt in ("-s", "--snow_bands"):
            snow_bands = True
        elif opt in ("-a", "--new_snow_albedo"):
            new_snow_alb_supplied = True
        elif opt in ("-f", "--fcanopy"):
            fcanopy = True
        elif opt in ("-r", "--irrigation"):
            irrigation = True
        elif opt in ("-o", "--outfile"):
            outfile = arg

    # Define fill values
    fill_value_mask = np.int32(0)
    fill_value_int = np.int32(-1)
    fill_value_str = ''
    fill_value_float = 9.96920996838687e+36

    # Read domain file
    ds = xr.open_dataset(domainfile)

    # Read relevant dims and vars
    lat = ds['lat']
    lon = ds['lon']
#    lats = ds['lats']
#    lons = ds['lons']
    mask = np.where(ds['mask']==1,1,0).astype(np.int32)
#    lats = np.where(mask==1,lats,fill_value_float)
#    lons = np.where(mask==1,lons,fill_value_float)
    nLat = len(lat)
    nLon = len(lon)
    minlat = np.asscalar(lat[0])
    maxlat = np.asscalar(lat[-1])
    minlon = np.asscalar(lon[0])
    maxlon = np.asscalar(lon[-1])
    cellsize = ( maxlat - minlat ) / (nLat - 1)
    minlat -= (0.5 * cellsize)
    minlon -= (0.5 * cellsize)
    maxlat += (0.5 * cellsize)
    maxlon += (0.5 * cellsize)
    llcorner = [minlat, minlon]

    # Read infile
    ds_in = xr.open_dataset(infile)

    # Read coordinate variables
    lat_in = ds_in['lat']
    lon_in = ds_in['lon']
    nLat_in = len(lat_in)
    nLon_in = len(lon_in)
    minlat_in = np.asscalar(lat_in[0])
    maxlat_in = np.asscalar(lat_in[-1])
    minlon_in = np.asscalar(lon_in[0])
    maxlon_in = np.asscalar(lon_in[-1])
    cellsize_in = ( maxlat_in - minlat_in ) / (nLat_in - 1)
    minlat_in -= (0.5 * cellsize_in)
    minlon_in -= (0.5 * cellsize_in)
    maxlat_in += (0.5 * cellsize_in)
    maxlon_in += (0.5 * cellsize_in)
    llcorner_in = [minlat_in, minlon_in]

    month = ds_in['month']
    layer = ds_in['layer']
    veg_class = ds_in['veg_class']
    veg_descr = ds_in['veg_descr']
    root_zone = ds_in['root_zone']
    nMonth = len(month)
    nLayer = len(layer)
    nClass = len(veg_class)
    nRoot = len(root_zone)
    if snow_bands:
        snow_band = ds_in['snow_band']
        nBand = len(snow_band)
    else:
        nBand = 1
        snow_band = [1]

    # Variables
    varnames_1d_soil = [
                        'layer',
                       ]
    varnames_1d_snow = [
                        'snow_band',
                       ]
    varnames_1d_veg = [
                       'veg_descr',
                      ]

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
    if (new_snow_alb_supplied):
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

    if (irrigation):
        varnames_3d_veg.append('irr_active')
        varnames_3d_veg.append('ithresh')
        varnames_3d_veg.append('itarget')

    varnames_4d_root = [
                        'root_depth',
                        'root_fract',
                       ]

    varnames_4d_clim = ['LAI']
    if (fcanopy):
      varnames_4d_clim.append('fcanopy')
    varnames_4d_clim.append('albedo')
    varnames_4d_clim.append('veg_rough')
    varnames_4d_clim.append('displacement')
    if (irrigation):
        varnames_4d_clim.append('fcrop')
        varnames_4d_clim.append('firr')
        varnames_4d_clim.append('irr_clim')

    varnames = []
#    for varname in varnames_1d_soil:
#        varnames.append(varname)
#    if snow_bands:
#        for varname in varnames_1d_snow:
#            varnames.append(varname)
#    for varname in varnames_1d_veg:
#        varnames.append(varname)
    for varname in varnames_2d_soil:
        varnames.append(varname)
    for varname in varnames_2d_veg:
        varnames.append(varname)
    for varname in varnames_3d_soil:
        varnames.append(varname)
    if snow_bands:
        for varname in varnames_3d_snow:
            varnames.append(varname)
    for varname in varnames_3d_veg:
        varnames.append(varname)
    for varname in varnames_4d_root:
        varnames.append(varname)
    for varname in varnames_4d_clim:
        varnames.append(varname)

    varnames_1d = []
    for varname in varnames_1d_soil:
        varnames_1d.append(varname)
    if snow_bands:
        for varname in varnames_1d_snow:
            varnames_1d.append(varname)
    for varname in varnames_1d_veg:
        varnames_1d.append(varname)

    varnames_2d = []
    for varname in varnames_2d_soil:
        varnames_2d.append(varname)
    for varname in varnames_2d_veg:
        varnames_2d.append(varname)

    varnames_int = [
                    'layer',
                    'run_cell',
                    'gridcell',
                    'fs_active',
                    'Nveg',
                    'overstory',
                   ]
    if (irrigation):
        varnames_int.append('irr_active')
        varnames_int.append('irr_clim')

    varnames_str = ['veg_descr']
    if (irrigation):
        varnames_str.append('ithresh')
        varnames_str.append('itarget')

    # Write output file
    ds_out = xr.Dataset(
                        {
#                         'lats': (['lat','lon'],lats),
#                         'lons': (['lat','lon'],lons),
                         'mask': (['lat','lon'],mask),
                        },
                        coords={
                                'lat': (['lat'],lat),
                                'lon': (['lon'],lon),
                                'month': (['month'],month),
                                'layer': (['nlayer'],layer),
                                'snow_band': (['snow_band'],snow_band),
                                'veg_class': (['veg_class'],veg_class),
                                'veg_descr': (['veg_class'],veg_descr),
                                'root_zone': (['root_zone'],root_zone),
                               }
                       )

    for varname in ds.coords:
        ds_out[varname].attrs = ds[varname].attrs
        ds_out[varname].encoding = ds[varname].encoding
    ds_out['mask'].attrs = ds['mask'].attrs
    if nodata_default:
        ds_out['mask'].encoding['_FillValue'] = fill_value_mask
    else:
        ds_out['mask'].encoding['_FillValue'] = ds['mask'].encoding['_FillValue']

    # Do the subsampling
    out_dict = {}
    for varname in varnames:

        print('varname',varname)
        if varname in varnames_1d:
            print('1d var:',varname)
            if (varname == 'layer'):
                shape = [nLayer]
                dimstr = ['nlayer']
            elif (varname == 'snow_band'):
                shape = [nBand]
                dimstr = ['snow_band']
            elif (varname == 'veg_descr'):
                shape = [nClass]
                dimstr = ['veg_class']
        elif varname in varnames_2d:
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
        elif varname in varnames_3d_veg:
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
        else:
            print('no case for',varname)

        if varname in varnames_int:
            dtype = np.int32
            fill_value = fill_value_int
        elif varname in varnames_str:
            dtype = np.str
            fill_value = fill_value_str
        else:
            dtype = np.single
            fill_value = fill_value_float

        if (varname in varnames_1d):
            out_dict[varname] = ds_in[varname]
        else:
            out_dict[varname] = subsample_variable(ds_in[varname],
                                                   shape,
                                                   dtype,
                                                   llcorner_in,
                                                   cellsize_in,
                                                   llcorner,
                                                   cellsize,
                                                   mask, fill_value)

        ds_out[varname] = (dimstr, out_dict[varname])
        ds_out[varname].attrs = ds_in[varname].attrs
        if nodata_default:
            if varname in varnames_int:
                ds_out[varname].encoding['_FillValue'] = fill_value_int
            elif varname not in varnames_str:
                ds_out[varname].encoding['_FillValue'] = fill_value_float
        else:
            ds_out[varname].encoding = ds_in[varname].encoding

    # HACK
    ds_out['Cv'] = np.around(ds_out['Cv'], 6)

    print('writing to',outfile)
    ds_out.to_netcdf(outfile)

    ds_in.close()
    ds_out.close()
    ds.close()

if __name__ == "__main__":
    main()    
