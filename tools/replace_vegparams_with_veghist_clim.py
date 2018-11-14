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

    # Overwrite values from data_in with fill_value where mask != 1
    data[..., mask != 1] = fill_value

    return data


def compute_depth_area_rel(area_max, nNodes, type):

    depth = np.empty(nNodes, dtype=np.single)
    area = np.empty(nNodes, dtype=np.single)

    depth_min = 0
    if type == 'lake':
        depth_max = 1000
        area_min = area_max * 0.99
    else:
        depth_max = 10
        area_min = area_max / nNodes

    for i in range(nNodes):
        depth[i] = depth_max * (1 - (i / nNodes))
        area[i] = area_max + (area_min - area_max) * i / (nNodes - 1)

    return [depth,area]


def main():
    domainfile = ''
    paramfile = ''
    libfile = ''
    rootfile = ''
    climfile = ''
    outfile = ''
    lakes = False
    open_water_class = -1
    fill_class = 0
    new_snow_alb_supplied = False
    snow_bands = False
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hd:p:b:r:c:l:w:f:o:as",["domainfile=","paramfile=","libfile=","rootfile=","climfile=","lake_class=","open_water_class=","fill_class=","outfile="])
    except getopt.GetoptError:
        print(sys.argv[0], ' -d <domainfile> -p <paramfile> -b <libfile> -r <rootfile> -c <climfile> [-l <lake_class>] [-w <open_water_class>] [-f <fill_class>] -o <outfile> [-a] [-s]')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0], ' -d <domainfile> -p <paramfile> -b <libfile> -r <rootfile> -c <climfile> [-l <lake_class>] [-w <open_water_class>] [-f <fill_class>] -o <outfile> [-a] [-s]')
            sys.exit()
        elif opt in ("-d", "--domainfile"):
            domainfile = arg
        elif opt in ("-p", "--paramfile"):
            paramfile = arg
        elif opt in ("-b", "--libfile"):
            libfile = arg
        elif opt in ("-r", "--rootfile"):
            rootfile = arg
        elif opt in ("-c", "--climfile"):
            climfile = arg
        elif opt in ("-l", "--lake_class"):
            lakes = True
            lake_class = int(arg)
        elif opt in ("-w", "--open_water_class"):
            open_water_class = int(arg)
        elif opt in ("-f", "--fill_class"):
            fill_class = int(arg)
        elif opt in ("-o", "--outfile"):
            outfile = arg
        elif opt in ("-a", "--new_snow_alb_supplied"):
            new_snow_alb_supplied = True
        elif opt in ("-s", "--snow_bands"):
            snow_bands = True


    # Define fill values
    fill_value_mask = np.int32(0)
    fill_value_int = np.int32(-9999)
    fill_value_str = ''
#NOTE: maybe should put this back?
#    fill_value_float = 9.96920996838687e+36
    fill_value_float = np.nan

    # Default soil/snow parameters for cells that exist in mask but not in soil parameter file
    default = {
        'run_cell': 0,
        'gridcell': 0,
        'infilt': 0.1,
        'Ds': 0.01,
        'Dsmax': 10.0,
        'Ws': 0.75,
        'c': 2.0,
        'avg_T': 20,
        'dp': 4.0,
        'rough': 0.01,
        'snow_rough': 0.03,
        'annual_prec': fill_value_float,
        'fs_active': 0,
        'new_snow_albedo': 0.85,
        'expt': [12,12,12,],
        'Ksat': [1200,1200,1200],
        'phi_s': [-99,-99,-99],
        'init_moist': [40,120,400],
        'depth': [0.1,0.3,1.0],
        'bubble': [14.0,14.0,14.0],
        'quartz': [0.69,0.69,0.69],
        'bulk_density': [1570,1570,1570],
        'soil_density': [2620,2620,2620],
        'Wcr_FRACT': [0.39,0.39,0.39],
        'Wpwp_FRACT': [0.23,0.23,0.23],
        'resid_moist': [0,0,0],
        'AreaFract': 1.0,
        'Pfactor': 1.0,
    }

    # Define groups of variables
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
# Getting elev from domain file
#                        'elev',
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
# Getting Nveg from internal computation
#                       'Nveg',
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
# Getting Cv from climfile
#                       'Cv',
                       'overstory',
                       'rarc',
                       'rmin',
                       'wind_h',
                       'RGL',
                       'rad_atten',
                       'wind_atten',
                       'trunk_ratio',
                      ]

    varnames_4d_root = [
                        'root_depth',
                        'root_fract',
                       ]

    varnames_4d_clim = [
# Getting these from climfile
#                        'LAI',
#                        'fcanopy',
#                        'albedo',
                        'veg_rough',
                        'displacement',
                       ]

    # Not actually in the climfile, but derived from Cv
    varnames_2d_climfile = ['Nveg']
    varnames_3d_climfile = ['Cv']
    varnames_4d_climfile = ['LAI','fcanopy','albedo']

    varnames = []
    # Need to have Cv and Nveg first in the loop
    for varname in varnames_3d_climfile:
        varnames.append(varname)
    for varname in varnames_2d_climfile:
        varnames.append(varname)
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
    for varname in varnames_4d_climfile:
        varnames.append(varname)

    varnames_soil = []
    for varname in varnames_1d_soil:
        varnames_soil.append(varname)
    if snow_bands:
        for varname in varnames_1d_snow:
            varnames_soil.append(varname)
    for varname in varnames_2d_soil:
        varnames_soil.append(varname)
    for varname in varnames_3d_soil:
        varnames_soil.append(varname)
    if snow_bands:
        for varname in varnames_3d_snow:
            varnames_soil.append(varname)

    varnames_veg = []
    for varname in varnames_1d_veg:
        varnames_veg.append(varname)
    for varname in varnames_2d_veg:
        varnames_veg.append(varname)
    for varname in varnames_3d_veg:
        varnames_veg.append(varname)
    for varname in varnames_4d_root:
        varnames_veg.append(varname)
    for varname in varnames_4d_clim:
        varnames_veg.append(varname)
    for varname in varnames_2d_climfile:
        varnames_veg.append(varname)
    for varname in varnames_3d_climfile:
        varnames_veg.append(varname)
    for varname in varnames_4d_climfile:
        varnames_veg.append(varname)

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
    for varname in varnames_2d_climfile:
        varnames_veg.append(varname)

    varnames_int = [
                    'layer',
                    'snow_band',
                    'run_cell',
                    'gridcell',
                    'fs_active',
                    'Nveg',
                    'overstory',
                   ]

    varnames_str = ['veg_descr']

    # Read domain file
    ds_dom = xr.open_dataset(domainfile)

    # Read relevant dims and vars
    lat = ds_dom['lat']
    lon = ds_dom['lon']
    # Ensure that mask uses preferred value for nodata
    ds_dom['mask'].encoding['_FillValue'] = fill_value_mask
    mask = np.where(ds_dom['mask']==1,1,fill_value_mask).astype(np.int32)
    frac = ds_dom['frac']
    elev = ds_dom['elev']
    elev = np.where(mask==1,elev,fill_value_float)
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

    # Open original parameter file
    ds_param = xr.open_dataset(paramfile)

    # Read coord vars
    layer = ds_param['layer']
    root_zone_old = ds_param['root_zone']
    month_old = ds_param['month']
    veg_class_old = ds_param['veg_class']
    nLayer = len(layer)
    nRoot_old = len(root_zone_old)
    nMonth_old = len(month_old)
    nClass_old = len(veg_class_old)
    if snow_bands:
        snow_band = ds_param['snow_band']
        nBand = len(snow_band)
    else:
        nBand = 1
        snow_band = [1]

    # Read geographic coords
    lat_param = ds_param['lat']
    lon_param = ds_param['lon']
    nLat_param = len(lat_param)
    nLon_param = len(lon_param)
    minlat_param = np.asscalar(lat_param[0])
    maxlat_param = np.asscalar(lat_param[-1])
    minlon_param = np.asscalar(lon_param[0])
    maxlon_param = np.asscalar(lon_param[-1])
    cellsize_param = ( maxlat_param - minlat_param ) / (nLat_param - 1)
    minlat_param -= (0.5 * cellsize_param)
    minlon_param -= (0.5 * cellsize_param)
    maxlat_param += (0.5 * cellsize_param)
    maxlon_param += (0.5 * cellsize_param)
    llcorner_param = [minlat_param, minlon_param]
    mask_param = ds_param['mask']

    # A few hacks
    for varname in ['annual_prec','fs_active']:
        ds_param[varname][:] = default[varname]
        ds_param[varname][:] *= mask_param

    # Read new veg library file
    lib_data = np.loadtxt(libfile, dtype=bytes, delimiter=',').astype(str)
    nClass = lib_data.shape[0]

    # Read new root zone file
    root_data = np.loadtxt(rootfile, dtype=bytes, delimiter=' ').astype(str)
    if (root_data.shape[0] != nClass):
        print('ERROR: number of classes in',rootfile,'(',root_data.shape[0],') not equal to number of classes in',libfile,'(',nClass,')')
    nRoot = int((root_data.shape[1]-1)/2)
    root_zone = np.arange(nRoot).astype(np.int32) + 1

    # Store veg lib params in dict
    veg_lib_params = {}
    veg_lib_params['veg_class'] = np.empty([nClass,nLat,nLon],np.str)
    veg_lib_params['overstory'] = np.empty([nClass,nLat,nLon],np.int32)
    veg_lib_params['rarc'] = np.empty([nClass,nLat,nLon],np.double)
    veg_lib_params['rmin'] = np.empty([nClass,nLat,nLon],np.double)
#    veg_lib_params['LAI'] = np.empty([nClass,12,nLat,nLon],np.double)
#    veg_lib_params['albedo'] = np.empty([nClass,12,nLat,nLon],np.double)
    veg_lib_params['veg_rough'] = np.empty([nClass,12,nLat,nLon],np.double)
    veg_lib_params['displacement'] = np.empty([nClass,12,nLat,nLon],np.double)
    veg_lib_params['wind_h'] = np.empty([nClass,nLat,nLon],np.double)
    veg_lib_params['RGL'] = np.empty([nClass,nLat,nLon],np.double)
    veg_lib_params['rad_atten'] = np.empty([nClass,nLat,nLon],np.double)
    veg_lib_params['wind_atten'] = np.empty([nClass,nLat,nLon],np.double)
    veg_lib_params['trunk_ratio'] = np.empty([nClass,nLat,nLon],np.double)
    veg_lib_params['veg_descr'] = np.empty([nClass],np.str)
    for i in range(nClass):
        veg_lib_params['veg_class'][i] = lib_data[i,0]
        veg_lib_params['overstory'][i] = lib_data[i,1].astype(np.int32)
        veg_lib_params['rarc'][i] = np.around(lib_data[i,2].astype(np.double),6)
        veg_lib_params['rmin'][i] = np.around(lib_data[i,3].astype(np.double),6)
        for j in range(12):
#            veg_lib_params['LAI'][i,j] = np.around(lib_data[i,4+j].astype(np.double),6)
#            veg_lib_params['albedo'][i,j] = np.around(lib_data[i,16+j].astype(np.double),6)
            veg_lib_params['veg_rough'][i,j] = np.around(lib_data[i,28+j].astype(np.double),6)
            veg_lib_params['displacement'][i,j] = np.around(lib_data[i,40+j].astype(np.double),6)
        veg_lib_params['wind_h'][i] = np.around(lib_data[i,52].astype(np.double),6)
        veg_lib_params['RGL'][i] = np.around(lib_data[i,53].astype(np.double),6)
        veg_lib_params['rad_atten'][i] = np.around(lib_data[i,54].astype(np.double),6)
        veg_lib_params['wind_atten'][i] = np.around(lib_data[i,55].astype(np.double),6)
        veg_lib_params['trunk_ratio'][i] = np.around(lib_data[i,56].astype(np.double),6)
        veg_lib_params['veg_descr'][i] = lib_data[i,57]

    veg_lib_params['root_depth'] = np.empty([nClass,nRoot,nLat,nLon],np.double)
    veg_lib_params['root_fract'] = np.empty([nClass,nRoot,nLat,nLon],np.double)
    for i in range(nClass):
        for j in range(nRoot):
            veg_lib_params['root_depth'][i,j] = np.around(root_data[i,1+2*j].astype(np.double),4)
            veg_lib_params['root_fract'][i,j] = np.around(root_data[i,2+2*j].astype(np.double),4)

    i = fill_class
    default['veg_class'] = lib_data[i,0]
    default['overstory'] = lib_data[i,1].astype(np.int32)
    default['rarc'] = np.around(lib_data[i,2].astype(np.double),6)
    default['rmin'] = np.around(lib_data[i,3].astype(np.double),6)
    default['LAI'] = np.empty(12)
    default['fcanopy'] = np.empty(12)
    default['albedo'] = np.empty(12)
    default['veg_rough'] = np.empty(12)
    default['displacement'] = np.empty(12)
    for j in range(12):
        default['LAI'][j] = np.around(lib_data[i,4+j].astype(np.double),6)
        default['fcanopy'][j] = 0.5
        default['albedo'][j] = np.around(lib_data[i,16+j].astype(np.double),6)
        default['veg_rough'][j] = np.around(lib_data[i,28+j].astype(np.double),6)
        default['displacement'][j] = np.around(lib_data[i,40+j].astype(np.double),6)
    default['wind_h'] = np.around(lib_data[i,52].astype(np.double),6)
    default['RGL'] = np.around(lib_data[i,53].astype(np.double),6)
    default['rad_atten'] = np.around(lib_data[i,54].astype(np.double),6)
    default['wind_atten'] = np.around(lib_data[i,55].astype(np.double),6)
    default['trunk_ratio'] = np.around(lib_data[i,56].astype(np.double),6)
    default['veg_descr'] = lib_data[i,57]
    default['root_depth'] = np.empty(nRoot)
    default['root_fract'] = np.empty(nRoot)
    for j in range(nRoot):
        default['root_depth'][j] = np.around(root_data[i,1+2*j].astype(np.double),4)
        default['root_fract'][j] = np.around(root_data[i,2+2*j].astype(np.double),4)

    # Open monthly climatological seasonal cycle
    ds_clim = xr.open_dataset(climfile)

    # Read coord vars
    veg_class = ds_clim['veg_class']
    class_name = ds_clim['class_name']
    nClass_clim = len(veg_class)
    nMonth = 12
    month = np.arange(1,nMonth+1)
    if (nClass_clim != nClass):
        print('ERROR: number of classes in',climfile,'(',nClass_clim,') not equal to number of classes in',libfile,'(',nClass,')')

    # Read geographic coords
    lat_clim = ds_clim['lat']
    lon_clim = ds_clim['lon']
    nLat_clim = len(lat_clim)
    nLon_clim = len(lon_clim)
    minlat_clim = np.asscalar(lat_clim[0])
    maxlat_clim = np.asscalar(lat_clim[-1])
    minlon_clim = np.asscalar(lon_clim[0])
    maxlon_clim = np.asscalar(lon_clim[-1])
    cellsize_clim = ( maxlat_clim - minlat_clim ) / (nLat_clim - 1)
    minlat_clim -= (0.5 * cellsize_clim)
    minlon_clim -= (0.5 * cellsize_clim)
    maxlat_clim += (0.5 * cellsize_clim)
    maxlon_clim += (0.5 * cellsize_clim)
    llcorner_clim = [minlat_clim, minlon_clim]

    # Prepare output ds
    ds_out = xr.Dataset(
                        {
                         'mask': (['lat', 'lon'], mask),
                         'elev': (['lat', 'lon'], elev),
                        },
                        coords={
                                'lat': (['lat'], lat),
                                'lon': (['lon'], lon),
                                'month': (['month'], month),
                                'layer': (['nlayer'], layer),
                                'snow_band': (['snow_band'], snow_band),
                                'veg_class': (['veg_class'], veg_class),
                                'veg_descr': (['veg_class'], class_name),
                                'root_zone': (['root_zone'], root_zone),
                               }
                       )

    for varname in ['lat', 'lon', 'mask', 'elev']:
        ds_out[varname].attrs = ds_dom[varname].attrs
        ds_out[varname].encoding = ds_dom[varname].encoding
    param_coordvarnames = ['layer', 'root_zone']
    if snow_bands:
        param_coordvarnames.append('snow_band')
    for varname in param_coordvarnames:
        ds_out[varname].attrs = ds_param[varname].attrs
        ds_out[varname].encoding = ds_param[varname].encoding
    for varname in ['veg_class']:
        ds_out[varname].attrs = ds_clim[varname].attrs
        ds_out[varname].encoding = ds_clim[varname].encoding
    ds_out['month'].attrs = ds_clim['month_of_year'].attrs
    ds_out['month'].encoding = ds_clim['month_of_year'].encoding
    ds_out['veg_descr'].attrs = ds_clim['class_name'].attrs
    if ~snow_bands:
        ds_out['snow_band'].attrs['long_name'] = 'Snow elevation band'
  

    # Do the subsampling
    out_dict = {}
    for varname in varnames:

        print('varname',varname)
        if varname in varnames_1d:
#            print('1d var:',varname)
#            if (varname == 'layer'):
#                shape = [nLayer]
#                dimstr = ['nlayer']
#            elif (varname == 'snow_band'):
#                shape = [nBand]
#                dimstr = ['snow_band']
#            elif (varname == 'veg_descr'):
#                shape = [nClass]
#                dimstr = ['veg_class']
            continue
        elif varname in varnames_2d or varname in varnames_2d_climfile:
            print('2d var:', varname)
            shape = [nLat, nLon]
            dimstr = ['lat', 'lon']
        elif varname in varnames_3d_soil:
            print('3d_soil var:', varname)
            shape = [nLayer, nLat, nLon]
            dimstr = ['nlayer', 'lat', 'lon']
        elif varname in varnames_3d_snow:
            print('3d_snow var:', varname)
            shape = [nBand, nLat, nLon]
            dimstr = ['snow_band', 'lat', 'lon']
        elif varname in varnames_3d_veg or varname in varnames_3d_climfile:
            print('3d_veg var:', varname)
            shape = [nClass, nLat, nLon]
            dimstr = ['veg_class', 'lat', 'lon']
        elif varname in varnames_4d_root:
            print('4d_root var:', varname)
            shape = [nClass, nRoot, nLat, nLon]
            dimstr = ['veg_class', 'root_zone', 'lat', 'lon']
        elif varname in varnames_4d_clim or varname in varnames_4d_climfile:
            print('4d_clim var:', varname)
            shape = [nClass, nMonth, nLat, nLon]
            dimstr = ['veg_class', 'month', 'lat', 'lon']
        else:
            print('no case for',varname)

        if varname in varnames_int:
            dtype = np.int32
            fill_value = fill_value_int
        elif varname in varnames_str:
            dtype = np.str
            fill_value = fill_value_str
        else:
            dtype = np.double
            fill_value = fill_value_float

        if varname in varnames_soil:
            fill_tmp = ds_param[varname].encoding['_FillValue']
            data_tmp = np.where(ds_param[varname] == fill_tmp, fill_value, ds_param[varname])
            ds_param[varname][:] = data_tmp[:]
            out_dict[varname] = subsample_variable(ds_param[varname],
                                                   shape,
                                                   dtype,
                                                   llcorner_param,
                                                   cellsize_param,
                                                   llcorner,
                                                   cellsize,
                                                   mask, fill_value)
        elif varname in varnames_veg:
            if varname in veg_lib_params.keys():
                fill_tmp = ds_param[varname].encoding['_FillValue']
                data_tmp = np.where(ds_param[varname] == fill_tmp, fill_value, ds_param[varname])
                ds_param[varname][:] = data_tmp[:]
                # Mask veg params to only exist where associated veg is present
                if varname in varnames_3d_veg:
                    veg_lib_params[varname][Cv_mask != 1] = fill_value
                elif varname in varnames_4d_clim:
                    for m in range(nMonth):
                        tmp = np.empty([nClass, nLat, nLon], dtype=np.double)
                        tmp[:] = veg_lib_params[varname][:,m,:,:]
                        tmp[Cv_mask != 1] = fill_value
                        veg_lib_params[varname][:,m,:,:] = tmp[:]
                elif varname in varnames_4d_root:
                    for i in range(nRoot):
                        tmp = np.empty([nClass, nLat, nLon], dtype=np.double)
                        tmp[:] = veg_lib_params[varname][:,i,:,:]
                        tmp[Cv_mask != 1] = fill_value
                        veg_lib_params[varname][:,i,:,:] = tmp[:]
                out_dict[varname] = subsample_variable(veg_lib_params[varname],
                                                       shape,
                                                       dtype,
                                                       llcorner,
                                                       cellsize,
                                                       llcorner,
                                                       cellsize,
                                                       mask, fill_value)
            elif varname in varnames_3d_climfile:
                fill_tmp = ds_clim[varname].encoding['_FillValue']
                data_tmp = np.where(ds_clim[varname] == fill_tmp, fill_value, ds_clim[varname])
                ds_clim[varname][:] = data_tmp[:]
                out_dict[varname] = subsample_variable(ds_clim[varname],
                                                       shape,
                                                       dtype,
                                                       llcorner_clim,
                                                       cellsize_clim,
                                                       llcorner,
                                                       cellsize,
                                                       mask, fill_value)
                if varname == 'Cv':
                    # Might be worth validating Cv between 0 and 1
                    if open_water_class >= 0:
                        # Remove open water where land fraction < 1.0
                        # This assumes that open water corresponding to ocean or
                        # large lakes caused the land fraction to be < 1.0, so
                        # that for consistency with frac, we must remove it.
                        # This means that the remaining classes present in the
                        # cell represent the land portion only.
                        # NOTE: This logic assumes that frac is nan where mask != 1
                        frac_tmp = np.where(np.isnan(frac), 100, frac)
                        a = np.where(frac_tmp < 1.0)[0]
                        b = np.where(frac_tmp < 1.0)[1]
                        for i,j in zip(a,b):
                            if not np.isnan(out_dict['Cv'][open_water_class,i,j]):
                                out_dict['Cv'][open_water_class,i,j] = 0

                    # Normalize Cv such that Cv sums exactly to 1.0
                    Cv_tmp = np.empty(out_dict['Cv'].shape)
                    Cv_tmp[:] = out_dict['Cv'][:]
                    a = np.where(mask != 1)[0]
                    b = np.where(mask != 1)[1]
                    for i,j in zip(a,b):
                        out_dict['Cv'][:,i,j] = np.nan
                        Cv_tmp[0,i,j] = 1.0
                    Cv_sum = np.nansum(Cv_tmp, axis=0)
                    a = np.where(Cv_sum != 1.0)[0]
                    b = np.where(Cv_sum != 1.0)[1]
                    for i,j in zip(a,b):
                        new_sum = 0
                        for v in range(nClass):
                            if not np.isnan(out_dict['Cv'][v,i,j]) and \
                                out_dict['Cv'][v,i,j] > 0:
                                out_dict['Cv'][v,i,j] /= Cv_sum[i,j]
                                new_sum += out_dict['Cv'][v,i,j]
                        error = new_sum - 1.0
                        if error != 0.0:
                            Cv_max = max(out_dict['Cv'][:,i,j])
                            for v in range(nClass):
                                if out_dict['Cv'][v,i,j] == Cv_max and \
                                    Cv_max > error:
                                    out_dict['Cv'][v,i,j] -= error
                                    break
                    Cv_tmp[:] = out_dict['Cv'][:]
                    a = np.where(mask != 1)[0]
                    b = np.where(mask != 1)[1]
                    for i,j in zip(a,b):
                        Cv_tmp[0,i,j] = 1.0
                    Cv_sum = np.nansum(Cv_tmp, axis=0)
                    a = np.where(Cv_sum != 1.0)[0]
                    b = np.where(Cv_sum != 1.0)[1]
                    for i,j in zip(a,b):
                        if Cv_sum[i,j] == 0.0:
                            out_dict['Cv'][:,i,j] = np.nan
                        print('WARNING: Normalization failed; Cv_sum is',Cv_sum[i,j],'at i,j',i,j)
                        
            elif varname == 'Nveg':
                Cv_mask = np.zeros([nClass,nLat,nLon],dtype=np.int32)
                Cv_mask[out_dict['Cv'] > 0] = 1
                Cv_mask[out_dict['Cv'] == fill_value_float] = 0
                Nveg = np.sum(Cv_mask, axis=0, dtype=np.int32).astype(np.int32)
                Nveg[Nveg == 0] = fill_value_int
                out_dict['Nveg'] = Nveg
            elif varname == 'root_fract':
                # Might be worth validating root_fract between 0 and 1
                # Normalize root_fract such that it sums exactly to 1.0
                root_tmp = np.empty(out_dict['root_fract'].shape)
                root_tmp[:] = out_dict['root_fract'][:]
                a = np.where(mask != 1)[0]
                b = np.where(mask != 1)[1]
                for i,j in zip(a,b):
                    out_dict['root_fract'][...,i,j] = np.nan
                    root_tmp[:,0,i,j] = 1.0
                root_sum = np.nansum(root_tmp, axis=1)
                a = np.where(root_sum != 1.0)[0]
                b = np.where(root_sum != 1.0)[1]
                c = np.where(root_sum != 1.0)[2]
                for v,i,j in zip(a,b,c):
                    new_sum = 0
                    for r in range(nRoot):
                        if not np.isnan(out_dict['root_fract'][v,r,i,j]) and \
                            out_dict['root_fract'][v,r,i,j] > 0:
                            out_dict['root_fract'][v,r,i,j] /= root_sum[v,i,j]
                            new_sum += out_dict['root_fract'][v,r,i,j]
                    error = new_sum - 1.0
                    if error != 0.0:
                        root_max = max(out_dict['root_fract'][v,:,i,j])
                        for r in range(nRoot):
                            if out_dict['root_fract'][v,r,i,j] == root_max and \
                                root_max > error:
                                out_dict['root_fract'][v,r,i,j] -= error
                                break
                root_tmp[:] = out_dict['root_fract'][:]
                a = np.where(mask != 1)[0]
                b = np.where(mask != 1)[1]
                for i,j in zip(a,b):
                    out_dict['root_fract'][...,i,j] = np.nan
                    root_tmp[:,0,i,j] = 1.0
                root_sum = np.nansum(root_tmp, axis=1)
                a = np.where(root_sum != 1.0)[0]
                b = np.where(root_sum != 1.0)[1]
                b = np.where(root_sum != 1.0)[2]
                for v,i,j in zip(a,b,c):
                    if root_sum[v,i,j] == 0.0:
                        out_dict['root_fract'][v,:,i,j] = np.nan
                    print('WARNING: Normalization failed; root_sum is',root_sum[v,i,j],'at v,i,j',v,i,j)
            elif varname in varnames_4d_climfile:
                tmpvar = varname + '_mean'
                tmp = ds_clim[tmpvar]
                tmp2 = np.empty([tmp.shape[1],tmp.shape[0],tmp.shape[2],tmp.shape[3]], dtype=dtype)
                for c in range(nClass):
                    for m in range(nMonth):
                        tmp2[c,m] = tmp[m,c]
                out_dict[varname] = subsample_variable(tmp2,
                                                       shape,
                                                       dtype,
                                                       llcorner_clim,
                                                       cellsize_clim,
                                                       llcorner,
                                                       cellsize,
                                                       mask, fill_value)

        ds_out[varname] = (dimstr, out_dict[varname])
        if (varname in varnames_soil):
            ds_out[varname].attrs = ds_param[varname].attrs
        elif varname in varnames_veg:
            if varname in veg_lib_params.keys():
                ds_out[varname].attrs = ds_param[varname].attrs
            elif varname == 'Nveg':
                ds_out[varname].attrs = ds_param[varname].attrs
            elif varname in varnames_3d_climfile:
                ds_out[varname].attrs = ds_clim[varname].attrs
            elif varname in varnames_4d_climfile:
                tmpvar = varname + '_mean'
                ds_out[varname].attrs = ds_clim[tmpvar].attrs
        ds_out[varname].encoding['_FillValue'] = fill_value
        ds_out[varname].encoding['zlib'] = True

        if varname == 'off_gmt':
            out_dict[varname][:] = np.around(out_dict[varname], 2)

    # Assign default values for mismatches
    mismatches = np.empty([nLat,nLon], dtype=np.int32)
    mismatches[:] = mask[:]
    mismatches[out_dict['Ksat'][0] > 0] = 0
    a = np.where(mismatches == 1)[0]
    b = np.where(mismatches == 1)[1]
    for varname in default.keys():
        if varname in varnames_2d_soil:
            for i,j in zip(a,b):
                ds_out[varname][i,j] = default[varname]
        elif varname in varnames_3d_soil:
            for k in range(nLayer):
                for i,j in zip(a,b):
                    ds_out[varname][k,i,j] = default[varname][k]
        elif snow_bands and varname in varnames_3d_snow:
            for k in range(nBand):
                if k == 0:
                    for i,j in zip(a,b):
                        ds_out[varname][k,i,j] = default[varname]
                else:
                    for i,j in zip(a,b):
                        ds_out[varname][k,i,j] = 0
    for i,j in zip(a,b):
        sign = lon[j] / abs(lon[j])
        tmp = abs(lon[j]) * 12 / 180
        ds_out['off_gmt'][i,j] = sign * int(tmp + 0.5)
        if snow_bands:
            ds_out['elevation'][0,i,j] = elev[i,j]
    mismatches[:] = mask[:]
    mismatches[Cv_sum > 0] = 0
    a = np.where(mismatches == 1)[0]
    b = np.where(mismatches == 1)[1]
    for v in range(nClass):
        for i,j in zip(a,b):
            if v == fill_class:
                ds_out['Cv'][v,i,j] = 1.0
            else:
                ds_out['Cv'][v,i,j] = 0.0
    v = fill_class
    for i,j in zip(a,b):
        ds_out['Nveg'][i,j] = 1
    for varname in veg_lib_params.keys():
        if varname in varnames_3d_veg:
            for i,j in zip(a,b):
                ds_out[varname][v,i,j] = default[varname]
        elif varname in varnames_4d_clim:
            for m in range(nMonth):
                for i,j in zip(a,b):
                    ds_out[varname][v,m,i,j] = default[varname][m]
        elif varname in varnames_4d_root:
            for r in range(nRoot):
                for i,j in zip(a,b):
                    ds_out[varname][v,r,i,j] = default[varname][r]
    for varname in varnames_4d_climfile:
        for m in range(nMonth):
            for i,j in zip(a,b):
                ds_out[varname][v,m,i,j] = default[varname][m]


    if lakes:

        # Lake parameters
        nNodes = 10
        lake_node = np.arange(nNodes).astype(np.int32).tolist()
        lake_idx = np.full([nLat,nLon],fill_value_int,dtype=np.int32)
        numnod = np.full([nLat,nLon],fill_value_int,dtype=np.int32)
        mindepth = np.full([nLat,nLon],fill_value_float,dtype=np.double)
        wfrac = np.full([nLat,nLon],fill_value_float,dtype=np.double)
        rpercent = np.full([nLat,nLon],fill_value_float,dtype=np.double)
        depth_in = np.full([nLat,nLon],fill_value_float,dtype=np.double)
        basin_depth = np.full([nNodes,nLat,nLon],fill_value_float,dtype=np.double)
        basin_area = np.full([nNodes,nLat,nLon],fill_value_float,dtype=np.double)
        area_max = np.full([nLat,nLon],0,dtype=np.double)
        area_max[:,:] = ds_out['Cv'][lake_class,:,:]
        area_max[area_max == fill_value_float] = -1
#        a = np.where(area_max>0.25)[0]
#        b = np.where(area_max>0.25)[1]
#        for i,j in zip(a,b):
#            ds_out['Dsmax'][i,j] = 0
        lake_idx[area_max==0] = -1
        lake_idx[area_max>0] = 0
        numnod[area_max>0] = nNodes
        mindepth[area_max>0] = 0
        wfrac[area_max>0] = 0
        rpercent[area_max>0] = 0
        for i in range(nLat):
            for j in range(nLon):
                if area_max[i,j] > 0:
                    [tmp_depth,tmp_area] = compute_depth_area_rel(area_max[i,j],numnod[i,j],'lake')
                    basin_depth[:,i,j] = tmp_depth
                    basin_area[:,i,j] = tmp_area
                    depth_in[i,j] = tmp_depth[0]
        area_max[area_max == -1] = fill_value_float

        ds_out['lake_node'] = (['lake_node'],lake_node)
        ds_out['lake_node'].attrs['long_name'] = 'lake basin node'
        ds_out.set_coords('lake_node',inplace=True)
        ds_out['lake_idx'] = (['lat','lon'],lake_idx)
        ds_out['lake_idx'].attrs['long_name'] = 'index of veg tile that contains the lake/wetland'
        ds_out['lake_idx'].encoding['_FillValue'] = fill_value_int
        ds_out['lake_idx'].encoding['zlib'] = True
        ds_out['numnod'] = (['lat','lon'],numnod)
        ds_out['numnod'].attrs['long_name'] = 'Maxium number of lake layers in the grid cell'
        ds_out['numnod'].encoding['_FillValue'] = fill_value_int
        ds_out['numnod'].encoding['zlib'] = True
        ds_out['mindepth'] = (['lat','lon'],mindepth)
        ds_out['mindepth'].attrs['long_name'] = 'Minimum lake water depth for channel runoff to occur'
        ds_out['mindepth'].encoding['_FillValue'] = fill_value_float
        ds_out['mindepth'].encoding['zlib'] = True
        ds_out['wfrac'] = (['lat','lon'],wfrac)
        ds_out['wfrac'].attrs['long_name'] = 'Channel outlet width (expressed as a fraction of lake perimeter)'
        ds_out['wfrac'].encoding['_FillValue'] = fill_value_float
        ds_out['wfrac'].encoding['zlib'] = True
        ds_out['rpercent'] = (['lat','lon'],rpercent)
        ds_out['rpercent'].attrs['long_name'] = 'Initial lake depth'
        ds_out['rpercent'].encoding['_FillValue'] = fill_value_float
        ds_out['rpercent'].encoding['zlib'] = True
        ds_out['depth_in'] = (['lat','lon'],depth_in)
        ds_out['depth_in'].attrs['long_name'] = 'Fraction of runoff from other veg tiles that flows into the lake'
        ds_out['depth_in'].encoding['_FillValue'] = fill_value_float
        ds_out['depth_in'].encoding['zlib'] = True
        ds_out['basin_depth'] = (['lake_node','lat','lon'],basin_depth)
        ds_out['basin_depth'].attrs['long_name'] = 'Elevation (above lake bottom) of points on lake depth-area curve'
        ds_out['basin_depth'].encoding['_FillValue'] = fill_value_float
        ds_out['basin_depth'].encoding['zlib'] = True
        ds_out['basin_area'] = (['lake_node','lat','lon'],basin_area)
        ds_out['basin_area'].attrs['long_name'] = 'Surface area (expressed as fraction of grid cell area) of points on lake depth-area curve'
        ds_out['basin_area'].encoding['_FillValue'] = fill_value_float
        ds_out['basin_area'].encoding['zlib'] = True

    # Write to output file
    ds_out.to_netcdf(outfile)

    ds_dom.close()
    ds_param.close()
    ds_clim.close()
    ds_out.close()

if __name__ == "__main__":
    main()

