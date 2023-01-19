# !/usr/bin/env python
# coding: utf-8
# pylint: disable=C0411, C0413
# ^^Disables import order warnings because i cannot get rid of them
# pylint: disable=W0105
# ^^ Disables "Pointless" String Statements Warning
# pylint: disable=C0103
# ^^ Disables "Does not conform to snake case" for variables ds, js

'''
regridder.py
v1.3
Developer: Maggie Trimpin
Updated 03/23/2022
'''
from scipy.fft import dst
import xarray as xr
import xesmf as xe
import numpy as np
import json
import sys
import os
import stat
import time
import re
from multiprocessing.dummy import freeze_support
import warnings
warnings.filterwarnings("ignore")
sys.dont_write_bytecode = True

EXT_N4 = ["nc", "nc3", "nc4", "nc4c", "n4", "NC", "NC4"]

def process_input(input_json):
    '''
    Input: path to a .json file containing information about input nc, regridding parameters, output key

    Extracts information about input file/files, calls regrid() to process each accordingly. 
    Then, calls output() to save .nc output files to the system
    '''
    with open(input_json, encoding="utf8") as file:
        js = json.load(file)  # stores info from json file in js dictionary

    # =============== Extract info from json =======================
    nc_filepath = js['in_dir']  # path to input netcdf file
    # list of allowed extensions, ie. nc4, nc
    ext = nc_filepath.split(".")[-1]

    if ext in EXT_N4:  # input filepath is netcdf file, process single file
        ds = xr.open_dataset(nc_filepath)
        regridded_ds = regrid(ds, js)
        output(regridded_ds, js)

    else:  # directory of netcdf files, process one by one
        print("multiple nc4's")
        for file in os.listdir(nc_filepath):
            ext = file.split(".")[-1]
            if ext in EXT_N4:  # check all files in directory to ensure nc4
                # extract information from each nc4 file in directory
                name = os.path.basename(file)
                print(f"======File: {name} ========")
                file_path = os.path.join(nc_filepath, file)
                # regrid each dataset
                ds = xr.open_dataset(file_path)
                regridded_ds = regrid(ds, js)
                # output regridded datasets
                output(regridded_ds, js)


def regrid(ds, js):
    '''
    Called in the process_input function
    Input: current xarray dataset for processing, json file with regridding parameters

    Interprets interpolation and resolution information from json file. 
    Creates a gridfile from resolution parameters found in pregrid_generate()
    
    If mask == True, converts each DataArray in ds.data_vars to a np.masked_array.
    Once conversion is complete, creates a new dataset of masked arrays, and regrids them accordingly. 

    If mask == False, regrids the input dataset without masking

    Returns: new regridded dataset
    '''
    # js['regridding'] has info in form ['interpolation', 'resolution']
    interpolation = js['regridding'][0]
    # interpolation options = bilinear, nearest neighbor
    # interperet interpolations from params.py
    if interpolation == 'remapbil':
        interpolation = 'bilinear'
    elif interpolation == 'remapnn':
        interpolation = 'nearest_s2d'
    elif interpolation == 'remapcon':
        interpolation = 'conservative'
    else:
        print(f"Interpolation {interpolation} \
            not supported by this regridder.")

    # generate timestamp for temp_dir
    # create temporary directory to hold gridfile
    if not os.path.exists("../tmp"):
        os.mkdir("../tmp")
    
    js["TMPDIR"] = f"../tmp/regridder_temporary_{time.time()}"
    os.mkdir(js["TMPDIR"])

    js["gridfile"] = pregrid_generate(js)  # Create resolution gridfile
    # Read gridfile and extract values
    with open(js["gridfile"], 'r', encoding="utf8") as gridfile:
        gridfile_list = gridfile.readlines()
    str_values = "".join(gridfile_list)
    values = re.findall(r'-?\d*\.\d+|\d+|-?\d+', str_values)
    lonll = float(values[2]) # longitude lower limit (x_start)
    lonul = float(values[2]) + float(values[0]) * float(values[3])
    # upper limit = lower limit + (# of values * step size)
    lon_res = float(values[3]) # longitude step size
    latll = float(values[4]) # latitude lower limit (y_start)
    latul = float(values[4]) + float(values[1]) * float(values[5])
    # upper limit = lower limit + (# of values * step size)
    lat_res = float(values[5]) # latitude step size

    lon_count = int(values[0]) # number of steps in longitude coords (x_size)
    lat_count = int(values[1]) # number of steps in latitude coords (y_size)

    #keep latitude within -90, 90 range
    print(f"latll = {latll}, latul = {latul}")
    if latll < -90:
        latll = latll + lat_res
    if latul > 90:
        latul = latul - lat_res

    #define output grid 
    if interpolation != "conservative":
        ds_out = xr.Dataset(
                data_vars=ds.data_vars,
                coords={'lon': (['lon'], np.arange(lonll, lonul, lon_res)),
                                'lat': (['lat'], np.arange(latll, latul, lat_res)),
                                },
                attrs=ds.attrs,
            )
    
    else:
        #if conservative interpolation, define bounds for input and output grid
        #collect resolution from input grid
        input_lon_res = ds.lon.values[1] - ds.lon.values[0]
        input_lat_res = ds.lat.values[1] - ds.lat.values[0]
        input_lon_count = len(ds.lon.values)
        input_lat_count = len(ds.lat.values)
        
        ds = ds.assign({
                        'lon_b':(['lon_b'], np.linspace(min(ds.lon)-input_lon_res/2, max(ds.lon)+input_lon_res/2, input_lon_count+1)),  
                        'lat_b': (['lat_b'], np.linspace(latll-input_lat_res/2, latul+input_lat_res/2, input_lat_count+1))
                        })
        
        ds_out = xr.Dataset(
                data_vars=ds.data_vars,
                coords={'lon': (['lon'], np.linspace(lonll, lonul, lon_count,endpoint=True)),
                        'lat': (['lat'], np.linspace(latll, latul, lat_count, endpoint=True)),
                        'lon_b': (['lon_b'], np.linspace(lonll-lon_res/2, lonul+lon_res/2, lon_count+1)),
                        'lat_b': (['lat_b'], np.linspace(latll-lat_res/2, latul+lat_res/2, lat_count+1))
                        },
                attrs=ds.attrs,
            )

    #create regridder object with input and target grids
    regridder = xe.Regridder(ds, ds_out, interpolation) #extrap_method="inverse_dist") 
    #process dataset with regridder
    ds_out = regridder(ds, keep_attrs=True)
 
    '''
    #tiny issues in western hemisphere, might need to manually rotate globe
    #re-adjust "new_ds" grid range (ie. -180:180 to 0:360)
    
    original = 180
    if original == 360:
        ds_out.coords['lon'] = (ds_out.coords['lon'] + 180) % 360 - 180 #convert from 0:360 to -180:180
        ds_out = ds_out.sortby(ds_out.lon)
        print(ds_out)
    elif original == 180:
        ds_out.coords['lon'] = (ds_out.coords['lon'] + 180) % 360 #convert from -180:180 to 0:360
        ds_out = ds_out.sortby(ds_out.lon)
        ds_out.assign_attrs(ds.attrs)
        print(ds_out)

    '''
    return ds_out


def output(ds, js):
    '''
    Outputs the processed dataset to .nc file on the system.

    Input: dataset, json file
    If output key is empty, outputs the file in the form "json file name + in_dir"
    to a temp directory. 

    If output key exists, outputs the file in the form "json name + in_dir + outkey + .nc4"
    '''
    debug = js['debug']  # if true, output differences. if false output nothing
    out_key = js['out_key']  # out key used for output identification
    path = "../tmp"  # output files to temporary directory
    nc_filename = os.path.basename(js['in_dir'])

    if out_key == "":  # if empty, output file name = json file name + in_dir
        output_path = os.path.join(path, f"regridder_XESMF_{os.path.splitext(nc_filename)[0]}")
        if debug:
            print(f"output to: {output_path}")
    else:  # otherwise the out filename = json name + in_dir + outkey + .nc4
        output_path = os.path.join(path, f"regridder_XESMF_{os.path.splitext(nc_filename)[0]}_{out_key}.nc4")
        if debug:
            print(f"output to: {output_path}")

    if os.path.isfile(output_path):
        # change file permissions so can be overwritten
        os.chmod(output_path, stat.S_IRWXU)
    ds.to_netcdf(output_path)  # , encoding={'_FillValue': None})


def pregrid_generate(config):
    '''
    Creates gridfile based on the given resolution for regridding
    '''
    # creates gridfile with grid/ resolution info, returns path to gridfile
    # Error string
    param = config["regridding"]
    OUTDIR = config["TMPDIR"]
    gridtype  = 'lonlat'

    # Set grid template values
    if param[1] == 'JRA-55':
        xsize     = 288
        ysize     = 145
        xfirst    = 0
        xinc      = 1.25
        yfirst    = 90
        yinc      = -1.25
    elif param[1] == '20cr2x2':
        # NOAA (ESRL) Version
        xsize     = 180
        ysize     = 91
        xfirst    = 0
        xinc      = 2
        yfirst    = 90
        yinc      = -2
    elif param[1] == 'MERRA0.5':
        # DISC Version
        xsize     = 540
        ysize     = 361
        xfirst    = -180
        xinc      = 0.666666666666
        yfirst    = -90
        yinc      = 0.5
    elif param[1] == 'MERRA1.25':
        # DISC Version
        xsize     = 288
        ysize     = 144
        xfirst    = -179.375
        xinc      = 1.25
        yfirst    = -89.375
        yinc      = 1.25
    elif param[1] == 'MERRA2':
        # DISC Version
        xsize     = 576
        ysize     = 361
        xfirst    = -180
        xinc      = 0.625
        yfirst    = -90
        yinc      = 0.5
    elif param[1] == 'gpcp2.5':
        xsize     = 144
        ysize     = 72
        xfirst    = 1.25
        xinc      = 2.5
        yfirst    = -88.75
        yinc      = 2.5
    elif param[1] == 'cfsr0.5a':
        xsize     = 720
        ysize     = 361
        xfirst    = 0
        xinc      = 0.5
        yfirst    = 90
        yinc      = -0.5
    elif param[1] == 'cfsr0.5b':
        xsize     = 720
        ysize     = 360
        xfirst    = 0.25
        xinc      = 0.5
        yfirst    = 89.75
        yinc      = -0.5
    elif param[1] == 'cfsr1.0':
        xsize     = 360
        ysize     = 181
        xfirst    = 0
        xinc      = 1.0
        yfirst    = 90
        yinc      = -1.0
    elif param[1] == 'cfsr2.5':
        xsize     = 144
        ysize     = 73
        xfirst    = 0
        xinc      = 2.5
        yfirst    = 90
        yinc      = -2.5
    elif param[1] == 'ncepncar2.5':
        # NOAA (ESRL) Version
        xsize     = 144
        ysize     = 73
        xfirst    = 0
        xinc      = 2.5
        yfirst    = 90
        yinc      = -2.5
    elif param[1] == 'geos1x125':
        # http://acmg.seas.harvard.edu/geos/doc/man/appendix_2.html#GMAO_4x5
        xsize     = 288
        ysize     = 181
        xfirst    = -180
        xinc      = 1.25
        yfirst    = -90
        yinc      = 1.0
    elif param[1] == 'geos1x1':
        # http://acmg.seas.harvard.edu/geos/doc/man/appendix_2.html#GMAO_4x5
        xsize     = 360
        ysize     = 181
        xfirst    = -180
        xinc      = 1.0
        yfirst    = -90
        yinc      = 1.0
    elif param[1] == 'geos4x5':
        # http://acmg.seas.harvard.edu/geos/doc/man/appendix_2.html#GMAO_4x5
        xsize     = 72
        ysize     = 46
        xfirst    = -180
        xinc      = 5
        yfirst    = -90
        yinc      = 4
    elif param[1] == 'geos2x25':
        # http://acmg.seas.harvard.edu/geos/doc/man/appendix_2.html#GMAO_4x5
        xsize     = 144
        ysize     = 91
        xfirst    = -180
        xinc      = 2.5
        yfirst    = -90
        yinc      = 2.0
    elif param[1] == 'geos0.25':
        # http://acmg.seas.harvard.edu/geos/doc/man/appendix_2.html#GMAO_4x5
        xsize     = 1152
        ysize     = 721
        xfirst    = -180
        xinc      = 0.3125
        yfirst    = -90
        yinc      = 0.25
    elif param[1] == 'geos0.5':
        # http://acmg.seas.harvard.edu/geos/doc/man/appendix_2.html#GMAO_4x5
        xsize     = 540
        ysize     = 361
        xfirst    = -180
        xinc      = 0.666666667
        yfirst    = -90
        yinc      = 0.5
    elif param[1] == 'fv1x125':
        # Extrapolated from geos1x125
        xsize     = 288
        ysize     = 181
        xfirst    = 0
        xinc      = 1.25
        yfirst    = -90
        yinc      = 1.0
    elif param[1] == 'fv2x25':
        # Extrapolated from geos2x25
        xsize     = 144
        ysize     = 91
        xfirst    = 0
        xinc      = 2.5
        yfirst    = -90
        yinc      = 2.0
    elif param[1] == 'fv4x5':
        # Extrapolated from geos4x5
        xsize     = 72
        ysize     = 46
        xfirst    = 0
        xinc      = 5.0
        yfirst    = -90
        yinc      = 4.0
    elif param[1] == 'ERA-40':
        # ECMWF Version
        xsize     = 144
        ysize     = 73
        xfirst    = 0
        xinc      = 2.5
        yfirst    = 90
        yinc      = -2.5
    elif param[1] == 'ERA2.5':
        # ECMWF Version
        xsize     = 144
        ysize     = 73
        xfirst    = 0
        xinc      = 2.5
        yfirst    = 90
        yinc      = -2.5
    elif param[1] == 'ERA-I':
        # ECMWF Version
        xsize     = 240
        ysize     = 121
        xfirst    = 0
        xinc      = 1.5
        yfirst    = 90
        yinc      = -1.5
    elif param[1] == 'ERA1.5':
        # ECMWF Version
        xsize     = 240
        ysize     = 121
        xfirst    = 0
        xinc      = 1.5
        yfirst    = 90
        yinc      = -1.5
    elif param[1] == 'ERA.75':
        # ECMWF Version
        xsize     = 480
        ysize     = 241
        xfirst    = 0
        xinc      = 0.75
        yfirst    = 90
        yinc      = -0.75
    elif param[1] == 'GPCC2.5':
        # NOAA Version
        xsize     = 144
        ysize     = 72
        xfirst    = 1.25
        xinc      = 2.5
        yfirst    = 88.75
        yinc      = -2.5
    elif param[1] == 'GPCC1.0':
        # NOAA Version
        xsize     = 360
        ysize     = 180
        xfirst    = 0.5
        xinc      = 1.0
        yfirst    = 89.5
        yinc      = -1.0
    elif param[1] == 'GPCC0.5':
        # NOAA Version
        xsize     = 720
        ysize     = 360
        xfirst    = 0.25
        xinc      = 0.5
        yfirst    = 89.75
        yinc      = -0.5
    elif param[1] == 'CMORPH0.25':
        xsize     = 1440
        ysize     = 480
        xfirst    = -180
        xinc      = 0.25
        yfirst    = -59.875
        yinc      = 0.25
    elif param[1] == 'NLDAS-2':
        # DISC Version
        xsize     = 464
        ysize     = 224
        xfirst    = -124.9375
        xinc      = 0.125
        yfirst    = 25.0625
        yinc      = 0.125
    elif param[1] == 'GLDAS-2_0.25':
        # DISC Version
        xsize     = 1440
        ysize     = 600
        xfirst    = -179.875
        xinc      = 0.25
        yfirst    = -59.875
        yinc      = 0.25
    elif param[1] == 'GLDAS-2_1':
        # DISC Version
        xsize     = 360
        ysize     = 150
        xfirst    = -179.5
        xinc      = 1.0
        yfirst    = -59.5
        yinc      = 1.0
    elif param[1] == 'NCA-LDAS':
        # DISC Version
        xsize     = 464
        ysize     = 224
        xfirst    = -124.9375
        xinc      = 0.125
        yfirst    = 25.0625
        yinc      = 0.125
    elif param[1] == 'GFDL':
        xsize     = 288
        ysize     = 180
        xfirst    = 0.625
        xinc      = 1.25
        yfirst    = -89.5
        yinc      = 1.
    elif param[1] == 'gpcp3':
        xsize     = 720
        ysize     = 360
        xfirst    = -179.75
        xinc      = 0.5
        yfirst    = 89.75
        yinc      = -0.5
    elif param[1] == 'TMPA':
        xsize     = 1440
        ysize     = 400
        xfirst    = -179.875
        xinc      = 0.25
        yfirst    = -49.875
        yinc      = 0.25

    else:
        print("grid '" + str(param[1]) + "' not defined")

    # Open grid template file for writing
    gridfile_name = OUTDIR + "/gridfile"

    try:
        with open(gridfile_name, 'w+', encoding="utf8") as gridfile:
            gridfile.write("gridtype  = "+gridtype+"\n")
            gridfile.write("xsize     = "+str(xsize)+"\n")
            gridfile.write("ysize     = "+str(ysize)+"\n")
            gridfile.write("xfirst    = "+str(xfirst)+"\n")
            gridfile.write("xinc      = "+str(xinc)+"\n")
            gridfile.write("yfirst    = "+str(yfirst)+"\n")
            gridfile.write("yinc      = "+str(yinc)+"\n")

    except (IOError, OSError):
        print("error opening gridfile "+gridfile_name+"\n")
    return gridfile_name


if __name__ == "__main__":
    freeze_support()
    if len(sys.argv) != 2:
        raise ValueError('Please input the path to a json file\n')
    process_input(sys.argv[1])
