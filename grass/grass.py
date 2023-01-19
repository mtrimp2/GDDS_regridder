#!/usr/bin/env python
# coding: utf-8
# pylint: disable=C0413
#^^Disables import order warnings because i cannot get rid of them
# pylint: disable=W0105
#^^ Disables "Pointless" String Statements Warning

'''
GRaSS: GESDISC Regridding and Subsetting System
Previously: LRS-8/ MEDSAR
Updated 02-23-2022 
Developer: Maggie Trimpin
'''

import warnings
warnings.filterwarnings("ignore")
import sys
import os
import json
import re
import xesmf
import xarray
import numpy as np
import time
import shutil
import stat
import netCDF4 as nc
sys.dont_write_bytecode = True
import params

def grass(json_filepath):
    '''
    function: grass
    inputs: json_filepath (String)
    output: If debug = True in the json file, the function outputs diagnostic information.
    Otherwise, outputs nothing.
    An output file is produced regardless of "debug" status,
    to a directory specified in output_file

    Description: Opens the json file passed to the function,
    and uses xarray to open the dataset(s)
    to be processed.
    It then calls perform_operation to perform subsetting/regridding on the specified dataset.
    '''
    file = open(json_filepath)
    data = json.load(file)

    ####################### Extract info from json ############################
    debug = data['debug'] #if true, output differences. if false, output nothing
    out_key = data['out_key']
    nc_filepath = data['in_dir']
    #list of allowed extensions, ie. nc4, nc, 
    ext = nc_filepath.split(".")[-1]

    ############# Pass nc4 file(s) to method to perform operations ############
    if nc_filepath.startswith("https://"):
        try:
           ds = xarray.open_dataset(nc_filepath+'#mode=bytes')
           opened = True
        except: # catch *all* exceptions
            opened = False
            print("File not found. Exiting")
        
        if opened:
            new_ds = perform_operation(ds, data)
            #output that file
            #Do not cast to float32 for remote files
            output_file(new_ds, nc_filepath, out_key,debug)
            #remember to pass debug to determine what to output


    elif ext in (params.EXT_N4):  #if the file is of type netcdf #from params?
        ds = xarray.open_dataset(nc_filepath)
        data['NaN_to_missing_val'] = False
        if data['NaN_to_missing_val'] == True:
            ds = nan_to_fillval(ds, nc_filepath)
        #run process on the one file
        new_ds = perform_operation(ds, data)
        #output that file
        ds_formatted = new_ds.astype('float32')#Cast for equality to lears datatypes
        output_file(ds_formatted, nc_filepath, out_key,debug)
        #remember to pass debug to determine what to output


    else: #dealing with a directory of .nc4 files
        for file in os.listdir(nc_filepath):
            ext = file.split(".")[-1]
            if ext in params.EXT_N4:
                name = os.path.basename(file)
                print(f"======File: {name} ========")
                file_path = os.path.join(nc_filepath,file)
                ds = xarray.open_dataset(file_path)
                data['NaN_to_missing_val'] = False
                if data['NaN_to_missing_val'] == True:
                    ds = nan_to_fillval(ds, nc_filepath)
                #run process on each file
                new_ds = perform_operation(ds, data)
                #output each file
                ds_formatted = new_ds.astype('float32')#Cast for equality to lears datatypes
                output_file(ds_formatted, nc_filepath, out_key,debug)
                #remember to pass debug to determine what to output
    

############ Perform Desired Operations, return working dataset ###################
def perform_operation(ds, data):
    '''
    function: perform_operation
    inputs: dataset to be operated on (xarray DataSet), data (json file dictonary)
    "data" is input from json file, can be indexed by key name (ie. data['debug'] = True or False)
    output: new xarray dataset, with all specified operations performed

    Description: perform_operation first iterates through the dataset's variable 
    names, and sets the coordinate labels to lat and lon. The method  
    then iterates through each key in the json file,
    and for each operation, if the key value is not empty,
    the function is called to carry out that particular operation.
    (ie. if data['subset_temporal'] is not empty, then subset the dataset
    based on the temporal bounds specified in the json file)
    '''

    #set coordinate names to lat and lon
    long_renamed = False
    lat_renamed = False
    lev_renamed = True
    if len(ds.dims) > 3:# dims = time, lat, lon, lev
        lev_renamed = False
    for name in list(ds.variables):
        try:
            current_var = eval(f"ds.{name}.units")
            if current_var == 'degrees_east':
                ds = eval(f"ds.rename({{'{name}': 'lon'}})")
                long_renamed == True
            elif current_var == 'degrees_north':
                ds = eval(f"ds.rename({{'{name}': 'lat'}})")
                lat_renamed == True
            elif current_var == 'hPa': #?? What else could it be
                ds = eval(f"ds.rename({{'{name}': 'lev'}})")
                lev_renamed == True
        except:
            if lat_renamed and long_renamed and lev_renamed:
                break

    working_ds = ds.transpose("time", "lon", "lat", ...)
    working_ds = ds

    debug = data['debug']
    if data['subset_variable'] != []:
        working_ds = variable_subset(working_ds, data)

    if data['regridding'] != []:
        working_ds = regrid(working_ds, data, debug)

    if data['subset_spatial'] != []:
        working_ds = spatial_subset(working_ds, data)

    if data['subset_temporal'] != []:
        working_ds = temporal_subset(working_ds, data)

    if data['subset_vertical'] != []:
        working_ds = vertical_subset(working_ds, data)
 
    return working_ds

def nan_to_fillval(ds, nc_filepath):
    #extract 'missing value' from ds, change all occurences of that value to NaN
    #nc_ds = nc.Dataset(nc_filepath)
    fillval = nc.default_fillvals['f8']
    #fillval = float(fillval)
    ds.fillna(fillval) 
    for var in list(ds.variables):
        ds[var].attrs['_FillValue'] = fillval
        #print(ds[var].attrs['_FillValue'])
    ds.attrs['_FillValue'] = fillval
    #print(fillval)
    return ds


############ Export the Completed ds in netCDF form ###############
def output_file(ds, nc_filepath, out_key, debug):
    '''
    function: output_file
    inputs: dataset to be outputted (xarray DataSet),
    nc_filepath (String), out_key (String), debug (boolean)
    output: nothing

    Description: Outputs file to a path based on the tool name & version +
    json name + netcdf name + out_key
    Example output file name: grass-ver2_sample_MERRA2_400.tavg1_2d_slv_Nx.20200101.nc4
    '''

    path = "../../../tmp"
    nc_filename = os.path.basename(nc_filepath)

    if out_key == "": #if empty, output file name = json file name + in_dir
        output_path = os.path.join(path, f"grass_{os.path.splitext(nc_filename)[0]}")
        if debug:
            print(f"output to: {output_path}")
    else: # otherwise the output file name = json file name + in_dir + outkey + .nc4
        output_path = os.path.join(path, \
        f"grass_{os.path.splitext(nc_filename)[0]}_{out_key}.nc4")
        if debug:
            print(f"output to: {output_path}")

    if os.path.isfile(output_path):
        #change file permissions so can be overwritten
        os.chmod(output_path, stat.S_IRWXU)
    ds.to_netcdf(output_path)#, encoding={'_FillValue': None})


############## Operation Function Definitions ######################
def variable_subset(ds, data):
    '''
    function: variable_subset
    inputs: dataset to be subsetted (xarray DataSet), data (json file dictionary)

    Variable_subset extracts the string or list of string corresponding to variables 
    to include, supplied by data['subset_variable'].

    Output: new netCDF, containing dataset information subsetted to include only variables passed to function
    '''
    subset_variable = data['subset_variable']
    new_ds = ds[subset_variable]
    return new_ds

def regrid(ds, data, debug):
    '''
    function: regrid
    inputs: dataset to be subsetted (xarray DataSet),
    data (json file dictionary), debug (boolean)

    The regrid function reads and interprets the interpolation method from data['regridding'].
    It then creates a temporary directory, in which it writes a gridfile containing information
    about grid resolution. This gridfile is generated by calling params.py.
    Regrid then reads from the gridfile and regrids the dataset accordingly. 

    Once regridding is complete, if debug == false, the temporary directory is deleted, 
    along with any files created during runtime.

    output: new dataset, regridded according to the
    inputted interpolation and resolution information
    '''
    #create a temporary directory for gridfile
    interpolation = data['regridding'][0]
    #interpolation options = [bilinear, nearest neighbor, inverse distance,
    # spline, binning, spectral and triangulation.]
    #interperet interpolations from params.py
    if interpolation == 'remapbil':
        interpolation = 'bilinear'
    elif interpolation == 'remapdis':
        interpolation = 'inverse distance'
    elif interpolation == 'remapnn':
        interpolation = 'nearest_s2d'
    elif interpolation == 'remapcon':
        interpolation = 'conservative'
    else:
        print(f"Interpolation {interpolation} not supported by GRaSS.")
    
    #generate timestamp for temp_dir
    time_stamp = time.time()
    temporary_dir = f"../../../tmp/grass_temporary_dir_{time_stamp}"
    os.mkdir(temporary_dir)
    data["TMPDIR"] = temporary_dir #temporary output dir created at beginning of perform_operation
    data["gridfile"] = params.pregrid_generate(data) #outputs a path to the gridfile location
    #Read gridfile and extract values
    gridfile = open(data["gridfile"], 'r')
    gridfile_list = gridfile.readlines()
    str_values = "".join(gridfile_list)
    values = re.findall(r'-?\d*\.\d+|\d+|-?\d+', str_values)
    lonll = float(values[2])
    lonul = float(values[2]) + float(values[0]) * float(values[3])
    lon_res = float(values[3])
    latll = float(values[4])
    latul = float(values[4]) + float(values[1]) * float(values[5])
    lat_res = float(values[5])

    if lat_res < 0:
        data["neg_lat_step"] = True
    else:
        data["neg_lat_step"] = False
    if lon_res < 0:
        data["neg_lon_step"] = True
    else:
        data["neg_lon_step"] = False
    

    #re-adjust grid range (ie. -180:180 to 0:360)
    if min(ds['lon']) != lonll:
        if min(ds['lon']) < lonll:
            difference = lonll - min(ds['lon'])
            ds = ds.assign_coords({"lon": (ds.lon + difference)})
        else:
            difference = min(ds['lon']) - lonll
            ds = ds.assign_coords({"lon": (ds.lon + difference)})

    if min(ds['lat']) != latll:
        if min(ds['lat']) < latll:
            difference = latll - min(ds['lat'])
            ds = ds.assign_coords({"lat": (ds.lat + difference)})
        else:
            difference = min(ds['lat']) - latll
            ds = ds.assign_coords({"lat": (ds.lat + difference)})
    
    #OPTIONAL: Masking?
    #must we regrid EVERY variable separately? Let's see if it works ig
    for var in list(ds.data_vars):
        ds["mask"] = xarray.where(~np.isnan(ds[var].isel(time=0)), 1, 0)
        ds_out = xarray.Dataset(
            data_vars=ds.data_vars,
            coords={'lon': (['lon'], np.arange(lonll, lonul, lon_res)),
                            'lat': (['lat'], np.arange(latll, latul, lat_res)),
                            },
            attrs=ds.attrs,
        )
        ds_out["mask"] = xarray.where(~np.isnan(ds_out[var].isel(time=0)), 1, 0)
        #create regridder object with input and target grids
        regridder = xesmf.Regridder(ds, ds_out, interpolation, extrap_method="nearest_s2d", periodic=True)#, unmapped_to_nan = True)#, ignore_degenerate=True)
        #process dataset with regridder
        ds = regridder(ds, keep_attrs=True)

    '''
    #Multi-var regrid
    #Create output grid, using values specified by input file
    output_grid = xarray.Dataset(
    data_vars=ds.data_vars,
    coords={'lon': (['lon'], np.arange(lonll, lonul, lon_res)),
                    'lat': (['lat'], np.arange(latll, latul, lat_res)),
                    },
    attrs=ds.attrs,
    )

    #create regridder object with input and target grids
    regridder = xesmf.Regridder(ds, output_grid, interpolation, extrap_method="nearest_s2d")
    #process dataset with regridder
    ds = regridder(ds, keep_attrs=True)
    '''
    
    if not debug: #delete temporary directory if created
        shutil.rmtree(temporary_dir, ignore_errors=True)

    return ds

def spatial_subset(ds, data):
    '''
    function: spatial_subset
    inputs: dataset to be subsetted (xarray DataSet), data (json file dictionary)
    Spatial_subset reads and interprets the lat/lon boundaries 
    from data['subset_spatial'] and then
    subsets the dataset based on the lat/lon values. 

    If regridding has occurred, longitude upper/lower limits may have changed. Make conversions for spatial subsetting as well

    Since xarray is unable to interpret subsetting across the International Date Line (ie. latLL = 160, latUL = -160),
    if that is the case, two separate datasets are processed, one for each side of the IDL, and concatenated, 
    then returned as one new dataset.

    output: new dataset, subsetted spatially in accordance with input lat/lon values
    '''
    coordinates = data['subset_spatial']
    lonll = int(coordinates[0])
    lonul = int(coordinates[1])
    latll = int(coordinates[2])
    latul = int(coordinates[3])

    if data['regridding'] != []: #if dataset is regridded, we need to pull in grid range info to make it coincide with our dataset
        #begin cases
        #we always assume that the coordinates entered by user are in the format [-180:180, -90:90]
        #assuming regridding has already happened => min and max (ds[lat] and ds[lon]) are the new lower and upper bounds of the dataset
        

        #Apparently not needed? Am I stupid
        '''if min(ds['lon']) > -180:
            difference = min(ds['lon']) + 180
            lonul = lonul + (difference)
            lonll = lonll + (difference)
            lonll = lonll.item(0)
            lonul = lonul.item(0)
        elif min(ds['lon']) < -180:
            difference = min(ds['lon']) - 180
            lonul = lonul + (difference)
            lonll = lonll + (difference)
            lonll = lonll.item(0)
            lonul = lonul.item(0)

        if min(ds['lat']) > -90:
            difference = min(ds['lat']) + 90
            latll = latll + (difference)
            latul = latul + (difference)
            latul = latul.item(0)
            latll = latll.item(0)
        elif min(ds['lat']) < -90:
            difference = min(ds['lat']) - 90
            latll = latll + (difference)
            latul = latul + (difference)
            latul = latul.item(0)
            latll = latll.item(0)'''
    else:
        data["neg_lat_step"] = False
        data["neg_lon_step"] = False
    
    if lonll < lonul: #normal order, doesn't cross anti-meridian
        if not data["neg_lon_step"] and not data["neg_lat_step"]: #if steps increase for both lat and lon, slice normally
            new_ds = ds.sel(time = ds.variables['time'],\
            lat=slice(latll,latul), lon = slice(lonll,lonul))
        elif data["neg_lon_step"]: #if step is negative for lon, slice lonul:lonll
            new_ds = ds.sel(time = ds.variables['time'],\
            lat=slice(latll,latul), lon = slice(lonul,lonll))
        elif data["neg_lat_step"]: #if step is negative for lat, slice latul:latll
            new_ds = ds.sel(time = ds.variables['time'],\
            lat=slice(latul,latll), lon = slice(lonll,lonul))
        #print(new_ds)
        return new_ds

    else: #if crosses antimeridian
        lat_index = slice(latll,latul)
        #split into two, where antimeridian is crossed
        lon_left_of_idl = slice(lonll, 180)
        lon_right_of_idl = slice(-180, lonul)
        lon_index = xarray.concat((ds.lon.sel(lon=lon_right_of_idl),
                                ds.lon.sel(lon=lon_left_of_idl)), dim='lon')
        indexers = {'lon': lon_index, 'lat': lat_index}
        ds = ds.sel(**indexers)
    reindexed_ds = ds.reindex_like(ds.lon)
    new_ds = reindexed_ds.where(((reindexed_ds.lon < lonul).any() and (reindexed_ds.lon > lonll).any()), drop=True)
    return new_ds

def temporal_subset(ds, data):
    '''
    function: temporal_subset
    Functionality includes ability to subset temporal indicies (1,3,5), slices (5-8), or both (1, 3-6, 9)
    inputs: dataset to be subsetted (xarray DataSet), data (json file dictionary)

    Temporal_subset reads the input indicies/slices and appends them to a list of indicies. These indicies are then 
    converted to a dictionary, and xarray then uses them to index the dataset, with each number corresponding to the 
    number of hours after datetime64 00:30:00.00
    
    output: new dataset, subsetted tempotally by selecting only the 
    temporal values included in the time dictionary.
    '''
    subset_temporal = data['subset_temporal']
    indicies = []
    #"subset_temporal":[1,3/5,7,9]
    for item in subset_temporal:
        if isinstance(item, int):
            #if current item is an integer, add it to index array
            indicies.append(item - 1)
        else: #current item is of type "start-end"
            #extract integers on right and left of "-"
            index_range = re.findall('\d+', item)
            start_indx = int(index_range[0])
            stop_indx = int(index_range[1])
            for i in range(start_indx - 1, stop_indx):
                indicies.append(i)
    dictionary = {"time": indicies}
    new_ds = ds.isel(dictionary)
    return new_ds


def vertical_subset(ds, data):
    '''
    function: vertical_subset

    output: new dataset, subsetted tempotally by selecting only the 
    temporal values included in the time dictionary.
    '''
    subset_vertical = data['subset_vertical']
    #"subset_vertical":[1,3/5,7,9]
    indicies = []
    for item in subset_vertical:
        if isinstance(item, int):
            #if current item is an integer, add it to index array
            indicies.append(item - 1)
        else: #current item is of type "start-end"
            #extract integers on right and left of "-"
            index_range = re.findall('\d+', item)
            start_indx = int(index_range[0])
            stop_indx = int(index_range[1])
            for i in range(start_indx - 1, stop_indx):
                indicies.append(i)
    dictionary = {"lev": indicies}
    new_ds = ds.isel(dictionary)

    return new_ds

if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise ValueError('Please input: the file path of a json file containing information about the procedures to be performed.\n')
    grass(sys.argv[1])