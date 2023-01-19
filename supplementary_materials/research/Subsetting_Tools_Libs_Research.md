# Subsetting Tools & Libraries Research Documentation
Author: Maggie Trimpin
ver 1.
Updated 9-29-2021

---
## What is it?

Test data subsetting is extracting a smaller sized – referential integer set of data from a ‘production’ database to a non-production environment.

**Subsetting** is the process of retrieving just the parts of large files which are of interest for a specific purpose. This occurs usually in a client—server setting, where the extraction of the parts of interest occurs on the server before the data is sent to the client over a network. The main purpose of subsetting is to save bandwidth on the network and storage space on the client computer.

## Why Subset Data?
* restrict or divide the time range
* select cross sections of data
* select particular kinds of time series
* exclude particular observations

## Methods of Subsetting 
1. **Variable Subsetting** - one or multiple variables in the dataset
2. **Spatial Subsetting** - block out an area based on coordinate lower/upper bounds
3. **Temporal Subsetting** - splicing the data into daily, monthly, hourly units

---

## Available Tools and Libraries
Note: Tools are listed alphabetically. Ordering of the tools does not reflect heirarchy of use/ priority.
+ **CDO** [documentation](https://code.mpimet.mpg.de/projects/cdo/wiki/Cdo#Documentation): 
The Climate Data Operators (CDO) is a large tool set for working on climate and Numerical Weather Prediction (NWP) model data. CDO can also be used to analyse any kind of gridded data not related to climate science.
+ **GDAL** [documentation](https://gdal.org/): 
GDAL is a translator library for raster and vector geospatial data formats that is released under an X/MIT style Open Source License by the Open Source Geospatial Foundation. As a library, it presents a single raster abstract data model and single vector abstract data model to the calling application for all supported formats. It also comes with a variety of useful command line utilities for data translation and processing. 
+ **NetCDF4** [documentation](https://unidata.github.io/netcdf4-python/):
netcdf4-python is a Python interface to the netCDF C library. This module can read, write, and create files in netCDF format. Features of netCDF4 include multiple unlimited dimensions, groups and zlib data compression. Subsetting is possible through netCDF4 for all types (variable, spatial, temporal).
+ **xarray** [documentation](http://xarray.pydata.org/en/stable/):
xarray (formerly xray)  is a Python package tailored to work with netCDF files. It introduces dataset labels in the form of dimensions, coordinates and attributes on top of raw NumPy-like arrays. These dimensional labels and data structures allow for all types of subsetting, and also allow for advanced analytics and visualization.


  ## Tool Comparison Table
| Tool        |   Subsetting Supported  | 
| :------- | :----------: | 
| **CDO**         |   <ul><li>spatial (lat/long)</li><li> temporal</li> <li>variable</li></ul> |
| **GDAL**         |   <ul><li>spatial (lat/long)</li><li> temporal</li> <li>variable</li></ul> |
|   **NCO**      |   <ul><li>spatial (lat/long)</li><li> temporal</li> <li>variable</li></ul>          |
|  **NetCDF4**       | <ul><li>spatial (lat/long)</li><li> temporal</li> <li>variable</li></ul>          |  
|    **xarray**         |  <ul><li>spatial (lat/long)</li><li> temporal</li> <li>variable</li></ul>               |


---
## Vertical Subsetting: Introductory Research

### **What is it?**
(If I'm not mistaken) it's quite similar to spatial subsetting, just vertical? Height-specific bounds?

### **Tool: hdfread- MATLAB Function**

| Parameter        |   Description  | 
| :------- | :----------: | 
| **Fields**         |   String naming the data set field to be read. You can specify only one field name for a Grid data set
| **Box** | Two-element cell array, {longitude,latitude}, specifying the longitude and latitude coordinates that define a region. longitude and latitude are each two-element vectors specifying longitude and latitude coordinates. 
| **Time** | Two-element cell array, [start stop], where start and stop are numbers that specify the start and end-point for a period of time.
| **Vertical** | Two-element cell array, {dimension, range}<ul><li>dimension -- String specifying the name of the data set field to be read from. You can specify only one field name for a Grid data set.</li><li> range -- Two-element array specifying the minimum and maximum range for the subset. If dimension is a dimension name, then range specifies the range of elements to extract. If dimension is a field name, then range specifies the range of values to extract.</li></ul> 'Vertical' subsetting may be used alone or in conjunction with 'Box' or 'Time'. To subset a region along multiple dimensions, vertical subsetting may be used up to eight times in one call to hdfread
