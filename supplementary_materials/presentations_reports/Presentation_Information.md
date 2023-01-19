# Internship Final Report
Updated: 12/09/2021\
Author: Maggie Trimpin

---

# GRaSS: The New Subsetter/Regridder For GES DISC

## Abstract

The goal of this research project has been to develop a subsetter and regridder to replace the current tool in use, LEARS (Level 3/4 Regridder and Subsetter). Though working, LEARS has limited functionality and is far from optimized; both factors that prompted this research. 

The tool that has been developed over the course of this project has been named GRaSS (the GESDISC Regridding and Subsetting System). This tool was developed in Python, using the Xarray package for netCDF dataset processing, and all forms of subsetting. xESMF was used for regridding, with bilinear, conservative, nearest neighbor and patch interpolations supported. The data used for experimentation were MERRA-2 netCDF files (mention exact ones?) obtained from the GESDISC website. 

At this point in time, GRaSS can perform all operations that LEARS can (spatial, temporal, and variable subsetting and regridding), and produces identical output. Additionally, through time and memory resource allocation monitors, GRaSS has proven to be more memory-efficient than LEARS on all operations, and more time efficient on many.

Special features of GRaSS include the ability to subset spatially across the anti-meridian, or International Date Line. This was not originally a capability of Xarray, but GRaSS is able to divide a subset on both sides of the anti-meridian and concatenate them, producing the expected result. Additionally, GRaSS is able to temporally subset not only by range(ie. 4-8), but also by selection of individual indicies (ie. 1,3,5), and a combination of the two (ie. 1-3, 5, 7-9, 11).

Though working, GRaSS is not fully ready for implementation. GRaSS still has limitations on global regridding- Any regridded plots must be spacially subsetted as well in order to produce accurate output. Additionally, vertical subsetting has not yet been implemented. 

## Introduction/ Background

The NASA Goddard Earth Sciences Data and Information Services Center (GES DISC) Level 3/4 Regridder and Subsetter (LEARS) has been the primary subset and regrid tool used for Earth Satellite Data since early 2018. It is one of the most popular HTTP services at the center in terms of number of files and volume of data moved. This tool provides the ability to subset spatially, temporally, vertically and by variable, and regrid based on a variety of different interpolations and resolutions. These processes are executed using the CDO Python package and command line tool. However, the LEARS tool is still limited in some functionality and is far from optimized. Thus, this project aims to develop a new, more optimal subsetter and regridder, for eventual replacement of the current LEARS tool.

## Methods/approach
The methods section should discuss original data collection procedures
as well as where the data is located and how the data was analyzed. Be sure to reference
the data sources in the Reference and Acknowledgement section.

--- 

GRaSS was developed using Xarray, as well as netCDF data from the GES DISC website for testing purposes.

### Why xarray?

When compared to other data processing and subsetting tools, xarray displayed far more functionality and customizability.

xarray benefits:
* Apply operations over dimensions by name: x.sum('time').
* Select values by label (or logical location) instead of integer location: x.loc['2014-01-01'] or x.sel(time='2014-01-01').
* Mathematical operations (e.g., x - y) vectorize across multiple dimensions (array broadcasting) based on dimension names, not shape.
* Easily use the split-apply-combine paradigm with groupby: x.groupby('time.dayofyear').mean().
* Database-like alignment based on coordinate labels that smoothly handles missing values: x, y = xr.align(x, y, join='outer').
* Keep track of arbitrary metadata in the form of a Python dictionary: x.attrs.
* The N-dimensional nature of xarray’s data structures makes it suitable for dealing with multi-dimensional scientific data,

Additionally, xarray datasets are comparable with xESMF, a fast, easy regridding package that uses ESMF/ESMPy as backend

xESMF benefits:
* Increased usability and simplicity
* Increased grid types (regridding ability for rectangular, curvilinear, quadrilateral)
* Increased functionality such as tracking metadata
* Compatible with xarray datasets, as well as numpy ndarrays

GRaSS takes input in  the form of a json file containing the following information:

```
{ 
	"in_dir": input file, if the path is a directory such as ../data/*.nc4, then deal with all the nc4 files,
	"debug": true or false. when it's false, dont output anything to display,
	"subset_spatial":[lon1, lon2, lat1, lat2],
	"subset_temporal": [Range (ie. "1/5"), individual indicies (ie. 1,2,4,6,7), or combination]
    "subset_variable": [One or more variables],
	"subset_vertical":[], #Yet to be implemented
	"regridding":[regridding interpolation, target grid resolution],
	"out_key": ""  when it = '' (empty), then the output file name = json file name + in_dir, otherwise the output file name = in_dir + outkey + .nc4
}
```

The json file is passed to GRaSS as a sole parameter, where each key value is analyzed. 

As shown in the json input file above, the "in_dir" is where the user should supply the path to the netcdf file(s) for processing. If the "in_dir" field is a directory path rather than a file path, the script runs the following processes on each .nc4 file in the directory. Otherwise, it simply performs the operation on the one input file. Xarray also includes support for OPeNDAP (via the netCDF4 library or Pydap), which lets us access datasets over HTTP. Thus, the "in_dir" field can also be a link to a remote dataset, which the code can open and process accordingly. 

Next, the "debug" field can be set to True or False, based on whether the user would like the program to output diagnostic updates to the compiler. Additionally, if "debug" is set to false, the program will delete any files created during processing upon completion (not including the processed dataset that is output). Otherwise, any files created during processing will remain. 

For the five fields reserved for subsetting and regridding, the program checks the values for each key. If a given key is not empty, a method is called to perform the subsetting or regridding operation that corresponds to the key title and value(s).  

GRaSS checks for each key in the following order:
Variable subsetting, regridding, spatial] subsetting, temporal subsetting, vertical subsetting (not yet implemented).

---

GRaSS reads the json file in the main function, simply called grass(). This function opens the json file passed to the function, and uses xarray to open the dataset(s) to be processed. It then passes the dataset to a function called perform_operation(), where the json file values are read, and the corresponding operations are performed on the dataset.

The first task performed by the function is to read the dataset units, identify the variables corresponding to "degrees_east" and "degrees_north", and change the value names to "lon" and "lat", accordingly. Many datasets have inconsistently named coordinate values (ie. X and Y, latitude and longitude, etc.), so the program must name them uniformly before performing any operations.

Then, the script sets "debug" to the value under data['debug'] from the json file, and reads each of the fields for subsetting and regridding. For each field, the program simply uses if statements to evaluate "if [field] != []", or if the key is not empty. If not, the program calls the function corresponding to the value. For example, 

```
if data['subset_spatial'] != []:
        new_ds = spatial_subset(working_ds, data)
        working_ds = new_ds.copy(deep=True)
```
In plain English, this code means that if the field for spatial subsetting is not empty, then the program should call the spatial subset function, pass the working dataset and json information, and then update the working dataset with the new subsetted dataset.

Once all of the operations are complete, a new nc4 file is outputted to the folder of the original input file, with a name identifiable by the name of the tool, the input json, the nc4 file for processing, and the out_key provided by the json key. If an out_key is not provided, then the output file name = json file name + input file name. Example: GRaSS_regridding_MERRA2_400.tavg1_2d_slv_Nx.20200101_key.nc4.

An additional parameter, "NaN_to_missing_val", is used to determine whether or not to convert any NaN values within a dataset to the "missing_value" specified in the attributes of each dataset. With a default value of "True," this parameter promts the system to read from the dataset attributes, obtain the default "missing value," and convert all NaN values within the dataset to that value. This functionality comes in handy when we use our comparison tool, TEA, to evaluate the results of multiple different tools' output files of the same dataset, to ensure they have their missing values set to the same number. 

## Results

In order to evaluate the accuracy and performance of the new tool, an additional program was created to compare the output files of GRaSS against those of LEARS. This comparison script is called TEA, the Tool Equivalence Analyzer.

### TEA Overview

TEA's purpose is to compare files subsetted and regridded by GRaSS to those subsetted and regridded by LEARS, using the same input values. The files will be compared using the [NCO](http://nco.sourceforge.net/) tool nccmp, in our TEA. 

TEA reads the input json file and searches for existing output files from each tool for that file. If the files do not exist, then TEA runs each tool on the input file, where an output netcdf file is generated. TEA then runs the ncdiff operation on the generated output files in the command line. If the output files already exist, TEA simply compares the existing fies. 

If the files are identical, the program reports as such. If differences in the comparison files are detected, they are outputted to a text file in a temp directory, with a name corresponding to the names of the two tools and json file under evaluation. The program will display whether the files are identical, the path to the output text file, and the memory and time report from runtime (only if TEA must run either of the tools; it is not necessary if TEA compares existing output files).

Upon comparison of the completed GRaSS tool with the currently implemented LEARS tool, GRaSS was observed to produce identical output files to those of LEARS. All functions were tested by TEA, with many different input varieties for each processs (ie. Spatial subsetting in many different quadrants, temporal subsetting using integers or slices, regridding using many different interpolations and resolutions). Additionally, a test was run executing GRaSS and LEARS on a json file containing every possible operation (Spatial, temporal, and variable subsetting and regridding), and both tools executed with no problem.


### Efficiency/ Functionality Evaluation against Lears-mini

As a whole, GRaSS provides more functionality than the existing LEARS tool. Thanks to its xarray sctructure, it provides a variety of advanced operations and more customizability than otherwise available with CDO (The library used by LEARS). Aspects that make GRaSS superior to LEARS include:
* Functionality for reading remote files into memory without the need to download them. 
* Ability to subset spatially across the prime meridian. This was not originally a capability of Xarray, but GRaSS is able to divide a subset on both sides of the prime meridian and concatenate them, producing the expected result. 
* While LEARS cannot process json files as input, GRaSS is able to take a json file, or even a directory of json files as an input, and process all of them accordingly. 

In addition to increased functionality, the time and memory resource allocation for each tool was analyzed. A script was developed to run each tool 25 times, and take the average amount of time and memory used for execution. These values were then compared, and the following conclusions were drawn:
* For all operations, GRaSS requires less memory than LEARS, with an average of 25.33 KiB of memory saved per execution
* In the case of spatial and multi-variable subsetting, GRaSS is faster than LEARS, with an average of 5.628 seconds saved per execution
* GRaSS is slower than LEARS in temporal subsetting by indicies and slices. However, in the special case of temporal subsetting using *both* indicies and slices (ie. [1, '3-5', 8, '10-11']), GRaSS is far more efficient- saving an average of 4.246 seconds and 20 KiB. 
* Though GRaSS is slower than LEARS in single-variable subsetting, regridding and some temporal subsetting cases, the average amount of time differene is only 3.266 seconds per execution. That value, in comparison to the amount of time saved with GRaSS for spatial, multi-variable and mixed temporal subsetting (5.628 seconds), poses a far less significant difference. 
* In regards to the conda environment dependencies required for the successful execution of each tool, the minimized LERARS environment has a size of 540M, while GRaSS's environment size is only 539M. A small difference, but a difference nonetheless. 

### Current Limitations

Though working, GRaSS is not fully ready for implementation. GRaSS still has limitations on global regridding- Any regridded plots must be spacially subsetted as well in order to produce accurate output. Additionally, vertical subsetting has not yet been implemented.

## Conclusions

Though not completely ready for deployment at the moment, GRaSS has proven to have strong potential to soon replace the LEARS tool for subsetting and regridding Earth Data for the DISC. GRaSS not only provides more extensive functionality and customizability, due to its Xarray backend. It also proves more efficient in memory and, in some cases, time resource allocation. All in all, GRaSS has proven itself to be more than capable to one day be a reliable and optimal resource for the GES DISC. 