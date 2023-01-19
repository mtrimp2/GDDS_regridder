# User guide for GRaSS - GESDISC Regridding and Subsetting System
Author & Developer: Maggie Trimpin
Updated: 07/19/2022

**Note**: This is not the most recent version of GRaSS (found at xds-1 repository noted in GDDS_regridder readme).

## Workflow:
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


### Current Limitations

Though working, GRaSS is not fully ready for implementation. GRaSS still has limitations on global regridding- Any regridded plots must be spacially subsetted as well in order to produce accurate output. Additionally, vertical subsetting has not yet been implemented.

## Use tea to compare output from grass to that of lears-mini

In tea dir:

python tea.py ../grass/sample.json ../lears-mini/lears-mini.py ../grass/grass.py