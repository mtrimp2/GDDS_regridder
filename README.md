# User guide for the GDDS_regridder repository

## Directories:
- **data**: contains datasets on which to test the regridders.
- **grass**: contains grass tool, usees xarray and xesmf to perform a variety of subsetting and regridding operations.
- **lears-mini**: contains minimized version of the L34RS regridder, which takes json input.
- **regridder**: contains regridder tool and sample json input file.
- **supplementary_materials**: contains reports/presentations and various markdown files with observations and research on various tools
- **tea**: contains TEA comparison tool, which takes one json file and two python tools to compare as input. Also contains a time/memory usage comparison tool, to evaluate tool efficiency.

GRaSS tool and documentation can be found at [This repo](https://git.earthdata.nasa.gov/projects/XDS/repos/xds-1/browse), you all should have access.