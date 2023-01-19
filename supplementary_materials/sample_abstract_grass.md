The NASA LDAS projects, including the North American Land Data Assimilation System (NLDAS), Global Land Data Assimilation System (GLDAS), Famine Early Warning Systems Network (FEWS NET) Land Data Assimilation System (FLDAS), and National Climate Assessment - Land Data Assimilation System (NCA-LDAS), are instances of the Land Information System (LIS) land surface modeling software that have been configured for specific domains and purposes. These LDAS datasets are archived by the NASA Goddard Earth Sciences (GES) Data and Information Services Center (DISC) and distributed to worldwide users for water resources applications, drought and wetness monitoring, numerical weather prediction studies, water and energy cycle scientific investigations, and interpretation of satellite and ground-based observations.

# GRaSS Abstract

The goal of this research project has been to develop a subsetter and regridder to replace the current tool in use, LEARS (Level 3/4 Regridder and Subsetter). Though working, LEARS has limited functionality and is far from optimized; both factors that prompted this research. 

The tool that has been developed over the course of this project has been named GRaSS (the GESDISC Regridding and Subsetting System). This tool was developed in Python, using the Xarray package for netCDF dataset processing, and all forms of subsetting. xESMF was used for regridding, with bilinear, conservative, nearest neighbor and patch interpolations supported. The data used for experimentation were MERRA-2 netCDF files (mention exact ones?) obtained from the GESDISC website. 

At this point in time, GRaSS can perform all operations that LEARS can (spatial, temporal, and variable subsetting and regridding), and produces identical output. Additionally, through time and memory resource allocation monitors, GRaSS has proven to be more memory-efficient than LEARS on all operations, and more time efficient on many.

Special features of GRaSS include the ability to subset spatially across the anti-meridian, or International Date Line. This was not originally a capability of Xarray, but GRaSS is able to divide a subset on both sides of the anti-meridian and concatenate them, producing the expected result. Additionally, GRaSS is able to temporally subset not only by range(ie. 4-8), but also by selection of individual indicies (ie. 1,3,5), and a combination of the two (ie. 1-3, 5, 7-9, 11).

Though working, GRaSS is not fully ready for implementation. GRaSS still has limitations on global regridding- Any regridded plots must be spacially subsetted as well in order to produce accurate output. Additionally, vertical subsetting has not yet been implemented. 


**NASA GES DISC provides various data services along with the data:**

* GES DISC Subsetter provides spatial subsetting and variable subsetting to save users both the bandwidth of the network and storage space. 
* If users like to compare LDAS data to other products with different latitude/longitude grids, GES DISC Subsetter can regrid the datasets with 7 interpolation methods and more than 30 grid options.
* Some LDAS products, such as NLDAS_FOR0125_H, are in the GRIB format. GES DISC Subsetter can convert them to the netCDF format. Conversely, products such as GLDAS_NOAH025, which are in netCDF, can be converted to GRIB.
* To facilitate GIS and other uses, the GeoTIFF and Cloud Optimized GeoTIFF (COG) output formats have been implemented and enabled for all the LDAS products. 
* “Data Rods for Hydrology” service provides users with fast access to long time series for single locations (for selected variables).
* Giovanni is an online (Web) environment for the display and analysis of geophysical parameters in which the provenance (data lineage) can easily be accessed.
* The GrADS Data Server (GDS) is a stable, secure data server that provides subsetting and analysis services across the internet. 
* The THREDDS Data Server (TDS) is a web server that provides metadata and data access for scientific datasets, using OPeNDAP, OGC WMS and WCS, HTTP, and other remote data access protocols.
