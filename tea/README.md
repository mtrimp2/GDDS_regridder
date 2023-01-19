# User guide for TEA- Tool Equivalence Analyzer 

Author: Maggie Trimpin
Updated: 7/26/2022

## Overview:
Compare netCDF4 outputs from 2 tools (eg, lears-mini and grass or regridder)

## Requirements:
- conda environment (conda env create -f tea_env.yml)
- input json file or directory of files for comparison 
    - Format: 
    ```
    { 
        "in_dir": input file, if the path is a directory such as ../data/*.nc4, then deal with all the nc4 files,
        "debug": true or false. when it's false, dont output anything to display,
        "regridding":[regridding interpolation, target grid resolution],
        "out_key": when it = "" (empty), then the output file name = json file name + in_dir, otherwise the output file name = in_dir + outkey + .nc4
    }
    ```
- input two tool paths (ie ../regridder/regridder.py ../lears-mini/lears-mini.py)

## Workflow

Pass input as: (Json_file_path, tool1_path, tool2_path)- Can take multiple jsons as input, will run both tools on each json file individually, output which produce identical results and which produce different 

Runs regridder on json input, outputs regridded datasets to two files, then uses  nccmp to compare 

At the end, outputs names of each json file that produced identical results, and then names of each json file that produced different results