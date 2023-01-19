# pylint: disable=W0105
#^^ Pointless String Statements

'''
Compare netCDF4 outputs from 2 tools (eg, lears-mini and grass)
with similar settings as our fibonnaci comparison

TEA: Tool Equivalence Analyzer
Previously: LRS-10
Updated 02-23-2022
Developer: Maggie Trimpin
'''

import os
import json
import sys
import tracemalloc
import time

def open_jsons(path):
    '''
    function: open_jsons
    inputs: path (String)
    Checks if the input path is one json file, or if it is a directory of json files, 
    opens all subdirectory jsons
    Returns a list of all json file paths for analysis
    '''
    json_list = []
    if os.path.isfile(path):#if the path is a file 
        if path.endswith(".json"):
            json_list.append(path)
    else: #path is a directory path
        for root, dirs, files in os.walk(path):
            for file in files:
                #append the file name to the list if it is a json
                if file.endswith(".json"):
                    json_list.append(os.path.join(root,file))
    return json_list


def brew_tea(json_filepath, tool1, tool2):
    '''
    function: brew_tea
    inputs: json_filepath (String), tool1 and tool2 (Python files)
    Description: Extracts json information about netcdf file for analysis,
    debug status, output key, etc.
    Then, checks if the output files for each tool exist.
    If the output file does not exist, then TEA runs each tool on the given nc4 file,
    and then compares the outputs. If the output files already exist,
    it simply compares the existing files.

    Once all processes are complete, if debug == false,
    TEA deletes any files created during runtime.
    '''
    file = open(json_filepath)
    data = json.load(file)
    debug = data["debug"]
    out_key = data['out_key']
    nc_filepath = data['in_dir']

    tool1_name = os.path.basename(os.path.splitext(tool1)[0])
    tool2_name = os.path.basename(os.path.splitext(tool2)[0])

    output_path1 = check_file_path(tool1_name, out_key, nc_filepath)
    output_path2 = check_file_path(tool2_name, out_key, nc_filepath)

    created_files = []

    #check if path to output file 1 exists
    if os.path.isfile(output_path1):
        nc_path1 = output_path1
    
    else: #If output file doesn't already exist, run process
        tracemalloc.start()
        start_time = time.time()

        run_tool(json_filepath, tool1)
        nc_path1 = output_path1
        created_files.append(nc_path1)

        exec_time = (time.time() - start_time)
        exec_mem = tracemalloc.get_traced_memory()[1]

        if debug:
            print('Time taken for tool 1 execution:', round( exec_time, 3 ), ' seconds')
            print('Memory used for tool 1 execution: ', round( exec_mem, 3 ), ' KiB')
        tracemalloc.stop()

    #check if path to output file 2 exists
    if os.path.isfile(output_path2):
        nc_path2 = output_path2

    else: #If output file doesn't already exist, run process
        tracemalloc.start()
        start_time = time.time()
        
        run_tool(json_filepath, tool2)
        nc_path2 = output_path2
        created_files.append(nc_path2)
        
        exec_time = (time.time() - start_time)
        exec_mem = tracemalloc.get_traced_memory()[1]

        if debug:
            print('Time taken for tool 2 execution:', round( exec_time, 3 ), ' seconds')
            print('Memory used for tool 2 execution: ', round( exec_mem, 3 ), ' KiB')
        tracemalloc.stop()

    #compare output files for the tools
    same = compare_tools(nc_path1, nc_path2, tool1_name, tool2_name, json_filepath, out_key)

    #After analysis/comparison, delete any files created if debug==false
    if not debug:
        if created_files != []:
            #delete them
            for file in created_files:
                os.remove(file)
    return same


def check_file_path(tool_name, out_key, nc_filepath):
    '''
    function: check_file_path
    inputs: tool_name (String), out_key (String), json_filepath (String), nc_filepath (String)
    Description: Checks if the output path for each particular tool/file combination exist.
    Returns the generated output path, so that the program can compare the output files
    '''

    path = "../../../tmp"
    nc_filename = os.path.basename(nc_filepath)

    if out_key == "": #if empty, output file name = json file name + in_dir
        output_path = os.path.join(path, f"{tool_name}_{os.path.splitext(nc_filename)[0]}")
    else: # otherwise the output file name = json file name + in_dir + outkey + .nc4
        output_path = os.path.join(path, f"{tool_name}_{os.path.splitext(nc_filename)[0]}_{out_key}.nc4")
    return output_path

def run_tool(json_filepath, tool):
    '''
    function: run_tool
    inputs: json_filepath (String) and tool (Python file path)
    Description: Runs the specified tool on the input from the json file,
    outputs new nc4 file to the path specified within the tool
    '''
    os.system(f'python {tool} {json_filepath}')
    #tool outputs new nc4 file to a specific path

def compare_tools(nc_path1, nc_path2, tool1_name, tool2_name, json_filepath, out_key):
    '''
    function: compare_tools
    inputs: nc_path1 and nc_path2 (Strings), tool1_name and tool2_name (Strings)
    Description: Calls nccmp, outputs if the nc output files from each tool are identical.
    If not, outputs statistics.

    nccmp flag descriptions:
        -d              : Compare data (variable values)
        -f              : Forcefully compare, do not stop after first difference.
        -S              : Report data difference statistics (file1 - file2).
        -q              : Do not print metadata or data differences.
        -N              : Allows NaN values to be equal.
        --tolerance=TOL : Compare if absolute difference > TOL.
        -M              : Ignore difference between values that have different missing_value/_FillValue.

    '''

    json_filename = os.path.basename(json_filepath)
    difference_filepath = f"{tool1_name}_{tool2_name}_{json_filename}_comparison_{out_key}"

    os.system(f'nccmp -d -f -Sq -N --tolerance=.0001 -M {nc_path1} {nc_path2} > ../../../tmp/{difference_filepath}.txt')
    #^^^^optionally, ignore variables with the NaN discrepancies? 
    filesize = os.path.getsize(f'../../../tmp/{difference_filepath}.txt')

    if filesize == 0:
        print("The files are identical; No output file.")
        #delete file
        os.remove(f'../../../tmp/{difference_filepath}.txt')
        same = True
    else:
        print(f"Output differences to {difference_filepath}.txt")
        same = False
    return same

#Pass as: (Json_file_path, tool1_path, tool2_path)
if __name__ == "__main__":
    #Captures input json file and tools from user, and runs TEA.
    if len(sys.argv) != 4:
        raise ValueError('Please input json file, Control tool, and at least one additional Tool for comparison')
    jsons_for_comparison = open_jsons(sys.argv[1])
    identical = []
    different = []
    for file in jsons_for_comparison:
        same = brew_tea(file, sys.argv[2], sys.argv[3])
        if same:
            identical.append(file)
        else:
            different.append(file)
    
    print()
    print(len(jsons_for_comparison), " json files processed.")
    print()
    print("Jsons producing identical output from tools: ")
    for file in identical:
        print(os.path.basename(file))
    print()
    print("Jsons producing different output from tools: ")
    for file in different:
        print(os.path.basename(file))

#python TEA.py ../../test/global_regrid.json ../lears-mini/lears-mini.py ../grass/grass.py

#chmod +rwx tmp
