'''
lears-mini
'''
#!/Users/jsu3/miniconda3/bin/python3.9 -u

import sys
import json
import time
import os
import glob
import shutil
from lears import params
from lears import netCDF_handler

def main():
    '''
    main
    '''
    with open(sys.argv[1]) as f_conf:
        config = json.load(f_conf)
    config["TMPDIR"] = params.TMPDIR + "t" + str(int(time.time())) + "_" + str(os.getpid())
    os.makedirs(config["TMPDIR"], mode=0o777)
    if config["regridding"]:
        config["gridfile"] = params.pregrid_generate(config)

    config["pwd"] = os.getcwd()

    for file in glob.glob(config['in_dir'], recursive=True):
        ext = file.split(".")[-1]
        config["ifile"] = file
        if ext in params.EXT_N4:
            # if config["debug"]:
            #     json_formatted_str = json.dumps(config, indent=4)
            #     print(json_formatted_str)
            netCDF_handler.netCDF_handler(config)

    if not config["debug"]:
        shutil.rmtree(config["TMPDIR"])

if __name__ == "__main__":
    main()
