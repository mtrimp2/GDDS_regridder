import sys
sys.dont_write_bytecode = True

EXT_N4 = ["nc", "nc3", "nc4", "nc4c", "n4", "NC", "NC4"] 
EXT_HD5 = ["hdf" , "h4" , "hdf4" , "he2" , "h5" , "hdf5" , "he5"] #can also be opened with xarray, same way as nc4 I believe? I'm gonna look into it I might be wrong

def pregrid_generate(config):
    "Generate grid file and cdo arguments"
    # Error string
    param = config["regridding"]
    OUTDIR = config["TMPDIR"]
    gridtype  = 'lonlat'

    # Set grid template values
    if (param[1] == 'JRA-55'):
        xsize     = 288
        ysize     = 145
        xfirst    = 0
        xinc      = 1.25
        yfirst    = 90
        yinc      = -1.25
    elif (param[1] == '20cr2x2'):
        # NOAA (ESRL) Version
        xsize     = 180
        ysize     = 91
        xfirst    = 0
        xinc      = 2
        yfirst    = 90
        yinc      = -2
    elif (param[1] == 'MERRA0.5'):
        # DISC Version
        xsize     = 540
        ysize     = 361
        xfirst    = -180
        xinc      = 0.666666666666
        yfirst    = -90
        yinc      = 0.5
    elif (param[1] == 'MERRA1.25'):
        # DISC Version
        xsize     = 288
        ysize     = 144
        xfirst    = -179.375
        xinc      = 1.25
        yfirst    = -89.375
        yinc      = 1.25
    elif (param[1] == 'MERRA2'):
        # DISC Version
        xsize     = 576
        ysize     = 361
        xfirst    = -180
        xinc      = 0.625
        yfirst    = -90
        yinc      = 0.5
    elif (param[1] == 'gpcp2.5'):
        xsize     = 144
        ysize     = 72
        xfirst    = 1.25
        xinc      = 2.5
        yfirst    = -88.75
        yinc      = 2.5   
    elif (param[1] == 'cfsr0.5a'):
        xsize     = 720
        ysize     = 361
        xfirst    = 0
        xinc      = 0.5
        yfirst    = 90
        yinc      = -0.5
    elif (param[1] == 'cfsr0.5b'):
        xsize     = 720
        ysize     = 360
        xfirst    = 0.25
        xinc      = 0.5
        yfirst    = 89.75
        yinc      = -0.5
    elif (param[1] == 'cfsr1.0'):
        xsize     = 360
        ysize     = 181
        xfirst    = 0
        xinc      = 1.0
        yfirst    = 90
        yinc      = -1.0
    elif (param[1] == 'cfsr2.5'):
        xsize     = 144
        ysize     = 73
        xfirst    = 0
        xinc      = 2.5
        yfirst    = 90
        yinc      = -2.5
    elif (param[1] == 'ncepncar2.5'):
        # NOAA (ESRL) Version
        xsize     = 144
        ysize     = 73
        xfirst    = 0
        xinc      = 2.5
        yfirst    = 90
        yinc      = -2.5
    elif (param[1] == 'geos1x125'):
        # From http://acmg.seas.harvard.edu/geos/doc/man/appendix_2.html#GMAO_4x5
        xsize     = 288
        ysize     = 181
        xfirst    = -180
        xinc      = 1.25
        yfirst    = -90
        yinc      = 1.0
    elif (param[1] == 'geos1x1'):
        # From http://acmg.seas.harvard.edu/geos/doc/man/appendix_2.html#GMAO_4x5
        
        xsize     = 360
        ysize     = 181
        xfirst    = -180
        xinc      = 1.0
        yfirst    = -90
        yinc      = 1.0
    elif (param[1] == 'geos4x5'):
        # From http://acmg.seas.harvard.edu/geos/doc/man/appendix_2.html#GMAO_4x5
        
        xsize     = 72
        ysize     = 46
        xfirst    = -180
        xinc      = 5
        yfirst    = -90
        yinc      = 4
    elif (param[1] == 'geos2x25'):
        # From http://acmg.seas.harvard.edu/geos/doc/man/appendix_2.html#GMAO_4x5
        
        xsize     = 144
        ysize     = 91
        xfirst    = -180
        xinc      = 2.5
        yfirst    = -90
        yinc      = 2.0
    elif (param[1] == 'geos0.25'):
        # From http://acmg.seas.harvard.edu/geos/doc/man/appendix_2.html#GMAO_4x5
        xsize     = 1152
        ysize     = 721
        xfirst    = -180
        xinc      = 0.3125
        yfirst    = -90
        yinc      = 0.25
    elif (param[1] == 'geos0.5'):
        # From http://acmg.seas.harvard.edu/geos/doc/man/appendix_2.html#GMAO_4x5
        xsize     = 540
        ysize     = 361
        xfirst    = -180
        xinc      = 0.666666667
        yfirst    = -90
        yinc      = 0.5
    elif (param[1] == 'fv1x125'):
        # Extrapolated from geos1x125
        xsize     = 288
        ysize     = 181
        xfirst    = 0
        xinc      = 1.25
        yfirst    = -90
        yinc      = 1.0
    elif (param[1] == 'fv2x25'):
        # Extrapolated from geos2x25
        xsize     = 144
        ysize     = 91
        xfirst    = 0
        xinc      = 2.5
        yfirst    = -90
        yinc      = 2.0
    elif (param[1] == 'fv4x5'):
        # Extrapolated from geos4x5
        xsize     = 72
        ysize     = 46
        xfirst    = 0
        xinc      = 5.0
        yfirst    = -90
        yinc      = 4.0
    elif (param[1] == 'ERA-40'):
        # ECMWF Version
        xsize     = 144
        ysize     = 73
        xfirst    = 0
        xinc      = 2.5
        yfirst    = 90
        yinc      = -2.5
    elif (param[1] == 'ERA2.5'):
        # ECMWF Version
        xsize     = 144
        ysize     = 73
        xfirst    = 0
        xinc      = 2.5
        yfirst    = 90
        yinc      = -2.5
    elif (param[1] == 'ERA-I'):
        # ECMWF Version
        xsize     = 240
        ysize     = 121
        xfirst    = 0
        xinc      = 1.5
        yfirst    = 90
        yinc      = -1.5
    elif (param[1] == 'ERA1.5'):
        # ECMWF Version
        xsize     = 240
        ysize     = 121
        xfirst    = 0
        xinc      = 1.5
        yfirst    = 90
        yinc      = -1.5
    elif (param[1] == 'ERA.75'):
        # ECMWF Version
        xsize     = 480
        ysize     = 241
        xfirst    = 0
        xinc      = 0.75
        yfirst    = 90
        yinc      = -0.75
    elif (param[1] == 'GPCC2.5'):
        # NOAA Version
        xsize     = 144
        ysize     = 72
        xfirst    = 1.25
        xinc      = 2.5
        yfirst    = 88.75
        yinc      = -2.5
    elif (param[1] == 'GPCC1.0'):
        # NOAA Version
        xsize     = 360
        ysize     = 180
        xfirst    = 0.5
        xinc      = 1.0
        yfirst    = 89.5
        yinc      = -1.0
    elif (param[1] == 'GPCC0.5'):
        # NOAA Version
        xsize     = 720
        ysize     = 360
        xfirst    = 0.25
        xinc      = 0.5
        yfirst    = 89.75
        yinc      = -0.5
    elif (param[1] == 'CMORPH0.25'): 
        xsize     = 1440
        ysize     = 480
        xfirst    = -180
        xinc      = 0.25
        yfirst    = -59.875
        yinc      = 0.25
    elif (param[1] == 'NLDAS-2'):
        # DISC Version
        xsize     = 464
        ysize     = 224
        xfirst    = -124.9375
        xinc      = 0.125
        yfirst    = 25.0625
        yinc      = 0.125
    elif (param[1] == 'GLDAS-2_0.25'):
        # DISC Version
        xsize     = 1440
        ysize     = 600
        xfirst    = -179.875
        xinc      = 0.25
        yfirst    = -59.875
        yinc      = 0.25
    elif (param[1] == 'GLDAS-2_1'):
        # DISC Version
        xsize     = 360
        ysize     = 150
        xfirst    = -179.5
        xinc      = 1.0
        yfirst    = -59.5
        yinc      = 1.0
    elif (param[1] == 'NCA-LDAS'):
        # DISC Version
        xsize     = 464 
        ysize     = 224
        xfirst    = -124.9375
        xinc      = 0.125
        yfirst    = 25.0625
        yinc      = 0.125
    elif (param[1] == 'GFDL'):
        xsize     = 288
        ysize     = 180
        xfirst    = 0.625
        xinc      = 1.25
        yfirst    = -89.5
        yinc      = 1.
    elif (param[1] == 'gpcp3'):
        xsize     = 720
        ysize     = 360
        xfirst    = -179.75
        xinc      = 0.5
        yfirst    = 89.75
        yinc      = -0.5
    elif (param[1] == 'TMPA'):
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
    cmd = ""

    try:
        gridfile = open(gridfile_name,'w+')
        gridfile.write("gridtype  = "+gridtype+"\n")
        gridfile.write("xsize     = "+str(xsize)+"\n")
        gridfile.write("ysize     = "+str(ysize)+"\n")
        gridfile.write("xfirst    = "+str(xfirst)+"\n")
        gridfile.write("xinc      = "+str(xinc)+"\n")
        gridfile.write("yfirst    = "+str(yfirst)+"\n")
        gridfile.write("yinc      = "+str(yinc)+"\n")
        gridfile.close()
        
        cmd = "-"+str(param[0])+","+gridfile_name 
    except (IOError, OSError) as err:
        print("error opening gridfile "+gridfile_name+"\n")
    return gridfile_name