#!/usr/bin/env python
# coding: utf-8

import getopt, sys

import os
from nco import Nco
nco = Nco()
import pyart
import datetime
import numpy as np
import shutil
import glob

import getpass, socket
import cftime
import netCDF4 as nc4

def trim_ncfile(filename,outpath,ray_start,ray_end):
    file_dtstr = filename.split("_")[1]

    print(filename);
    print(file_dtstr);
    DS = nc4.Dataset(filename);

    outfile = filename.replace(file_dtstr,file_dtstr+'-rev1')
    outfile = os.path.join(outpath,outfile);
    print(outfile);

    optstr = "-d time,{},{}".format(ray_start,ray_end);
    nco.ncks(input=filename, output=outfile, options=[optstr]);
    
    DS.close()
    
try:
    opts, args = getopt.getopt(sys.argv[1:], "f:p:s:e:", ["file=","path=","start=","end="])
except getopt.GetoptError as err:
    # print help information and exit:                                                                              
    print(err) # will print something like "option -a not recognized                                                
    sys.exit(2)

for o, a in opts:
    if o == "-f":
        infile = a;
    elif o == "-p":
        outpath = a;
    elif o == "-s":
        ray_start = int(a);
    elif o == "-e":
        ray_end = int(a);
    else:
        assert False, "unhandled option"

trim_ncfile(infile,outpath,ray_start,ray_end);
    

