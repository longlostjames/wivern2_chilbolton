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
import fnmatch

import getpass, socket
import cftime
import netCDF4 as nc4

import itertools
import operator

 
    
def trim_ncfile(infile,outfile,ray_start,ray_end):
    file_dtstr = infile.split("_")[1]

    print(infile);
    print(file_dtstr);
    DS = nc4.Dataset(infile);

    optstr = "-d time,{},{}".format(ray_start,ray_end);
    nco.ncks(input=infile, output=outfile, options=[optstr]);
    
    DS.close()
    
try:
    opts, args = getopt.getopt(sys.argv[1:], "i:o:s:e:", ["in=","out=","start=","end="])
except getopt.GetoptError as err:
    # print help information and exit:                                                                              
    print(err) # will print something like "option -a not recognized                                                
    sys.exit(2)

for o, a in opts:
    if o == "-i":
        infile = a;
    elif o == "-o":
        outfile = a;
    elif o == "-s":
        ray_start=int(a);
    elif o == "-e":
        ray_end=int(a)
    else:
        assert False, "unhandled option"

trim_ncfile(infile,outfile,ray_start,ray_end);
    

