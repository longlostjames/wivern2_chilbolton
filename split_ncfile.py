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

def split_ncfile(filename,outpath,nray):
    file_dtstr = filename.split("_")[1]

    print(filename);
    print(file_dtstr);
    DS = nc4.Dataset(filename);

    nraytotal = int(DS['time'][:].shape[0]);
    print(nraytotal);

    
    for i in np.arange(0,nraytotal,nray,dtype=int):
        dtime  = cftime.num2date(DS['time'][i],DS['time'].units)
        dtime0 = dtime.replace(microsecond=0);
        imin = i;
        imax = np.min([imin+nray-1,nraytotal-1]);
        dtime  = cftime.num2date(DS['time'][imax],DS['time'].units)
        dtime1 = dtime.replace(microsecond=0)+datetime.timedelta(seconds=1)
        print(dtime0,dtime1)
        dtstr0 = cftime.datetime.strftime(dtime0,'%Y%m%d%H%M%S')
        dtstr1 = cftime.datetime.strftime(dtime1,'%Y%m%d%H%M%S')
        print(dtstr0)
        outfile = filename.replace(file_dtstr,dtstr0+'-'+dtstr1)
        outfile = outfile.replace('nc','nc4');
        outfile = os.path.join(outpath,outfile);
        print(outfile);
        
        optstr = "-4 -L 1 -d time,{},{}".format(imin,imax);
        nco.ncks(input=filename, output=outfile, options=[optstr]);
    DS.close()
    
try:
    opts, args = getopt.getopt(sys.argv[1:], "f:p:n:", ["file=","path=","nray="])
except getopt.GetoptError as err:
    # print help information and exit:                                                                              
    print(err) # will print something like "option -a not recognized                                                
    sys.exit(2)

for o, a in opts:
    if o == "-f":
        infile = a;
    elif o == "-p":
        outpath = a;
    elif o == "-n":
        nray = int(a);
    else:
        assert False, "unhandled option"

split_ncfile(infile,outpath,nray);
    

