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

def trim_camra(filename,outpath):
    DS = nc4.Dataset(filename);
    #bad_rays = np.where(np.all(DS['IRX'][:,:,-1]>=4095,axis=1));
    
    temp = np.squeeze(DS['IRX'][:,6000,-1]);

    good = temp*0;
    good[np.where(temp<4095)] = 1;

    print(good);

    x = np.nditer(good);
    
    good_runs = [[i for i,value in it] for key,it in itertools.groupby(enumerate(x), key=operator.itemgetter(1))  if key==1]
    
    print(good_runs);

    good_runs = [run for run in good_runs if len(run)>1];

    print(good_runs);

    dt = cftime.num2pydate(DS['time'][:],DS['time'].units);

    timesteps = dt[1:]-dt[0:-1];
    timesteps = np.insert(timesteps[0],0,timesteps);

    run_starts = [cftime.num2pydate(DS['time'][run[0]],DS['time'].units) for run in good_runs];

    print(run_starts);
    run_starts = [run - datetime.timedelta(seconds=10) for run in run_starts];
    start_strings = [run_start.strftime('%Y%m%d%H%M%S') for run_start in run_starts];

    DS.close()
    
#    print(run_starts);
    print(start_strings);
    
    irun = 0
 #   revfile = filename.replace("_orig.nc","-rev{:02d}.nc".format(irev))
 #   os.rename(filename,revfile);

    for run in good_runs:

        if run[0]==0:
            revfile = filename.replace("_orig.nc",".nc")                                                      
#        filename = revfile;
#        revfile = filename.replace("-rev{:02d}.nc".format(irev-1),"-rev{:02d}.nc".format(irev))
#        revfile = filename.replace("_orig.nc","-rev{:02d}.nc".format(irev))
        else:
            revfile = 'radar-camra_{}_fix-ts.nc'.format(start_strings[irun]);
        
        ray_start = run[0];
        ray_end   = run[-1];

        print(run[0],run[-1]);

        outfile = os.path.join(outpath,revfile);
        print(outfile);

        trim_ncfile(filename,outfile,ray_start,ray_end);

        irun=irun+1;

#    outfile = filename.replace("_orig.nc","-out.nc");
#    filestem = filename.replace("_orig.nc","");
#    pattern = "{}-rev??.nc".format(filestem); 
#    files2cat = [f for f in fnmatch.filter(os.listdir('.'), pattern)];
#    print(files2cat);

#    nco.ncrcat(input=files2cat,output=outfile);
 
    
def trim_ncfile(infile,outfile,ray_start,ray_end):
    file_dtstr = infile.split("_")[1]

    print(infile);
    print(file_dtstr);
    DS = nc4.Dataset(infile);

    optstr = "-d time,{},{}".format(ray_start,ray_end);
    nco.ncks(input=infile, output=outfile, options=[optstr]);
    
    DS.close()
    
try:
    opts, args = getopt.getopt(sys.argv[1:], "f:p:", ["file=","path="])
except getopt.GetoptError as err:
    # print help information and exit:                                                                              
    print(err) # will print something like "option -a not recognized                                                
    sys.exit(2)

for o, a in opts:
    if o == "-f":
        infile = a;
    elif o == "-p":
        outpath = a;
    else:
        assert False, "unhandled option"

trim_camra(infile,outpath);
    

