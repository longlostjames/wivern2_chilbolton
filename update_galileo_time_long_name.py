#!/usr/bin/env python
# coding: utf-8

import sys, getopt

import wivern_chilbolton_utils as wivern
import os, getpass, glob
import numpy as np
import fnmatch

user = getpass.getuser()


def main(argv):

   path = '.';

   try:
      opts, args = getopt.getopt(argv,"hp:d:",["path==", "date="])
   except getopt.GetoptError:
      print ('update_galileo_time_long_name.py -p <path> -d <date>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print ('update_galileo_time_long_name.py -p <path> -d <date>')
         sys.exit()
      elif opt == '-p':
         path = arg;
      elif opt == '-d':
         datestr = arg;

   pattern = '*_fix-ts_l1_v1.0.nc'

   datepath = os.path.join(path,datestr);

   print(datepath);

   l1dirs = [];

   for root,dirs,files in os.walk(datepath):
      l1dirs += [os.path.join(root,d) for d in dirs];

   if l1dirs==[]:
      l1dirs = datepath;
      l1files = [os.path.join(root,f) for f in fnmatch.filter(files, pattern)];
      for f in l1files:
         wivern.update_galileo_time_long_name(f);
   else:
      for d in l1dirs:
         print(d);
         l1files = []; 
         for root,dirs,files in os.walk(d):
            l1files += [os.path.join(root,f) for f in fnmatch.filter(files, pattern)];
            event = os.path.split(d)[-1];
         if len(l1files)>0:
            for f in l1files:
               wivern.update_galileo_time_long_name(f);


if __name__ == "__main__":
   main(sys.argv[1:])



