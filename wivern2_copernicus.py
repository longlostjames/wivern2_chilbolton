#!/usr/bin/env python
# coding: utf-8

import netCDF4 as nc4
import glob
import os
import fnmatch
import numpy as np;
import matplotlib as mpl
import matplotlib.pyplot as plt
import cftime
import numpy.ma as ma
import pyart;
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
from matplotlib import colors
import getopt, sys
import cmocean
import datetime


try:
    opts, args = getopt.getopt(sys.argv[1:], "d:e:", ["date=","event="])
except getopt.GetoptError as err:
    # print help information and exit:                                                                                            
    print(err) # will print something like "option -a not recognized                                                              
    sys.exit(2)

eventstr = "";

for o, a in opts:
    if o == "-d":
        datestr = a;
    elif o == "-e":
        eventstr = a;

mpl.use('Agg')


def load_35GHz(infile):
    header_dt = np.dtype([('ngates','<i4'),('npulses','<i4')]);
    time_dt   = np.dtype([('year', '<i4'), ('month', '<i4'),('day','<i4'),('hour','<i4'),('minute','<i4'),('second','<i4'),('centisecond','<i4')]);

    A = np.fromfile(infile, dtype=header_dt, count=1)[0];
    ngate = A['ngates'];
    npulse = A['npulses'];
    iq_dt = np.dtype([('IQ',"({},{})<i2".format(npulse,ngate))]);   
    record_dt = np.dtype([('time',time_dt),('I',iq_dt),('Q',iq_dt)]);
    records = np.fromfile(infile,dtype=record_dt, offset=8);
    times = [datetime.datetime(time['year'],time['month'],time['day'],time['hour'],time['minute'],time['second'],time['centisecond']*10000) for time in records['time']];
    nray = len(times);
    print(nray);
    Itmp = [record['I'][0].reshape(npulse,ngate) for record in records];
    Qtmp = [record['Q'][0].reshape(npulse,ngate) for record in records];
    I = np.stack(Itmp);
    Q = np.stack(Qtmp);
    print(I.shape);
    
    Imean = np.mean(I,axis=1);
    Qmean = np.mean(Q,axis=1);

    for iray in np.arange(nray):
        for isample in np.arange(ngate):
            I[iray,:,isample] = I[iray,:,isample]-Imean[iray,isample];
            Q[iray,:,isample] = Q[iray,:,isample]-Qmean[iray,isample];
   
    gate_len=0.0599585;
    range_offset=-0.660;#Range offset (Default:-630.0m; 17/04/2014: -660.0m)
    freq=34.96e9; # Hz
    PRF=5000; # Hz     
    Zcal=-146.8-1.0+115; #dB (17/04/2014: -1dB)
      
    blindrange=15; # Blind range in 60-m range gates 
    #noise_nom=5e-11; % Nominal noise level
    range_km=np.arange(ngate)*gate_len+range_offset;
    c=299792458; #m/s
    wavelength=c/freq;
    vf=wavelength*PRF/4;
    Zcal=Zcal+np.log10(PRF/512)*10

    I=I*np.sqrt(10**(Zcal/10.0));
    Q=Q*np.sqrt(10**(Zcal/10.0));
    
    Zpp   = np.empty([nray,ngate]);
    velpp = np.empty([nray,ngate]);
        
    for iray in np.arange(nray):
        for igate in np.arange(ngate):
            Zpp[iray,igate]= np.mean(I[iray,:,igate]**2+Q[iray,:,igate]**2);
            
            I0 = I[iray,1:npulse,igate];
            Q0 = Q[iray,1:npulse,igate];
            Im1 = I[iray,0:npulse-1,igate];
            Qm1 = Q[iray,0:npulse-1,igate];
              
            real_vh_arr = I0*Im1+Q0*Qm1;
            imag_vh_arr = Q0*Im1-I0*Qm1;
            
            real_vh = np.sum(real_vh_arr);
            imag_vh = np.sum(imag_vh_arr);
            
            velpp[iray,igate]=ma.arctan2(imag_vh,real_vh)*vf/np.pi;
 
    Zpp=10*np.log10(Zpp)+10*np.log10((range_km**2)*np.ones([nray,1]));
    
    blindzone = Zpp*0;
    blindzone[:,np.arange(15)] = 1;
    
    # Also mask first ray
    blindzone[0,:] = 1;
    
    Zpp = ma.masked_where(blindzone==1,Zpp);
    velpp = ma.masked_where(blindzone==1,velpp);
    
    DS_out = dict();
    DS_out['I'] = I;
    DS_out['Q'] = Q;
    DS_out['Zpp'] = Zpp;
    DS_out['velpp'] = velpp;
    DS_out['range'] = range_km;
    DS_out['datetime'] = times;
    DS_out['folding_velocity'] = vf;
    
            
    return DS_out

#datestr = '20211005';

# 3GHz files to determine start and end times
raw_ts_root_3ghz = '/gws/nopw/j04/ncas_obs/cao/raw_data/ncas-radar-camra-1/data/campaign/wivern2/ts/ESA/L0'
raw_ts_path_3ghz = os.path.join(raw_ts_root_3ghz,datestr,eventstr);
#raw_ts_path_3ghz = os.path.join(raw_ts_root_3ghz,datestr);

pattern = '*{}*_fix-ts.nc'.format(datestr);

for root,dirs,files in os.walk(raw_ts_path_3ghz):
    tsfiles_3ghz = [os.path.join(root,f) for f in fnmatch.filter(files, pattern)];
    
istartfile = 0;
iendfile = -1;

DS0 = nc4.Dataset(tsfiles_3ghz[istartfile]);
DS1 = nc4.Dataset(tsfiles_3ghz[iendfile]);

print(tsfiles_3ghz[istartfile:iendfile]);
nfiles = len(tsfiles_3ghz[istartfile:iendfile])+1;

dt0 = cftime.num2pydate(DS0['time'][:],DS0['time'].units);
dt1 = cftime.num2pydate(DS1['time'][:],DS1['time'].units);

dt_min_3ghz = dt0[0];
dt_max_3ghz = dt1[-1];
if dt1[-1]<dt1[-2]:
        dt_max_3ghz=dt_max_3ghz+datetime.timedelta(days=1)

print(dt_min_3ghz, dt_max_3ghz);


DS0.close()
DS1.close()

raw_ts_root_35ghz = '/gws/nopw/j04/ncas_obs/cao/raw_data/ncas-radar-ka-band-1/data/campaign/wivern2/ts/ESA/L0/'
raw_ts_path_35ghz = os.path.join(raw_ts_root_35ghz,datestr,eventstr);

pattern = '{}*_iqdata.bin'.format(datestr);

for root,dirs,files in os.walk(raw_ts_path_35ghz):
    tsfiles_35ghz_all = [os.path.join(root,f) for f in fnmatch.filter(files, pattern)];
    
filetimes = np.array([datetime.datetime.strptime(f.split('/')[-1],'%Y%m%d%H%M%S%f_iqdata.bin') for f in tsfiles_35ghz_all]);

#print(filetimes);


indices2keep = np.squeeze(np.where((filetimes>=dt_min_3ghz) & (filetimes<=dt_max_3ghz)));

print(indices2keep);

tsfiles_35ghz = [tsfiles_35ghz_all[ifile] for ifile in indices2keep];
filetimes_35ghz = [filetimes[ifile] for ifile in indices2keep];

#filetimes = [datetime.datetime.strptime(f.split('/')[-1],'%Y%m%d%H%M%S%f_iqdata.bin') for f in tsfiles_35ghz];
for ifile in range(len(tsfiles_35ghz)):
    print(tsfiles_35ghz[ifile].split('/')[-1],filetimes_35ghz[ifile]);
                   
infile = tsfiles_35ghz[0];

print(infile);

DS_out = load_35GHz(infile);

istartfile = 0;
iendfile = -1;

DS0 = load_35GHz(tsfiles_35ghz[istartfile]);
DS1 = load_35GHz(tsfiles_35ghz[iendfile]);

nfiles = len(tsfiles_35ghz);

dt0 = DS0['datetime'][:];
dt1 = DS1['datetime'][:];

dt_min = dt0[0];
dt_max = dt1[-1];

print(dt_min_3ghz, dt_max_3ghz);
print(dt_min, dt_max);


infile = tsfiles_35ghz[istartfile];

ifile = 1;
print('{}/{} {}'.format(ifile,nfiles,infile));

#dBZh_offset  = 9.0; # dB
#ZDR_offset   = 0.6; # dB
#range_offset = -24.0; # metres

DS0 = load_35GHz(infile);

vf = DS0['folding_velocity'];

fig, ax = plt.subplots(2,1,figsize=(12,12),constrained_layout=True)
fig.set_constrained_layout_pads(w_pad=2 / 72, h_pad=2 / 72, hspace=0.2,wspace=0.2)

fig.patch.set_facecolor('white')

myFmt = mdates.DateFormatter('%H:%M')

ax[0].xaxis.set_major_formatter(myFmt);
ax[1].xaxis.set_major_formatter(myFmt);

hmax = 12;

ax[0].set_xlim(dt_min_3ghz,dt_max_3ghz);
ax[0].set_ylim(0,hmax);
ax[1].set_xlim(dt_min_3ghz,dt_max_3ghz);
ax[1].set_ylim(0,hmax);

ax[0].set_xlabel('Time (UTC)');
ax[1].set_xlabel('Time (UTC)');

ax[0].set_ylabel('Height above radar (km)');
ax[1].set_ylabel('Height above radar (km)');

gate_width = DS0['range'][1]-DS0['range'][0];
gate_starts = DS0['range'][:] - gate_width/2.0;

vel_cmap = plt.cm.get_cmap('RdBu')
vel_cmap = vel_cmap.reversed()


h1 = ax[0].pcolormesh(DS0['datetime'],gate_starts,DS0['Zpp'][1:,1:].transpose(),vmin=-20,vmax=50,cmap='pyart_HomeyerRainbow');
cb1=plt.colorbar(h1,ax=ax[0],orientation='vertical');
cb1.ax.set_ylabel("Reflectivity Factor (dBZ)");
titlestr = "Chilbolton 35GHz Copernicus Radar";
ax[0].set_title(titlestr,loc='left')
ax[0].set_title('{}'.format(dt0[0].date()),loc='right')

ax[1].xaxis.set_major_formatter(myFmt);
h2 = ax[1].pcolormesh(DS0['datetime'],gate_starts,DS0['velpp'][1:,1:].transpose(),vmin=-vf,vmax=vf,cmap=vel_cmap);
cb2=plt.colorbar(h2,ax=ax[1],orientation='vertical');
cb2.ax.set_ylabel("Doppler velocity (m/s)");
titlestr = "Chilbolton 35GHz Copernicus Radar";
ax[1].set_title(titlestr,loc='left')
ax[1].set_title('{}'.format(dt0[0].date()),loc='right')

for ifile in np.arange(istartfile+1,istartfile+nfiles):
    
    file = tsfiles_35ghz[ifile];
        
    print('{}/{} {}'.format(ifile+1,nfiles,file));

    DS0 = load_35GHz(file);
    gate_width = DS0['range'][1]-DS0['range'][0];
    gate_starts = DS0['range'][:]-gate_width/2.0;
    
    ax[0].pcolormesh(DS0['datetime'],gate_starts,DS0['Zpp'][1:,1:].transpose(),vmin=-20,vmax=50,cmap='pyart_HomeyerRainbow');
    ax[1].pcolormesh(DS0['datetime'],gate_starts,DS0['velpp'][1:,1:].transpose(),vmin=-vf,vmax=vf,cmap=vel_cmap);
    

ax[0].grid(True);
ax[1].grid(True);



figfile = "esa-wivern2-radar-{}-{}_35ghz.png".format(datestr,eventstr);

figpath = '/home/users/cjwalden/jupyter/wivern2'

plt.savefig(os.path.join(figpath,figfile),dpi=200)
