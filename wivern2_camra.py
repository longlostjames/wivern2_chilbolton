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
import cmocean
import getopt, sys
import datetime

try:
    opts, args = getopt.getopt(sys.argv[1:], "d:", ["date="])
except getopt.GetoptError as err:
    # print help information and exit:                                                                                            
    print(err) # will print something like "option -a not recognized                                                              
    sys.exit(2)

for o, a in opts:
    if o == "-d":
        datestr = a;
    

mpl.use('Agg')


def load_3GHz(infile,dBZh_offset,ZDR_offset,range_offset):
    
    #infile:  e.g.  radar-camra_20210621112731_fix-ts.nc'

    DS = nc4.Dataset(infile);
    dt = cftime.num2pydate(DS['time'][:],DS['time'].units)
    
    if dt[-1]<dt[-2]:
        dt[-1]=dt[-1]+datetime.timedelta(days=1)
        
    cable_losses = DS.cable_losses;   # 4.8
    radar_const  = DS.radar_constant; # 64.7 
    rec_gain     = DS.receiver_gain;  # 45.5
    freq         = DS['frequency'][:];
    prf          = DS['prf'][:];
    

    adcmax       = 4096;
    zlothresh    = 3840;
    zlomin       =  -70;
    zhimin       =  -38;
    zcxmin       =  -77;
    zloscale     =    0.015625;
    zhiscale     =    0.015625;
    
    zlotable     = np.power(10,((zlomin+np.arange(adcmax)*zloscale)/10.));
    zhitable     = np.power(10,((zhimin+np.arange(adcmax)*zhiscale)/10.));
    
    dBZv_offset = dBZh_offset-ZDR_offset;
    
    dBZcal = radar_const-rec_gain+cable_losses;
    
    Zh_cal = 10**((dBZcal+dBZh_offset)/10);
    Zv_cal = 10**((dBZcal+dBZv_offset)/10);
   
    temp = DS['ZHI'][:,:];
    
    nray,npulse,ngate = temp.shape;

    print(nray,npulse,ngate);
    
    rng      = DS['range'][:ngate];
    range_km = (rng+range_offset)/1000.; # range in km
    
    ITX = DS['ITX'][:,:,:]; 
    QTX = DS['QTX'][:,:,:];
    IRX = DS['IRX'][:,:,:];
    QRX = DS['QRX'][:,:,:];

    ITXmean = np.mean(ITX,axis=1);
    QTXmean = np.mean(QTX,axis=1);
    IRXmean = np.mean(IRX,axis=1);
    QRXmean = np.mean(QRX,axis=1);

    for iray in np.arange(nray):
        for isample in np.arange(ngate):
            ITX[iray,:,isample] = ITX[iray,:,isample]-ITXmean[iray,isample];
            QTX[iray,:,isample] = QTX[iray,:,isample]-QTXmean[iray,isample];
            IRX[iray,:,isample] = IRX[iray,:,isample]-IRXmean[iray,isample];
            QRX[iray,:,isample] = QRX[iray,:,isample]-QRXmean[iray,isample];
    
    Vtx = ITX - 1j*QTX;
    Vrx = IRX + 1j*QRX;
 
    V = np.multiply(Vrx, Vtx);
    
    V = ma.masked_where(V==0,V);

    V = ma.divide(np.conj(V),np.abs(V));
    V[:,1::2,:]=V[:,1::2,:]*-1.;
    
    zlo_int  = DS['ZLO'][:,:,:].astype(int);
    zhi_int  = DS['ZHI'][:,:,:].astype(int);
    
    ZED          = zlo_int.astype(float);
    index        = zlo_int<=zlothresh;  
    ZED[index]   = zlotable[zlo_int[index]];
    index        = zlo_int>zlothresh;
    ZED[index]   = zhitable[zhi_int[index]];
    ZED[:, ::2,:] = ZED[:, ::2,:]*Zh_cal;
    ZED[:,1::2,:] = ZED[:,1::2,:]*Zv_cal;

    #ZED = ZED*Zscalefact;
    V = V*np.sqrt(ZED);
    

    
    I = np.real(V);
    Q = np.imag(V);
    
    I=np.squeeze(I);
    Q=np.squeeze(Q);
    
    Zpp   = np.empty([nray,ngate]);
    velpp = np.empty([nray,ngate]);
    

    PRF=305; #Hz
    freq=3e9; #Hz

    c=299792458; #m/s
    lamda=c/freq;
    vf=lamda*PRF/4;
    
    for iray in np.arange(nray):
        for igate in np.arange(ngate):
            Zpp[iray,igate]= np.mean(I[iray,:,igate]**2+Q[iray,:,igate]**2);
            
            I0 = I[iray,2:npulse:2,igate];
            Q0 = Q[iray,2:npulse:2,igate];
            Im1 = I[iray,1:npulse-1:2,igate];
            Qm1 = Q[iray,1:npulse-1:2,igate];
            Im2 = I[iray,0:npulse-2:2,igate];
            Qm2 = Q[iray,0:npulse-2:2,igate];
            Ip1 = I[iray,3:npulse:2,igate];
            Qp1 = Q[iray,3:npulse:2,igate];
  
            
            real_vh_arr = I0*Im2+Q0*Qm2;
            imag_vh_arr = Q0*Im2-I0*Qm2;
            real_vv_arr = Ip1*Im1+Qp1*Qm1;
            imag_vv_arr = Qp1*Im1-Ip1*Qm1;
            
            real_vh = np.sum(real_vh_arr);
            imag_vh = np.sum(imag_vh_arr);
            real_vv = np.sum(real_vv_arr);
            imag_vv = np.sum(imag_vv_arr);
            
            velpp[iray,igate]=ma.arctan2(imag_vh+imag_vv,real_vh+real_vv)*vf/np.pi;

            
            #real_vh=0; imag_vh=0;
            #real_vv=0; imag_vv=0;

            #for ip in np.arange(2,np.int(npulse/2),2):
            #   real_vh=real_vh+I[iray,ip,igate]*I[iray,ip-2,igate]+Q[iray,ip,igate]*Q[iray,ip-2,igate]; 
            #   imag_vh=imag_vh+Q[iray,ip,igate]*I[iray,ip-2,igate]-I[iray,ip,igate]*Q[iray,ip-2,igate];
         
            #   real_vv=real_vv+I[iray,ip+1,igate]*I[iray,ip-1,igate]+Q[iray,ip+1,igate]*Q[iray,ip-1,igate]; 
            #   imag_vv=imag_vv+Q[iray,ip+1,igate]*I[iray,ip-1,igate]-I[iray,ip+1,igate]*Q[iray,ip-1,igate];              
            #velpp[iray,igate]=ma.arctan2(imag_vh+imag_vv,real_vh+real_vv)*vf/np.pi;
        
    Zpp=10*np.log10(Zpp)+10*np.log10((range_km**2)*np.ones([nray,1]));
    
    blindzone = Zpp*0;
    blindzone[:,np.arange(15)] = 1;
    
    # Als0 mask first ray which is corrupted
    blindzone[0,:] = 1;
    
    Zpp = ma.masked_where(blindzone==1,Zpp);
    velpp = ma.masked_where(blindzone==1,velpp);
    

    PRF=prf/2.0; #Hz - effective PRF for each polarization
    #freq=3e9; #Hz

    c=299792458; #m/s
    wavelength=c/freq;
    vf=wavelength*PRF/4; # Folding velocity
    
    DS_out = dict();
    DS_out['I'] = I;
    DS_out['Q'] = Q;
    DS_out['Zpp'] = Zpp;
    DS_out['velpp'] = velpp;
    DS_out['range'] = range_km;
    DS_out['datetime'] = dt;
    DS_out['folding_velocity'] = vf;
    
    DS.close();
    
    return DS_out


#datestr = '20201114';

raw_ts_root_3ghz = '/gws/nopw/j04/ncas_obs/cao/raw_data/ncas-radar-camra-1/data/campaign/wivern2/ts/ESA/'
#raw_ts_path_3ghz = os.path.join(raw_ts_root_3ghz,datestr,'EventB');
raw_ts_path_3ghz = os.path.join(raw_ts_root_3ghz,datestr);


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

dt_min = dt0[0];
dt_max = dt1[-1];
if dt1[-1]<dt1[-2]:
        dt_max=dt_max+datetime.timedelta(days=1)

print(dt_min, dt_max);


DS0.close()
DS1.close()

infile = tsfiles_3ghz[istartfile];

ifile = 1;
print('{}/{} {}'.format(ifile,nfiles,infile));

dBZh_offset  = 9.0; # dB
ZDR_offset   = 0.6; # dB
range_offset = -24.0; # metres

DS0 = load_3GHz(infile,dBZh_offset,ZDR_offset,range_offset);

vf = DS0['folding_velocity'];


fig, ax = plt.subplots(2,1,figsize=(12,12),constrained_layout=True)
fig.set_constrained_layout_pads(w_pad=2 / 72, h_pad=2 / 72, hspace=0.2,wspace=0.2)

fig.patch.set_facecolor('white')

myFmt = mdates.DateFormatter('%H:%M')

ax[0].xaxis.set_major_formatter(myFmt);
ax[1].xaxis.set_major_formatter(myFmt);

hmax = 12;

ax[0].set_xlim(dt_min,dt_max);
ax[0].set_ylim(0,hmax);
ax[1].set_xlim(dt_min,dt_max);
ax[1].set_ylim(0,hmax);

ax[0].set_xlabel('Time (UTC)');
ax[1].set_xlabel('Time (UTC)');

ax[0].set_ylabel('Height above radar (km)');
ax[1].set_ylabel('Height above radar (km)');

gate_width = DS0['range'][1]-DS0['range'][0];
gate_starts = DS0['range'][:] - gate_width/2.0;

vel_cmap = plt.cm.get_cmap('RdBu')
vel_cmap = vel_cmap.reversed()


h1 = ax[0].pcolor(DS0['datetime'],gate_starts,DS0['Zpp'][1:,1:].transpose(),vmin=-20,vmax=50,cmap='pyart_HomeyerRainbow');
cb1=plt.colorbar(h1,ax=ax[0],orientation='vertical');
cb1.ax.set_ylabel("Reflectivity Factor (dBZ)");
titlestr = "Chilbolton 3GHz CAMRa Radar";
ax[0].set_title(titlestr,loc='left')
ax[0].set_title('{}'.format(dt0[0].date()),loc='right')

ax[1].xaxis.set_major_formatter(myFmt);
h2 = ax[1].pcolor(DS0['datetime'],gate_starts,DS0['velpp'][1:,1:].transpose(),vmin=-vf,vmax=vf,cmap=vel_cmap);
cb2=plt.colorbar(h2,ax=ax[1],orientation='vertical');
cb2.ax.set_ylabel("Doppler velocity (m/s)");
titlestr = "Chilbolton 3GHz CAMRa Radar";
ax[1].set_title(titlestr,loc='left')
ax[1].set_title('{}'.format(dt0[0].date()),loc='right')


for file in tsfiles_3ghz[istartfile+1:]:
    
    ifile = ifile+1;
    
    print('{}/{} {}'.format(ifile,nfiles,file));

    DS0 = load_3GHz(file,dBZh_offset,ZDR_offset,range_offset);
    gate_width = DS0['range'][1]-DS0['range'][0];
    gate_starts = DS0['range'][:]-gate_width/2.0;
    
    ax[0].pcolor(DS0['datetime'],gate_starts,DS0['Zpp'][1:,1:].transpose(),vmin=-20,vmax=50,cmap='pyart_HomeyerRainbow');
    ax[1].pcolor(DS0['datetime'],gate_starts,DS0['velpp'][1:,1:].transpose(),vmin=-vf,vmax=vf,cmap=vel_cmap);
    

ax[0].grid(True);
ax[1].grid(True);

figfile = "esa-wivern2-radar-3ghz_{}.png".format(datestr);

figpath = '/home/users/cjwalden/jupyter/wivern2'

plt.savefig(os.path.join(figpath,figfile),dpi=200)

plt.close();


