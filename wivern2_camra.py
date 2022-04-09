#!/usr/bin/env python
# coding: utf-8

# In[1]:


import netCDF4 as nc4
import glob
import os
import fnmatch
import numpy as np;
import matplotlib.pyplot as plt
import cftime
import numpy.ma as ma
import pyart;
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
from matplotlib import colors
import cmocean
import getpass
import socket

from datetime import tzinfo, datetime, time


# In[2]:


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


# In[4]:


module_version = 1.0;


def convert_camra_l0a2l0b(infile,outfile):
    
    """This routine converts raw (Level 0a) time series data from the Chilbolton Advanced Meteorological Radar (CAMRa) to Level 0b data.
    Processing involves removing redundant dimensions, and removing bias from the ADC samples of transmit and receive I and Q.
    Single estimates of I and Q for each transmit pulse are produced and stored in the output file.  
    Metadata are added specific to the WIVERN-2 project ground-based observations collected in 2020-2021.
    
    :param infile: Full path of NetCDF Level 0a raw data file, e.g. `<path-to-file>/radar-camra_<YYYYddmmHHMMSS_fix-ts.nc`
    :type infile: str
    
    :param outfile: Full path of NetCDF Level 0b output file, e.g. `<path-to-file>/ncas-radar-camra-1_cao_20201210-212823_fix-ts_l0b_v1.0.nc`
    """
    
    DSin = nc4.Dataset(infile);    
    
    dt_start = cftime.num2pydate(DSin['time'][0],DSin['time'].units)
    dt_end   = cftime.num2pydate(DSin['time'][-1],DSin['time'].units)

    
    try: outfile.close()  # just to be safe, make sure dataset is not already open.
    except: pass
    DSout = nc4.Dataset(outfile,mode='w',format='NETCDF4') 
    print(outfile)

    # ------------------------
    # Set up global attributes
    # ------------------------
    DSout.product_version = "v1.0" ;
    DSout.processing_level = "0b" ;
    
    DSout.licence = "Data usage licence - UK Open Government Licence agreement: \n http://www.nationalarchives.gov.uk/doc/open-government-licence" ;
    DSout.acknowledgement = "Acknowledgement is required of UK Research and Innovation as the data provider (in partnership with the National Centre for Atmospheric Science) whenever and wherever these data are used." ;
    DSout.platform = "Chilbolton Atmospheric Observatory" ;
    DSout.platform_type = "stationary_platform" ;
    DSout.title = "Time series from 3 GHz CAMRa radar collected for ESA WIVERN-2 campaign at Chilbolton Observatory";
    DSout.creator_name = "Chris Walden" ;
    DSout.creator_email = "chris.walden@ncas.ac.uk" ;
    DSout.creator_url = "https://orcid.org/0000-0002-5718-466X" ;
    DSout.institution = "National Centre for Atmospheric Science (NCAS)";
    DSout.instrument_name = "ncas-radar-camra-1";
    DSout.instrument_software = "radar-camra-rec" ;
    DSout.instrument_software_version = "1.4 Rev 58" ;

    DSout.references = "https://doi.org/10.1049/ecej:19940205; http://purl.org/net/epubs/work/63318; https://doi.org/10.1109/IGARSS.2006.429; https://doi.org/10.3390/atmos10110714";
    DSout.source = "3GHz Chilbolton Advanced Meteorological Radar (CAMRa)";
    DSout.comment = "";
    DSout.project = "WIVERN-2 Doppler Wind Radar Science Performance Study";
    DSout.project_principal_investigator = "Anthony Illingworth";
    DSout.project_principal_investigator_email = "a.j.illingworth@reading.ac.uk";
    DSout.project_principal_investigator_url = "https://orcid.org/0000-0002-5774-8410";
    
    DSout.processing_software_url = "https://github.com/longlostjames/wivern_chilbolton_utils.git";
    DSout.processing_software_version = "1.0"; 
        
    DSout.scantype = "vertical_pointing";
    
    DSout.time_coverage_start = datetime.strftime(dt_start,'%Y-%m-%dT%H:%M:%SZ');
    DSout.time_coverage_end = datetime.strftime(dt_end,'%Y-%m-%dT%H:%M:%SZ');
    DSout.geospatial_bounds = "51.1450N -1.4384E";
    
    DSout.pulse_compression = "false";
    
    DSout.ADC_bits_per_sample = DSin.ADC_bits_per_sample;
    DSout.ADC_channels        = DSin.ADC_channels;
    
    # ----------------
    # Scalar variables
    # ----------------
    varin = DSin['latitude'];
    varout = DSout.createVariable('latitude',varin.datatype);    
    varout.standard_name = 'latitude';
    varout.long_name = 'latitude of the antenna';
    varout.units = 'degree_north';
    varout[:]=51.1450;
    
    varin = DSin['longitude'];
    varout = DSout.createVariable('longitude',varin.datatype);    
    varout.standard_name = 'longitude';
    varout.long_name = 'longitude of the antenna';
    varout.units = 'degree_east';
    varout[:]=-1.4384;

    varin = DSin['height'];
    varout = DSout.createVariable('altitude',varin.datatype);    
    varout.standard_name = 'altitude';
    varout.long_name = 'altitude of the elevation axis above the geoid (WGS84)';
    varout.units = 'm';
    varout[:]=146.3;
    
    varin = DSin['frequency'];
    varout = DSout.createVariable('frequency',varin.datatype);    
    varout.standard_name = 'radiation_frequency';
    varout.long_name = 'frequency of transmitted radiation';
    varout.units = 'GHz';
    varout[:]=varin[:];
  
    varin = DSin['prf'];
    varout = DSout.createVariable('prf',varin.datatype);    
    varout.long_name = 'pulse repetition frequency';
    varout.units = 'Hz';
    varout[:]=varin[:];
    
    varin = DSin['beamwidthH'];
    varout = DSout.createVariable('beamwidthH',varin.datatype);    
    varout.long_name = 'horizontal angular beamwidth';
    varout.units = 'degree';
    varout[:]=varin[:];
    
    varin = DSin['beamwidthV'];
    varout = DSout.createVariable('beamwidthV',varin.datatype);    
    varout.long_name = 'vertical angular beamwidth';
    varout.units = 'degree';
    varout[:]=varin[:];
    
    varin = DSin['antenna_diameter'];
    varout = DSout.createVariable('antenna_diameter',varin.datatype);    
    varout.long_name = 'antenna diameter';
    varout.units = 'm';
    varout[:]=varin[:];
    
    varout = DSout.createVariable('antenna_focal_length',varin.datatype);    
    varout.long_name = 'focal length of antenna';
    varout.units = 'm';
    varout[:]=np.float32(9.0);
    
    varout = DSout.createVariable('antenna_focus_radial_location');
    varout.long_name = 'distance along boresight from elevation axis to antenna focus';
    varout.units = 'm';
    varout = np.float32(14.34);
    
    varin = DSin['pulse_period'];
    varout = DSout.createVariable('pulse_width',varin.datatype);    
    varout.long_name = 'pulse width';
    varout.units = 'us';
    varout[:]=varin[:];
    
    varin = DSin['transmit_power'];
    varout = DSout.createVariable('transmit_power',varin.datatype);    
    varout.long_name = 'peak transmitted power';
    varout.units = 'W';
    varout[:]=varin[:];
    
    varin = DSin['clock'];
    varout = DSout.createVariable('clock',varin.datatype);    
    varout.long_name = 'clock input to timer card';
    varout.units = 'Hz';
    varout[:]=varin[:];
    
    varout = DSout.createVariable('clock_divfactor','i4');    
    varout.long_name = 'clock divide factor';
    varout.units = '1';
    varout[:]=DSin['clock'].clock_divfactor; 
    
    varout = DSout.createVariable('delay_clocks','i4');    
    varout.long_name = 'clock cycles before sampling is initiated';
    varout.units = '1';
    varout[:]=DSin.delay_clocks; 
    
    varout = DSout.createVariable('samples_per_pulse','i4');    
    varout.long_name = 'number of samples per pulse';
    varout.units = '1';
    varout[:]=DSin.samples_per_pulse;    

    varout = DSout.createVariable('pulses_per_daq_cycle','i4');    
    varout.long_name = 'number of pulses per data acquisition cycle';
    varout.units = '1';
    varout[:]=DSin.pulses_per_daq_cycle; 
    
    varout = DSout.createVariable('pulses_per_ray','i4');    
    varout.long_name = 'number of pulses per ray';
    varout.units = '1';
    varout[:]=DSin.pulses_per_ray;    

    varout = DSout.createVariable('radar_constant','f4');    
    varout.long_name = 'radar constant';
    varout.units = 'dB';
    varout[:]=DSin.radar_constant; 
    
    varout = DSout.createVariable('receiver_gain','f4');    
    varout.long_name = 'receiver gain';
    varout.units = 'dB';
    varout[:]=DSin.receiver_gain; 

    varout = DSout.createVariable('cable_losses','f4');    
    varout.long_name = 'cable losses';
    varout.units = 'dB';
    varout[:]=DSin.cable_losses; 

    varout = DSout.createVariable('extra_attenuation','f4');    
    varout.long_name = 'extra attenuation';
    varout.units = 'dB';
    varout[:]=DSin.extra_attenuation; 

    # ---------------
    # Copy dimensions
    # ---------------
    the_dim = DSin.dimensions['time'];
    DSout.createDimension('time', len(the_dim) if not the_dim.isunlimited() else None)

    the_dim = DSin.dimensions['pulses'];
    DSout.createDimension('pulse', len(the_dim) if not the_dim.isunlimited() else None)

    the_dim = DSin.dimensions['samples'];
    DSout.createDimension('range', len(the_dim) if not the_dim.isunlimited() else None)
    
    # --------------------
    # Coordinate variables
    # --------------------
    varin = DSin['time'];
    varout = DSout.createVariable('time',varin.datatype,('time'));    
    varout.standard_name = 'time';
    varout.long_name = 'time at the end of each recorded ray';
    varout.units = varin.units;
    varout[:]=varin[:];

    varin = DSin['range'];
    varout = DSout.createVariable('range',varin.datatype,('range'));    
    varout.long_name = varin.long_name;
    varout.range_offset_applied = varin.range_offset;
    varout.units = varin.units;
    varout[:]=varin[:]; 
    
    # --------------------------
    # Antenna pointing variables
    # --------------------------
    varin = DSin['elevation'];
    varout = DSout.createVariable('elevation',varin.datatype,'time'));
    varout.long_name = 'elevation angle of the antenna boresight above the horizon';
    varout.units = 'degree';
    varout.elevation_offset_applied = np.float32(0.);
    varout[:] = varin[:];
    
    varin = DSin['azimuth'];
    varout = DSout.createVariable('azimuth',varin.datatype,'time'));
    varout.long_name = 'azimuth angle of the antenna boresight clockwise from the grid north';
    varout.comment = 'More generally this is the azimuth angle of the plane perpendicular to the elevation axis, which remains defined even when the elevation is 90 degree';
    varout.units = 'degree';
    varout.azimuth_offset_applied = np.float32(0.);
    varout[:] = varin[:];
        
    # ---------------
    # Field variables
    # ---------------
    varin = DSin['ZLO'];
    varout = DSout.createVariable('ZLO',varin.datatype,('time','pulse','range'),zlib=True);    
    varout.long_name = 'ZLO log amplifier output (channel with +12dB gain)';
    varout.units = 'dB';
    varout[:]=varin[:];
    comment_string  = "This is an estimator for co-polar radar equivalent reflectivity factor.\n"
    comment_string += "It does not take into account range correction, the radar constant, receiver gain or cable losses.\n"
    comment_string += "The data are packed and only values in the range [0,3840] (equivalent to the actual_range when the data are unpacked) should be used."
    varout.comment = comment_string;
    varout.scale_factor = np.float32(0.015625);
    varout.add_offset = np.float32(-70.);
    varout.actual_range = [np.float32(-70.),np.float32(-10.)];
    
    varin = DSin['ZHI'];
    varout = DSout.createVariable('ZHI',varin.datatype,('time','pulse','range'),zlib=True);    
    varout.long_name = 'ZHI log amplifier output (channel with -20dB attenuation)';
    varout.units = 'dB';
    varout[:] = varin[:];
    comment_string  = "This is an estimator for co-polar radar equivalent reflectivity factor.\n"
    comment_string += "It does not take into account range correction, the radar constant, receiver gain or cable losses.\n"
    comment_string += "The data are packed and only values in the range [1793,4095] (equivalent to the actual_range when the data are unpacked) should be used."
    varout.comment = comment_string;
    varout.scale_factor = np.float32(0.015625);
    varout.add_offset = np.float32(-38.);
    varout.actual_range = [np.float32(-9.984375), np.float32(25.984375)];
    
    varin = DSin['ZCX'];
    varout = DSout.createVariable('ZCX',varin.datatype,('time','pulse','range'),zlib=True);    
    varout.long_name = 'cross-polar log amplifier output';
    varout.units = 'dB';
    varout[:]=varin[:];
    comment_string  = "This is an estimator for cross-polar radar equivalent reflectivity factor.\n"
    comment_string += "It does not take into account range correction, the radar constant, receiver gain or cable losses."
    varout.comment = comment_string;
    varout.scale_factor = np.float32(0.03125);
    varout.add_offset = np.float32(-77.);
    
    
    # -------------------------------------------------------
    # Determine bias-corrected I and Q of each transmit pulse
    # -------------------------------------------------------
    delay_clocks = DSout['delay_clocks'][:];
    clock_divfactor = DSout['clock_divfactor'][:];

    pre_tx    = (18   - delay_clocks) // clock_divfactor;
    tx_pulse  = (24   - delay_clocks) // clock_divfactor;
    post_tx   = (68   - delay_clocks) // clock_divfactor;
    hold_end  = (4708 - delay_clocks) // clock_divfactor;
    post_hold = (4748 - delay_clocks) // clock_divfactor;

    ITXin = DSin['ITX'][:,:,:];
    QTXin = DSin['QTX'][:,:,:];

    ITX_bias_by_gate = np.mean(ITXin,axis=1);
    QTX_bias_by_gate = np.mean(QTXin,axis=1);
    
    ITXnew = ITXin - ITX_bias_by_gate[:,None,:];
    QTXnew = QTXin - QTX_bias_by_gate[:,None,:];

    # Only use data while sample and hold is active
    # ---------------------------------------------
    sample_end = min([hold_end,len(DSin.dimensions['samples'])]);
    
    ITXout = np.mean(ITXnew[:,:,post_tx:sample_end],axis=2);
    QTXout = np.mean(QTXnew[:,:,post_tx:sample_end],axis=2);

    varout = DSout.createVariable('ITX','f4',('time','pulse'),zlib=True);    
    varout.long_name = 'bias-corrected samples of in-phase video signal for each transmitted pulse';
    varout.units = '1';
    varout[:]=ITXout;
    
    varout = DSout.createVariable('QTX','f4',('time','pulse'),zlib=True);    
    varout.long_name = 'bias-corrected samples of quadrature video signal for each transmitted pulse';
    varout.units = '1';
    varout[:]=QTXout;
    
    # ----------------------------------------
    # Determine bias-corrected receive I and Q
    # ----------------------------------------
    IRXin = DSin['IRX'][:,:,:];
    QRXin = DSin['QRX'][:,:,:];
    IRX_bias_by_gate = np.mean(IRXin,axis=1);
    QRX_bias_by_gate = np.mean(QRXin,axis=1);
    IRXout = IRXin - IRX_bias_by_gate[:,None,:];
    QRXout = QRXin - QRX_bias_by_gate[:,None,:];
    
    varout = DSout.createVariable('IRX','f4',('time','pulse','range'),zlib=True);    
    varout.long_name = 'bias-corrected samples of received in-phase video signal at output of IF limiting amplifier';
    varout.units = '1';
    varout[:]=IRXout;
    
    varin = DSin['QRX'];
    varout = DSout.createVariable('QRX','f4',('time','pulse','range'),zlib=True);    
    varout.long_name = 'bias-corrected samples of received quadrature video signal at output of IF limiting amplifier';
    varout.units = '1';
    varout[:]=QRXout;
    
    # -----------------------
    # Update history metadata
    # -----------------------
    user = getpass.getuser()
    
    updttime = datetime.utcnow()
    updttimestr = updttime.ctime()

    history = updttimestr + (" - user:" + user
    + " machine: " + socket.gethostname()
    + " program: wivern_chilbolton_utils.py convert_camra_l0a2l0b"
    + " version:" + str(module_version));
    
    DSout.history = DSin.history + "\n" + history;

    DSout.last_revised_date = datetime.strftime(updttime,'%Y-%m-%dT%H:%M:%SZ')

    DSin.close();
    DSout.close();


# In[3]:


def convert_camra_l0b2l1(infile,dBZh_offset,ZDR_offset,range_offset,outfile,version):
     
    """This routine converts Level 0a time series data from the Chilbolton Advanced Meteorological Radar (CAMRa) to Level 0b data.
    Processing is applied to rmenove redundant dimensions, and to

    It adds metadata specific to the WIVERN-2 project ground-based observations collected in 2020-2021.
    
    :param infile: Full path of NetCDF Level 0a raw data file, e.g. `<path-to-file>/radar-camra_<YYYYddmmHHMMSS_fix-ts.nc`
    :type infile: str
    
    :param outfile: Full path of NetCDF Level 0b output file, e.g. `<path-to-file>/ncas-radar-camra-1_cao_20201210-212823_fix-ts_l0b_v1.0.nc`
    """
    toexclude = ['ZLO', 'ZHI', 'ZCX', 'ITX', 'QTX', 'IRX', 'QRX'];

    with nc4.Dataset(infile) as src, nc4.Dataset(outfile,mode='w',format='NETCDF4') as dst:
        
        # copy global attributes all at once via dictionary
        dst.setncatts(src.__dict__)
        
        # copy dimensions
        for name, dimension in src.dimensions.items():
            dst.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))
        
        # copy all file data except for the excluded
        for name, variable in src.variables.items():
            if name not in toexclude:
                x = dst.createVariable(name, variable.datatype, variable.dimensions)
                dst[name].setncatts(src[name].__dict__)
                dst[name][:] = src[name][:]
                # copy variable attributes all at once via dictionary
    
    try: outfile.close()  # just to be safe, make sure dataset is not already open.
    except: pass
    DSout = nc4.Dataset(outfile,mode='r+',format='NETCDF4') 
    print(outfile)

    DSin = nc4.Dataset(infile);
    dt = cftime.num2pydate(DSin['time'][:],DSin['time'].units);
    
    DSout.product_version = "v{:.1f}".format(version) ;
    DSout.processing_level = "1" ;

    DSout.title = "Calibrated time series from 3 GHz CAMRa radar collected for ESA WIVERN-2 campaign at Chilbolton Observatory";
    
    DSout.comment = DSin.comment + "\n Correction to account for inverse square power loss with range has not been applied";
    
    varout = DSout.createVariable('I','f4',('time','pulse','range'),zlib=True);    
    varout.long_name = 'in-phase video signal';
    varout.units = '1';
    varout.comment = 'Values are derived from I/Q on unit circle mutiplied by square root of linear reflectivity factor';
    
    varout = DSout.createVariable('Q','f4',('time','pulse','range'),zlib=True);    
    varout.long_name = 'quadrature video signal';
    varout.units = '1';
    varout.comment = 'Values are derived from I/Q on unit circle mutiplied by square root of linear reflectivity factor';
    
    varout = DSout.createVariable('ZCX','f4',('time','pulse','range'),zlib=True);    
    varout.long_name = 'cross-polar radar equivalent reflectivity factor';
    varout.units = 'dBZ';
    varout.comment = '';
    
    cable_losses = DSout['cable_losses'][:];
    radar_const  = DSout['radar_constant'][:]; 
    rec_gain     = DSout['receiver_gain'][:];
    freq         = DSout['frequency'][:];
    prf          = DSout['prf'][:];
    
    dBZv_offset = dBZh_offset-ZDR_offset;
    
    dBZcal = radar_const-rec_gain+cable_losses;
    
    Zh_cal = 10**((dBZcal+dBZh_offset)/10);
    Zv_cal = 10**((dBZcal+dBZv_offset)/10);
   
    range_m  = DSin['range'][:];
    range_km = (range_m+range_offset)/1000.; # range in km
    
    ITX = DSin['ITX'][:,:]; 
    QTX = DSin['QTX'][:,:];
    IRX = DSin['IRX'][:,:,:];
    QRX = DSin['QRX'][:,:,:];

    Vtx = ITX - 1j*QTX;
    Vrx = IRX + 1j*QRX;
    
    V = np.multiply(Vrx[:,:,:], Vtx[:,:,None]);
    V = ma.masked_where(V==0,V);
    #V = ma.divide(np.conj(V),np.abs(V));
    V = ma.divide(np.conj(V),np.abs(V));

    V[:,1::2,:] = V[:,1::2,:]*-1.;
    
    ZLO  = DSin['ZLO'][:,:,:];
    ZHI  = DSin['ZHI'][:,:,:];
    
    threshold          = DSin['ZLO'].actual_range[1];
    
    ZED                = ZLO.copy();
    ZED[ZLO>threshold] = ZHI[ZLO>threshold];
        
    # Convert to linear ZED
    ZED = np.power(10,ZED/10.);
    
    ZED[:, ::2,:] = ZED[:, ::2,:]*Zh_cal;
    ZED[:,1::2,:] = ZED[:,1::2,:]*Zv_cal;

    ZCX = DSin['ZCX'][:,:,:];
    ZCX[:, ::2,:] = ZCX[:, ::2,:] + dBZcal + dBZh_offset;
    ZCX[:,1::2,:] = ZCX[:,1::2,:] + dBZcal + dBZv_offset;
    DSout['ZCX'][:] = ZCX;
    
    V = V*np.sqrt(ZED);
    
    I = np.real(V);
    Q = np.imag(V);

    DSout['I'][:] = I;
    DSout['Q'][:] = Q;
    DSout['range'][:] = range_m;
    DSout['range'].range_offset_applied += range_offset;

    
    user = getpass.getuser()
    
    updttime = datetime.utcnow()
    updttimestr = updttime.ctime()
    module_version = 1.0;

    history = updttimestr + (" - user:" + user
    + " machine: " + socket.gethostname()
    + " program: wivern_chilbolton_utils.py convert_camra_l0b2l1"
    + " version:" + str(module_version));
    
    print(history);

    DSout.history = DSin.history + "\n" + history;
    
    DSin.close();
    DSout.close();
    
    return


# In[ ]:


now = datetime.utcnow();
datestr = datetime.strftime(now,'%Y-%m-%dT%H:%M:%SZ')
print(datestr)


# In[5]:


datestr = '20201210';

raw_ts_root_3ghz = '/gws/nopw/j04/ncas_obs/cao/raw_data/ncas-radar-camra-1/data/campaign/wivern2/ts/ESA/L0'
raw_ts_path_3ghz = os.path.join(raw_ts_root_3ghz,datestr);


pattern = '*{}*_fix-ts.nc'.format(datestr);

for root,dirs,files in os.walk(raw_ts_path_3ghz):
    tsfiles_3ghz = [os.path.join(root,f) for f in fnmatch.filter(files, pattern)];

infile = tsfiles_3ghz[2];
data_version = 1.0;
outfile_string = os.path.split(infile)[1].replace('.nc','_l0b.nc');
splits = outfile_string.split('_');
instrument_name =splits[0].replace('radar-camra','ncas-radar-camra-1');
platform = 'cao';
datestr = splits[1][0:8];
timestr = splits[1][8:];
level = splits[3].split('.')[0];
outfile = '{}_{}_{}-{}_{}_{}_v{:.1f}.nc'.format(instrument_name,platform,datestr,timestr,splits[2],level,data_version)

print(outfile);

outfile.replace('radar-camra','ncas-radar-camra-1_cao');

outfile = os.path.join('/home/users/cjwalden/git',outfile);

convert_camra_L0a2L0b(infile,outfile);

l0bfile = outfile;
l1file = outfile.replace('l0b','l1');

convert_camra_l0b2l1(l0bfile,9.0,0.6,0.0,l1file,1.0);


# In[ ]:





# In[ ]:


DSin = nc4.Dataset(infile);
#print(DSin);
DS = nc4.Dataset(outfile);
#print(DS['ZLO']);
print(DS);
print(DS['latitude'][:]);
print(DS['longitude'][:]);
print(DS['altitude'][:]);
print(DS['ZHI']);

dtime = cftime.num2pydate(DS['time'][:],DS['time'].units)
print(DSin);
print(DSin['range']);
print(DSin['time']);

DSnew = load_3GHz_l0b(outfile,9.0,0.6,-24.0);


# In[ ]:


DS = nc4.Dataset(l1file);

dtime = cftime.num2pydate(DS['time'][:],DS['time'].units)

ZED = 5.*np.log10(DS['I'][...]**2+DS['Q'][...]**2);

LDR  = np.mean(np.power(10,DS['ZCX'][:,:,:]/10.0),axis=1);
LDR /= np.mean(np.power(10,ZED[:,:,:]/10.0),axis=1);

LDR = 10.*np.log10(LDR);




fig = plt.figure(figsize=(12,30))
fig.set_constrained_layout(True)
fig.set_constrained_layout_pads(w_pad=2 / 72, h_pad=2 / 72, hspace=0.2,wspace=0.2)

gs0 = fig.add_gridspec(6, 2)
ax1 = fig.add_subplot(gs0[0, :2]);
ax2 = fig.add_subplot(gs0[2, 0]);
ax3 = fig.add_subplot(gs0[2, 1]);
ax4 = fig.add_subplot(gs0[3, 0]);
ax5 = fig.add_subplot(gs0[3, 1]);
ax6 = fig.add_subplot(gs0[4, 0]);
ax7 = fig.add_subplot(gs0[4, 1]);
ax9 = fig.add_subplot(gs0[5, 0]);
ax10 = fig.add_subplot(gs0[5, 1]);
ax8 = fig.add_subplot(gs0[1, :2]);

h1=ax1.pcolormesh(dtime,DS['range'],ZED[:,1,:].T,vmin=-20,vmax=50,cmap='pyart_HomeyerRainbow');
cb1=plt.colorbar(h1,ax=ax1,orientation='vertical');
cb1.ax.set_ylabel("Reflectivity Factor (dBZ)");
myFmt = mdates.DateFormatter('%H:%M')
ax1.xaxis.set_major_formatter(myFmt);

h8=ax8.pcolormesh(dtime,DS['range'],LDR[:,:].T,vmin=-35,vmax=5,cmap='pyart_HomeyerRainbow');
cb8=plt.colorbar(h8,ax=ax8,orientation='vertical');
cb8.ax.set_ylabel("Cross-polar reflectivity factor (dBZ)");
myFmt = mdates.DateFormatter('%H:%M')
ax8.xaxis.set_major_formatter(myFmt);

irange = 13;
ax2.plot(DS['I'][:,::2,irange],DS['Q'][:,::2,irange],'.r');
ax2.plot(DS['I'][:,1::2,irange],DS['Q'][:,1::2,irange],'.b');
ax2.axis('equal');
ax2.set_xlabel('I');
ax2.set_ylabel('Q');
ax2.grid(True);
ax2.set_title('Range = {:5.2f}'.format(DS['range'][irange]));

irange = 14;
ax3.plot(DS['I'][:,::2,irange],DS['Q'][:,::2,irange],'.r');
ax3.plot(DS['I'][:,1::2,irange],DS['Q'][:,1::2,irange],'.b');
ax3.axis('equal');
ax3.set_xlabel('I');
ax3.set_ylabel('Q');
ax3.grid(True);
ax3.set_title('Range = {:5.2f}'.format(DS['range'][irange]));

irange = 15;
ax4.plot(DS['I'][:,::2,irange],DS['Q'][:,::2,irange],'.r');
ax4.plot(DS['I'][:,1::2,irange],DS['Q'][:,1::2,irange],'.b');
ax4.axis('equal');
ax4.set_xlabel('I');
ax4.set_ylabel('Q');
ax4.grid(True);
ax4.set_title('Range = {:5.2f}'.format(DS['range'][irange]));

irange = 16;
ax5.plot(DS['I'][:,::2,irange],DS['Q'][:,::2,irange],'.r');
ax5.plot(DS['I'][:,1::2,irange],DS['Q'][:,1::2,irange],'.b');
ax5.axis('equal');
ax5.set_xlabel('I');
ax5.set_ylabel('Q');
ax5.grid(True);
ax5.set_title('Range = {:5.2f}'.format(DS['range'][irange]));

irange = 17;
ax6.plot(DS['I'][:,::2,irange],DS['Q'][:,::2,irange],'.r');
ax6.plot(DS['I'][:,1::2,irange],DS['Q'][:,1::2,irange],'.b');
ax6.axis('equal');
ax6.set_xlabel('I');
ax6.set_ylabel('Q');
ax6.grid(True);
ax6.set_title('Range = {:5.2f}'.format(DS['range'][irange]));

irange = 18;
ax7.plot(DS['I'][:,::2,irange],DS['Q'][:,::2,irange],'.r');
ax7.plot(DS['I'][:,1::2,irange],DS['Q'][:,1::2,irange],'.b');
ax7.axis('equal');
ax7.set_xlabel('I');
ax7.set_ylabel('Q');
ax7.grid(True);
ax7.set_title('Range = {:5.2f}'.format(DS['range'][irange]));

irange = 19;
ax9.plot(DS['I'][:,::2,irange],DS['Q'][:,::2,irange],'.r');
ax9.plot(DS['I'][:,1::2,irange],DS['Q'][:,1::2,irange],'.b');
ax9.axis('equal');
ax9.set_xlabel('I');
ax9.set_ylabel('Q');
ax9.grid(True);
ax9.set_title('Range = {:5.2f}'.format(DS['range'][irange]));

irange = 20;
ax10.plot(DS['I'][:,::2,irange],DS['Q'][:,::2,irange],'.r');
ax10.plot(DS['I'][:,1::2,irange],DS['Q'][:,1::2,irange],'.b');
ax10.axis('equal');
ax10.set_xlabel('I');
ax10.set_ylabel('Q');
ax10.grid(True);
ax10.set_title('Range = {:5.2f}'.format(DS['range'][irange]));


DS.close()


# In[ ]:


ZLO = DS['ZLO'][:,:,:];
ZHI = DS['ZHI'][:,:,:];


ITX = DS['ITX'];
QTX = DS['QTX'];

IRX = DS['IRX'];
QRX = DS['QRX'];

ZED = ZLO.copy();
ZED[ZLO>-10] = ZHI[ZLO>-10];

dBZcal = DS['radar_constant'][:]-DS['receiver_gain'][:]+DS['cable_losses'][:];

ZED += dBZcal;


#inan = np.argwhere(np.isnan(ZLO));
#ZED[inan] = ZHI[inan];

IQ_TX_amp = np.sqrt(ITX[:,:]**2+QTX[:,:]**2);

ITXnorm = ITX[:,:]/IQ_TX_amp;
QTXnorm = QTX[:,:]/IQ_TX_amp;

IQ_RX_amp = np.sqrt(IRX[:,:,:]**2+QRX[:,:,:]**2);

IRXnorm = IRX[:,:,:]/IQ_RX_amp;
QRXnorm = QRX[:,:,:]/IQ_RX_amp;

fig = plt.figure(figsize=(12,30))
fig.set_constrained_layout(True)
fig.set_constrained_layout_pads(w_pad=2 / 72, h_pad=2 / 72, hspace=0.2,wspace=0.2)

gs0 = fig.add_gridspec(6, 2)
ax1 = fig.add_subplot(gs0[0, :2]);
ax2 = fig.add_subplot(gs0[2, 0]);
ax3 = fig.add_subplot(gs0[2, 1]);
ax4 = fig.add_subplot(gs0[3, 0]);
ax5 = fig.add_subplot(gs0[3, 1]);
ax6 = fig.add_subplot(gs0[4, 0]);
ax7 = fig.add_subplot(gs0[5, 1]);
ax8 = fig.add_subplot(gs0[1, :2]);

h1=ax1.pcolormesh(dtime,DS['range'],ZLO[:,1,:].T,vmin=-70,vmax=50,cmap='pyart_HomeyerRainbow');
cb1=plt.colorbar(h1,ax=ax1,orientation='vertical');
cb1.ax.set_ylabel("Reflectivity Factor (dBZ)");
myFmt = mdates.DateFormatter('%H:%M')
ax1.xaxis.set_major_formatter(myFmt);

h8=ax8.pcolormesh(dtime,DS['range'],ZED[:,1,:].T,vmin=-70,vmax=50,cmap='pyart_HomeyerRainbow');
cb8=plt.colorbar(h8,ax=ax8,orientation='vertical');
cb8.ax.set_ylabel("Reflectivity Factor (dBZ)");
myFmt = mdates.DateFormatter('%H:%M')
ax8.xaxis.set_major_formatter(myFmt);

ax2.plot(ITX,QTX,'.');
ax2.axis('equal');
ax2.set_xlabel('ITX');
ax2.set_ylabel('QTX');
ax2.grid(True);

irange = 5;
ax3.plot(IRX[...,irange],QRX[...,irange],'.');
ax3.axis('equal');
ax3.set_xlabel('IRX');
ax3.set_ylabel('QRX');
ax3.grid(True);
ax3.set_title('Range = {:5.2f}'.format(DS['range'][irange]));

irange = 10;
ax4.plot(IRX[...,irange],QRX[...,irange],'.');
ax4.axis('equal');
ax4.set_xlabel('IRX');
ax4.set_ylabel('QRX');
ax4.grid(True);
ax4.set_title('Range = {:5.2f}'.format(DS['range'][irange]));

irange = 55;
ax5.plot(IRX[...,irange],QRX[...,irange],'.');
ax5.axis('equal');
ax5.set_xlabel('IRX');
ax5.set_ylabel('QRX');
ax5.grid(True);
ax5.set_title('Range = {:5.2f}'.format(DS['range'][irange]));


# In[ ]:


#ax[1].plot(np.arange(pre_tx,tx_pulse),DS['QTX'][1,ipulse,pre_tx:tx_pulse]);
#ax[1].plot(np.arange(tx_pulse,post_tx),DS['QTX'][1,ipulse,tx_pulse:post_tx]);
ax[1].pcolormesh(DS['QTX'][1,:,0:8].T);

h2= ax[2].pcolormesh(ITX_bias[:,:]);
cb2=plt.colorbar(h2,ax=ax[2],orientation='vertical')


fig, ax = plt.subplots(3,3,figsize=(24,24),constrained_layout=True)
fig.set_constrained_layout_pads(w_pad=2 / 72, h_pad=2 / 72, hspace=0.2,wspace=0.2)

IQ_TX_amp = np.sqrt(ITXpulse**2+QTXpulse**2);

ITXnorm = ITXpulse/IQ_TX_amp;
QTXnorm = QTXpulse/IQ_TX_amp;


TXphase1 = 180./np.pi*np.arctan2(QTX[1,:,:]-QTX_mean_bias,ITX[1,:,:]-ITX_mean_bias);
TXphase2 = 180./np.pi*np.arctan2(QTX[1,:,:]-2047.,ITX[1,:,:]-2047.);


h0=ax[0,0].pcolormesh(IRX[1,:,:].T);
cb0 = plt.colorbar(h0,ax=ax[0],orientation='vertical')

ax[0,1].plot(ITXpulse[1,:],QTXpulse[1,:],'.');

print(ITXpulse.shape);
print(QTXpulse.shape);

#ax[1].plot(np.arange(tx_pulse,post_tx),DS['QTX'][1,ipulse,tx_pulse:post_tx]);
for pulse in np.arange(6100): 
ax[2,1].plot(ITX[1,pulse,7]-ITX_mean_bias,QTX[1,pulse,7]-QTX_mean_bias,'.y');
ax[2,1].plot(ITX[1,pulse,7]-ITX_bias_by_gate[1,7],QTX[1,pulse,7]-QTX_bias_by_gate[1,7],'.m');

ax[1,2].plot(ITX[1,pulse,:8]-ITX_bias_by_gate[1,:8],QTX[1,pulse,:8]-QTX_bias_by_gate[1,:8],'.r');
ax[1,2].plot(ITX[1,pulse,7]-ITX_mean_bias,QTX[1,pulse,7]-QTX_mean_bias,'.y');
ax[1,2].plot(ITX[1,pulse,7]-ITX_bias_by_gate[1,7],QTX[1,pulse,7]-QTX_bias_by_gate[1,7],'.m');

ax[1,1].plot(ITX_new[1,pulse,:],QTX_new[1,pulse,:],'.b');
ax[2,2].plot(ITXnorm[1,pulse],QTXnorm[1,pulse],'.b');

#ax[1,1].plot(ITX[1,pulse,8:]-ITX_mean_bias,QTX[1,pulse,8:]-QTX_mean_bias,'.g');
#ax[1,2].plot(ITX[1,pulse,8]-ITX_bias_by_gate[1,8],QTX[1,pulse,8]-QTX_bias_by_gate[1,8],'.y');

#ax[2].plot(TXphase1[pulse,30:],TXphase2[pulse,30:],'.');
#ax[2].legend(['TXphase1','TXphase2']);



ax[0,1].grid(True);
ax[0,1].axis('equal');
ax[0,2].grid(True);
ax[0,2].axis('equal');
ax[1,0].grid(True);
ax[1,0].axis('equal');
ax[1,1].grid(True);
ax[1,1].axis('equal');
ax[1,2].grid(True);
ax[1,2].axis('equal');
ax[2,2].grid(True);
ax[2,2].axis('equal');

#ax[0].set_xlim(0,30);
#ax[1].set_xlim(0,30);
#ax[0].set_ylim(0,4095);
#ax[1].set_ylim(0,4095);


# In[ ]:


x = np.array([0,1,2,3,4,5,6,7,8,9]);

print(x[:8]);
print(x[8:]);


# In[ ]:


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
    


# In[ ]:





# In[ ]:


testfile = '/home/users/cjwalden/git/test3ghz_L0b.nc';
DStest = nc4.Dataset(testfile);
print(DStest);


# In[ ]:


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

h1 = ax[0].pcolormesh(DS0['datetime'],gate_starts,DS0['Zpp'][1:,1:].transpose(),vmin=-20,vmax=50,cmap='pyart_HomeyerRainbow');
cb1=plt.colorbar(h1,ax=ax[0],orientation='vertical');
cb1.ax.set_ylabel("Z (dBZ)");
titlestr = "Chilbolton 3GHz CAMRa Radar";
ax[0].set_title(titlestr,loc='left')
ax[0].set_title('{}'.format(dt0[0].date()),loc='right')

ax[1].xaxis.set_major_formatter(myFmt);
h2 = ax[1].pcolormesh(DS0['datetime'],gate_starts,DS0['velpp'][1:,1:].transpose(),vmin=-vf,vmax=vf,cmap=cmocean.cm.balance);
cb2=plt.colorbar(h2,ax=ax[1],orientation='vertical');
cb2.ax.set_ylabel("VEL (m/s)");


for file in tsfiles_3ghz[istartfile+1:]:
    
    ifile = ifile+1;
    
    print('{}/{} {}'.format(ifile,nfiles,file));

    DS0 = load_3GHz(file,dBZh_offset,ZDR_offset,range_offset);
    gate_width = DS0['range'][1]-DS0['range'][0];
    gate_starts = DS0['range'][:]-gate_width/2.0;
    
    ax[0].pcolormesh(DS0['datetime'],gate_starts,DS0['Zpp'][1:,1:].transpose(),vmin=-20,vmax=50,cmap='pyart_HomeyerRainbow');
    ax[1].pcolormesh(DS0['datetime'],gate_starts,DS0['velpp'][1:,1:].transpose(),vmin=-vf,vmax=vf,cmap=cmocean.cm.balance);
    

ax[0].grid(True);
ax[1].grid(True);

figfile = "wivern2-radar-3ghz_{}.png".format(datestr);

figpath = '/home/users/cjwalden/jupyter/wivern2'

plt.savefig(os.path.join(figpath,figfile),dpi=200)


# In[ ]:


datestr = '20201114';

raw_ts_root_3ghz = '/gws/nopw/j04/ncas_obs/cao/raw_data/ncas-radar-camra-1/data/campaign/wivern2/ts/'
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


for file in tsfiles_3ghz[istartfile+1:iendfile]:
    
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


# In[ ]:


testfile = 'radar-camra_20200924023440_fix-ts-b.nc'
testpath = '/gws/nopw/j04/ncas_obs/cao/raw_data/ncas-radar-camra-1/data/campaign/wivern2/ts/ESA/20200924'
testfile=os.path.join(testpath,testfile);
DStest = nc4.Dataset(testfile);

dt = cftime.num2pydate(DStest['time'][:],DStest['time'].units);

print(dt);

print(DStest)
fig, ax = plt.subplots(2,1,figsize=(12,12),constrained_layout=True)
fig.set_constrained_layout_pads(w_pad=2 / 72, h_pad=2 / 72, hspace=0.2,wspace=0.2)
#myFmt = mdates.DateFormatter('%H:%M')

ax[0].plot(DStest['ZLO'][-5,3049,:],'.');
#ax[0].plot(DS2['ZHI'][68,:,-1]);


# In[ ]:


fig, ax = plt.subplots(4,1,figsize=(12,12),constrained_layout=True)
fig.set_constrained_layout_pads(w_pad=2 / 72, h_pad=2 / 72, hspace=0.2,wspace=0.2)
myFmt = mdates.DateFormatter('%H:%M')

ax[0].plot(ts_out['datetime'],ts_out['Zpp'][:,100]);
ax[0].xaxis.set_major_formatter(myFmt);

ax[1].xaxis.set_major_formatter(myFmt);

h1 = ax[1].pcolor(ts_out['datetime'],ts_out['range'],ts_out['Zpp'][:,:].transpose(),cmap='pyart_HomeyerRainbow');
cb1=plt.colorbar(h1,ax=ax[1],orientation='vertical');
cb1.ax.set_ylabel("ZED_HC (dB)");

ax[2].xaxis.set_major_formatter(myFmt);

h2 = ax[2].pcolor(ts_out['datetime'],ts_out['range'],ts_out['velpp'][:,:].transpose(),cmap='pyart_HomeyerRainbow');
cb2=plt.colorbar(h2,ax=ax[2],orientation='vertical');
cb2.ax.set_ylabel("VEL (m/s)");


# In[ ]:


time=double(ncread(infile,'time'));
num_rays=length(time);
I=zeros(ngates,npulses,num_rays); Q=zeros(ngates,npulses,num_rays);
  i1=zeros(ngates,npulses); q1=zeros(ngates,npulses);
  i2=zeros(ngates,npulses); q2=zeros(ngates,npulses);    
  zlo=zeros(ngates,npulses); zhi=zeros(ngates,npulses);
  
  
  ZLO=double(ncread(infile,'ZLO')); ZHI=double(ncread(infile,'ZHI'));
  I1=double(ncread(infile,'ITX')); Q1=double(ncread(infile,'QTX'));
  I2=double(ncread(infile,'IRX')); Q2=double(ncread(infile,'QRX'));

--


---
for raynum=1:num_rays 

  
    zlo=ZLO(:,:,raynum); zhi=ZHI(:,:,raynum);
    i1=I1(:,:,raynum); i2=I2(:,:,raynum);
    q1=Q1(:,:,raynum); q2=Q2(:,:,raynum);
    
  for ir=1:size(I,1),
    i1(ir,:)=i1(ir,:)-mean(i1(ir,:)); i2(ir,:)=i2(ir,:)-mean(i2(ir,:));
    q1(ir,:)=q1(ir,:)-mean(q1(ir,:)); q2(ir,:)=q2(ir,:)-mean(q2(ir,:));             
  end
  Vtx=complex(i1,-q1); Vrx=complex(i2,q2);
  V=Vrx.*Vtx; V=V./abs(V);
  V=conj(V);
  V(:,2:2:end)=V(:,2:2:end).*-1;
  
  get_ipython().run_line_magic('Merge', 'power measurements from the two amplifiers (ZLO preferred')
  get_ipython().run_line_magic('unless', 'it is close to saturation, then use ZHI)')
  Z=nan(size(zlo));
  index=find(zlo<=zlothresh);  Z(index)=zlotable(zlo(index)+1);
  index=find(zlo>zlothresh); Z(index)=zhitable(zhi(index)+1);
  
  get_ipython().run_line_magic('Alternating', 'polarisation (odd pulses H, even pulses V)')
  Z(:,1:2:end)=Z(:,1:2:end).*Zh_cal;  
  Z(:,2:2:end)=Z(:,2:2:end).*Zv_cal;
  
  V=V.*sqrt(Z);
  I(:,:,raynum)=real(V); 
  Q(:,:,raynum)=imag(V);
  
end

I=squeeze(I(:,:,1:raynum));
Q=squeeze(Q(:,:,1:raynum));
  
Zpp=nan(size(I,1),size(I,3));
velpp=nan(size(I,1),size(I,3));


PRF=305; %Hz
freq=3e9; %Hz

c=299792458; %m/s
lambda=c/freq;
vf=lambda*PRF/4;

for it=1:raynum
    for ir=1:size(I,1)
        Zpp(ir,it)=mean(I(ir,:,it).^2+Q(ir,:,it).^2);
        real_vh=0; imag_vh=0;
        real_vv=0; imag_vv=0;
        for ip=3:2:round(size(I,2)./2)
           real_vh=real_vh+I(ir,ip,it).*I(ir,ip-2,it)+Q(ir,ip,it).*Q(ir,ip-2,it); 
           imag_vh=imag_vh+Q(ir,ip,it).*I(ir,ip-2,it)-I(ir,ip,it).*Q(ir,ip-2,it);
     
           real_vv=real_vv+I(ir,ip+1,it).*I(ir,ip-1,it)+Q(ir,ip+1,it).*Q(ir,ip-1,it); 
           imag_vv=imag_vv+Q(ir,ip+1,it).*I(ir,ip-1,it)-I(ir,ip+1,it).*Q(ir,ip-1,it);              
        end
        velpp(ir,it)=atan2(imag_vh+imag_vv,real_vh+real_vv).*vf./pi;
    end
end
Zpp=log10(Zpp).*10+log10((range.^2)*ones(1,raynum)).*10;

blindrange=15;    

Zpp(1:blindrange,:)=NaN;
velpp(1:blindrange,:)=NaN;

figure
pcolor(time./(24*3600),range,Zpp)
shading flat
colormap(jet)
colorbar 
caxis([-20 50])
ylabel('range (km)')
xlabel('Time (UTC)')
set(gca,'FontSize',16)
set(gca,'LineWidth',2)
datetick('x','keeplimits')

figure
pcolor(time./(24*3600),range,velpp)
shading flat
colormap(hsv)
colorbar
ylim([0 10])
caxis([-vf vf])
ylabel('range (km)')
xlabel('Time (UTC)')
set(gca,'FontSize',16)
set(gca,'LineWidth',2)
datetick('x','keeplimits')    
   
return Zpp,velpp,I,Q,time,range


# In[ ]:


def load_camra_ts(infile):

    # Author: Chris Walden (chris.walden@ncas.ac.uk)
    # Based on MATLAB code written by John Nicol and Chris Westbrook

    Zcal         = -60; # dBZ calibration of stored IQ data  (CJW: ?because range in km)
    Zscalefact   = np.power(10,(-Zcal/10.));
    
    cable_losses = 4.8;
    radar_const  = 64.7;
    rec_gain     = 45.5;
    range_offset = -24.0; # metres
    adcmax       = 4096;
    zlothresh    = 3840;
    zlomin       = -70;
    zhimin       = -38;
    zcxmin       = -77;
    zloscale     = 0.015625;
    zhiscale     = 0.015625;
    zlotable     = np.power(10,((zlomin+np.arange(adcmax)*zloscale)/10.));
    zhitable     = np.power(10,((zhimin+np.arange(adcmax)*zhiscale)/10.));
    
    dBZh_offset = 7.3;
    ZDR_offset  = 0.62;
    dBZv_offset = dBZh_offset-ZDR_offset;
    
    dBZcal      = radar_const-rec_gain+cable_losses;
    Zh_cal      = np.power(10,((dBZcal+dBZh_offset)/10.));
    Zv_cal      = np.power(10,((dBZcal+dBZv_offset)/10.));

    tsfile = os.path.join(tspath,infile);

    print(tsfile);
   
    ncfile = tsfile.replace("ts","raw");
    print(ncfile);

    ncfile = tsfile.replace("ts","raw");
    if ncfile.endswith('.dat'):
        ncfile = ncfile.replace("dat","nc");
    print(ncfile);

    DS = nc4.Dataset(ncfile)
    
    camra_range = dataset['range'][:]+range_offset)/1000.;
    print(dataset.variables.keys())
    #print(camra_range[:]);

    # ---------------------------------------------------------------------
    # Read in I-Q time series
    # ---------------------------------------------------------------------
    if tsfile.endswith('.dat'):
        ts_dict = read_camra_ts_ascii(tsfile);
    else:
        ts_dict = read_camra_ts_netcdf(tsfile);
    
    npulse   = ts_dict['npulse'];
    nsample  = ts_dict['nsample'];
    nray     = ts_dict['nray'];
    dt       = ts_dict['datetime'];
    unixtime = ts_dict['unixtime'];
    
    ITX = ts_dict['ITX']; QTX = ts_dict['QTX'];
    IRX = ts_dict['IRX']; QRX = ts_dict['QRX'];
    print(npulse); print(nsample); print(nray); 
    
    ngate = camra_range.shape[0];
    npulse_new = int(np.floor(npulse*nsample/ngate));
    newsize = int(npulse_new*ngate);
    
    tmp0 = ITX.reshape(nray,-1);
    tmp1 = tmp0[:,0:newsize];
    ITXnew = tmp1.reshape(nray,npulse_new,-1);
    
    tmp0 = QTX.reshape(nray,-1);
    tmp1 = tmp0[:,0:newsize];
    QTXnew = tmp1.reshape(nray,npulse_new,-1);
    
    tmp0 = IRX.reshape(nray,-1);
    tmp1 = tmp0[:,0:newsize];
    IRXnew = tmp1.reshape(nray,npulse_new,-1);
    
    tmp0 = QRX.reshape(nray,-1);
    tmp1 = tmp0[:,0:newsize];
    QRXnew = tmp1.reshape(nray,npulse_new,-1);    
        
    ITXmean = np.mean(ITXnew,axis=1);
    QTXmean = np.mean(QTXnew,axis=1);
    IRXmean = np.mean(IRXnew,axis=1);
    QRXmean = np.mean(QRXnew,axis=1);
    
    print(ITXmean);
    print(ITXnew[:,:,0]);
    #for iray in np.arange(nray):
    #    for isample in np.arange(nsample):
    #        ITX[iray,:,isample] = ITX[iray,:,isample]-ITXmean[iray,isample];
    #        QTX[iray,:,isample] = QTX[iray,:,isample]-QTXmean[iray,isample];
    #        IRX[iray,:,isample] = IRX[iray,:,isample]-IRXmean[iray,isample];
    #        QRX[iray,:,isample] = QRX[iray,:,isample]-QRXmean[iray,isample];
    
    for iray in np.arange(nray):
        for isample in np.arange(nsample):
            ITXnew[iray,:,isample] = ITXnew[iray,:,isample]-ITXmean[iray,isample];
            QTXnew[iray,:,isample] = QTXnew[iray,:,isample]-QTXmean[iray,isample];
            IRXnew[iray,:,isample] = IRXnew[iray,:,isample]-IRXmean[iray,isample];
            QRXnew[iray,:,isample] = QRXnew[iray,:,isample]-QRXmean[iray,isample];
    
    Vtx = ITXnew - 1j*QTXnew;
    Vrx = IRXnew + 1j*QRXnew;
 

    V = np.multiply(Vrx, Vtx);
    V = np.divide(np.conj(V),np.abs(V));
    V[:,1::2,:]=V[:,1::2,:]*-1.;
    
    zlo_int  = ts_dict['ZLO'].astype(int);
    tmp0 = zlo_int.reshape(nray,-1);
    tmp1 = tmp0[:,0:newsize];
    zlo_new  = tmp1.reshape(nray,npulse_new,-1);
    zhi_int  = ts_dict['ZHI'].astype(int);

    tmp0 = zhi_int.reshape(nray,-1);
    tmp1 = tmp0[:,0:newsize];
    zhi_new  = tmp1.reshape(nray,npulse_new,-1);
    
    #ZED          = ts_dict['ZLO'].astype(float);
    ZED          = zlo_new;
    index        = zlo_new<=zlothresh;  
    ZED[index]   = zlotable[zlo_new[index]];
    index        = zlo_new>zlothresh;
    ZED[index]   = zhitable[zhi_new[index]];
    ZED[:, ::2,:] = ZED[:, ::2,:]*Zh_cal;
    ZED[:,1::2,:] = ZED[:,1::2,:]*Zv_cal;

    ZED = ZED*Zscalefact;
    V = V*np.sqrt(ZED);
    
    ts_out = dict();
    ts_out['Itx'] = np.real(Vtx);
    ts_out['Qtx'] = np.imag(Vtx);
    ts_out['I'] = np.real(V);
    ts_out['Q'] = np.imag(V);
    ts_out['ZED'] = ZED;
    ts_out['range'] = camra_range;
    ts_out['datetime'] = dt;
    ts_out['unixtime'] = unixtime;
    ts_out['Zcal'] = Zcal;
    
    return ts_out

