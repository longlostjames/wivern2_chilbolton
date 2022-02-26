#!/usr/bin/env python
# coding: utf-8

# ==========================================================================
# Module for processing raw radar data from CAMRa (3GHz), Copernicus (35GHz)
# and Galileo (94GHz) radars at Chilbolton, collected as part of the
# ESA WIVERN-2 project in 2020-2021.
# Author: Chris Walden, UK Research & Innovation and
#                       National Centre for Atmospheric Science
# Last modified: 04-02-2022
# ==========================================================================

"""Module for processing raw radar data from CAMRa (3GHz), Copernicus (35GHz)
and Galieo (94GHz) radars at Chilbolton, collected as part of the
ESA WIVERN-2 project in 2020-2021."""

module_version = 1.0;

import glob
import os
import getpass
import socket

import fnmatch

import numpy as np;
import numpy.ma as ma

import netCDF4 as nc4
import cftime
from datetime import tzinfo, datetime, time

import pyart
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
from matplotlib import colors
import cmocean




# ===================
# CONVERSION ROUTINES
# ===================

def convert_camra_ts_l0a2l0b(infile,outfile):

    """This routine converts raw (Level 0a) time series data from the Chilbolton Advanced Meteorological Radar (CAMRa) to Level 0b data.
    Processing involves removing redundant dimensions, and removing bias from the ADC samples of transmit and receive I and Q.
    Single estimates of I and Q for each transmit pulse are produced and stored in the output file.
    Metadata are added specific to the WIVERN-2 project ground-based observations collected in 2020-2021.

    :param infile: Full path of NetCDF Level 0a raw data file, e.g. `<path-to-file>/radar-camra_20201210212823_fix-ts.nc`
    :type infile: str

    :param outfile: Full path of NetCDF Level 0b output file, e.g. `<path-to-file>/ncas-radar-camra-1_cao_20201210-212823_fix-ts_l0b_v1.0.nc`
    :type outfile: str
    """

    DSin = nc4.Dataset(infile);

    dt_start = cftime.num2pydate(DSin['time'][0],DSin['time'].units)
    dt_end   = cftime.num2pydate(DSin['time'][-1],DSin['time'].units)

    try: outfile.close()  # just to be safe, make sure dataset is not already open.
    except: pass
    DSout = nc4.Dataset(outfile,mode='w',format='NETCDF4')
    print("Creating {}".format(outfile));

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

    DSout.ADC_bits_per_sample = np.int(DSin.ADC_bits_per_sample);
    DSout.ADC_channels        = np.int(DSin.ADC_channels);

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
    varout[:]=146.7;

    varout = DSout.createVariable('altitude_agl',varin.datatype);
    varout.standard_name = 'altitude';
    varout.long_name = 'altitude of the elevation axis above ground';
    varout.units = 'm';
    varout[:]=16.0;

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

    varout = DSout.createVariable('antenna_focal_length','f4');
    varout.long_name = 'focal length of antenna';
    varout.units = 'm';
    varout[:] = 9.0;

    varout = DSout.createVariable('antenna_focus_radial_location','f4');
    varout.long_name = 'distance along boresight from elevation axis to antenna focus';
    varout.units = 'm';
    varout[:] = 14.34;

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
    varout.range_offset_applied = np.float32(varin.range_offset);
    varout.units = varin.units;
    varout[:]=varin[:];

    # --------------------------
    # Antenna pointing variables
    # --------------------------
    varin = DSin['elevation'];
    varout = DSout.createVariable('elevation',varin.datatype,'time');
    varout.long_name = 'elevation angle of the antenna boresight above the horizon';
    varout.units = 'degree';
    varout.elevation_offset_applied = np.float32(0.);
    varout[:] = varin[:];

    varin = DSin['azimuth'];
    varout = DSout.createVariable('azimuth',varin.datatype,'time');
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
#    add_offset = np.min(ITXout[:]);
#    scale_factor = (np.max(ITXout[:])-np.min(ITXout[:])) / (2**16-1)
#    packed_data = np.int16(np.rint((ITXout[:,:] - add_offset)/scale_factor));
#    varout.scale_factor = np.float32(scale_factor);
#    varout.add_offset = np.float32(add_offset);
    varout.long_name = 'bias-corrected samples of in-phase video signal for each transmitted pulse';
    varout.units = '1';
#    varout[:] = packed_data;
    varout[:]=ITXout;

    varout = DSout.createVariable('QTX','f4',('time','pulse'),zlib=True);
#    add_offset = np.min(QTXout[:]);
#    scale_factor = (np.max(QTXout[:])-np.min(QTXout[:])) / (2**16-1)
#    packed_data = np.int16(np.rint((QTXout[:,:] - add_offset)/scale_factor));
#    varout.scale_factor = np.float32(scale_factor);
#    varout.add_offset = np.float32(add_offset);
    varout.long_name = 'bias-corrected samples of quadrature video signal for each transmitted pulse';
    varout.units = '1';
#    varout[:] = packed_data;
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
#    add_offset = np.min(IRXout[:]);
#    scale_factor = (np.max(IRXout[:])-np.min(IRXout[:])) / (2**16-1)
#    packed_data = np.int16(np.rint((IRXout[:,:,:] - add_offset)/scale_factor));
#    varout.scale_factor = np.float32(scale_factor);
#    varout.add_offset = np.float32(add_offset);
    varout.long_name = 'bias-corrected samples of received in-phase video signal at output of IF limiting amplifier';
    varout.units = '1';
#    varout[:] = packed_data;
    varout[:]=IRXout;

    varin = DSin['QRX'];
    varout = DSout.createVariable('QRX','f4',('time','pulse','range'),zlib=True);
#    add_offset = np.min(QRXout[:]);
#    scale_factor = (np.max(QRXout[:])-np.min(QRXout[:])) / (2**16-1)
#    packed_data = np.int16(np.rint((QRXout[:,:,:] - add_offset)/scale_factor));
#    varout.scale_factor = np.float32(scale_factor);
#    varout.add_offset = np.float32(add_offset);
    varout.long_name = 'bias-corrected samples of received quadrature video signal at output of IF limiting amplifier';
    varout.units = '1';
#    varout[:] = packed_data;
    varout[:]=QRXout;

    # -----------------------
    # Update history metadata
    # -----------------------
    user = getpass.getuser()

    updttime = datetime.utcnow()
    updttimestr = updttime.ctime()

    history = updttimestr + (" - user:" + user
    + " machine: " + socket.gethostname()
    + " program: wivern_chilbolton_utils.py convert_camra_ts_l0a2l0b"
    + " version:" + str(module_version));

    DSout.history = history + "\n" + DSin.history;

    DSout.last_revised_date = datetime.strftime(updttime,'%Y-%m-%dT%H:%M:%SZ')

    DSin.close();
    DSout.close();

    return

def convert_camra_ts_l0b2l1(infile,outfile,dBZh_offset,ZDR_offset,range_offset,version):

    """This routine converts Level 0b time series data from the Chilbolton Advanced Meteorological Radar (CAMRa) to Level 1 data.

    Processing is applied to produce I and Q time series for received data at each range gate.
    I and Q values from the limiting amplifiers in the receiver are compared with transmit I and Q for each pulse.
    These are then scaled to the unit circle and multiplied the square root of the sampled linear reflectivity.
    Processing includes separate calibration offsets for even (H-polarized) and odd (V-polarized) pulses.
    Metadata are added specific to the WIVERN-2 project ground-based observations collected in 2020-2021.

    :param infile: Full path of NetCDF Level 0b input file, e.g. `<path-to-file>/ncas-radar-camra-1_cao_20201210-212823_fix-ts_l0b_v1.0.nc`
    :type infile: str
    :param outfile: Full path of NetCDF Level 1 output file, e.g. `<path-to-file>/ncas-radar-camra-1_cao_20201210-212823_fix-ts_l1_v1.0.nc`
    :type outfile: str
    :param dBZh_offset: Calibration offset to be applied to H polarized reflectivity factor
    :type dBZh_offset: float
    :param ZDR_offset: Calibration offset that would be applied to differential reflectivity (used to calculate the calibration offset to apply to V polarized reflectivity factor)
    :type ZDR_offset: float
    :param range_offset: Additional range offset in metres to be applied
    :type range_offset: float
    :param version: Version of data product in the format `N.m`, where `N` denotes the major verion and `m` a minor revision.
    :type version: str
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
                # copy variable attributes all at once via dictionary
                dst[name].setncatts(src[name].__dict__)

                dst[name][:] = src[name][:]


    try: outfile.close()  # just to be safe, make sure dataset is not already open.
    except: pass
    DSout = nc4.Dataset(outfile,mode='r+',format='NETCDF4')
    print(outfile)

    DSin = nc4.Dataset(infile);
    dt = cftime.num2pydate(DSin['time'][:],DSin['time'].units);

    DSout.product_version = "v{}".format(version) ;
    DSout.processing_level = "1" ;

    DSout.title = "Calibrated time series from 3 GHz CAMRa radar collected for ESA WIVERN-2 campaign at Chilbolton Observatory";

    comment_string = "Correction to account for inverse square power loss with range has not been applied";
    if len(DSin.comment)>0:
        DSout.comment = DSin.comment + "\n " + comment_string;
    else:
        DSout.comment = comment_string;

    varout = DSout.createVariable('dBZ_offsets_applied','f4',('pulse'),zlib=True);
    varout.long_name = 'dBZ calibration offsets applied for even and odd pulses';
    varout.units = 'dB';
    varout.comment = 'dBZ offsets for even pulses (H-polarized) and odd pulses (V-polarized)';

    varout = DSout.createVariable('I','f4',('time','pulse','range'),zlib=True);
    varout.ancillary_variables = 'qc_flag';
    varout.long_name = 'co-polar in-phase video signal';
    varout.units = '1';
    varout.comment = 'Values are derived from I/Q on unit circle multiplied by square root of linear reflectivity factor';

    varout = DSout.createVariable('Q','f4',('time','pulse','range'),zlib=True);
    varout.ancillary_variables = 'qc_flag';
    varout.long_name = 'co-polar quadrature video signal';
    varout.units = '1';
    varout.comment = 'Values are derived from I/Q on unit circle multiplied by square root of linear reflectivity factor';

    varout = DSout.createVariable('ZCX','f4',('time','pulse','range'),zlib=True);
    varout.ancillary_variables = 'qc_flag';
    varout.long_name = 'cross-polar radar equivalent reflectivity factor';
    varout.units = 'dBZ';
    varout.comment = '';

    varout = DSout.createVariable('qc_flag','u1',('time','pulse','range'),fill_value=255);
    varout.is_quality = 'true';
    varout.qualified_variables = 'I Q ZCX';
    varout.long_name = 'Quality control flag';
    varout.flag_values = np.uint8(0),np.uint8(1), np.uint8(2), np.uint8(3), np.uint8(4), np.uint8(255);
    varout.flag_meanings = 'not_used good_data probably_good_data bad_data data_in_blind_range no_qc_performed';
    varout[:] = 2;

    cable_losses = DSout['cable_losses'][:];
    radar_const  = DSout['radar_constant'][:];
    rec_gain     = DSout['receiver_gain'][:];
    freq         = DSout['frequency'][:];
    prf          = DSout['prf'][:];

    dBZv_offset = dBZh_offset-ZDR_offset;

    dBZcal = radar_const-rec_gain+cable_losses;

    DSout['dBZ_offsets_applied'][::2]  = dBZh_offset;
    DSout['dBZ_offsets_applied'][1::2] = dBZv_offset;

    Zh_cal = 10**((dBZcal+dBZh_offset)/10);
    Zv_cal = 10**((dBZcal+dBZv_offset)/10);

    range_m  = DSin['range'][:];
    range_km = (range_m+range_offset)/1000.; # range in km

    ITX = DSin['ITX'][:,:]; #*DSin['ITX'].scale_factor+DSin['ITX'].add_offset;
    QTX = DSin['QTX'][:,:]; #*DSin['QTX'].scale_factor+DSin['QTX'].add_offset;
    IRX = DSin['IRX'][:,:,:]; #*DSin['IRX'].scale_factor+DSin['IRX'].add_offset;
    QRX = DSin['QRX'][:,:,:]; #*DSin['QRX'].scale_factor+DSin['QRX'].add_offset;

    Vtx = ITX - 1j*QTX;
    Vrx = IRX + 1j*QRX;

    V = np.multiply(Vrx[:,:,:], Vtx[:,:,None]);
    V = ma.masked_where(V==0,V);
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
 #   add_offset = np.min(ZCX[:]);
 #   scale_factor = (np.max(ZCX[:])-np.min(ZCX[:])) / (2**16-1)
 #   packed_data = np.int16(np.rint((ZCX[:,:,:] - add_offset)/scale_factor));
 #   DSout['ZCX'].scale_factor = np.float32(scale_factor);
 #   DSout['ZCX'].add_offset = np.float32(add_offset);
 #   DSout['ZCX'][:] = packed_data;
    DSout['ZCX'][:] = ZCX;

    V = V*np.sqrt(ZED);

    I = V.real;
    Q = V.imag;

#    add_offset = np.min(I[:]);
#    scale_factor = (np.max(I[:])-np.min(I[:])) / (2**16-1)
#    packed_data = np.int16(np.rint((I[:,:,:] - add_offset)/scale_factor));
#    DSout['I'].scale_factor = np.float32(scale_factor);
#    DSout['I'].add_offset = np.float32(add_offset);
    DSout['I'][:] = I;

#    add_offset = np.min(Q[:]);
#    scale_factor = (np.max(Q[:])-np.min(Q[:])) / (2**16-1)
#    packed_data = np.int16(np.rint((Q[:,:,:] - add_offset)/scale_factor));
#    DSout['Q'].scale_factor = np.float32(scale_factor);
#    DSout['Q'].add_offset = np.float32(add_offset);
    DSout['Q'][:] = Q

    blind_range = np.arange(15);
    DSout['qc_flag'][:,:,blind_range] = 4;

    DSout['range'][:] = range_m;
    DSout['range'].range_offset_applied += range_offset;
    DSout['range'].comment = "range_offset_applied includes offset applied by the data acquisition program";

    user = getpass.getuser();

    updttime = datetime.utcnow();
    updttimestr = updttime.ctime();

    history = updttimestr + (" - user:" + user
    + " machine: " + socket.gethostname()
    + " program: wivern_chilbolton_utils.py convert_camra_ts_l0b2l1"
    + " version:" + str(module_version));

    print(history);

    DSout.history = history + "\n" + DSin.history;

    DSin.close();
    DSout.close();

    return

def convert_copernicus_ts_l0a2l1(infile,outfile,dBZ_offset,range_offset,data_version):

    """This routine converts raw (Level 0a) time series data from the Chilbolton 35GHz Cloud Radar (Copernicus) to Level 1 data.
    Processing involves removing redundant dimensions, and removing bias from the ADC samples of the received I and Q.
    Metadata are added specific to the WIVERN-2 project ground-based observations collected in 2020-2021.

    :param infile: Full path of binary Level 0a raw data file, e.g. `<path-to-file>/2021061813445611_iqdata.bin`
    :type infile: str

    :param outfile: Full path of NetCDF Level 0b output file, e.g. `<path-to-file>/ncas-radar-ka-band-1_cao_20210618-134456_fix-ts_l1_v1.0.nc`
    :type outfile: str

    :param dBZ_offset: dB calibration offset to apply to reflectivity.  This is converted to linear units and I and Q values are multiplied by the square root of this value.
    :type dBZoffset: float

    :param range_offset: Range offset to apply in m.
    :type range_offset: float

    :param data_version: Version of data product in the format `N.m`, where `N` denotes the major verion and `m` a minor revision
    :type data_version: str
    """

    # -------------------------------
    # Read binary IQ time-series file
    # -------------------------------
    header_dt = np.dtype([('ngates','<i4'),('npulses','<i4')]);
    time_dt   = np.dtype([('year', '<i4'), ('month', '<i4'),('day','<i4'),('hour','<i4'),('minute','<i4'),('second','<i4'),('centisecond','<i4')]);

    A         = np.fromfile(infile, dtype=header_dt, count=1)[0];
    ngate     = A['ngates'];
    npulse    = A['npulses'];
    iq_dt     = np.dtype([('IQ',"({},{})<i2".format(npulse,ngate))]);
    record_dt = np.dtype([('time',time_dt),('I',iq_dt),('Q',iq_dt)]);

    records = np.fromfile(infile,dtype=record_dt, offset=8);
    times = np.array([datetime(time['year'],time['month'],time['day'],time['hour'],time['minute'],time['second'],time['centisecond']*10000) for time in records['time']]);
    nray = len(times);

    dt_start = times[0];
    dt_end   = times[-1];
    time_units_out = datetime.strftime(dt_start,'seconds since %Y-%m-%d 00:00:00Z');

    # -----------------------------
    # Set up NetCDF file for output
    # -----------------------------
    DSout = nc4.Dataset(outfile,mode='w',format='NETCDF4')

    # ------------------------
    # Set up global attributes
    # ------------------------
    DSout.product_version = "v{}".format(data_version);
    DSout.processing_level = "1" ;

    DSout.licence = "Data usage licence - UK Open Government Licence agreement: \n http://www.nationalarchives.gov.uk/doc/open-government-licence" ;
    DSout.acknowledgement = "Acknowledgement is required of UK Research and Innovation as the data provider (in partnership with the National Centre for Atmospheric Science) whenever and wherever these data are used." ;
    DSout.platform = "Chilbolton Atmospheric Observatory" ;
    DSout.platform_type = "stationary_platform" ;
    DSout.title = "Time series from 35 GHz Copernicus radar collected for ESA WIVERN-2 campaign at Chilbolton Observatory";
    DSout.creator_name = "Chris Walden" ;
    DSout.creator_email = "chris.walden@ncas.ac.uk" ;
    DSout.creator_url = "https://orcid.org/0000-0002-5718-466X" ;
    DSout.institution = "National Centre for Atmospheric Science (NCAS)";
    DSout.instrument_name = "ncas-radar-ka-band-1";
    DSout.instrument_software = "radar-copernicus-iq-rec" ;
    DSout.instrument_software_version = "0.1" ;

    DSout.references = "";
    DSout.source = "35GHz Copernicus Radar";
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

    DSout.ADC_bits_per_sample = np.int(12);
    DSout.ADC_channels        = np.int(8);

    # ----------------
    # Scalar variables
    # ----------------

    varout = DSout.createVariable('latitude','f4');
    varout.standard_name = 'latitude';
    varout.long_name = 'latitude of the antenna';
    varout.units = 'degree_north';
    varout[:]=51.1450;

    varout = DSout.createVariable('longitude','f4');
    varout.standard_name = 'longitude';
    varout.long_name = 'longitude of the antenna';
    varout.units = 'degree_east';
    varout[:]=-1.4384;

    varout = DSout.createVariable('altitude','f4');
    varout.standard_name = 'altitude';
    varout.long_name = 'altitude of the antenna above the geoid (WGS84)';
    varout.units = 'm';

    varout = DSout.createVariable('altitude_agl','f4');
    varout.long_name = 'altitude of the antenna above ground';
    varout.units = 'm';
    varout[:] = 1.9;

    varout = DSout.createVariable('frequency','f4');
    varout.standard_name = 'radiation_frequency';
    varout.long_name = 'frequency of transmitted radiation';
    varout.units = 'GHz';
    varout[:]=34.96;

    varout = DSout.createVariable('prf','f4');
    varout.long_name = 'pulse repetition frequency';
    varout.units = 'Hz';
    varout[:]=5000;

    varout = DSout.createVariable('beamwidthH','f4');
    varout.long_name = 'horizontal angular beamwidth';
    varout.units = 'degree';
    varout[:]=0.25;

    varout = DSout.createVariable('beamwidthV','f4');
    varout.long_name = 'vertical angular beamwidth';
    varout.units = 'degree';
    varout[:]=0.25;

    varout = DSout.createVariable('antenna_diameter','f4');
    varout.long_name = 'antenna diameter';
    varout.units = 'm';
    varout[:]=2.4;

    varout = DSout.createVariable('antenna_focal_length','f4');
    varout.long_name = 'focal length of antenna';
    varout.units = 'm';
    varout[:] = 0.75;

    varout = DSout.createVariable('pulse_width','f4');
    varout.long_name = 'pulse width';
    varout.units = 'us';
    varout[:]=0.4;

    varout = DSout.createVariable('transmit_power','f4');
    varout.long_name = 'peak transmitted power';
    varout.units = 'W';
    varout[:]=1500;

    varout = DSout.createVariable('clock','f4');
    varout.long_name = 'clock input to timer card';
    varout.units = 'Hz';
    varout[:]=10000000;

    varout = DSout.createVariable('clock_divfactor','i4');
    varout.long_name = 'clock divide factor';
    varout.units = '1';
    varout[:]=2;

    varout = DSout.createVariable('delay_clocks','i4');
    varout.long_name = 'clock cycles before sampling is initiated';
    varout.units = '1';
    varout[:]=3;

    varout = DSout.createVariable('samples_per_pulse','i4');
    varout.long_name = 'number of samples per pulse';
    varout.units = '1';
    varout[:] = ngate;

    varout = DSout.createVariable('pulses_per_daq_cycle','i4');
    varout.long_name = 'number of pulses per data acquisition cycle';
    varout.units = '1';

    varout = DSout.createVariable('pulses_per_ray','i4');
    varout.long_name = 'number of pulses per ray';
    varout.units = '1';
    varout[:] = npulse;

    varout = DSout.createVariable('radar_constant','f4');
    varout.long_name = 'radar constant';
    varout.units = 'dB';

    varout = DSout.createVariable('receiver_gain','f4');
    varout.long_name = 'receiver gain';
    varout.units = 'dB';

    varout = DSout.createVariable('cable_losses','f4');
    varout.long_name = 'cable losses';
    varout.units = 'dB';

    varout = DSout.createVariable('extra_attenuation','f4');
    varout.long_name = 'extra attenuation';
    varout.units = 'dB';
    varout[:] = 0.0;

    varout = DSout.createVariable('dBZ_offset','f4');
    varout.long_name = 'dBZ offset applied';
    varout.comment = 'When converted to linear units this is included in the linear reflectivity factor.  The square root of the latter is used to scale the I and Q values.';

    lightspeed = 299792458; # m s-1
    pulse_width_s = DSout['pulse_width'][:]*1e-6;
    gate_width = lightspeed*pulse_width_s/2.0;
    range_m    = np.arange(ngate)*gate_width;
    range_m += range_offset;

    # -----------------
    # Create dimensions
    # -----------------
    DSout.createDimension('time', None);
    DSout.createDimension('pulse', npulse);
    DSout.createDimension('range', ngate);

    # --------------------
    # Coordinate variables
    # --------------------
    varout = DSout.createVariable('time','f4',('time'));
    varout.standard_name = 'time';
    varout.long_name = 'time at the start of each recorded ray';
    varout.units = time_units_out;
    varout[:]=nc4.date2num(times,varout.units);

    varout = DSout.createVariable('range','f4',('range'));
    varout.long_name = 'distance from the antenna to the middle of each range gate';
    varout.range_offset_applied = np.float32(range_offset);
    varout.units = 'm';
    varout[:]=range_m;

    # --------------------------
    # Antenna pointing variables
    # --------------------------
    varout = DSout.createVariable('elevation','f4','time');
    varout.long_name = 'elevation angle of the antenna boresight above the horizon';
    varout.units = 'degree';
    varout.elevation_offset_applied = np.float32(0.);

    varout = DSout.createVariable('azimuth','f4','time');
    varout.long_name = 'azimuth angle clockwise from grid north of the plane perpendicular to axis of the elevation tilting mechanism';
    varout.units = 'degree';
    varout.azimuth_offset_applied = np.float32(0.);

    # ---------------
    # Field variables
    # ---------------
    varout = DSout.createVariable('I','f4',('time','pulse','range'),zlib=True);
    varout.ancillary_variables = 'qc_flag';
    varout.long_name = 'co-polar in-phase video signal';
    varout.units = '1';
    varout.comment = 'Scaled to account for calibration for square root of linear reflectivity factor';

    varout = DSout.createVariable('Q','f4',('time','pulse','range'),zlib=True);
    varout.ancillary_variables = 'qc_flag';
    varout.long_name = 'co-polar quadrature video signal';
    varout.units = '1';
    varout.comment = 'Scaled to account for calibraton for square root of linear reflectivity factor';

    varout = DSout.createVariable('qc_flag','u1',('time','pulse','range'));
    varout.is_quality = 'true';
    varout.qualified_variables = 'I Q';
    varout.long_name = 'Quality control flag';
    varout.flag_values = np.uint8(0),np.uint8(1), np.uint8(2), np.uint8(3), np.uint8(4), np.uint8(255);
    varout.flag_meanings = 'not_used good_data probably_good_data bad_data data_in_blind_range no_qc_performed';
    varout[:] = 2;

    Itmp = [record['I'][0].reshape(npulse,ngate) for record in records];
    Qtmp = [record['Q'][0].reshape(npulse,ngate) for record in records];

    I = np.stack(Itmp);
    Q = np.stack(Qtmp);

    Imean = np.mean(I,axis=1);
    Qmean = np.mean(Q,axis=1);

    I[:,:,:] = I[:,:,:]-Imean[:,None,:];
    Q[:,:,:] = Q[:,:,:]-Qmean[:,None,:];

#    Zcal=-146.8-1.0+115;

    Zcal=dBZ_offset+10*np.log10(DSout['prf'][:]/512);

    DSout['dBZ_offset'][:] = Zcal;

    I=I*np.sqrt(10**(Zcal/10.0));
    Q=Q*np.sqrt(10**(Zcal/10.0));

    DSout['I'][:] = I;
    DSout['Q'][:] = Q;

    blind_range = np.arange(15);
    DSout['qc_flag'][:,:,blind_range] = 4;

    return

def convert_galileo_ts_l0b2l1(infile,outfile,dBZ_offset,range_offset,data_version):

    """This routine converts raw (Level 0b) time series data from the Chilbolton 94GHz Cloud Radar (Galileo) to Level 1 data.
    Level 0b data have been generated from as-recorded Level 0a data by splitting the files along the time dimension and converting to NetCDF4.
    Processing by this routine involves removing redundant dimensions, and removing bias from the ADC samples of the received I and Q.
    Metadata are added specific to the WIVERN-2 project ground-based observations collected in 2020-2021.

    :param infile: Full path of NetCDF Level 0b data file, e.g. `<path-to-file>/radar-galileo_<YYYYddmmHHMMSS>-<YYYYddmmHHMMSS>_fix-ts.nc4`
    :type infile: str

    :param outfile: Full path of NetCDF Level 0b output file, e.g. `<path-to-file>/ncas-radar-w-band-1_cao_20201210-212823_fix-ts_l1_v1.0.nc`
    :type outfile: str

    :param dBZ_offset: dB calibration offset to apply to reflectivity.  This is converted to linear units and I and Q values are multiplied by the square root of this value.
    :type dBZoffset: float

    :param range_offset: Range offset to apply in metres.
    :type range_offset: float

    :param data_version: Version of data product in the format `N.m`, where `N` denotes the major verion and `m` a minor revision.
    :type data_version: str
    """

    # -----------------------------------------------------
    # Read NetCDF spectra file with embedded IQ time series
    # -----------------------------------------------------
    DSin = nc4.Dataset(infile);

    dt_start = cftime.num2pydate(DSin['time'][0],DSin['time'].units)
    dt_end   = cftime.num2pydate(DSin['time'][-1],DSin['time'].units)

    # -----------------------------
    # Set up NetCDF file for output
    # -----------------------------
    try: outfile.close()  # just to be safe, make sure dataset is not already open.
    except: pass
    DSout = nc4.Dataset(outfile,mode='w',format='NETCDF4')
    print(outfile)

    # ------------------------
    # Set up global attributes
    # ------------------------
    DSout.product_version = "v{}".format(data_version);
    DSout.processing_level = "1" ;

    DSout.licence = "Data usage licence - UK Open Government Licence agreement: \n http://www.nationalarchives.gov.uk/doc/open-government-licence" ;
    DSout.acknowledgement = "Acknowledgement is required of UK Research and Innovation as the data provider (in partnership with the National Centre for Atmospheric Science) whenever and wherever these data are used." ;
    DSout.platform = "Chilbolton Atmospheric Observatory" ;
    DSout.platform_type = "stationary_platform" ;
    DSout.title = "Time series from 94 GHz Galileo radar collected for ESA WIVERN-2 campaign at Chilbolton Observatory";
    DSout.creator_name = "Chris Walden" ;
    DSout.creator_email = "chris.walden@ncas.ac.uk" ;
    DSout.creator_url = "https://orcid.org/0000-0002-5718-466X" ;
    DSout.institution = "National Centre for Atmospheric Science (NCAS)";
    DSout.instrument_name = "ncas-radar-w-band-1";
    DSout.instrument_software = "radar-galileo-rec" ;
    DSout.instrument_software_version = "" ;

    DSout.references = "";
    DSout.source = "94GHz Galileo Radar";
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

    DSout.ADC_bits_per_sample = np.int(12);
    DSout.ADC_channels        = np.int(8);

    # ----------------
    # Scalar variables
    # ----------------

    varout = DSout.createVariable('latitude','f4');
    varout.standard_name = 'latitude';
    varout.long_name = 'latitude of the antenna';
    varout.units = 'degree_north';
    varout[:]=51.1447;

    varout = DSout.createVariable('longitude','f4');
    varout.standard_name = 'longitude';
    varout.long_name = 'longitude of the antenna';
    varout.units = 'degree_east';
    varout[:]=-1.4385;

    varout = DSout.createVariable('altitude','f4');
    varout.standard_name = 'altitude';
    varout.long_name = 'altitude of the antenna above the geoid (WGS84)';
    varout.units = 'm';
    varout[:] = 131.6;

    varout = DSout.createVariable('altitude_agl','f4');
    varout.long_name = 'altitude of the antenna above ground';
    varout.units = 'm';
    varout[:] = 0.9;

    varout = DSout.createVariable('frequency','f4');
    varout.standard_name = 'radiation_frequency';
    varout.long_name = 'frequency of transmitted radiation';
    varout.units = 'GHz';
    varout[:]=94.008;

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
    varout.comment = "Refers to diameter of each antenna in bistatic pair.  Separation of antennae is 0.66m."
    varout[:] = varin[:];

    varout = DSout.createVariable('antenna_focal_length','f4');
    varout.long_name = 'focal length of antenna';
    varout.units = 'm';
    varout[:] = 0.18;

    varin = DSin['pulse_period'];
    varout = DSout.createVariable('pulse_width',varin.datatype);
    varout.long_name = 'pulse width';
    varout.units = 'us';
    varout[:] = varin[:];

    varin = DSin['transmit_power'];
    varout = DSout.createVariable('transmit_power',varin.datatype);
    varout.long_name = 'peak transmitted power';
    varout.units = 'W';
    varout[:] = varin[:];

    varin = DSin['clock'];
    varout = DSout.createVariable('clock',varin.datatype);
    varout.long_name = 'clock input to timer card';
    varout.units = 'Hz';
    varout[:] = varin[:];

    clock_divfactor = DSin['clock'].clock_divfactor;
    varout = DSout.createVariable('clock_divfactor','i4');
    varout.long_name = 'clock divide factor';
    varout.units = '1';
    varout[:]=clock_divfactor;

    delay_clocks = DSin.getncattr('delay_clocks');
    varout = DSout.createVariable('delay_clocks','i4');
    varout.long_name = 'clock cycles before sampling is initiated';
    varout.units = '1';
    varout[:] = delay_clocks;

    samples_per_pulse = DSin.getncattr('samples_per_pulse');
    varout = DSout.createVariable('samples_per_pulse','i4');
    varout.long_name = 'number of samples per pulse';
    varout.units = '1';
    varout[:] = samples_per_pulse;

    pulses_per_daq_cycle = DSin.getncattr('pulses_per_daq_cycle');
    varout = DSout.createVariable('pulses_per_daq_cycle','i4');
    varout.long_name = 'number of pulses per data acquisition cycle';
    varout.units = '1';
    varout[:] = pulses_per_daq_cycle;

    pulses_per_ray = DSin.getncattr('pulses_per_ray');
    varout = DSout.createVariable('pulses_per_ray','i4');
    varout.long_name = 'number of pulses per ray';
    varout.units = '1';
    varout[:] = pulses_per_ray;

    varout = DSout.createVariable('radar_constant','f4');
    varout.long_name = 'radar constant';
    varout.units = 'dB';

    varout = DSout.createVariable('receiver_gain','f4');
    varout.long_name = 'receiver gain';
    varout.units = 'dB';

    varout = DSout.createVariable('cable_losses','f4');
    varout.long_name = 'cable losses';
    varout.units = 'dB';

    varout = DSout.createVariable('extra_attenuation','f4');
    varout.long_name = 'extra attenuation';
    varout.units = 'dB';
    varout[:] = 0.0;


    # ---------------
    # Copy dimensions
    # ---------------
    the_dim = DSin.dimensions['time'];
    DSout.createDimension('time', len(the_dim) if not the_dim.isunlimited() else None)

    fft_bin_dim = DSin.dimensions['fft_bin_dim'];
    spectra_number_dim = DSin.dimensions['spectra_number_dim'];
    DSout.createDimension('pulse', len(fft_bin_dim)*len(spectra_number_dim));

    the_dim = DSin.dimensions['range'];
    DSout.createDimension('range', len(the_dim) if not the_dim.isunlimited() else None)

    print(DSout.dimensions['pulse']);

    # --------------------
    # Coordinate variables
    # --------------------
    varin = DSin['time'];
    varout = DSout.createVariable('time',varin.datatype,('time'));
    varout.standard_name = 'time';
    varout.long_name = 'time at the start of each recorded ray';
    varout.units = varin.units;
    varout[:]=varin[:];

    varin = DSin['range'];
    varout = DSout.createVariable('range',varin.datatype,('range'));
    varout.long_name = 'distance from the antenna to the middle of each range gate';
    varout.range_offset_applied = np.float32(range_offset);
    varout.units = varin.units;
    varout[:]=varin[:]+range_offset-varin.range_offset;

    # --------------------------
    # Antenna pointing variables
    # --------------------------
    varout = DSout.createVariable('elevation','f4','time');
    varout.long_name = 'elevation angle above the horizon of the plane containing the transmit and receive antenna boresights';
    varout.units = 'degree';
    varout.elevation_offset_applied = np.float32(0.);

    varout = DSout.createVariable('azimuth','f4','time');
    varout.long_name = 'azimuth angle clockwise from grid north of the plane in which the antenna boresights are tilted away from zenith';
    varout.units = 'degree';
    varout.azimuth_offset_applied = np.float32(0.);

    # --------------------------------
    # Determine bias-corrected I and Q
    # --------------------------------
    # Input Level 0b file has I and Q dependent on the following dimensions (time,spectra,pulses,sample).
    # We rearrange this into a three-dimensional array dependent on (time,pulse,range);

    Ico_in = DSin['IPF_HH'][:,:,:,:];
    Qco_in = DSin['QPF_HH'][:,:,:,:];

    Icx_in = DSin['IPF_HV'][:,:,:,:];
    Qcx_in = DSin['QPF_HV'][:,:,:,:];

    nray     = Ico_in.shape[0];
    nspectra = Ico_in.shape[1];
    npulse   = Ico_in.shape[2];
    ngate    = Ico_in.shape[3];

    Ico = Ico_in.reshape(nray,npulse*nspectra,ngate);
    Qco = Qco_in.reshape(nray,npulse*nspectra,ngate);

    Ico_mean = np.mean(Ico,axis=1);
    Qco_mean = np.mean(Qco,axis=1);

    Ico[:,:,:] = Ico[:,:,:]-Ico_mean[:,None,:];
    Qco[:,:,:] = Qco[:,:,:]-Qco_mean[:,None,:];

    Zcal=dBZ_offset+10*np.log10(DSout['prf'][:]/512);

    Ico=Ico*np.sqrt(10**(Zcal/10.0));
    Qco=Qco*np.sqrt(10**(Zcal/10.0));

    Icx = Icx_in.reshape(nray,npulse*nspectra,ngate);
    Qcx = Qcx_in.reshape(nray,npulse*nspectra,ngate);

    Icx_mean = np.mean(Icx,axis=1);
    Qcx_mean = np.mean(Qcx,axis=1);

    Icx[:,:,:] = Icx[:,:,:]-Icx_mean[:,None,:];
    Qcx[:,:,:] = Qcx[:,:,:]-Qcx_mean[:,None,:];

    Icx=Icx*np.sqrt(10**(Zcal/10.0));
    Qcx=Qcx*np.sqrt(10**(Zcal/10.0));

    # ---------------
    # Field variables
    # ---------------

    varout = DSout.createVariable('I','f4',('time','pulse','range'),zlib=True);
    varout.ancillary_variables = 'qc_flag';
    varout.long_name = 'co-polar in-phase video signal';
    varout.units = '1';
    varout.comment = 'Scaled to account for calibration for square root of linear reflectivity factor';
#    add_offset = np.min(Ico[:]);
#    scale_factor = (np.max(Ico[:])-np.min(Ico[:])) / (2**16-1)
#    packed_data = np.int16(np.rint((Ico[:,:,:] - add_offset)/scale_factor));
    #varout.scale_factor = np.float32(scale_factor);
    #varout.add_offset = np.float32(add_offset);
    #varout[:] = packed_data;
    varout[:] = Ico;

    varout = DSout.createVariable('Q','f4',('time','pulse','range'),zlib=True);
    varout.ancillary_variables = 'qc_flag';
    varout.long_name = 'cross-polar quadrature video signal';
    varout.units = '1';
    varout.comment = 'Scaled to account for calibration for square root of linear reflectivity factor';
#    add_offset = np.min(Qco);
#    scale_factor = (np.max(Qco)-np.min(Qco)) / (2**16-1)
#    packed_data = np.int16(np.rint((Qco - add_offset)/scale_factor));
    #varout.scale_factor = np.float32(scale_factor);
    #varout.add_offset = np.float32(add_offset);
    #varout[:] = packed_data;
    varout[:] = Qco;

    varout = DSout.createVariable('Icx','f4',('time','pulse','range'),zlib=True);
    varout.ancillary_variables = 'qc_flag';
    varout.long_name = 'cross-polar in-phase video signal';
    varout.units = '1';
    varout.comment = 'Scaled to account for calibration for square root of linear reflectivity factor';
#    add_offset = np.min(Icx);
#    scale_factor = (np.max(Icx)-np.min(Icx)) / (2**16-1);
#    packed_data = np.int16(np.rint((Icx - add_offset)/scale_factor));
#    varout.scale_factor = np.float32(scale_factor);
#    varout.add_offset = np.float32(add_offset);
#    varout[:] = packed_data;
    varout[:] = Icx;

    varout = DSout.createVariable('Qcx','f4',('time','pulse','range'),zlib=True);
    varout.ancillary_variables = 'qc_flag';
    varout.long_name = 'cross-polar quadrature video signal';
    varout.units = '1';
    varout.comment = 'Scaled to account for calibration for square root of linear reflectivity factor';
#    add_offset = np.min(Qcx);
#    scale_factor = (np.max(Qcx)-np.min(Qcx)) / (2**16-1);
#    packed_data = np.int16(np.rint((Qcx - add_offset)/scale_factor));
#    varout.scale_factor = np.float32(scale_factor);
#    varout.add_offset = np.float32(add_offset);
#    varout[:] = packed_data;
    varout[:] = Qcx;

    varout = DSout.createVariable('qc_flag','u1',('time','pulse','range'));
    varout.is_quality = 'true';
    varout.qualified_variables = 'I Q Icx Qcx';
    varout.long_name = 'Quality control flag';
    varout.flag_values = np.uint8(0),np.uint8(1), np.uint8(2), np.uint8(3), np.uint8(4), np.uint8(255);
    varout.flag_meanings = 'not_used good_data probably_good_data bad_data data_in_blind_range no_qc_performed';
    varout[:] = 2;

    blind_range = np.arange(14);
    DSout['qc_flag'][:,:,blind_range] = 4;

    return

def process_wivern2_camra_ts(datestr,inpath,outpath):

    pattern = '*{}*_fix-ts.nc'.format(datestr);

    print(datestr);
    print(inpath);
    datepath = os.path.join(inpath,datestr);

    tsfiles = [];
    tsdirs = [];

    for root,dirs,files in os.walk(datepath):
        tsfiles += [os.path.join(root,f) for f in fnmatch.filter(files, pattern)];
        tsdirs += dirs;

    data_version = "1.0";

    dBZ_offset = 9.0;
    ZDR_offset = 0.6;
    range_offset = -865.56+864.0;

    l0bpath = os.path.join(outpath,'L0b',datestr);
    l1path = os.path.join(outpath,'L1',datestr);

    os.makedirs(l0bpath,exist_ok=True);
    os.makedirs(l1path,exist_ok=True);

    print(tsdirs);
    for dir in tsdirs:
        print("I am Here!");
        os.makedirs(os.path.join(l0bpath,dir),exist_ok=True);
        os.makedirs(os.path.join(l1path,dir),exist_ok=True);

    for f in tsfiles:
        outfile_splits = os.path.split(f);

        outfile_string = outfile_splits[1].replace('.nc','_l0b.nc');
        splits = outfile_string.split('_');
        instrument_name =splits[0].replace('radar-camra','ncas-radar-camra-1');
        platform = 'cao';
        datestr = splits[1][0:8];
        timestr = splits[1][8:];
        level = splits[3].split('.')[0];
        l0bfile = '{}_{}_{}-{}_{}_{}_v{}.nc'.format(instrument_name,platform,datestr,timestr,splits[2],level,data_version)

        l0bfile.replace('radar-camra','ncas-radar-camra-1_cao');

        event = outfile_splits[0].split('/')[-1];

        if not event in tsdirs:
            l0bfile = os.path.join(outpath,'L0b',datestr,l0bfile);
        else:
            l0bfile = os.path.join(outpath,'L0b',datestr,event,l0bfile);

        convert_camra_ts_l0a2l0b(f,l0bfile);

        l1file = l0bfile.replace('l0b','l1');
        l1file = l1file.replace('L0b','L1');

        print(l0bfile);
        print(l1file);

        convert_camra_ts_l0b2l1(l0bfile,l1file,dBZ_offset,ZDR_offset,range_offset,data_version);

    return

def process_wivern2_copernicus_ts(datestr,inpath,outpath):

    pattern = '{}*_iqdata.bin'.format(datestr);

    print(datestr);
    print(inpath);
    datepath = os.path.join(inpath,datestr);

    tsfiles = [];
    tsdirs = [];

    for root,dirs,files in os.walk(datepath):
        tsfiles += [os.path.join(root,f) for f in fnmatch.filter(files, pattern)];
        tsdirs += dirs;

    data_version = "1.0";

    dBZ_offset = -147.8+115;
    range_offset = -645.0;

    l1path = os.path.join(outpath,'L1',datestr);

    os.makedirs(l1path,exist_ok=True);

    print(tsdirs);
    for dir in tsdirs:
        print("I am Here!");
        os.makedirs(os.path.join(l1path,dir),exist_ok=True);

    for f in tsfiles:
        outfile_splits = os.path.split(f);

        outfile_string = outfile_splits[1].replace('iqdata.bin','fix-ts_l1.nc');
        splits = outfile_string.split('_');
        instrument_name ='ncas-radar-ka-band-1';
        platform = 'cao';
        datestr = splits[0][0:8];
        timestr = splits[0][8:];
        scantype = splits[1];
        level = splits[2].split('.')[0];
        l1file = '{}_{}_{}-{}_{}_{}_v{}.nc'.format(instrument_name,platform,datestr,timestr,scantype,level,data_version)

        event = outfile_splits[0].split('/')[-1];

        if not event in tsdirs:
            l1file = os.path.join(outpath,'L1',datestr,l1file);
        else:
            l1file = os.path.join(outpath,'L1',datestr,event,l1file);

        print(l1file);

        convert_copernicus_ts_l0a2l1(f,l1file,dBZ_offset,range_offset,data_version);

    return

def process_wivern2_galileo_ts(datestr,inpath,outpath):

    pattern = '*{}*_fix-fft-raw.nc4'.format(datestr);

    print(datestr);
    print(inpath);
    datepath = os.path.join(inpath,datestr);

    tsfiles = [];
    tsdirs = [];

    for root,dirs,files in os.walk(datepath):
        tsfiles += [os.path.join(root,f) for f in fnmatch.filter(files, pattern)];
        tsdirs += dirs;

    print(tsfiles);

    data_version = "1.0";

    dBZ_offset = -167.0+111.5;
    range_offset = -330.0;

    l1path = os.path.join(outpath,'L1',datestr);

    os.makedirs(l1path,exist_ok=True);

    print(tsdirs);
    for dir in tsdirs:
        print("I am Here!");
        os.makedirs(os.path.join(l1path,dir),exist_ok=True);

    for f in tsfiles:
        outfile_splits = os.path.split(f);

        outfile_string = outfile_splits[1].replace('.nc4','_l1.nc');
        splits = outfile_string.split('_');
        instrument_name =splits[0].replace('radar-galileo','ncas-radar-w-band-1');
        platform = 'cao';
        datestr0 = splits[1].split('-')[0];
        datestr = datestr0[0:8];
        timestr = datestr0[8:];
        level = splits[3].split('.')[0];
        l1file = '{}_{}_{}-{}_{}_{}_v{}.nc'.format(instrument_name,platform,datestr,timestr,'fix-ts',level,data_version)

        event = outfile_splits[0].split('/')[-1];

        if not event in tsdirs:
            l1file = os.path.join(outpath,'L1',datestr,l1file);
        else:
            l1file = os.path.join(outpath,'L1',datestr,event,l1file);

        print(l1file);

        convert_galileo_ts_l0b2l1(f,l1file,dBZ_offset,range_offset,data_version);

    return
