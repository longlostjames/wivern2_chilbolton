

Filename conventions
====================

Each radar file follows the following naming convention:

`<instrument_name>_<platform_name>_<date>-<time>_<scan_type>_<processing_level>_v<version>.nc`,

where the platform name is set to `cao`, to denote Chilbolton Atmospheric
Observatory, and the date is given in the format `<YYYYmmdd>`.

The scan_type is `fix-ts`, as all files relate to time series from fixed vertical
pointing dwells.


There are three instruments:

`ncas-radar-camra-1` is the 3 GHz Chilbolton Advanced Meteorological Radar
(CAMRa), which uses the 25 m Chilbolton antenna.

`ncas-radar-ka-band-1` is the 35 GHz Copernicus cloud radar, which has a 2.4 m
fixed vertically pointing antenna.

`ncas-radar-w-band-1` is the 94 GHz Galileo radar, which is a bistatic system
with two 0.46 m antennae separated by 0.66 m.


.. _file-format:

File format description
=======================

Level 1 files are provided, and follow the description below.  They have been
derived by processing Level 0a (as-recorded) data, with the processing steps
indicated in the history attribute.  See :ref:`raw-data-conversion` for a
summary of processing software.


Level 1 files
-------------

These files are in NetCDF-4 format with the following content:

Dimensions
..........

The following dimensions are present for all three radars: `time`, `pulse` and `range`.


Common global attributes
........................

The left column of the following table lists global attributes that are present
for all three radars. The right column provides example attribute values
relevant to the 3GHz CAMRa radar.

.. tabularcolumns:: |>{\raggedright\arraybackslash}\X{3}{8}|>{\raggedright\arraybackslash}\X{5}{8}|

.. table::
   :widths: auto
   :class: longtable


   +-------------------------------------+----------------------------------------------------------------------------------+
   |Name                                 |Example                                                                           |
   +=====================================+==================================================================================+
   |product_version                      |v1.0                                                                              |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |licence                              |This dataset is released for use under a Creative Commons Attribution 4.0         |
   |                                     |International (CC-BY 4.0) license                                                 |
   |                                     |(see https://creativecommons.org/licenses/by/4.0/ for terms and conditions)       |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |acknowledgement                      |This dataset was developed as part of the activity                                |
   |                                     |"Doppler Wind Radar Science Performance Study (WIVERN-2)", funded by the          |
   |                                     |European Space Agency under Contract no. 4000130864/20/NL/CT.  Users should       |
   |                                     |acknowledge UK Research and Innovation as the data provider (in partnership       |
   |                                     |with the National Centre for Atmospheric Science)                                 |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |platform                             |Chilbolton Atmospheric Observatory                                                |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |platform_type                        |stationary_platform.                                                              |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |creator_name                         |Chris Walden                                                                      |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |creator_email                        |chris.walden@ncas.ac.uk                                                           |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |creator_url                          |https://orcid.org/0000-0002-5718-466X                                             |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |institution                          |National Centre for Atmospheric Science (NCAS)                                    |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |instrument_name                      |ncas-radar-camra-1                                                                |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |instrument_software                  |radar-camra-rec                                                                   |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |instrument_software_version          |1.4 Rev 58                                                                        |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |references                           |https://doi.org/10.1049/ecej:19940205; http://purl.org/net/epubs/work/63318;      |
   |                                     |https://doi.org/10.1109/IGARSS.2006.429; https://doi.org/10.3390/atmos10110714;   |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |source                               |3GHz Chilbolton Advanced Meteorological Radar (CAMRa)                             |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |project                              |WIVERN-2 Doppler Wind Radar Science Performance Study                             |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |project_principal_investigator       |Anthony Illingworth                                                               |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |project_principal_investigator_email |a.j.illingworth@reading.ac.uk                                                     |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |project_principal_investigator_url   |https://orcid.org/0000-0002-5774-8410                                             |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |processing_software_url              |https://github.com/longlostjames/wivern_chilbolton_utils.git                      |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |processing_software_version          |1.0                                                                               |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |processing_level                     |1                                                                                 |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |scantype                             |vertical_pointing                                                                 |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |time_coverage_start                  |2020-10-01T18:09:28Z                                                              |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |time_coverage_end                    |2020-10-01T18:21:22Z                                                              |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |geospatial_bounds                    |51.1450N -1.4384E                                                                 |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |pulse_compression                    |false                                                                             |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |ADC_bits_per_sample                  |12                                                                                |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |ADC_channels                         |8                                                                                 |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |last_revised_date                    |2022-03-28T14:19:36Z                                                              |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |title                                |Calibrated time series from 3 GHz CAMRa radar collected for ESA WIVERN-2          |
   |                                     |campaign at Chilbolton Observatory                                                |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |comment                              |Correction to account for inverse square power loss with range has not been       |
   |                                     |applied                                                                           |
   +-------------------------------------+----------------------------------------------------------------------------------+
   |history                              |Mon Feb 28 17:20:39 2022 - user:cjwalden machine: host293.jc.rl.ac.uk program:    |
   |                                     |wivern_chilbolton_utils.py convert_camra_ts_l0b2l1 version:1.0\n                  |
   |                                     |Mon Feb 28 17:20:00 2022 - user:cjwalden machine: host293.jc.rl.ac.uk program:    |
   |                                     |wivern_chilbolton_utils.py convert_camra_ts_l0a2l0b version:1.0\n                 |
   |                                     |Mon Dec  6 16:58:37 2021: /home/users/cjwalden/anaconda3/envs/cao_3_8/bin/ncks    |
   |                                     |-d time,6,37 --output=radar-camra_20201001180645_fix-ts.nc                        |
   |                                     |radar-camra_20201001180645_fix-ts_orig.nc\n                                       |
   |                                     |Thu Oct 01 18:07:00 2020 - /usr/local/bin/radar-camra-rec -fix 3600 115 90        |
   |                                     |-gates 5 201 -cellsize 1 -pulse_pairs 3050 -op rad -id 0 -file 8030 -scan 7827    |
   |                                     |-date 20201001180645 -tsdump -tssamples 200                                       |
   +-------------------------------------+----------------------------------------------------------------------------------+


Scalar variables
................

The following scalar variables are present for all three radars:

.. tabularcolumns:: |>{\raggedright\arraybackslash}\X{3}{10}|>{\raggedright\arraybackslash}\X{1}{10}|>{\raggedright\arraybackslash}\X{4}{10}|>{\raggedright\arraybackslash}\X{2}{10}|

.. table::
   :widths: auto
   :class: longtable

   +------------------------------+---------------+-----------------------------------------------------------------------------------+----------------------------------------+
   |Name                          |Data type      |Long name                                                                          |Units                                   |
   +==============================+===============+===================================================================================+========================================+
   |latitude                      |float32        |latitude of the antenna                                                            |degree_north                            |
   +------------------------------+---------------+-----------------------------------------------------------------------------------+----------------------------------------+
   |longitude                     |float32        |longitude of the antenna                                                           |degree_east                             |
   +------------------------------+---------------+-----------------------------------------------------------------------------------+----------------------------------------+
   |frequency                     |float32        |frequency of transmitted radiation                                                 |GHz                                     |
   +------------------------------+---------------+-----------------------------------------------------------------------------------+----------------------------------------+
   |prf                           |float32        |pulse repetition frequency                                                         |Hz                                      |
   +------------------------------+---------------+-----------------------------------------------------------------------------------+----------------------------------------+
   |beamwidthH                    |float32        |horizontal angular beamwidth                                                       |degree                                  |
   +------------------------------+---------------+-----------------------------------------------------------------------------------+----------------------------------------+
   |beamwidthV                    |float32        |vertical angular beamwidth                                                         |degree                                  |
   +------------------------------+---------------+-----------------------------------------------------------------------------------+----------------------------------------+
   |antenna_diameter              |float32        |antenna diameter                                                                   |m                                       |
   +------------------------------+---------------+-----------------------------------------------------------------------------------+----------------------------------------+
   |antenna_focal_length          |float32        |focal length of antenna                                                            |m                                       |
   +------------------------------+---------------+-----------------------------------------------------------------------------------+----------------------------------------+
   |pulse_width                   |float32        |pulse width                                                                        |us                                      |
   +------------------------------+---------------+-----------------------------------------------------------------------------------+----------------------------------------+
   |transmit_power                |float32        |peak transmitted power                                                             |W                                       |
   +------------------------------+---------------+-----------------------------------------------------------------------------------+----------------------------------------+
   |clock                         |float32        |clock input to timer card                                                          |Hz                                      |
   +------------------------------+---------------+-----------------------------------------------------------------------------------+----------------------------------------+
   |clock_divide_factor           |float32        |clock divide factor                                                                |1                                       |
   +------------------------------+---------------+-----------------------------------------------------------------------------------+----------------------------------------+
   |delay_clocks                  |float32        |clock cycles before sampling is initiated                                          |1                                       |
   +------------------------------+---------------+-----------------------------------------------------------------------------------+----------------------------------------+
   |samples_per_pulse             |float32        |number of samples per pulse                                                        |1                                       |
   +------------------------------+---------------+-----------------------------------------------------------------------------------+----------------------------------------+
   |pulses_per_daq_cycle          |float32        |number of pulses per data acquisition cycle                                        |1                                       |
   +------------------------------+---------------+-----------------------------------------------------------------------------------+----------------------------------------+
   |pulses_per_ray                |float32        |number of pulses per ray                                                           |1                                       |
   +------------------------------+---------------+-----------------------------------------------------------------------------------+----------------------------------------+
   |radar_constant                |float32        |radar constant                                                                     |dB                                      |
   +------------------------------+---------------+-----------------------------------------------------------------------------------+----------------------------------------+
   |receiver_gain                 |float32        |receiver gain                                                                      |dB                                      |
   +------------------------------+---------------+-----------------------------------------------------------------------------------+----------------------------------------+
   |cable_losses                  |float32        |cable losses                                                                       |dB                                      |
   +------------------------------+---------------+-----------------------------------------------------------------------------------+----------------------------------------+
   |extra_attenuation             |float32        |extra attenuation                                                                  |dB                                      |
   +------------------------------+---------------+-----------------------------------------------------------------------------------+----------------------------------------+

The following scalar variables are instrument specific, or have an instrument-specific meaning.

**3GHz CAMRa radar**

.. tabularcolumns:: |>{\raggedright\arraybackslash}\X{3}{10}|>{\raggedright\arraybackslash}\X{1}{10}|>{\raggedright\arraybackslash}\X{4}{10}|>{\raggedright\arraybackslash}\X{2}{10}|

.. table::
   :widths: auto
   :class: longtable


   +------------------------------+---------------+-----------------------------------------------------------------------------------+----------------------------------------+
   |Name                          |Data type      |Long name                                                                          |Units                                   |
   +==============================+===============+===================================================================================+========================================+
   |altitude                      |float32        |altitude of the elevation axis above the geoid (WGS84)                             |m                                       |
   +------------------------------+---------------+-----------------------------------------------------------------------------------+----------------------------------------+
   |altitude_agl                  |float32        |altitude of the elevation axis above ground                                        |m                                       |
   +------------------------------+---------------+-----------------------------------------------------------------------------------+----------------------------------------+
   |antenna_focus_radial_location |float32        |distance along boresight from elevation axis to focus                              |m                                       |
   +------------------------------+---------------+-----------------------------------------------------------------------------------+----------------------------------------+

**35GHz Copernicus and 94GHz Galileo radars**

.. tabularcolumns:: |>{\raggedright\arraybackslash}\X{3}{10}|>{\raggedright\arraybackslash}\X{1}{10}|>{\raggedright\arraybackslash}\X{4}{10}|>{\raggedright\arraybackslash}\X{2}{10}|

.. table::
  :widths: auto
  :class: longtable

  +------------------------------+---------------+-----------------------------------------------------------------------------------+----------------------------------------+
  |Name                          |Data type      |Long name                                                                          |Units                                   |
  +==============================+===============+===================================================================================+========================================+
  |altitude                      |float32        |altitude of the antenna above the geoid (WGS84)                                    |m                                       |
  +------------------------------+---------------+-----------------------------------------------------------------------------------+----------------------------------------+
  |altitude_agl                  |float32        |altitude of the antenna above ground                                               |m                                       |
  +------------------------------+---------------+-----------------------------------------------------------------------------------+----------------------------------------+
  |dBZ_offset                    |float32        |dBZ offset applied                                                                 |dB                                      |
  +------------------------------+---------------+-----------------------------------------------------------------------------------+----------------------------------------+


Coordinate variables
....................

These have slightly different interpretation depending on the particular radar.

**3GHz CAMRa radar**

.. tabularcolumns:: |>{\raggedright\arraybackslash}\X{2}{12}|>{\raggedright\arraybackslash}\X{1}{12}|>{\raggedright\arraybackslash}\X{2}{12}|>{\raggedright\arraybackslash}\X{4}{12}|>{\raggedright\arraybackslash}\X{3}{12}|

.. table::
  :widths: auto
  :class: longtable

  +------------------------------+---------------+-----------------+-------------------------------------------------------------------------------------+----------------------------------------+
  |Name                          |Data type      |Dimension        |Long name                                                                            |Units                                   |
  +==============================+===============+=================+=====================================================================================+========================================+
  |time                          |float32        |time             |time at the end of each recorded ray                                                 |seconds since 2020-09-22 00:00:00 +00:00|
  +------------------------------+---------------+-----------------+-------------------------------------------------------------------------------------+----------------------------------------+
  |range                         |float32        |range            |distance from the antenna to the middle of each range gate                           |m                                       |
  +------------------------------+---------------+-----------------+-------------------------------------------------------------------------------------+----------------------------------------+
  |elevation                     |float32        |time             |elevation angle above the horizon of the antenna boresight                           |degree                                  |
  +------------------------------+---------------+-----------------+-------------------------------------------------------------------------------------+----------------------------------------+
  |azimuth                       |float32        |time             |azimuth angle of the antenna boresight clockwise from grid north.                    |degree                                  |
  +------------------------------+---------------+-----------------+-------------------------------------------------------------------------------------+----------------------------------------+

**35GHz Copernicus radar**

.. tabularcolumns:: |>{\raggedright\arraybackslash}\X{2}{12}|>{\raggedright\arraybackslash}\X{1}{12}|>{\raggedright\arraybackslash}\X{2}{12}|>{\raggedright\arraybackslash}\X{4}{12}|>{\raggedright\arraybackslash}\X{3}{12}|

.. table::
  :widths: auto
  :class: longtable

  +------------------------------+---------------+-----------------+-------------------------------------------------------------------------------------+----------------------------------------+
  |Name                          |Data type      |Dimension        |Long name                                                                            |Units                                   |
  +==============================+===============+=================+=====================================================================================+========================================+
  |time                          |float32        |time             |time at the start of each recorded ray                                               |seconds since 2020-09-22 00:00:00 +00:00|
  +------------------------------+---------------+-----------------+-------------------------------------------------------------------------------------+----------------------------------------+
  |range                         |float32        |range            |distance from the antenna to the middle of each range gate                           |m                                       |
  +------------------------------+---------------+-----------------+-------------------------------------------------------------------------------------+----------------------------------------+
  |elevation                     |float32        |time             |elevation angle above the horizon of the antenna boresight                           |degree                                  |
  +------------------------------+---------------+-----------------+-------------------------------------------------------------------------------------+----------------------------------------+
  |azimuth                       |float32        |time             |azimuth angle from grid north of the plane containing the antenna boresight and      |degree                                  |
  |                              |               |                 |zenith vectors                                                                       |                                        |
  +------------------------------+---------------+-----------------+-------------------------------------------------------------------------------------+----------------------------------------+

**94GHz Galileo radar**

.. tabularcolumns:: |>{\raggedright\arraybackslash}\X{2}{12}|>{\raggedright\arraybackslash}\X{1}{12}|>{\raggedright\arraybackslash}\X{2}{12}|>{\raggedright\arraybackslash}\X{4}{12}|>{\raggedright\arraybackslash}\X{3}{12}|

.. table::
  :widths: auto
  :class: longtable

  +------------------------------+---------------+-----------------+-------------------------------------------------------------------------------------+----------------------------------------+
  |Name                          |Data type      |Dimension        |Long name                                                                            |Units                                   |
  +==============================+===============+=================+=====================================================================================+========================================+
  |time                          |float32        |time             |time at the end of each recorded ray                                                 |seconds since 2020-09-22 00:00:00 +00:00|
  +------------------------------+---------------+-----------------+-------------------------------------------------------------------------------------+----------------------------------------+
  |range                         |float32        |range            |distance from the antenna to the middle of each range gate                           |m                                       |
  +------------------------------+---------------+-----------------+-------------------------------------------------------------------------------------+----------------------------------------+
  |elevation                     |float32        |time             |elevation angle above the horizon of the antenna boresight                           |degree                                  |
  +------------------------------+---------------+-----------------+-------------------------------------------------------------------------------------+----------------------------------------+
  |azimuth                       |float32        |time             |azimuth angle from grid north of the plane containing the antenna boresight and      |degree                                  |
  |                              |               |                 |zenith vectors                                                                       |                                        |
  +------------------------------+---------------+-----------------+-------------------------------------------------------------------------------------+----------------------------------------+

Field variables
...............

**3GHz CAMRa radar**

.. tabularcolumns:: |>{\raggedright\arraybackslash}\X{2}{12}|>{\raggedright\arraybackslash}\X{1}{12}|>{\raggedright\arraybackslash}\X{2}{12}|>{\raggedright\arraybackslash}\X{4}{12}|>{\raggedright\arraybackslash}\X{3}{12}|

.. table::
  :widths: auto
  :class: longtable

  +------------------------------+---------------+-------------------------+-----------------------------------------------------------------------------+----------------------------------------+
  |Name                          |Date type      |Dimensions               |Long name                                                                    |Units                                   |
  +==============================+===============+=========================+=============================================================================+========================================+
  |I                             |float32        |time, pulse, range       |co-polar in-phase video signal                                               |1                                       |
  +------------------------------+---------------+-------------------------+-----------------------------------------------------------------------------+----------------------------------------+
  |Q                             |float32        |time, pulse, range       |co-polar quadrature video signal                                             |1                                       |
  +------------------------------+---------------+-------------------------+-----------------------------------------------------------------------------+----------------------------------------+
  |ZCX                           |float32        |time, pulse, range       |cross-polar radar equivalent reflectivity factor                             |dBZ                                     |
  +------------------------------+---------------+-------------------------+-----------------------------------------------------------------------------+----------------------------------------+

**35GHz Copernicus radar**

.. tabularcolumns:: |>{\raggedright\arraybackslash}\X{2}{12}|>{\raggedright\arraybackslash}\X{1}{12}|>{\raggedright\arraybackslash}\X{2}{12}|>{\raggedright\arraybackslash}\X{4}{12}|>{\raggedright\arraybackslash}\X{3}{12}|

.. table::
  :widths: auto
  :class: longtable

  +------------------------------+---------------+-------------------------+-----------------------------------------------------------------------------+----------------------------------------+
  |Name                          |Date type      |Dimensions               |Long name                                                                    |Units                                   |
  +==============================+===============+=========================+=============================================================================+========================================+
  |I                             |float32        |time, pulse, range       |co-polar in-phase video signal                                               |1                                       |
  +------------------------------+---------------+-------------------------+-----------------------------------------------------------------------------+----------------------------------------+
  |Q                             |float32        |time, pulse, range       |co-polar quadrature video signal                                             |1                                       |
  +------------------------------+---------------+-------------------------+-----------------------------------------------------------------------------+----------------------------------------+

**94GHz Galileo radar**

.. tabularcolumns:: |>{\raggedright\arraybackslash}\X{2}{12}|>{\raggedright\arraybackslash}\X{1}{12}|>{\raggedright\arraybackslash}\X{2}{12}|>{\raggedright\arraybackslash}\X{4}{12}|>{\raggedright\arraybackslash}\X{3}{12}|

.. table::
  :widths: auto
  :class: longtable

  +------------------------------+---------------+-------------------------+-----------------------------------------------------------------------------+----------------------------------------+
  |Name                          |Date type      |Dimensions               |Long name                                                                    |Units                                   |
  +==============================+===============+=========================+=============================================================================+========================================+
  |I                             |float32        |time, pulse, range       |co-polar in-phase video signal                                               |1                                       |
  +------------------------------+---------------+-------------------------+-----------------------------------------------------------------------------+----------------------------------------+
  |Q                             |float32        |time, pulse, range       |co-polar quadrature video signal                                             |1                                       |
  +------------------------------+---------------+-------------------------+-----------------------------------------------------------------------------+----------------------------------------+
  |Icx                           |float32        |time, pulse, range       |cross-polar in-phase video signal                                            |1                                       |
  +------------------------------+---------------+-------------------------+-----------------------------------------------------------------------------+----------------------------------------+
  |Qcx                           |float32        |time, pulse, range       |cross-polar quadrature video signal                                          |1                                       |
  +------------------------------+---------------+-------------------------+-----------------------------------------------------------------------------+----------------------------------------+


Quality control variables
.........................

All radars have the following quality-control flag:

.. tabularcolumns:: |>{\raggedright\arraybackslash}\X{2}{12}|>{\raggedright\arraybackslash}\X{1}{12}|>{\raggedright\arraybackslash}\X{2}{12}|>{\raggedright\arraybackslash}\X{4}{12}|>{\raggedright\arraybackslash}\X{3}{12}|

.. table::
  :widths: auto
  :class: longtable

  +------------------------------+---------------+-------------------------+-----------------------------------------------------------------------------+----------------------------------------+
  |Name                          |Date type      |Dimensions               |Long name                                                                    |Units                                   |
  +==============================+===============+=========================+=============================================================================+========================================+
  |qc_flag                       |uint8          |time, pulse, range       |quality control flag                                                         |                                        |
  +------------------------------+---------------+-------------------------+-----------------------------------------------------------------------------+----------------------------------------+

**3GHz CAMRa radar**

Alternating H- and V-polarised pulses are transmitted by the 3GHz radar.
Different dBZ offsets are applied for the two polarisations, and this is captured
in the following variable:

.. tabularcolumns:: |>{\raggedright\arraybackslash}\X{2}{12}|>{\raggedright\arraybackslash}\X{1}{12}|>{\raggedright\arraybackslash}\X{2}{12}|>{\raggedright\arraybackslash}\X{4}{12}|>{\raggedright\arraybackslash}\X{3}{12}|

.. table::
  :widths: auto
  :class: longtable

  +------------------------------+---------------+-------------------------+-----------------------------------------------------------------------------+----------------------------------------+
  |Name                          |Date type      |Dimensions               |Long name                                                                    |Units                                   |
  +==============================+===============+=========================+=============================================================================+========================================+
  |dBZ_offsets_applied           |float32        |pulse                    |dBZ calibration offset applied for even and odd pulses                       |dB                                      |
  +------------------------------+---------------+-------------------------+-----------------------------------------------------------------------------+----------------------------------------+
