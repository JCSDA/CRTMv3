===============================================================================
TELSEM2: a Tool to Estimate Land Surface Emissivities at Microwave to
         Millimetre frequencies
===============================================================================

TELSEM2 consists of the atlas datasets and an associated interpolator which
calculates emissivities at the specified spatial location and MW frequency.

The TELSEM2 atlas can be used with RTTOV: see the RTTOV user guide for more
information. The Fortran module included here is based on the implementation
of TELSEM2 in RTTOV, but this module may be compiled and run independently of
RTTOV.

The subroutine interfaces for this version of the TELSEM2 module are described
below.

------------------------------------------------------------------------------

The following files are present alongside the atlas data:
- readme.txt         : this readme file.
- README_TELSEM2.pdf : the documentation for the TELSEM2 atlas data and the
                       original Fortran module.
- mod_mwatlas_m2.F90 : the RTTOV TELSEM2 Fortran module modified for
                       stand-alone compilation.
- test_telsem2.F90   : example code which calls the TELSEM2 atlas.
- Makefile           : the Makefile for compiling the example.


To compile and run the example, first edit the Makefile for your compiler
with suitable flags. Then do:

$ make
$ ./test_telsem2

------------------------------------------------------------------------------

The module defines three Fortran KINDs which are used for the arguments to the
module subroutines:

jpim - INTEGER KIND
jprb - REAL KIND
jplm - LOGICAL KIND


The atlas accessed by the following subroutines:

rttov_readmw_atlas
  - read and initialise the atlas for a specified month

emis_interp_ind_sing
  - calculate emissivity at the native atlas resolution for a single frequency

emis_interp_ind_mult
  - calculate emissivities at the native atlas resolution for multiple frequencies

emis_interp_int_sing
  - calculate spatially-averaged emissivities for a single frequency

emis_interp_int_mult
  - calculate spatially-averaged emissivities for multiple frequencies

rttov_closemw_atlas
  - deallocate memory used by the atlas


---------------------------------------

To read and initialise the atlas:

First declare a variable, say "atlas" of derived type telsem2_atlas_data.
This will contain the loaded atlas data and is passed to all subroutines.

CALL rttov_readmw_atlas(dir, month, verbose, atlas, err, lat1, lat2, lon1, lon2)

  dir        : IN,    CHARACTER, location of TELSEM atlas files
  month      : IN,    INTEGER(KIND=jpim), month 1-12 of data to read
  verbose    : IN,    LOGICAL(KIND=jplm), if .TRUE. additional information is printed
  atlas      : INOUT, TYPE(telsem2_atlas_data), contains the loaded atlas data
  err        : OUT,   INTEGER(KIND=jpim), return status, non-zero implies an error condition

  By default data for the whole globe is read.
  Optionally, data for a limited area may be read in by specifying:

  lat1, lon1 : IN,  INTEGER(KIND=jpim), OPTIONAL, south-west corner of region to read
  lat2, lon2 : IN,  INTEGER(KIND=jpim), OPTIONAL, north-east corner of region to read

  NB The longitudes lon1 and lon2 must conform to 0 <= lon1 < lon2 < 360

---------------------------------------

To deallocate the atlas:

CALL rttov_closemw_atlas(atlas)

  atlas      : INOUT, TYPE(telsem2_atlas_data), atlas data to deallocate

---------------------------------------

To return V- and H-pol emissivities and optionally standard deviations and/or covariances:

CALL emis_interp_ind_sing(lat, lon, theta, freq, atlas, ev, eh, stdv, stdh, covvh, verb)

  lat        : IN,  REAL(KIND=jprb), latitude, degrees, -90 <= lat <= 90
  lon        : IN,  REAL(KIND=jprb), longitude, degrees, 0 <= lon < 360
  theta      : IN,  REAL(KIND=jprb), zenith angle, degrees, 0 <= theta <= 60
  freq       : IN,  REAL(KIND=jprb), central frequency of channel, GHz
  atlas      : IN,  TYPE(telsem2_atlas_data), atlas data to use
  ev, eh     : OUT, REAL(KIND=jprb), output V- and H-pol emissivities
  stdv, stdh : OUT, REAL(KIND=jprb), OPTIONAL, output standard deviations of V- and H-pol emissivities
  covvh      : OUT, REAL(KIND=jprb), OPTIONAL, covariance of V- and H-pol emissivities
  verb       : IN,  INTEGER(KIND=jpim), if 1 additional information is printed


CALL emis_interp_ind_mult(lat, lon, theta, freq, nchan, atlas, ev, eh, std, verb)

  lat                  : IN,  REAL(KIND=jprb), latitude, degrees, -90 <= lat <= 90
  lon                  : IN,  REAL(KIND=jprb), longitude, degrees, 0 <= lon < 360
  theta                : IN,  REAL(KIND=jprb), zenith angle, degrees, 0 <= theta <= 60
  freq(nchan)          : IN,  REAL(KIND=jprb), central frequencies of channel, GHz
  nchan                : IN,  INTEGER(KIND=jpim), number of frequencies
  atlas                : IN,  TYPE(telsem2_atlas_data), atlas data to use
  ev(nchan), eh(nchan) : OUT, REAL(KIND=jprb), output V- and H-pol emissivities
  std(2*nchan,2*nchan) : OUT, REAL(KIND=jprb), OPTIONAL, output covariance matrix
  verb                 : IN,  INTEGER(KIND=jpim), if 1 additional information is printed


CALL emis_interp_int_sing(lat, lon, resol, theta, freq, atlas, ev, eh, stdv, stdh, covvh, verb)

  lat        : IN,  REAL(KIND=jprb), latitude, degrees, -90 <= lat <= 90
  lon        : IN,  REAL(KIND=jprb), longitude, degrees, 0 <= lon < 360
  resol      : IN,  REAL(KIND=jprb), horizontal resolution, degrees, resol>=0.25
  theta      : IN,  REAL(KIND=jprb), zenith angle, degrees, 0 <= theta <= 60
  freq       : IN,  REAL(KIND=jprb), central frequency of channel, GHz
  atlas      : IN,  TYPE(telsem2_atlas_data), atlas data to use
  ev, eh     : OUT, REAL(KIND=jprb), output V- and H-pol emissivities
  stdv, stdh : OUT, REAL(KIND=jprb), OPTIONAL, output standard deviations of V- and H-pol emissivities
  covvh      : OUT, REAL(KIND=jprb), OPTIONAL, covariance of V- and H-pol emissivities
  verb       : IN,  INTEGER(KIND=jpim), if 1 additional information is printed


CALL emis_interp_int_mult(lat, lon, resol, theta, freq, nchan, atlas, ev, eh, std, verb)

  lat                  : IN,  REAL(KIND=jprb), latitude, degrees, -90 <= lat <= 90
  lon                  : IN,  REAL(KIND=jprb), longitude, degrees, 0 <= lon < 360
  resol                : IN,  REAL(KIND=jprb), horizontal resolution, degrees, resol>=0.25
  theta                : IN,  REAL(KIND=jprb), zenith angle, degrees, 0 <= theta <= 60
  freq(nchan)          : IN,  REAL(KIND=jprb), central frequencies of channel, GHz
  nchan                : IN,  INTEGER(KIND=jpim), number of frequencies
  atlas                : IN,  TYPE(telsem2_atlas_data), atlas data to use
  ev(nchan), eh(nchan) : OUT, REAL(KIND=jprb), output V- and H-pol emissivities
  std(2*nchan,2*nchan) : OUT, REAL(KIND=jprb), OPTIONAL, output covariance matrix
  verb                 : IN,  INTEGER(KIND=jpim), if 1 additional information is printed


Notes:
  - resol : the native atlas resolution is 0.25 degrees in latitude and
            longitude. By specifying resol with a value greater than 0.25, the
            atlas emissivities are averaged over adjacent cells.

  - std(:,:) : the covariance matrices are ordered with indices 1:2*nchan
               corresponding to 1:nchan V-pol frequencies followed by 1:nchan
               H-pol frequencies.

------------------------------------------------------------------------------
