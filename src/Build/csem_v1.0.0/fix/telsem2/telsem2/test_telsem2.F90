PROGRAM test_telsem

!
! Demonstration program showing how to call the TELSEM2 atlas and interpolator.
!

  USE mod_mwatlas_m2

  IMPLICIT NONE

  INTEGER(jpim) :: error_status

  LOGICAL(jplm) :: verbose   ! For atlas reading subroutine
  INTEGER(jpim) :: verb      !=1 for TRUE and 0 for FALSE - for emissivity routines

  INTEGER(jpim) :: month     ! (1->12)
  CHARACTER(LEN=256) :: dir  ! directory of emis database

  REAL(jprb)    :: resol     ! horizontal resolution for the user
  REAL(jprb)    :: lat       ! (-90->90)
  REAL(jprb)    :: lon       ! (0->360)
  REAL(jprb)    :: theta     ! (0->60°)

  ! For individual freq interpolations
  REAL(jprb)    :: ev, eh, stdv, stdh, covvh
  REAL(jprb)    :: freq      ! (19->85GHz)
  INTEGER(jpim) :: i

  ! For multiple freq interpolations
  INTEGER(jpim), PARAMETER :: nchan = 5
  REAL(jprb)    :: ev2(nchan), eh2(nchan), std(2*nchan,2*nchan)
  REAL(jprb)    :: freq2(nchan)

  ! Structure containing atlas data
  TYPE(telsem2_atlas_data) :: atlas

  !--- End of header ----------------------------------


  verbose = .TRUE.  ! Verbose output for reading subroutine
  verb = 0          ! No verbose output for emissivity subroutines

  !====================================================
  ! Read the atlas
  !====================================================

  dir = '../'
  month = 9

  WRITE(0,'(a,i3)') 'Reading atlas for month ',month
  CALL rttov_readmw_atlas(TRIM(dir), month, atlas, verbose, error_status)

  IF (error_status /= 0) THEN
    WRITE(0,'(a)') 'Error reading atlas'
    STOP 1
  ENDIF


  !====================================================
  ! Calculate emissivities
  !====================================================

  resol    = 0.25
  theta    = 15.

  lat      = -30.
  lon      = 302.
  freq     = 30.                          ! For single frequency experiment
  freq2(:) = (/30., 25., 38., 60., 90. /) ! For multiple frequency experiment
  i        = 1                            ! Index of freq in freq2 array

  WRITE(0,'(a)') ' '
  WRITE(0,'(a)') 'Inputs:'
  WRITE(0,'(a,f8.2)') 'lat   = ', lat
  WRITE(0,'(a,f8.2)') 'lon   = ', lon
  WRITE(0,'(a,f8.2)') 'theta = ', theta
  WRITE(0,'(a,f8.2)') 'freq  = ', freq
  WRITE(0,'(a)') ' '


  ! This example first makes calls to the atlas at it's native resolution (resol=0.25 degrees)
  ! and compares the outputs at 30GHz. These should all be identical for all 4 subroutines.

  WRITE(0,'(a)') 'The first four sets of output are identical:'
  WRITE(0,'(a)') ' '

  !----------------------------------------------------------------------------
  ! First call with no spatial interpolation
  !----------------------------------------------------------------------------

  ! Single frequency, no spatial averaging
  CALL emis_interp_ind_sing(lat, lon, theta, freq, atlas, ev, eh, stdv, stdh, covvh, verb)

  ! For the single-frequency call, stdv and stdh contain the standard deviations of the
  ! two polarisations and covvh contains the cross-polarisation covariance.

  CALL print_emis(ev, eh, stdv, stdh, covvh, 'Single freq, no spatial averaging')


  ! Multiple frequencies, no spatial averaging
  CALL emis_interp_ind_mult(lat, lon, theta, freq2, nchan, atlas, ev2, eh2, std, verb)

  ! For the multiple-frequency call, note the index used to access the H-pol standard
  ! deviation for the 30GHz channel: all V-pol emissivites are represented in std(:,:)
  ! followed by the H-pol emissivities. The std(:,:) array is the covariance matrix so
  ! we must take the SQRT of the diagonal elements to obtain the standard deviation.
  ! The off-diagonal elements of std(:,:) contain the cross-channel and cross-polarisation
  ! covariances.
                                                ! Recall i is the index of 30GHz in freq2(:) so:
  CALL print_emis(ev2(i),                     & ! 30GHz V-pol
                  eh2(i),                     & ! 30GHz H-pol
                  SQRT(std(i,i)),             & ! 30GHz V-pol standard deviation
                  SQRT(std(i+nchan,i+nchan)), & ! 30GHz H-pol standard deviation
                  std(i,i+nchan),             & ! 30GHz V-/H-pol covariance
                  'Multiple freq, no spatial averaging')



  !----------------------------------------------------------------------------
  ! Now call with spatial interpolation, but with resol = 0.25 degrees
  !----------------------------------------------------------------------------

  ! Single frequency, spatially averaged
  CALL emis_interp_int_sing(lat, lon, resol, theta, freq, atlas, ev, eh, stdv, stdh, covvh, verb)

  CALL print_emis(ev, eh, stdv, stdh, covvh, 'Single freq, with spatial averaging, native resol.')


  ! Multiple frequencies, spatially averaged
  CALL emis_interp_int_mult(lat, lon, resol, theta, freq2, nchan, atlas, ev2, eh2, std, verb)

  CALL print_emis(ev2(i), eh2(i), SQRT(std(i,i)), SQRT(std(i+nchan,i+nchan)), &
                  std(i,i+nchan), 'Multiple freq, with spatial averaging, native resol.')



  !----------------------------------------------------------------------------
  ! Now make another call with spatial averaging active
  !----------------------------------------------------------------------------

  WRITE(0,'(a)') ' '
  WRITE(0,'(a)') 'Now the spatial averaging is active and results are different to those above:'
  WRITE(0,'(a)') ' '

  resol = 0.8

  ! Single frequency, spatially averaged
  CALL emis_interp_int_sing(lat, lon, resol, theta, freq, atlas, ev, eh, stdv, stdh, covvh, verb)

  CALL print_emis(ev, eh, stdv, stdh, covvh, 'Single freq, with spatial averaging, non-native resol.')


  ! Multiple frequencies, spatially averaged
  CALL emis_interp_int_mult(lat, lon, resol, theta, freq2, nchan, atlas, ev2, eh2, std, verb)

  CALL print_emis(ev2(i), eh2(i), SQRT(std(i,i)), SQRT(std(i+nchan,i+nchan)), &
                  std(i,i+nchan), 'Multiple freq, with spatial averaging, non-native resol.')



  !----------------------------------------------------------------------------
  ! We can make a call only emissivities and display results for all frequencies
  !----------------------------------------------------------------------------

  WRITE(0,'(a)') ' '
  WRITE(0,'(a)') 'Now only return emissivities and print them for all frequencies:'
  WRITE(0,'(a)') ' '

  resol = 0.8

  ! Multiple frequencies, spatially averaged, return emissivities only
  CALL emis_interp_int_mult(lat, lon, resol, theta, freq2, nchan, atlas, ev2, eh2, verb = verb)

  WRITE(0,'(a)') 'Multiple freq, with spatial averaging, non-native resol.'
  WRITE(0,'(a,10f10.3)')  'Freq (GHz) : ', freq2
  WRITE(0,'(a,10f10.6)') 'Emis V-pol = ', ev2
  WRITE(0,'(a,10f10.6)') 'Emis H-pol = ', eh2
  WRITE(0,'(a)') ' '


  !====================================================
  ! Close the atlas
  !====================================================
  CALL rttov_closemw_atlas(atlas)

CONTAINS

  SUBROUTINE print_emis(ev, eh, stdv, stdh, cov, title)
    REAL(jprb),       INTENT(IN) :: ev, eh, stdv, stdh, cov
    CHARACTER(LEN=*), INTENT(IN) :: title

    WRITE(0,'(a)') title
    WRITE(0,'(a,2f10.6)') 'Emis V-pol, H-pol   = ', ev, eh
    WRITE(0,'(a,2f10.6)') 'Stddev V-pol, H-pol = ', stdv, stdh
    WRITE(0,'(a,f10.6)')  'Covariance V-/H-pol = ', covvh
    WRITE(0,'(a)') ' '
  END SUBROUTINE print_emis

END PROGRAM test_telsem

