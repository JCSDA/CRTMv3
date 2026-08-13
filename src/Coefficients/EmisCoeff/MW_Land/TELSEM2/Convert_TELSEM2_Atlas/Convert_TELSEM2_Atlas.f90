! Convert_TELSEM2_Atlas
!
! One-time tool: convert the 12 monthly TELSEM2 ASCII atlas files
! (ssmi_mean_emis_climato_MM_cov_interpol_M2) plus the 'correlations' file into
! a single CRTM-native netCDF TELSEM2 atlas coefficient file (Release 2)
! holding all twelve months and the uncertainty content.
!
! Usage:
!   Convert_TELSEM2_Atlas <atlas_dir> <output_netcdf_file>
!
! <atlas_dir> must contain the 12 monthly ASCII files and, beside them, the
! 'correlations' file (the standard RTTOV TELSEM2 distribution layout). Each
! ASCII record holds cellnum, 7 emissivities, 7 emissivity error VARIANCES and
! the two class indices; the variances are stored as standard deviations
! (SQRT applied here, matching RTTOV's reader).
!
PROGRAM Convert_TELSEM2_Atlas

  USE Type_Kinds              , ONLY: fp, Long
  USE Message_Handler         , ONLY: SUCCESS
  USE TELSEM2Atlas_Define     , ONLY: TELSEM2Atlas_type, TELSEM2Atlas_Create, &
                                      TELSEM2Atlas_Create_Uncertainty, &
                                      TELSEM2Atlas_Inspect, N_CLASS1
  USE TELSEM2Atlas_netCDF_IO  , ONLY: TELSEM2Atlas_netCDF_WriteFile

  IMPLICIT NONE

  INTEGER,      PARAMETER :: N_MONTHS   = 12
  INTEGER,      PARAMETER :: N_CHANNELS = 7
  INTEGER,      PARAMETER :: N_BANDS    = 720      ! 180 / 0.25
  INTEGER,      PARAMETER :: N_CELLS    = 660066
  REAL(fp),     PARAMETER :: RESOLUTION = 0.25_fp
  CHARACTER(*), PARAMETER :: BASENAME   = 'ssmi_mean_emis_climato_'
  CHARACTER(*), PARAMETER :: SUFFIX     = '_cov_interpol_M2'
  CHARACTER(*), PARAMETER :: CORRELNAME = 'correlations'

  TYPE(TELSEM2Atlas_type) :: atlas
  CHARACTER(512) :: atlas_dir, out_file, fname
  CHARACTER(2)   :: mm
  INTEGER        :: month_count(N_MONTHS), month_offset(N_MONTHS)
  INTEGER        :: m, j, i, ndat, ipos, p, unit, ios, cellnum, c1, c2
  INTEGER(Long)  :: n_total
  REAL(fp)       :: ssmi(14)
  INTEGER        :: nargs, err
  INTEGER        :: n_var_clips

  nargs = COMMAND_ARGUMENT_COUNT()
  IF ( nargs >= 1 ) THEN
    CALL GET_COMMAND_ARGUMENT(1, atlas_dir)
  ELSE
    WRITE(*,'(a)') 'Usage: Convert_TELSEM2_Atlas <atlas_dir> [<output_netcdf_file>]'
    STOP 1
  END IF
  IF ( nargs >= 2 ) THEN
    CALL GET_COMMAND_ARGUMENT(2, out_file)
  ELSE
    out_file = 'TELSEM2.MWland.EmisCoeff.nc'
  END IF
  ! Ensure trailing slash on the directory
  IF ( LEN_TRIM(atlas_dir) > 0 ) THEN
    IF ( atlas_dir(LEN_TRIM(atlas_dir):LEN_TRIM(atlas_dir)) /= '/' ) &
      atlas_dir = TRIM(atlas_dir)//'/'
  END IF

  ! ---- Pass 1: count kept cells per month ----
  WRITE(*,'(a)') 'Pass 1: counting populated cells per month...'
  DO m = 1, N_MONTHS
    WRITE(mm,'(I2.2)') m
    fname = TRIM(atlas_dir)//BASENAME//mm//SUFFIX
    OPEN(NEWUNIT=unit, FILE=TRIM(fname), STATUS='old', FORM='formatted', IOSTAT=ios)
    IF ( ios /= 0 ) THEN
      WRITE(*,'(a)') 'ERROR opening '//TRIM(fname); STOP 1
    END IF
    READ(unit,*) ndat
    ipos = 0
    DO j = 1, ndat
      READ(unit,*) cellnum, (ssmi(i),i=1,14), c1, c2
      IF ( c1 > 0 .AND. c2 > 0 ) ipos = ipos + 1
    END DO
    CLOSE(unit)
    month_count(m) = ipos
    WRITE(*,'(3x,"Month ",i2.2,": raw=",i8," kept=",i8)') m, ndat, ipos
  END DO

  ! Offsets and total
  month_offset(1) = 0
  DO m = 2, N_MONTHS
    month_offset(m) = month_offset(m-1) + month_count(m-1)
  END DO
  n_total = SUM(month_count)
  WRITE(*,'(a,i0)') 'Total stacked cells = ', n_total

  ! ---- Allocate ----
  CALL TELSEM2Atlas_Create( atlas, INT(N_CHANNELS,Long), INT(N_BANDS,Long), &
                            INT(N_CELLS,Long), INT(N_MONTHS,Long), n_total )
  atlas%Resolution      = RESOLUTION
  atlas%Month_Data_Count = INT(month_count, Long)
  atlas%Month_Offset     = INT(month_offset, Long)
  CALL TELSEM2Atlas_Create_Uncertainty( atlas )
  IF ( .NOT. atlas%Has_Uncertainty ) THEN
    WRITE(*,'(a)') 'ERROR allocating uncertainty arrays'; STOP 1
  END IF

  ! ---- Pass 2: fill ----
  WRITE(*,'(a)') 'Pass 2: reading emissivities and error variances...'
  n_var_clips = 0
  DO m = 1, N_MONTHS
    WRITE(mm,'(I2.2)') m
    fname = TRIM(atlas_dir)//BASENAME//mm//SUFFIX
    OPEN(NEWUNIT=unit, FILE=TRIM(fname), STATUS='old', FORM='formatted', IOSTAT=ios)
    IF ( ios /= 0 ) THEN
      WRITE(*,'(a)') 'ERROR opening '//TRIM(fname); STOP 1
    END IF
    READ(unit,*) ndat
    p = month_offset(m)
    DO j = 1, ndat
      READ(unit,*) cellnum, (ssmi(i),i=1,14), c1, c2
      IF ( c1 > 0 .AND. c2 > 0 ) THEN
        p = p + 1
        DO i = 1, N_CHANNELS
          atlas%Emissivity(p,i) = ssmi(i)
          ! Columns 8-14 are error variances; store the std, as RTTOV does
          ! at read time (mod_mwatlas_m2.F90: emis_err = sqrt(variance)).
          IF ( ssmi(N_CHANNELS+i) < 0.0_fp ) n_var_clips = n_var_clips + 1
          atlas%Emis_Err(p,i) = SQRT( MAX( ssmi(N_CHANNELS+i), 0.0_fp ) )
        END DO
        atlas%Cell_Number(p) = cellnum
        atlas%Class1(p)      = c1
        atlas%Class2(p)      = c2
      END IF
    END DO
    CLOSE(unit)
  END DO
  IF ( n_var_clips > 0 ) &
    WRITE(*,'(a,i0,a)') 'WARNING: ', n_var_clips, ' negative error variances clipped to zero'

  ! ---- Correlations: N_CLASS1 blocks of [1 header line + 7 rows of 7F5.2] ----
  WRITE(*,'(a)') 'Reading inter-channel correlations...'
  fname = TRIM(atlas_dir)//CORRELNAME
  OPEN(NEWUNIT=unit, FILE=TRIM(fname), STATUS='old', FORM='formatted', IOSTAT=ios)
  IF ( ios /= 0 ) THEN
    WRITE(*,'(a)') 'ERROR opening '//TRIM(fname); STOP 1
  END IF
  DO i = 1, N_CLASS1
    READ(unit,*,IOSTAT=ios)                       ! "class N" header line
    IF ( ios /= 0 ) THEN
      WRITE(*,'(a,i0)') 'ERROR reading correlations header, class ', i; STOP 1
    END IF
    DO j = 1, N_CHANNELS
      READ(unit,'(7F5.2)',IOSTAT=ios) atlas%Correlation(i,j,:)
      IF ( ios /= 0 ) THEN
        WRITE(*,'(a,i0,a,i0)') 'ERROR reading correlations, class ', i, ' row ', j; STOP 1
      END IF
    END DO
  END DO
  CLOSE(unit)

  CALL TELSEM2Atlas_Inspect( atlas )

  ! ---- Write netCDF ----
  WRITE(*,'(a)') 'Writing '//TRIM(out_file)//' ...'
  err = TELSEM2Atlas_netCDF_WriteFile( TRIM(out_file), atlas )
  IF ( err /= SUCCESS ) THEN
    WRITE(*,'(a)') 'ERROR writing netCDF file'; STOP 1
  END IF
  WRITE(*,'(a)') 'Done.'

END PROGRAM Convert_TELSEM2_Atlas
