! Convert_TELSEM2_Atlas
!
! One-time tool: convert the 12 monthly TELSEM2 ASCII atlas files
! (ssmi_mean_emis_climato_MM_cov_interpol_M2) into a single CRTM-native netCDF
! TELSEM2 atlas coefficient file holding all twelve months.
!
! Usage:
!   Convert_TELSEM2_Atlas <atlas_dir> <output_netcdf_file>
!
! Only the emissivity and the two surface-class indices are retained; the atlas
! uncertainty (variances / correlations) is not used by the CRTM surface optics.
!
PROGRAM Convert_TELSEM2_Atlas

  USE Type_Kinds              , ONLY: fp, Long
  USE Message_Handler         , ONLY: SUCCESS
  USE TELSEM2Atlas_Define     , ONLY: TELSEM2Atlas_type, TELSEM2Atlas_Create, &
                                      TELSEM2Atlas_Inspect
  USE TELSEM2Atlas_netCDF_IO  , ONLY: TELSEM2Atlas_netCDF_WriteFile

  IMPLICIT NONE

  INTEGER,      PARAMETER :: N_MONTHS   = 12
  INTEGER,      PARAMETER :: N_CHANNELS = 7
  INTEGER,      PARAMETER :: N_BANDS    = 720      ! 180 / 0.25
  INTEGER,      PARAMETER :: N_CELLS    = 660066
  REAL(fp),     PARAMETER :: RESOLUTION = 0.25_fp
  CHARACTER(*), PARAMETER :: BASENAME   = 'ssmi_mean_emis_climato_'
  CHARACTER(*), PARAMETER :: SUFFIX     = '_cov_interpol_M2'

  TYPE(TELSEM2Atlas_type) :: atlas
  CHARACTER(512) :: atlas_dir, out_file, fname
  CHARACTER(2)   :: mm
  INTEGER        :: month_count(N_MONTHS), month_offset(N_MONTHS)
  INTEGER        :: m, j, i, ndat, ipos, p, unit, ios, cellnum, c1, c2
  INTEGER(Long)  :: n_total
  REAL(fp)       :: ssmi(14)
  INTEGER        :: nargs, err

  nargs = COMMAND_ARGUMENT_COUNT()
  IF ( nargs >= 1 ) THEN
    CALL GET_COMMAND_ARGUMENT(1, atlas_dir)
  ELSE
    atlas_dir = '/home/jbenjam/CRTM/telsem2_atlas/'
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

  ! ---- Pass 2: fill ----
  WRITE(*,'(a)') 'Pass 2: reading emissivities...'
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
        END DO
        atlas%Cell_Number(p) = cellnum
        atlas%Class1(p)      = c1
        atlas%Class2(p)      = c2
      END IF
    END DO
    CLOSE(unit)
  END DO

  CALL TELSEM2Atlas_Inspect( atlas )

  ! ---- Write netCDF ----
  WRITE(*,'(a)') 'Writing '//TRIM(out_file)//' ...'
  err = TELSEM2Atlas_netCDF_WriteFile( TRIM(out_file), atlas )
  IF ( err /= SUCCESS ) THEN
    WRITE(*,'(a)') 'ERROR writing netCDF file'; STOP 1
  END IF
  WRITE(*,'(a)') 'Done.'

END PROGRAM Convert_TELSEM2_Atlas
