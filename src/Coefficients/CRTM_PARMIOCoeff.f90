!
! CRTM_PARMIOCoeff
!
! Lifecycle wrapper for the shared PARMIOCoeff lookup table. Mirrors
! CRTM_MWwaterCoeff (which holds MWwaterC for FASTEM): exposes a single
! module-level SAVE record `PARMIOC` that the SfcOptics dispatcher hands to
! Compute_PARMIO. Loaded via CRTM_PARMIOCoeff_Load and freed via
! CRTM_PARMIOCoeff_Destroy.

MODULE CRTM_PARMIOCoeff

  USE Message_Handler,        ONLY: SUCCESS, FAILURE, Display_Message
  USE PARMIOCoeff_Define,     ONLY: PARMIOCoeff_type, &
                                    PARMIOCoeff_Associated, &
                                    PARMIOCoeff_Destroy
  USE PARMIOCoeff_netCDF_IO,  ONLY: PARMIOCoeff_netCDF_ReadFile

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: PARMIOC
  PUBLIC :: CRTM_PARMIOCoeff_Load
  PUBLIC :: CRTM_PARMIOCoeff_Destroy
  PUBLIC :: CRTM_PARMIOCoeff_IsLoaded

  INTEGER, PARAMETER :: ML = 512

  ! Shared PARMIO LUT data, populated by CRTM_PARMIOCoeff_Load.
  TYPE(PARMIOCoeff_type), TARGET, SAVE :: PARMIOC

CONTAINS

  FUNCTION CRTM_PARMIOCoeff_Load(Filename, Quiet) RESULT(err_stat)
    CHARACTER(*),      INTENT(IN) :: Filename
    LOGICAL, OPTIONAL, INTENT(IN) :: Quiet
    INTEGER :: err_stat
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_PARMIOCoeff_Load'
    err_stat = PARMIOCoeff_netCDF_ReadFile(PARMIOC, Filename, Quiet=Quiet)
    IF (err_stat /= SUCCESS) THEN
      CALL Display_Message(ROUTINE_NAME, &
        'Failed to load PARMIOCoeff from '//TRIM(Filename), FAILURE)
    END IF
  END FUNCTION CRTM_PARMIOCoeff_Load


  SUBROUTINE CRTM_PARMIOCoeff_Destroy()
    CALL PARMIOCoeff_Destroy(PARMIOC)
  END SUBROUTINE CRTM_PARMIOCoeff_Destroy


  PURE FUNCTION CRTM_PARMIOCoeff_IsLoaded() RESULT(is_loaded)
    LOGICAL :: is_loaded
    is_loaded = PARMIOCoeff_Associated(PARMIOC)
  END FUNCTION CRTM_PARMIOCoeff_IsLoaded

END MODULE CRTM_PARMIOCoeff
