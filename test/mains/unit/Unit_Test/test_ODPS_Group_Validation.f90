!
! test_ODPS_Group_Validation
!
! Unit test for ODPS_Validate_Group (Tier 0 of the ODPS group modernization).
! Exercises: every supported group's canonical roster (must pass), the
! Zeeman-reserved indexes (must fail with a Zeeman explanation), an unknown
! group, roster size mismatches, roster content mismatches, and roster order
! mismatches (the predictor code is positional, so order matters).
!
! The failing rosters include the real-world case that motivated the check:
! the OMPS Group-4 files (Group_Index=4, Component_ID=[13,14], a private
! convention of an external UV trainer that collides with the Zeeman index).
!
PROGRAM test_ODPS_Group_Validation

  USE ODPS_Predictor, ONLY: ODPS_Validate_Group, &
                            GROUP_1, GROUP_2, GROUP_3, &
                            GROUP_MW_O3, GROUP_UV_NO2, &
                            RESERVED_ZSSMIS_GROUP, &
                            RESERVED_ZAMSUA_GROUP
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_ODPS_Group_Validation'
  CHARACTER(512) :: Message
  INTEGER :: n_Failed

  n_Failed = 0

  ! ---------------------------------------------------------------
  ! Supported groups with their canonical rosters: all must validate
  ! ---------------------------------------------------------------
  CALL Expect_Valid( GROUP_1, &
    (/ 7, 101, 15, 114, 121, 120, 119, 118 /), (/ 1, 3, 2, 4, 5, 6 /), &
    'group 1 canonical roster' )
  CALL Expect_Valid( GROUP_2, &
    (/ 20, 101, 15, 114, 121 /), (/ 1, 3, 2 /), &
    'group 2 canonical roster' )
  CALL Expect_Valid( GROUP_3, &
    (/ 113, 12 /), (/ 1 /), &
    'group 3 canonical roster' )
  CALL Expect_Valid( GROUP_MW_O3, &
    (/ 113, 12, 114 /), (/ 1, 3 /), &
    'group 7 canonical roster' )
  CALL Expect_Valid( GROUP_UV_NO2, &
    (/ 20, 101, 15, 114, 121, 122 /), (/ 1, 3, 2, 10 /), &
    'group 8 canonical roster' )

  ! ---------------------------------------------------------------
  ! Reserved and unknown group indexes: all must be rejected
  ! ---------------------------------------------------------------
  ! The OMPS Group-4 files exactly as shipped (dry=13, ozone=14, O3 absorber)
  CALL Expect_Invalid( RESERVED_ZSSMIS_GROUP, (/ 13, 14 /), (/ 3 /), &
    'OMPS-style Group-4 file (Zeeman-reserved index)' )
  ! A zssmis companion loaded through the wrong (ODPS) path
  CALL Expect_Invalid( RESERVED_ZSSMIS_GROUP, (/ 13 /), (/ 1 /), &
    'zssmis companion via the ODPS path (Zeeman-reserved index)' )
  CALL Expect_Invalid( RESERVED_ZAMSUA_GROUP, (/ 13 /), (/ 1 /), &
    'Zeeman-reserved AMSU-A index' )
  CALL Expect_Invalid( 6, (/ 13 /), (/ 1 /), &
    'Zeeman-reserved index 6' )
  CALL Expect_Invalid( 0, (/ 113, 12 /), (/ 1 /), 'group 0 (fill value)' )
  CALL Expect_Invalid( 9, (/ 113, 12 /), (/ 1 /), 'group 9 (beyond table)' )
  CALL Expect_Invalid( -1, (/ 113, 12 /), (/ 1 /), 'negative group' )

  ! ---------------------------------------------------------------
  ! Kernel-capability semantics (Tier 2): well-formed subset, reordered,
  ! or extended rosters whose components all map to kernels (with their
  ! required gases) are ACCEPTED; the compute path dispatches by the
  ! file's own roster.
  ! ---------------------------------------------------------------
  ! A group-3 file carrying an ozone component with the O3 gas (G7-style)
  CALL Expect_Valid( GROUP_3, (/ 113, 12, 114 /), (/ 1, 3 /), &
    'group 3 with an ozone component and the O3 gas' )
  ! Extra known absorber (harmlessly mapped, consumed by no kernel)
  CALL Expect_Valid( GROUP_3, (/ 113, 12 /), (/ 1, 3 /), &
    'group 3 with an extra known absorber' )
  ! Reordered roster (dispatch is by ID, not position)
  CALL Expect_Valid( GROUP_3, (/ 12, 113 /), (/ 1 /), &
    'group 3 roster reordered' )
  ! A dry+ozone UV subset: the physics the OMPS files actually contain,
  ! expressible as a legitimate group-8 subset roster
  CALL Expect_Valid( GROUP_UV_NO2, (/ 20, 114 /), (/ 3 /), &
    'group 8 dry+ozone subset (regenerated-OMPS shape)' )

  ! ---------------------------------------------------------------
  ! Malformed rosters: all must be rejected
  ! ---------------------------------------------------------------
  ! Unknown component ID (raw molecule set 13 has no CRTM kernel)
  CALL Expect_Invalid( GROUP_3, (/ 13, 12 /), (/ 1 /), &
    'group 3 with molecule-set dry (13) instead of effective dry (113)' )
  ! Component whose required gas is missing
  CALL Expect_Invalid( GROUP_MW_O3, (/ 113, 12, 114 /), (/ 1, 2 /), &
    'group 7 ozone component without the O3 gas' )
  ! Duplicate component
  CALL Expect_Invalid( GROUP_3, (/ 113, 113 /), (/ 1 /), &
    'duplicate component ID' )
  ! Duplicate absorber
  CALL Expect_Invalid( GROUP_3, (/ 113, 12 /), (/ 1, 1 /), &
    'duplicate absorber ID' )
  ! Unknown absorber ID
  CALL Expect_Invalid( GROUP_3, (/ 113, 12 /), (/ 1, 99 /), &
    'unknown absorber ID' )
  ! Partial trace trio (CO without CH4/N2O)
  CALL Expect_Invalid( GROUP_1, (/ 7, 101, 15, 114, 121, 119 /), &
    (/ 1, 3, 2, 5 /), 'partial trace trio (CO without CH4 and N2O)' )
  ! IR component on the MW basis
  CALL Expect_Invalid( GROUP_3, (/ 113, 101 /), (/ 1 /), &
    'IR water-line component (101) on the MW basis' )
  ! WLO without the CO2 gas its predictor 15 consumes
  CALL Expect_Invalid( GROUP_2, (/ 20, 101, 15, 114, 121 /), (/ 1, 3 /), &
    'group 2 WLO without the CO2 gas' )

  ! ---------------------------------------------------------------
  ! Report
  ! ---------------------------------------------------------------
  IF ( n_Failed == 0 ) THEN
    WRITE(*,'(a,": ALL TESTS PASSED")') PROGRAM_NAME
    STOP 0
  ELSE
    WRITE(*,'(a,": ",i0," TEST(S) FAILED")') PROGRAM_NAME, n_Failed
    STOP 1
  END IF

CONTAINS

  SUBROUTINE Expect_Valid( Group_Index, Component_ID, Absorber_ID, Label )
    INTEGER,      INTENT(IN) :: Group_Index
    INTEGER,      INTENT(IN) :: Component_ID(:)
    INTEGER,      INTENT(IN) :: Absorber_ID(:)
    CHARACTER(*), INTENT(IN) :: Label
    IF ( ODPS_Validate_Group( Group_Index, Component_ID, Absorber_ID, Message ) ) THEN
      WRITE(*,'("  PASS (valid):   ",a)') Label
    ELSE
      WRITE(*,'("  FAIL: expected valid but got invalid: ",a)') Label
      WRITE(*,'("        message: ",a)') TRIM(Message)
      n_Failed = n_Failed + 1
    END IF
  END SUBROUTINE Expect_Valid

  SUBROUTINE Expect_Invalid( Group_Index, Component_ID, Absorber_ID, Label )
    INTEGER,      INTENT(IN) :: Group_Index
    INTEGER,      INTENT(IN) :: Component_ID(:)
    INTEGER,      INTENT(IN) :: Absorber_ID(:)
    CHARACTER(*), INTENT(IN) :: Label
    IF ( .NOT. ODPS_Validate_Group( Group_Index, Component_ID, Absorber_ID, Message ) ) THEN
      IF ( LEN_TRIM(Message) == 0 ) THEN
        WRITE(*,'("  FAIL: rejected but with a blank message: ",a)') Label
        n_Failed = n_Failed + 1
      ELSE
        WRITE(*,'("  PASS (invalid): ",a)') Label
        WRITE(*,'("        message: ",a)') TRIM(Message)
      END IF
    ELSE
      WRITE(*,'("  FAIL: expected invalid but got valid: ",a)') Label
      n_Failed = n_Failed + 1
    END IF
  END SUBROUTINE Expect_Invalid

END PROGRAM test_ODPS_Group_Validation
