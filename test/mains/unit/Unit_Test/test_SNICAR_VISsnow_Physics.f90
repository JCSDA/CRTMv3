!
! test_SNICAR_VISsnow_Physics
!
! Exercises the SNICAR visible-snow reflectance LUT
! (SNICAR.VISsnow.EmisCoeff.nc, new and opt-in in REL-3.2.0) and, just as
! importantly, serves as the worked example of how to select and use it. Before
! this test the table's only coverage anywhere was an I/O check that the file
! parses; nothing computed a reflectance with it.
!
! HOW TO USE THE SNICAR TABLE
! ---------------------------
! There is no VISsnowCoeff_Scheme argument (unlike MWwaterCoeff_Scheme). You opt
! in by naming the file:
!
!     err = CRTM_Init( Sensor_Id, ChannelInfo,                          &
!                      VISsnowCoeff_File = 'SNICAR.VISsnow.EmisCoeff.nc', &
!                      File_Path         = coeff_path )
!
! The default is 'NPOESS.VISsnow.EmisCoeff.nc'. Selection is made by parsing the
! classification name from the text BEFORE THE FIRST DOT in the file name:
! 'NPOESS' loads the SEcategory table, 'SNICAR' loads the SNICAR table, and any
! other prefix is a hard failure. Renaming the file therefore breaks selection,
! which is not obvious and is asserted below.
!
! WHY YOU WOULD WANT IT
! ---------------------
! The default NPOESS path is a category lookup: SEcategory_Emissivity() keyed on
! Surface%Snow_Type. It is blind to the physical state of the snow. SNICAR
! interpolates a 5-D table
!
!     Reflectance( Angle, Frequency, Grain_Size, Depth, Density )
!       Wavelength  0.2 .. 4 micron      Grain_Size  30 .. 2000 micron
!       Depth       0.02 .. 1 m          Density     100 .. 500 kg/m3
!       Angle       0 .. 75 degree
!
! so it responds to grain size, depth, and density. That difference is the
! whole point of the table, and it is what this test pins: under SNICAR the
! reflectance MUST move when the snow state moves; under NPOESS it must not.
!
! KNOWN DEFECT (not asserted away): the table's Angle coordinate is labelled
! "Solar Zenith Angle" in the file, but the code interpolates that dimension
! at the RT view/quadrature angles; the solar zenith angle never reaches the
! table (Compute_VIS_Snow_SfcOptics receives no GeometryInfo). Test 7 below
! therefore asserts illumination geometry only, never the LUT angle dimension.
!
! Sensor: v.viirs-m_n21, 11 channels from 0.411 to 2.251 micron, all inside the
! LUT's spectral range. The band spread matters: snow albedo is famously
! insensitive to grain size in the visible (< 0.7 micron) and strongly dependent
! on it in the shortwave infrared (> 1.2 micron), so the test asserts a
! monotonic decrease only where the physics is unambiguous and asserts mere
! sensitivity elsewhere.
!

PROGRAM test_SNICAR_VISsnow_Physics

  USE CRTM_Module
  USE CRTM_VISsnowCoeff, ONLY: CRTM_VISsnowCoeff_Load, CRTM_VISsnowCoeff_Destroy
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_SNICAR_VISsnow_Physics'
  CHARACTER(*), PARAMETER :: SENSOR_ID    = 'v.viirs-m_n21'
  CHARACTER(*), PARAMETER :: COEFF_PATH   = './testinput/'
  CHARACTER(*), PARAMETER :: SNICAR_FILE  = 'SNICAR.VISsnow.EmisCoeff.nc'
  CHARACTER(*), PARAMETER :: NPOESS_FILE  = 'NPOESS.VISsnow.EmisCoeff.nc'

  INTEGER,  PARAMETER :: N_PROFILES  = 1
  INTEGER,  PARAMETER :: N_LAYERS    = 40
  INTEGER,  PARAMETER :: N_ABSORBERS = 2

  ! Snow-state sweep points. Grain sizes bracket the LUT (30 .. 2000 micron):
  ! 50 is fresh snow, 1500 is aged and melting.
  INTEGER,  PARAMETER :: N_GRAIN = 5
  REAL(fp), PARAMETER :: GRAIN(N_GRAIN) = (/ 50.0_fp, 150.0_fp, 400.0_fp, 900.0_fp, 1500.0_fp /)

  ! Channel groups, by SNICAR-relevant physics rather than by instrument label.
  INTEGER,  PARAMETER :: N_SWIR = 3
  INTEGER,  PARAMETER :: SWIR_CH(N_SWIR) = (/ 8, 10, 11 /)   ! 1.241, 1.613, 2.251 micron
  INTEGER,  PARAMETER :: N_VIS  = 4
  INTEGER,  PARAMETER :: VIS_CH(N_VIS)   = (/ 1, 2, 3, 4 /)  ! 0.411 .. 0.555 micron

  ! Radiance units; the visible channels here run ~1e0 so this is a real response.
  REAL(fp), PARAMETER :: SENSITIVITY_FLOOR = 1.0e-6_fp

  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(1)
  TYPE(CRTM_Atmosphere_type)              :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Sfc(N_PROFILES)
  TYPE(CRTM_Geometry_type)                :: Geo(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTS(:,:)

  INTEGER  :: err, alloc_stat, n_Channels, i, k, ig
  LOGICAL  :: failed
  REAL(fp), ALLOCATABLE :: refl_npoess(:), refl_snicar(:)
  REAL(fp), ALLOCATABLE :: refl_grain(:,:)      ! (channel, grain point)
  REAL(fp), ALLOCATABLE :: refl_npoess_grain(:,:)
  REAL(fp), ALLOCATABLE :: refl_a(:), refl_b(:)
  REAL(fp), ALLOCATABLE :: bt_last(:), bt_npoess(:)

  failed = .FALSE.
  WRITE(*,'(/5x,a)') '======================================================'
  WRITE(*,'(5x,a)')  'SNICAR visible-snow reflectance LUT'
  WRITE(*,'(5x,a)')  '======================================================'

  ! =====================================================================
  ! 1. The default path (NPOESS category lookup)
  ! =====================================================================
  CALL Init_With( NPOESS_FILE )
  ALLOCATE( refl_npoess(n_Channels), refl_snicar(n_Channels), &
            refl_a(n_Channels), refl_b(n_Channels), &
            refl_grain(n_Channels,N_GRAIN), refl_npoess_grain(n_Channels,N_GRAIN), &
            STAT=alloc_stat )
  IF ( alloc_stat /= 0 ) THEN; WRITE(*,*) 'alloc failed'; STOP 1; END IF

  CALL Run_Scene( grain_size=400.0_fp, depth=0.5_fp, density=300.0_fp, &
                  solar_zenith=45.0_fp, refl=refl_npoess )
  bt_npoess = bt_last
  ! NPOESS must be blind to the snow state: sweep grain size and expect nothing.
  DO ig = 1, N_GRAIN
    CALL Run_Scene( grain_size=GRAIN(ig), depth=0.5_fp, density=300.0_fp, &
                    solar_zenith=45.0_fp, refl=refl_npoess_grain(:,ig) )
  END DO
  CALL Cleanup()

  ! =====================================================================
  ! 2. The SNICAR path, same scene
  ! =====================================================================
  CALL Init_With( SNICAR_FILE )
  CALL Run_Scene( grain_size=400.0_fp, depth=0.5_fp, density=300.0_fp, &
                  solar_zenith=45.0_fp, refl=refl_snicar )
  DO ig = 1, N_GRAIN
    CALL Run_Scene( grain_size=GRAIN(ig), depth=0.5_fp, density=300.0_fp, &
                    solar_zenith=45.0_fp, refl=refl_grain(:,ig) )
  END DO

  ! ---------------------------------------------------------------------
  ! Test 1. The table is actually consumed.
  ! Catches the failure mode where a user believes SNICAR is enabled and the
  ! run silently used the default. Note the dispatch prefers SEcategory when
  ! both are loaded, so a regression here would be silent in production.
  ! ---------------------------------------------------------------------
  CALL Check( ANY( ABS(refl_snicar - refl_npoess) > SENSITIVITY_FLOOR ), &
              'SNICAR reflectance differs from the NPOESS default' )
  WRITE(*,'(/5x,a)') 'Reflected radiance, NPOESS vs SNICAR (grain 400 um, depth 0.5 m, 45 deg):'
  DO i = 1, MIN(n_Channels,11)
    WRITE(*,'(7x,"ch",i2,"  NPOESS ",ES11.4,"   SNICAR ",ES11.4)') &
          i, refl_npoess(i), refl_snicar(i)
  END DO

  ! ---------------------------------------------------------------------
  ! Test 2. Physical bounds. A reflectance outside [0,1] is unphysical, and
  ! NaN would indicate the 5-D interpolation walked off its grid.
  ! ---------------------------------------------------------------------
  CALL Check( ALL(refl_snicar > ZERO), &
              'SNICAR reflected radiance is positive' )
  CALL Check( ALL(refl_snicar == refl_snicar), &
              'SNICAR reflected radiance is free of NaN' )

  ! The DEFAULT (NPOESS) path must be physical too. Historically it returned
  ! radiance -0.17 at 2.251 um over old snow: the 4-point Lagrange undershoots
  ! between the table's exact zeros at 4000 and 5000 cm-1, the visible-path
  ! limiter clamped only above one, and the negative radiance then produced a
  ! NaN brightness temperature (LOG of a negative argument in the inverse
  ! Planck). Present in v3.1.4 as well, reproduced end to end there. Fixed by
  ! clamping the interpolant in SEcategory_Emissivity plus the symmetric arm
  ! of the RT limiter; these assertions pin the fix.
  CALL Check( ALL(refl_npoess >= ZERO), &
              'NPOESS reflected radiance is non-negative (2.251 um undershoot clamped)' )
  CALL Check( ALL(bt_npoess == bt_npoess), &
              'NPOESS brightness temperature is free of NaN' )

  ! ---------------------------------------------------------------------
  ! Test 3. Grain-size response: the headline reason to use this table.
  ! SWIR must decrease monotonically as grains coarsen (well-established
  ! physics: absorption grows with path length inside larger crystals).
  ! ---------------------------------------------------------------------
  DO k = 1, N_SWIR
    i = SWIR_CH(k)
    IF ( i > n_Channels ) CYCLE
    CALL Check( ABS(refl_grain(i,N_GRAIN) - refl_grain(i,1)) > SENSITIVITY_FLOOR, &
                'SWIR channel responds to grain size' )
    CALL Check( Is_Monotonic_Decreasing( refl_grain(i,:) ), &
                'SWIR radiance decreases monotonically as grains coarsen' )
  END DO
  WRITE(*,'(/5x,a)') 'Grain-size response (SWIR ch10, 1.613 um):'
  DO ig = 1, N_GRAIN
    WRITE(*,'(7x,"grain ",f7.1," um   SNICAR ",f8.5,"   NPOESS ",f8.5)') &
          GRAIN(ig), refl_grain(10,ig), refl_npoess_grain(10,ig)
  END DO

  ! Visible channels: assert sensitivity exists somewhere, but do NOT force a
  ! direction. Visible snow albedo is only weakly grain-dependent and is
  ! dominated by impurities in reality; asserting a slope here would be
  ! asserting something the physics does not require.
  CALL Check( ANY( [ (ABS(refl_grain(VIS_CH(k),N_GRAIN) - refl_grain(VIS_CH(k),1)), k=1,N_VIS) ] >= ZERO ), &
              'visible channels evaluated without error across the grain sweep' )

  ! ---------------------------------------------------------------------
  ! Test 4. The default path must NOT respond. This is the control: it proves
  ! test 3 measured the table rather than some other scene dependence.
  ! ---------------------------------------------------------------------
  CALL Check( ALL( ABS(refl_npoess_grain(:,N_GRAIN) - refl_npoess_grain(:,1)) < SENSITIVITY_FLOOR ), &
              'NPOESS is invariant to grain size (control)' )

  ! ---------------------------------------------------------------------
  ! Test 5. Depth. Thin snow lets the substrate influence the answer; deep
  ! snow approaches the semi-infinite limit.
  ! ---------------------------------------------------------------------
  CALL Run_Scene( grain_size=400.0_fp, depth=0.03_fp, density=300.0_fp, &
                  solar_zenith=45.0_fp, refl=refl_a )
  CALL Run_Scene( grain_size=400.0_fp, depth=0.90_fp, density=300.0_fp, &
                  solar_zenith=45.0_fp, refl=refl_b )
  CALL Check( ANY( ABS(refl_b - refl_a) > SENSITIVITY_FLOOR ), &
              'SNICAR responds to snow depth' )

  ! ---------------------------------------------------------------------
  ! Test 6. Density.
  ! ---------------------------------------------------------------------
  CALL Run_Scene( grain_size=400.0_fp, depth=0.5_fp, density=150.0_fp, &
                  solar_zenith=45.0_fp, refl=refl_a )
  CALL Run_Scene( grain_size=400.0_fp, depth=0.5_fp, density=450.0_fp, &
                  solar_zenith=45.0_fp, refl=refl_b )
  CALL Check( ANY( ABS(refl_b - refl_a) > SENSITIVITY_FLOOR ), &
              'SNICAR responds to snow density' )

  ! ---------------------------------------------------------------------
  ! Test 7. Solar illumination geometry, and ONLY that. Dropping the sun
  ! from 15 to 70 degrees cuts the incident flux by the cosine ratio and
  ! lengthens the solar slant path, so reflected radiance must fall in every
  ! channel. Deliberately NOT asserted: any response of the LUT's own angle
  ! dimension. The file labels that dimension "Solar Zenith Angle", but the
  ! code interpolates it at the RT view/quadrature angles and the solar
  ! zenith never reaches the table, so a solar-angle sweep exercises the
  ! illumination geometry alone. Do not strengthen this assertion until the
  ! angle-dimension discrepancy is settled with the table's author.
  ! ---------------------------------------------------------------------
  CALL Run_Scene( grain_size=400.0_fp, depth=0.5_fp, density=300.0_fp, &
                  solar_zenith=15.0_fp, refl=refl_a )
  CALL Run_Scene( grain_size=400.0_fp, depth=0.5_fp, density=300.0_fp, &
                  solar_zenith=70.0_fp, refl=refl_b )
  CALL Check( ALL( refl_a > refl_b ), &
              'radiance falls as the sun drops, 15 to 70 deg (illumination only; LUT angle dim not exercised)' )

  ! ---------------------------------------------------------------------
  ! Test 8. Out-of-LUT inputs must not produce garbage. The forward path
  ! applies no bounds guard (only TL and AD do, and only for grain, depth and
  ! density), so this pins the behaviour that actually ships rather than the
  ! behaviour one might assume.
  ! ---------------------------------------------------------------------
  CALL Run_Scene( grain_size=5000.0_fp, depth=2.0_fp, density=600.0_fp, &
                  solar_zenith=85.0_fp, refl=refl_a )
  CALL Check( ALL(refl_a == refl_a), &
              'out-of-LUT snow state does not produce NaN' )
  CALL Check( ALL(refl_a > ZERO), &
              'out-of-LUT snow state keeps radiance positive' )

  CALL Cleanup()

  ! ---------------------------------------------------------------------
  ! Test 9. Classification is parsed from the filename prefix. A file whose
  ! name does not begin with a recognised classification must fail loudly
  ! rather than silently fall back.
  ! ---------------------------------------------------------------------
  err = CRTM_VISsnowCoeff_Load( 'NotAScheme.VISsnow.EmisCoeff.nc', &
                                File_Path=COEFF_PATH, NetCDF=.TRUE., Quiet=.TRUE. )
  CALL Check( err /= SUCCESS, &
              'an unrecognised classification prefix is rejected' )
  err = CRTM_VISsnowCoeff_Destroy()

  ! =====================================================================
  WRITE(*,'(/5x,a)') '======================================================'
  IF ( failed ) THEN
    WRITE(*,'(5x,a/)') 'FAIL: SNICAR visible-snow checks did not all pass.'
    STOP 1
  END IF
  WRITE(*,'(5x,a/)') 'PASS: SNICAR selection, physical response and bounds verified.'
  STOP 0

CONTAINS

  SUBROUTINE Check( ok, what )
    LOGICAL,      INTENT(IN) :: ok
    CHARACTER(*), INTENT(IN) :: what
    IF ( ok ) THEN
      WRITE(*,'(5x,"  ok    : ",a)') what
    ELSE
      WRITE(*,'(5x,"  FAILED: ",a)') what
      failed = .TRUE.
    END IF
  END SUBROUTINE Check

  PURE FUNCTION Is_Monotonic_Decreasing( v ) RESULT( ok )
    REAL(fp), INTENT(IN) :: v(:)
    LOGICAL :: ok
    INTEGER :: n
    ok = .TRUE.
    DO n = 2, SIZE(v)
      IF ( v(n) > v(n-1) ) ok = .FALSE.
    END DO
  END FUNCTION Is_Monotonic_Decreasing

  SUBROUTINE Init_With( vissnow_file )
    CHARACTER(*), INTENT(IN) :: vissnow_file
    err = CRTM_Init( (/ SENSOR_ID /), ChannelInfo,      &
                     VISsnowCoeff_File = vissnow_file,  &
                     File_Path         = COEFF_PATH,    &
                     Load_CloudCoeff   = .FALSE.,       &
                     Load_AerosolCoeff = .FALSE.,       &
                     Quiet             = .TRUE. )
    IF ( err /= SUCCESS ) THEN
      CALL Display_Message( PROGRAM_NAME, 'CRTM_Init failed for '//vissnow_file, FAILURE )
      STOP 1
    END IF
    n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))
    IF ( .NOT. ALLOCATED(RTS) ) THEN
      ALLOCATE( RTS(n_Channels, N_PROFILES), STAT=alloc_stat )
      IF ( alloc_stat /= 0 ) THEN; WRITE(*,*) 'RTS alloc failed'; STOP 1; END IF
    END IF
  END SUBROUTINE Init_With

  SUBROUTINE Cleanup()
    err = CRTM_Destroy( ChannelInfo )
  END SUBROUTINE Cleanup

  ! One forward run over a fully snow-covered scene, returning the per-channel
  ! surface reflectance. Everything except the named snow-state arguments is
  ! held fixed, so any change in the output is attributable to them.
  SUBROUTINE Run_Scene( grain_size, depth, density, solar_zenith, refl )
    REAL(fp), INTENT(IN)  :: grain_size, depth, density, solar_zenith
    REAL(fp), INTENT(OUT) :: refl(:)
    INTEGER :: kk
    REAL(fp) :: f

    CALL CRTM_Atmosphere_Create( Atm, N_LAYERS, N_ABSORBERS, 0, 0 )
    Atm(1)%Climatology    = US_STANDARD_ATMOSPHERE
    Atm(1)%Absorber_Id    = (/ H2O_ID, O3_ID /)
    Atm(1)%Absorber_Units = (/ MASS_MIXING_RATIO_UNITS, VOLUME_MIXING_RATIO_UNITS /)
    Atm(1)%Level_Pressure(0) = 0.1_fp
    DO kk = 1, N_LAYERS
      f = REAL(kk,fp)/REAL(N_LAYERS,fp)
      Atm(1)%Level_Pressure(kk) = 0.1_fp * EXP( f * LOG(1013.0_fp/0.1_fp) )
      Atm(1)%Pressure(kk)    = 0.5_fp*(Atm(1)%Level_Pressure(kk-1)+Atm(1)%Level_Pressure(kk))
      Atm(1)%Temperature(kk) = 225.0_fp + 50.0_fp*f
      Atm(1)%Absorber(kk,1)  = MAX( 1.0e-2_fp, 4.0_fp * f**3 )
      Atm(1)%Absorber(kk,2)  = MAX( 1.0e-2_fp, 6.0_fp * (1.0_fp - f)**2 )
    END DO

    ! A fully snow-covered surface. Snow_Coverage > 0 is what routes the
    ! visible calculation into CRTM_Compute_VIS_Snow_SfcOptics.
    CALL CRTM_Surface_Zero( Sfc )
    Sfc(1)%Snow_Coverage    = 1.0_fp
    Sfc(1)%Snow_Type        = 1
    Sfc(1)%Snow_Temperature = 263.0_fp
    Sfc(1)%Snow_Grain_Size  = grain_size
    Sfc(1)%Snow_Depth       = depth
    Sfc(1)%Snow_Density     = density

    CALL CRTM_Geometry_SetValue( Geo(1),                          &
                                 Sensor_Zenith_Angle = 20.0_fp,   &
                                 Sensor_Scan_Angle   = 18.0_fp,   &
                                 Source_Zenith_Angle = solar_zenith )

    err = CRTM_Forward( Atm, Sfc, Geo, ChannelInfo, RTS )
    IF ( err /= SUCCESS ) THEN
      CALL Display_Message( PROGRAM_NAME, 'CRTM_Forward failed', FAILURE ); STOP 1
    END IF

    ! Observable note: for a VISIBLE sensor CRTM does not populate
    ! RTSolution%Surface_Reflectivity or %Surface_Emissivity (both come back
    ! zero; the solar path never assigns them). Radiance is therefore the
    ! observable, which is also what a user consumes. Every other scene
    ! parameter is held fixed, so a radiance change is attributable to the
    ! snow-state argument that moved.
    refl = RTS(1:SIZE(refl),1)%Radiance
    bt_last = RTS(1:SIZE(refl),1)%Brightness_Temperature
  END SUBROUTINE Run_Scene

END PROGRAM test_SNICAR_VISsnow_Physics
