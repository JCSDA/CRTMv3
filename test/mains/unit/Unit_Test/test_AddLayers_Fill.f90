!
! test_AddLayers_Fill
!
! Regression test for JCSDA/CRTMv3#403: CRTM_Atmosphere_AddLayers fills the
! atmosphere between a user profile top and the CRTM top of atmosphere
! (0.005 hPa) from the climatology model profile. A user top within rounding
! of an internal level used to add a zero-thickness layer, and the shift of
! the whole fill became meaningless (an 8.7 K change in the fill for a one-bit
! change in the top pressure).
!
! The user profiles are one climatology model (plus 3 K) paired with a
! different climatology, on thin layers (the model levels) and thick layers
! (every fourth model level), so the shift is not zero for temperature or for
! either absorber. For every case:
!
!   - the number of added layers is the expected one, and the forward,
!     tangent-linear and adjoint routines agree on it
!   - no added layer is thinner than the tolerance
!   - the fill is the climatology, shifted by one constant so that its lowest
!     layer equals the user's top two layers extrapolated linearly in log(p)
!     (computed here independently of the library)
!   - the user layers and the user top pressure come through unchanged
!
! Cases:
!
!   A. top on an internal level, tops from 0.016 to 20 hPa
!   B. top between internal levels
!   C. top one bit, a few parts in 1e9, and just inside the 1e-6 tolerance,
!      either side of an internal level: the same layers and the same fill as
!      case A. Just outside the tolerance: one more (thin but real) layer
!   D. top one bit above the top of atmosphere: nothing to add
!   E. top below every internal level: an error from all three routines
!
! No coefficient files are needed.
!
PROGRAM test_AddLayers_Fill

  USE CRTM_Module
  USE CRTM_Atmosphere    , ONLY: CRTM_Atmosphere_AddLayers, &
                                 CRTM_Atmosphere_AddLayers_TL, &
                                 CRTM_Atmosphere_AddLayers_AD
  USE CRTM_Model_Profiles, ONLY: CRTM_Get_Model_Profile
  USE CRTM_Parameters    , ONLY: MINIMUM_ABSORBER_AMOUNT
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_AddLayers_Fill'
  INTEGER , PARAMETER :: N_MODELS    = 6
  INTEGER , PARAMETER :: N_ABSORBERS = 2
  INTEGER , PARAMETER :: KMAX        = 100     ! model layers
  INTEGER , PARAMETER :: K_BOT       = 60      ! lowest level of the user profiles
  INTEGER , PARAMETER :: K_CUT_MAX   = 30      ! lowest user top tested (about 20 hPa)
  REAL(fp), PARAMETER :: T_TOL       = 1.0e-07_fp   ! K, fill against the rule
  REAL(fp), PARAMETER :: A_RTOL      = 1.0e-07_fp   ! relative, fill against the rule
  REAL(fp), PARAMETER :: MERGED_TOL  = 1.0e-04_fp   ! both, for a top merged with a level
  REAL(fp), PARAMETER :: T_SAME      = 1.0e-03_fp   ! K, fill for tops that differ by rounding
  REAL(fp), PARAMETER :: OFFSET      = 3.0_fp       ! K added to the user temperatures
  REAL(fp), PARAMETER :: NUDGE(6) = (/ 1.0e-09_fp, -1.0e-09_fp, 4.0e-09_fp, -4.0e-09_fp, &
                                       0.5e-06_fp,  0.9e-06_fp /)

  ! pl, tl, al: the climatology the fill is taken from; tu, au: the user profile
  REAL(fp) :: pl(0:KMAX), tl(0:KMAX), al(0:KMAX,N_ABSORBERS)
  REAL(fp) :: pu(0:KMAX), tu(0:KMAX), au(0:KMAX,N_ABSORBERS)
  REAL(fp) :: t_ref
  INTEGER  :: absorber_id(N_ABSORBERS)
  INTEGER  :: n_checks, n_failed
  INTEGER  :: model, user_model, k, i, step

  absorber_id = (/ H2O_ID, O3_ID /)
  n_checks = 0; n_failed = 0

  DO model = 1, N_MODELS
    user_model = MOD(model, N_MODELS) + 1
    CALL CRTM_Get_Model_Profile( absorber_id, pl, tl, al, Model=model )
    CALL CRTM_Get_Model_Profile( absorber_id, pu, tu, au, Model=user_model )
    tu = tu + OFFSET

    DO step = 1, 4, 3      ! thin and thick user layers
      DO k = 1, K_CUT_MAX

        ! A. top on an internal level. Its top fill layer is the reference for case C
        CALL Run_Case( 'A on-level', pl(k), k, step, .FALSE. )

        ! C. top within rounding, or within the tolerance, of that level:
        !    k layers, none of them a sliver, and the same fill
        CALL Run_Case( 'C one bit above', NEAREST(pl(k),  ONE), k, step, .TRUE. )
        CALL Run_Case( 'C one bit below', NEAREST(pl(k), -ONE), k, step, .TRUE. )
        DO i = 1, SIZE(NUDGE)
          CALL Run_Case( 'C inside tolerance', pl(k)*(ONE + NUDGE(i)), k, step, .TRUE. )
        END DO
        ! ...just outside the tolerance a thin but real layer is added
        CALL Run_Case( 'C outside tolerance', pl(k)*(ONE + 1.1e-06_fp), k+1, step, .FALSE. )
        CALL Run_Case( 'C outside tolerance', pl(k)*(ONE + 2.0e-06_fp), k+1, step, .FALSE. )
        CALL Run_Case( 'C outside tolerance', pl(k)*(ONE - 2.0e-06_fp), k  , step, .FALSE. )

        ! B. top between internal levels (the last added layer is a partial one)
        CALL Run_Case( 'B off-level', SQRT(pl(k-1)*pl(k)), k, step, .FALSE. )
      END DO
    END DO

    ! D. top at, and one bit above, the top of atmosphere: no layers to add
    CALL Run_Case( 'D TOA exact', pl(0), 0, 1, .FALSE. )
    CALL Run_Case( 'D TOA one bit above', NEAREST(pl(0), ONE), 0, 1, .FALSE. )

    ! E. top below every internal level: an error, not a profile
    CALL Run_Case( 'E below the model', pl(KMAX)*1.01_fp, -1, 1, .FALSE. )
  END DO

  WRITE(*,'(/,A,I0,A,I0,A)') PROGRAM_NAME//': ', n_checks, ' checks, ', n_failed, ' failed'
  IF ( n_failed > 0 ) THEN
    CALL Display_Message( PROGRAM_NAME, 'FAILED', FAILURE )
    STOP 1
  END IF
  CALL Display_Message( PROGRAM_NAME, 'PASSED', SUCCESS )

CONTAINS

  ! Profile value at pressure p, linear in log(p) between the model levels
  FUNCTION Level_Value( p, x ) RESULT( v )
    REAL(fp), INTENT(IN) :: p, x(0:)
    REAL(fp) :: v, w
    INTEGER  :: j
    DO j = 1, KMAX
      IF ( pl(j) >= p ) EXIT
    END DO
    j = MIN(j, KMAX)
    w = (LOG(p) - LOG(pl(j-1))) / (LOG(pl(j)) - LOG(pl(j-1)))
    v = x(j-1) + w*(x(j) - x(j-1))
  END FUNCTION Level_Value

  ! Mean pressure of a layer
  FUNCTION Layer_Pressure( p0, p1 ) RESULT( p )
    REAL(fp), INTENT(IN) :: p0, p1
    REAL(fp) :: p
    p = (p1 - p0)/LOG(p1/p0)
  END FUNCTION Layer_Pressure

  ! One case: a user profile with its top at ptop_in and levels every
  ! "step" model levels below it. n_expected is the number of layers the fill
  ! must add (-1: the call must fail). With same_as_ref, the top fill layer
  ! must also equal the one of the last case run without it.
  SUBROUTINE Run_Case( label, ptop_in, n_expected, step, same_as_ref )
    CHARACTER(*), INTENT(IN) :: label
    REAL(fp)    , INTENT(IN) :: ptop_in
    INTEGER     , INTENT(IN) :: n_expected, step
    LOGICAL     , INTENT(IN) :: same_as_ref
    TYPE(CRTM_Atmosphere_type) :: atm, atm_out, atm_TL, atm_out_TL, atm_AD, atm_out_AD
    REAL(fp) :: lev(0:KMAX), tlev(0:KMAX), alev(0:KMAX,N_ABSORBERS)
    REAL(fp) :: x_model(KMAX), p0, p1, p_join, w, x_join, x_rule, worst, thin, t_tolerance, a_tolerance
    INTEGER  :: err, n_user, n_added, j, l, k_first
    LOGICAL  :: ok
    CHARACTER(256) :: why

    ok = .TRUE.; why = ''
    ! A top merged with a level sits up to 1e-6 (relative) off it, and the
    ! library and this test then interpolate the climatology in different
    ! intervals: allow for that
    t_tolerance = T_TOL; a_tolerance = A_RTOL
    IF ( same_as_ref ) THEN
      t_tolerance = MERGED_TOL; a_tolerance = MERGED_TOL
    END IF

    ! User levels: ptop_in, then model levels below it
    IF ( n_expected < 0 ) THEN
      n_user = 2
      lev(0:2) = ptop_in*(/ ONE, 1.01_fp, 1.02_fp /)
      tlev(0:2) = tu(KMAX); alev(0:2,1) = au(KMAX,1); alev(0:2,2) = au(KMAX,2)
    ELSE
      DO k_first = 0, KMAX
        IF ( pl(k_first) > ptop_in*(ONE + 1.0e-03_fp) ) EXIT
      END DO
      n_user = (K_BOT - k_first)/step + 1
      lev(0)  = ptop_in
      tlev(0) = Level_Value( ptop_in, tu )
      DO j = 1, N_ABSORBERS
        alev(0,j) = Level_Value( ptop_in, au(:,j) )
      END DO
      DO l = 1, n_user
        lev(l)    = pl(k_first + (l-1)*step)
        tlev(l)   = tu(k_first + (l-1)*step)
        alev(l,:) = au(k_first + (l-1)*step,:)
      END DO
    END IF

    CALL CRTM_Atmosphere_Create( atm, n_user, N_ABSORBERS, 0, 0 )
    atm%Climatology    = model
    atm%Absorber_Id    = absorber_id
    atm%Absorber_Units = (/ MASS_MIXING_RATIO_UNITS, VOLUME_MIXING_RATIO_UNITS /)
    atm%Level_Pressure(0:n_user) = lev(0:n_user)
    DO l = 1, n_user
      atm%Pressure(l)    = Layer_Pressure( lev(l-1), lev(l) )
      atm%Temperature(l) = POINT_5*(tlev(l-1) + tlev(l))
      atm%Absorber(l,:)  = POINT_5*(alev(l-1,:) + alev(l,:))
    END DO
    CALL CRTM_Atmosphere_Create( atm_TL, n_user, N_ABSORBERS, 0, 0 )
    CALL CRTM_Atmosphere_Zero( atm_TL )
    atm_TL%Temperature = ONE

    err = CRTM_Atmosphere_AddLayers( atm, atm_out )

    IF ( n_expected < 0 ) THEN

      IF ( err == SUCCESS ) THEN
        ok = .FALSE.; why = 'no error from the forward routine'
      END IF
      err = CRTM_Atmosphere_AddLayers_TL( atm, atm_TL, atm_out_TL )
      IF ( err == SUCCESS ) THEN
        ok = .FALSE.; why = 'no error from the tangent-linear routine'
      END IF
      CALL CRTM_Atmosphere_Create( atm_AD, n_user, N_ABSORBERS, 0, 0 )
      CALL CRTM_Atmosphere_Create( atm_out_AD, n_user, N_ABSORBERS, 0, 0 )
      err = CRTM_Atmosphere_AddLayers_AD( atm, atm_out_AD, atm_AD )
      IF ( err == SUCCESS ) THEN
        ok = .FALSE.; why = 'no error from the adjoint routine'
      END IF

    ELSE IF ( err /= SUCCESS ) THEN

      ok = .FALSE.; why = 'AddLayers returned an error'

    ELSE

      n_added = atm_out%n_Layers - n_user
      IF ( n_added /= n_expected ) THEN
        ok = .FALSE.; WRITE(why,'(A,I0,A,I0)') 'layers added: ', n_added, ', expected ', n_expected
      ELSE IF ( n_added > 0 ) THEN

        thin = HUGE(ONE)
        DO l = 1, n_added
          thin = MIN(thin, (atm_out%Level_Pressure(l) - atm_out%Level_Pressure(l-1))/atm_out%Level_Pressure(l))
        END DO
        IF ( thin < 1.0e-06_fp ) THEN
          ok = .FALSE.; WRITE(why,'(A,ES10.3)') 'added layer with relative thickness ', thin
        END IF
        IF ( atm_out%Level_Pressure(n_added) /= lev(0) ) THEN
          ok = .FALSE.; why = 'user top level pressure not preserved'
        END IF

        ! The rule: climatology layer values (mean of the two level values),
        ! all shifted by one constant so that the lowest added layer equals
        ! the user's top two layers extrapolated in log(p) to its pressure
        p_join = Layer_Pressure( atm_out%Level_Pressure(n_added-1), atm_out%Level_Pressure(n_added) )
        w = (LOG(p_join) - LOG(atm%Pressure(1))) / (LOG(atm%Pressure(2)) - LOG(atm%Pressure(1)))
        DO j = 0, N_ABSORBERS
          DO l = 1, n_added
            p0 = atm_out%Level_Pressure(l-1); p1 = atm_out%Level_Pressure(l)
            IF ( j == 0 ) THEN
              x_model(l) = POINT_5*(Level_Value(p0, tl) + Level_Value(p1, tl))
            ELSE
              x_model(l) = POINT_5*(Level_Value(p0, al(:,j)) + Level_Value(p1, al(:,j)))
            END IF
          END DO
          worst = ZERO
          IF ( j == 0 ) THEN
            x_join = atm%Temperature(1) + w*(atm%Temperature(2) - atm%Temperature(1))
            DO l = 1, n_added
              x_rule = x_model(l) + (x_join - x_model(n_added))
              worst = MAX(worst, ABS(atm_out%Temperature(l) - x_rule))
            END DO
            IF ( worst > t_tolerance ) THEN
              ok = .FALSE.; WRITE(why,'(A,ES10.3,A)') 'fill temperature off the rule by ', worst, ' K'
            END IF
          ELSE
            x_join = atm%Absorber(1,j) + w*(atm%Absorber(2,j) - atm%Absorber(1,j))
            DO l = 1, n_added
              x_rule = x_model(l) + (x_join - x_model(n_added))
              IF ( x_rule < ZERO ) x_rule = MINIMUM_ABSORBER_AMOUNT
              worst = MAX(worst, ABS(atm_out%Absorber(l,j) - x_rule)/MAX(x_rule, x_model(l)))
            END DO
            IF ( worst > a_tolerance ) THEN
              ok = .FALSE.; WRITE(why,'(A,I0,A,ES10.3)') 'fill absorber ', j, ' off the rule by (relative) ', worst
            END IF
          END IF
        END DO

        ! A top that differs from an internal level by rounding gets the fill
        ! of the top on that level
        IF ( same_as_ref ) THEN
          IF ( ABS(atm_out%Temperature(1) - t_ref) > T_SAME ) THEN
            ok = .FALSE.
            WRITE(why,'(A,ES10.3,A)') 'fill differs from the on-level one by ', atm_out%Temperature(1) - t_ref, ' K'
          END IF
        ELSE
          t_ref = atm_out%Temperature(1)
        END IF

      END IF

      ! The user layers must come through unchanged
      IF ( n_added >= 0 ) THEN
        IF ( ANY(atm_out%Temperature(n_added+1:) /= atm%Temperature) .OR. &
             ANY(atm_out%Absorber(n_added+1:,:) /= atm%Absorber) ) THEN
          ok = .FALSE.; why = 'user layers changed'
        END IF
      END IF

      ! TL and AD must use the same layer count as the forward routine
      err = CRTM_Atmosphere_AddLayers_TL( atm, atm_TL, atm_out_TL )
      IF ( err /= SUCCESS ) THEN
        ok = .FALSE.; why = 'tangent-linear returned an error'
      ELSE IF ( atm_out_TL%n_Layers /= atm_out%n_Layers ) THEN
        ok = .FALSE.; why = 'tangent-linear layer count differs from forward'
      ELSE IF ( ANY(atm_out_TL%Temperature(n_added+1:) /= ONE) ) THEN
        ok = .FALSE.; why = 'tangent-linear did not map the user layers one to one'
      ELSE
        CALL CRTM_Atmosphere_Create( atm_AD, n_user, N_ABSORBERS, 0, 0 )
        CALL CRTM_Atmosphere_Zero( atm_AD )
        CALL CRTM_Atmosphere_Create( atm_out_AD, atm_out%n_Layers, N_ABSORBERS, 0, 0 )
        CALL CRTM_Atmosphere_Zero( atm_out_AD )
        atm_out_AD%Temperature = ONE
        err = CRTM_Atmosphere_AddLayers_AD( atm, atm_out_AD, atm_AD )
        IF ( err /= SUCCESS ) THEN
          ok = .FALSE.; why = 'adjoint returned an error'
        ELSE IF ( ANY(atm_AD%Temperature /= ONE) ) THEN
          ok = .FALSE.; why = 'adjoint did not map the user layers one to one'
        END IF
      END IF

    END IF

    n_checks = n_checks + 1
    IF ( .NOT. ok ) THEN
      n_failed = n_failed + 1
      WRITE(*,'(A,I0,A,I0,A,ES23.16,A)') 'FAIL  '//label//', model ', model, ', step ', step, &
        ', top ', ptop_in, ' hPa: '//TRIM(why)
    END IF

    CALL CRTM_Atmosphere_Destroy( atm ); CALL CRTM_Atmosphere_Destroy( atm_out )
    CALL CRTM_Atmosphere_Destroy( atm_TL ); CALL CRTM_Atmosphere_Destroy( atm_out_TL )
    CALL CRTM_Atmosphere_Destroy( atm_AD ); CALL CRTM_Atmosphere_Destroy( atm_out_AD )
  END SUBROUTINE Run_Case

END PROGRAM test_AddLayers_Fill
