!
! CRTM_RTSolution_Diff
!
! Diagnostic helper used by the regression test suite when an
! RTSolution(channel, profile) array fails an exact-match comparison
! against a saved reference. Reports max / mean / RMS absolute
! differences for Brightness_Temperature and Radiance and lists the
! top offending (channel, profile) entries so the developer can judge
! whether a failure is a tiny floating-point divergence or a real
! algorithmic regression.
!
MODULE CRTM_RTSolution_Diff

  USE Type_Kinds                 , ONLY: fp
  USE CRTM_RTSolution_Define     , ONLY: CRTM_RTSolution_type

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Report_RTSolution_Diff

CONTAINS

  SUBROUTINE Report_RTSolution_Diff( actual, expected, label, top_n )
    TYPE(CRTM_RTSolution_type), INTENT(IN) :: actual(:,:)
    TYPE(CRTM_RTSolution_type), INTENT(IN) :: expected(:,:)
    CHARACTER(*), OPTIONAL,     INTENT(IN) :: label
    INTEGER,      OPTIONAL,     INTENT(IN) :: top_n

    INTEGER :: nL, nM, l, m, i, n_show, n_total
    REAL(fp) :: dBT, dRad, absBT, absRad
    REAL(fp) :: max_dBT, sum_dBT, sumsq_dBT
    REAL(fp) :: max_dRad, sum_dRad, sumsq_dRad
    INTEGER  :: max_dBT_l, max_dBT_m
    INTEGER  :: max_dRad_l, max_dRad_m
    INTEGER  :: n_BT_gt_1mK, n_BT_gt_10mK, n_BT_gt_100mK, n_BT_gt_1K
    REAL(fp), ALLOCATABLE :: rank_dBT(:)
    INTEGER,  ALLOCATABLE :: rank_l(:), rank_m(:)
    CHARACTER(64) :: tag

    tag = 'RTSolution_Diff'
    IF ( PRESENT(label) ) tag = label

    nL = SIZE(actual, DIM=1)
    nM = SIZE(actual, DIM=2)
    IF ( SIZE(expected,1) /= nL .OR. SIZE(expected,2) /= nM ) THEN
      WRITE(*,'(/5x,a,": shape mismatch (",i0,"x",i0,") vs (",i0,"x",i0,")")') &
        TRIM(tag), nL, nM, SIZE(expected,1), SIZE(expected,2)
      RETURN
    END IF

    n_total = nL * nM
    max_dBT  = 0.0_fp;  sum_dBT  = 0.0_fp;  sumsq_dBT  = 0.0_fp
    max_dRad = 0.0_fp;  sum_dRad = 0.0_fp;  sumsq_dRad = 0.0_fp
    max_dBT_l  = 0;  max_dBT_m  = 0
    max_dRad_l = 0;  max_dRad_m = 0
    n_BT_gt_1mK   = 0
    n_BT_gt_10mK  = 0
    n_BT_gt_100mK = 0
    n_BT_gt_1K    = 0

    ALLOCATE( rank_dBT(n_total), rank_l(n_total), rank_m(n_total) )
    i = 0

    DO m = 1, nM
      DO l = 1, nL
        dBT  = actual(l,m)%Brightness_Temperature - expected(l,m)%Brightness_Temperature
        dRad = actual(l,m)%Radiance               - expected(l,m)%Radiance
        absBT  = ABS(dBT)
        absRad = ABS(dRad)

        sum_dBT    = sum_dBT    + absBT
        sumsq_dBT  = sumsq_dBT  + absBT*absBT
        sum_dRad   = sum_dRad   + absRad
        sumsq_dRad = sumsq_dRad + absRad*absRad

        IF ( absBT  > max_dBT  ) THEN; max_dBT  = absBT;  max_dBT_l  = l; max_dBT_m  = m; END IF
        IF ( absRad > max_dRad ) THEN; max_dRad = absRad; max_dRad_l = l; max_dRad_m = m; END IF

        IF ( absBT > 1.0e-3_fp ) n_BT_gt_1mK   = n_BT_gt_1mK   + 1
        IF ( absBT > 1.0e-2_fp ) n_BT_gt_10mK  = n_BT_gt_10mK  + 1
        IF ( absBT > 1.0e-1_fp ) n_BT_gt_100mK = n_BT_gt_100mK + 1
        IF ( absBT > 1.0_fp    ) n_BT_gt_1K    = n_BT_gt_1K    + 1

        i = i + 1
        rank_dBT(i) = absBT
        rank_l(i)   = l
        rank_m(i)   = m
      END DO
    END DO

    n_show = 5
    IF ( PRESENT(top_n) ) n_show = top_n
    n_show = MIN(n_show, n_total)

    WRITE(*,'(/5x,a)') '==================== RTSolution diagnostic ===================='
    WRITE(*,'( 5x,a,": ",i0," channel(s) x ",i0," profile(s) = ",i0," entries")') &
      TRIM(tag), nL, nM, n_total

    WRITE(*,'(/5x,"Brightness_Temperature [K]")')
    WRITE(*,'( 7x,"max |diff| = ",es12.4,"  at (chan=",i0,", prof=",i0,")")') &
      max_dBT, max_dBT_l, max_dBT_m
    WRITE(*,'( 7x,"mean|diff| = ",es12.4,"   rms|diff| = ",es12.4)') &
      sum_dBT/REAL(n_total,fp), SQRT(sumsq_dBT/REAL(n_total,fp))
    WRITE(*,'( 7x,"channels with |dBT| > 1mK / 10mK / 100mK / 1K : ",i0," / ",i0," / ",i0," / ",i0)') &
      n_BT_gt_1mK, n_BT_gt_10mK, n_BT_gt_100mK, n_BT_gt_1K

    WRITE(*,'(/5x,"Radiance")')
    WRITE(*,'( 7x,"max |diff| = ",es12.4,"  at (chan=",i0,", prof=",i0,")")') &
      max_dRad, max_dRad_l, max_dRad_m
    WRITE(*,'( 7x,"mean|diff| = ",es12.4,"   rms|diff| = ",es12.4)') &
      sum_dRad/REAL(n_total,fp), SQRT(sumsq_dRad/REAL(n_total,fp))

    IF ( n_show > 0 .AND. max_dBT > 0.0_fp ) THEN
      CALL Print_Top_N_BT( actual, expected, rank_dBT, rank_l, rank_m, n_show )
    END IF
    WRITE(*,'(5x,a,/)') '================================================================'

    DEALLOCATE( rank_dBT, rank_l, rank_m )
  END SUBROUTINE Report_RTSolution_Diff


  SUBROUTINE Print_Top_N_BT( actual, expected, abs_dBT, l_idx, m_idx, n_show )
    TYPE(CRTM_RTSolution_type), INTENT(IN)    :: actual(:,:)
    TYPE(CRTM_RTSolution_type), INTENT(IN)    :: expected(:,:)
    REAL(fp),                   INTENT(INOUT) :: abs_dBT(:)
    INTEGER,                    INTENT(INOUT) :: l_idx(:), m_idx(:)
    INTEGER,                    INTENT(IN)    :: n_show

    INTEGER  :: i, k, kmax, n
    REAL(fp) :: vmax, vtmp
    INTEGER  :: ltmp, mtmp

    n = SIZE(abs_dBT)

    ! Selection sort the top n_show entries (n_show is small, n may be large).
    DO i = 1, n_show
      vmax = abs_dBT(i); kmax = i
      DO k = i+1, n
        IF ( abs_dBT(k) > vmax ) THEN
          vmax = abs_dBT(k); kmax = k
        END IF
      END DO
      IF ( kmax /= i ) THEN
        vtmp = abs_dBT(i); abs_dBT(i) = abs_dBT(kmax); abs_dBT(kmax) = vtmp
        ltmp = l_idx(i);   l_idx(i)   = l_idx(kmax);   l_idx(kmax)   = ltmp
        mtmp = m_idx(i);   m_idx(i)   = m_idx(kmax);   m_idx(kmax)   = mtmp
      END IF
    END DO

    WRITE(*,'(/5x,"Top ",i0," |dBT| offenders:")') n_show
    WRITE(*,'( 7x,"chan",2x,"prof",4x,"actual_BT",6x,"expected_BT",6x,"|dBT|")')
    DO i = 1, n_show
      IF ( abs_dBT(i) <= 0.0_fp ) EXIT
      WRITE(*,'(7x,i4,2x,i4,3(2x,es14.6))') &
        l_idx(i), m_idx(i), &
        actual  (l_idx(i), m_idx(i))%Brightness_Temperature, &
        expected(l_idx(i), m_idx(i))%Brightness_Temperature, &
        abs_dBT(i)
    END DO
  END SUBROUTINE Print_Top_N_BT

END MODULE CRTM_RTSolution_Diff
