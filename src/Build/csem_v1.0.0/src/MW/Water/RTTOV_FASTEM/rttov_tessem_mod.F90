! Description:
!> @file
!!   Subroutines for TESSEM2 MW sea surface emissivity model
!
!> @brief
!!   Subroutines for TESSEM2 MW sea surface emissivity model
!!
!! @details
!!   This contains the code which implements TESSEM2 for the direct, TL, and
!!   AD/K models.
!!
!!   TESSEM2 is a neural network-based emissivity model applicable between 10
!!   and 700GHz.
!!
!!   It is recommmended to use TESSEM2 for channels above 200GHz.
!!
!!   For frequencies below 200GHz TESSEM2 is based on FASTEM-6, but there is no
!!   azimuthal dependence.
!!
!!   Reference:
!!   Prigent, C., Aires, F., Wang, D., Fox, S. and Harlow, C. (2016)
!!   Sea surface emissivity parameterization from microwaves to millimeter
!!   waves. Q.J.R. Meteorol. Soc. Accepted Author Manuscript.
!!   doi:10.1002/qj.2953
!
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 25 November 1998, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, KNMI and MeteoFrance.
!
!    Copyright 2016, EUMETSAT, All Rights Reserved.
!
MODULE rttov_tessem_mod

  USE CSEM_Type_Kinds, ONLY: jprb => CSEM_fp, &
                             jpim => CSEM_ip
  ! Disable implicit typing
  IMPLICIT NONE

  INTEGER(jpim), PARAMETER :: tessem_nin    = 5   ! Number of input variables (freq, theta, windspeed, tskin, salinity)
  INTEGER(jpim), PARAMETER :: tessem_nout   = 1   ! Number of output variables (emissivity H-pol or V-pol)
  INTEGER(jpim), PARAMETER :: tessem_ncache = 15  ! Size of neural network

  TYPE tessem_net
    REAL(jprb) :: b1(tessem_ncache)
    REAL(jprb) :: b2(tessem_nout)
    REAL(jprb) :: w1(tessem_ncache,tessem_nin)
    REAL(jprb) :: w2(tessem_nout,tessem_ncache)
    REAL(jprb) :: x_min(tessem_nin)
    REAL(jprb) :: x_max(tessem_nin)
    REAL(jprb) :: y_min(tessem_nout)
    REAL(jprb) :: y_max(tessem_nout)
  ENDTYPE tessem_net

  TYPE(tessem_net) :: net_h, net_v

  DATA net_h%b1 / -0.492534_jprb,  3.549479_jprb,  1.209267_jprb, -1.201308_jprb,  0.657236_jprb, &
                   0.522424_jprb,  0.439041_jprb,  1.070375_jprb,  0.917859_jprb,  2.999629_jprb, &
                   0.150532_jprb, -0.225454_jprb, -1.418973_jprb, -5.616032_jprb,  3.077098_jprb /
  DATA net_h%b2 / -3.191476_jprb /
  DATA net_h%w1 /  0.914705_jprb, -0.244028_jprb,  0.217031_jprb,  0.022331_jprb, -0.757337_jprb, &
                  -0.514908_jprb, -0.682990_jprb,  1.319635_jprb,  1.052532_jprb,  0.327665_jprb, &
                  -0.197895_jprb,  0.529253_jprb, -1.846079_jprb, -5.075406_jprb,  3.362876_jprb, &
                  -0.907370_jprb, -1.694778_jprb, -0.628096_jprb,  2.230164_jprb,  0.003549_jprb, &
                  -0.290958_jprb, -0.013354_jprb,  0.122967_jprb,  0.001593_jprb, -0.061660_jprb, &
                   1.609950_jprb, -0.056635_jprb,  1.059906_jprb, -0.224954_jprb,  0.942301_jprb, &
                   1.806560_jprb,  0.821685_jprb, -1.057938_jprb,  0.082932_jprb,  0.045666_jprb, &
                   0.119414_jprb,  0.032902_jprb, -0.178454_jprb, -0.058426_jprb, -0.098223_jprb, &
                   0.291156_jprb, -0.036004_jprb, -0.524861_jprb,  0.089880_jprb,  0.268074_jprb, &
                  -0.237687_jprb,  0.040081_jprb, -0.019350_jprb,  0.028382_jprb,  0.831685_jprb, &
                  -1.168739_jprb,  0.715660_jprb,  0.498660_jprb, -0.716106_jprb,  0.916730_jprb, &
                   0.060666_jprb, -0.605959_jprb,  0.214584_jprb, -0.233417_jprb, -0.207224_jprb, &
                   0.013675_jprb, -0.000319_jprb,  0.008815_jprb, -0.004765_jprb,  0.823811_jprb, &
                  -0.244888_jprb,  0.265740_jprb,  1.596743_jprb, -2.614260_jprb, -2.705295_jprb, &
                  -0.003804_jprb,  0.105391_jprb, -0.008504_jprb,  0.003713_jprb, -0.018812_jprb /
  DATA net_h%w2 /  0.019943_jprb,  2.713515_jprb, -0.147582_jprb, -0.305444_jprb, -0.275328_jprb, &
                  -0.054163_jprb,  0.628523_jprb,  0.028561_jprb,  0.044241_jprb,  0.034885_jprb, &
                  -0.168618_jprb,  0.580301_jprb, -0.157557_jprb, -0.586725_jprb,  0.179877_jprb /
  DATA net_h%x_min / 10.000000_jprb, 0.000000_jprb, 0.000000_jprb, 273.149990_jprb, 0.000000_jprb /
  DATA net_h%x_max / 700.000000_jprb, 88.000000_jprb, 24.000000_jprb, 309.149994_jprb, 40.000000_jprb /
  DATA net_h%y_min / 0.018036_jprb /
  DATA net_h%y_max / 0.884545_jprb /

  DATA net_v%b1 / -4.425637_jprb, -9.833792_jprb,  1.338131_jprb, -1.688713_jprb,  0.121844_jprb, &
                  -0.020075_jprb,  0.983048_jprb, -5.365536_jprb,  0.759206_jprb,  1.705090_jprb, &
                  -0.861439_jprb,  2.639012_jprb,  4.438483_jprb, -1.376479_jprb,  0.920075_jprb /
  DATA net_v%b2 / -3.126669_jprb /
  DATA net_v%w1 /  0.988378_jprb, -0.053936_jprb, -1.166637_jprb,  0.213664_jprb, -0.264162_jprb, &
                  -0.946456_jprb, -0.084722_jprb, -2.169090_jprb,  1.038926_jprb,  1.814769_jprb, &
                  -0.588050_jprb,  1.553885_jprb,  3.672554_jprb,  0.198168_jprb,  0.520655_jprb, &
                  -0.865033_jprb,  6.091850_jprb,  0.514930_jprb,  1.203841_jprb,  1.349278_jprb, &
                  -0.875449_jprb, -1.582926_jprb,  1.683341_jprb,  1.679280_jprb,  0.533539_jprb, &
                   1.047984_jprb,  0.534832_jprb,  0.429764_jprb,  1.978477_jprb, -0.982343_jprb, &
                   0.024583_jprb, -2.832310_jprb, -0.249698_jprb, -0.264623_jprb, -0.385522_jprb, &
                   0.197084_jprb,  0.071463_jprb, -1.025112_jprb, -0.365062_jprb,  0.099878_jprb, &
                   0.308190_jprb, -0.034068_jprb,  0.011647_jprb,  0.158408_jprb, -0.512171_jprb, &
                  -1.243063_jprb, -0.008161_jprb,  0.323026_jprb, -0.040724_jprb,  0.526443_jprb, &
                  -1.044800_jprb,  0.035568_jprb,  0.057517_jprb, -0.559605_jprb, -0.086835_jprb, &
                   0.058220_jprb,  0.284328_jprb,  0.235683_jprb, -0.050183_jprb, -0.074034_jprb, &
                   2.514267_jprb,  0.002472_jprb, -0.249759_jprb, -0.017623_jprb,  2.265044_jprb, &
                   1.006211_jprb, -0.000011_jprb, -0.010390_jprb,  0.285759_jprb, -0.024043_jprb, &
                  -0.010873_jprb, -0.063544_jprb, -0.043210_jprb, -0.000275_jprb,  0.017342_jprb /
  DATA net_v%w2 / -0.105095_jprb, -5.132794_jprb, -0.068984_jprb, -1.133355_jprb, -0.016834_jprb, &
                   0.021879_jprb, -1.791876_jprb,  2.727342_jprb,  0.070811_jprb,  0.556869_jprb, &
                  -1.172738_jprb, -1.725109_jprb,  1.650515_jprb, -1.123763_jprb, -0.929365_jprb /
  DATA net_v%x_min / 10.000000_jprb, 0.000000_jprb, 0.000000_jprb, 273.149990_jprb, 0.000000_jprb /
  DATA net_v%x_max / 700.000000_jprb, 88.000000_jprb, 24.000000_jprb, 309.149994_jprb, 40.000000_jprb /
  DATA net_v%y_min / 0.272829_jprb /
  DATA net_v%y_max / 0.997013_jprb /

CONTAINS

  ! -------------
  ! Direct
  ! -------------

  SUBROUTINE prop_neuralnet(net, x, y)

    TYPE(tessem_net), INTENT(IN)    :: net
    REAL(jprb),       INTENT(IN)    :: x(tessem_nin)
    REAL(jprb),       INTENT(INOUT) :: y(tessem_nout)
    INTEGER(jpim) :: i, j
    REAL(jprb)    :: trans(tessem_ncache)
    REAL(jprb)    :: new_x(tessem_nin), new_y(tessem_nout)

    ! preprocessing
    DO i = 1, tessem_nin
      new_x(i) = -1._jprb + 2._jprb * (x(i) - net%x_min(i)) / (net%x_max(i) - net%x_min(i))
    ENDDO
    ! propagation
    DO i = 1, tessem_ncache
      trans(i) = net%b1(i)
      DO j = 1, tessem_nin
        trans(i) = trans(i) + net%w1(i,j) * new_x(j)
      ENDDO
      trans(i) = 2._jprb / (1._jprb + EXP(-2._jprb * trans(i))) - 1._jprb
    ENDDO
    DO i = 1, tessem_nout
      new_y(i) = net%b2(i)
      DO j = 1, tessem_ncache
        new_y(i) = new_y(i) + net%w2(i,j) * trans(j)
      ENDDO
    ENDDO
    ! postprocessing
    DO i = 1, tessem_nout
      y(i) = net%y_min(i) + (new_y(i) + 1._jprb) * 0.5_jprb * (net%y_max(i) - net%y_min(i))
    ENDDO

  END SUBROUTINE prop_neuralnet

  SUBROUTINE rttov_tessem(freq, theta, windspeed, tskin, salinity, emis_h, emis_v)

    REAL(jprb), INTENT(IN)  :: freq, theta, windspeed, tskin, salinity
    REAL(jprb), INTENT(OUT) :: emis_h, emis_v

    REAL(jprb) :: x(tessem_nin), y(tessem_nout)

    x(:) = (/ freq, theta, windspeed, tskin, salinity /)

    CALL prop_neuralnet(net_h, x, y)
    emis_h = y(1)

    CALL prop_neuralnet(net_v, x, y)
    emis_v = y(1)

  END SUBROUTINE rttov_tessem

  ! -------------
  ! TL
  ! -------------

  SUBROUTINE prop_neuralnet_tl(net, x, x_tl, y_tl)

    TYPE(tessem_net), INTENT(IN)    :: net
    REAL(jprb),       INTENT(IN)    :: x(tessem_nin)
    REAL(jprb),       INTENT(IN)    :: x_tl(tessem_nin)
    REAL(jprb),       INTENT(INOUT) :: y_tl(tessem_nout)
    INTEGER(jpim) :: i, j
    REAL(jprb)    :: trans(tessem_ncache), trans_tmp(tessem_ncache)
    REAL(jprb)    :: trans_tl(tessem_ncache), trans_tmp_tl(tessem_ncache)
    REAL(jprb)    :: new_x(tessem_nin)
    REAL(jprb)    :: new_x_tl(tessem_nin), new_y_tl(tessem_nout)

    ! preprocessing
    DO i = 1, tessem_nin
      new_x(i)    = -1._jprb + 2._jprb * (x(i) - net%x_min(i)) / (net%x_max(i) - net%x_min(i))
      new_x_tl(i) = 2._jprb * x_tl(i) / (net%x_max(i) - net%x_min(i))
    ENDDO
    ! propagation
    DO i = 1, tessem_ncache
      trans(i)    = net%b1(i)
      trans_tl(i) = 0._jprb
      DO j = 1, tessem_nin
        trans(i)    = trans(i) + net%w1(i,j) * new_x(j)
        trans_tl(i) = trans_tl(i) + net%w1(i,j) * new_x_tl(j)
      ENDDO
      trans_tmp(i)    = 2._jprb / (1._jprb + EXP(-2._jprb * trans(i)))
      trans_tmp_tl(i) = trans_tl(i) * EXP(-2._jprb * trans(i)) * trans_tmp(i)**2
    ENDDO
    DO i = 1, tessem_nout
      new_y_tl(i) = 0._jprb
      DO j = 1, tessem_ncache
        new_y_tl(i) = new_y_tl(i) + net%w2(i,j) * trans_tmp_tl(j)
      ENDDO
    ENDDO
    ! postprocessing
    DO i = 1, tessem_nout
      y_tl(i) = new_y_tl(i) * 0.5_jprb * (net%y_max(i) - net%y_min(i))
    ENDDO

  END SUBROUTINE prop_neuralnet_tl

  SUBROUTINE rttov_tessem_tl(freq, theta, windspeed, tskin, salinity, &
                             windspeed_tl, tskin_tl, salinity_tl, &
                             emis_h_tl, emis_v_tl)

    REAL(jprb), INTENT(IN)  :: freq, theta, windspeed, tskin, salinity
    REAL(jprb), INTENT(IN)  :: windspeed_tl, tskin_tl, salinity_tl
    REAL(jprb), INTENT(OUT) :: emis_h_tl, emis_v_tl

    REAL(jprb) :: x(tessem_nin), x_tl(tessem_nin), y_tl(tessem_nout)

    x(:)    = (/ freq, theta, windspeed, tskin, salinity /)
    x_tl(:) = (/ 0._jprb, 0._jprb, windspeed_tl, tskin_tl, salinity_tl /)

    CALL prop_neuralnet_tl(net_h, x, x_tl, y_tl)
    emis_h_tl = y_tl(1)

    CALL prop_neuralnet_tl(net_v, x, x_tl, y_tl)
    emis_v_tl = y_tl(1)

  END SUBROUTINE rttov_tessem_tl

  ! -------------
  ! AD
  ! -------------

  SUBROUTINE prop_neuralnet_ad(net, x, x_ad, y_ad)

    TYPE(tessem_net), INTENT(IN)    :: net
    REAL(jprb),       INTENT(IN)    :: x(tessem_nin)
    REAL(jprb),       INTENT(INOUT) :: x_ad(tessem_nin)
    REAL(jprb),       INTENT(IN)    :: y_ad(tessem_nout)
    INTEGER(jpim) :: i, j
    REAL(jprb)    :: trans(tessem_ncache), trans_tmp(tessem_ncache)
    REAL(jprb)    :: trans_ad(tessem_ncache), trans_tmp_ad(tessem_ncache)
    REAL(jprb)    :: new_x(tessem_nin)
    REAL(jprb)    :: new_x_ad(tessem_nin), new_y_ad(tessem_nout)

    ! direct

    ! preprocessing
    DO i = 1, tessem_nin
      new_x(i) = -1._jprb + 2._jprb * (x(i) - net%x_min(i)) / (net%x_max(i) - net%x_min(i))
    ENDDO
    ! propagation
    DO i = 1, tessem_ncache
      trans(i) = net%b1(i)
      DO j = 1, tessem_nin
        trans(i) = trans(i) + net%w1(i,j) * new_x(j)
      ENDDO
      trans_tmp(i) = 2._jprb / (1._jprb + EXP(-2._jprb * trans(i)))
    ENDDO

    ! ad

    new_y_ad = 0._jprb
    new_x_ad = 0._jprb
    trans_ad = 0._jprb
    trans_tmp_ad = 0._jprb

    ! postprocessing
    DO i = 1, tessem_nout
      new_y_ad(i) = new_y_ad(i) + y_ad(i) * 0.5_jprb * (net%y_max(i) - net%y_min(i))
    ENDDO
    ! propagation
    DO i = 1, tessem_nout
      DO j = 1, tessem_ncache
        trans_tmp_ad(j) = trans_tmp_ad(j) + net%w2(i,j) * new_y_ad(i)
      ENDDO
      new_y_ad(i) = 0._jprb
    ENDDO
    DO i = 1, tessem_ncache
      trans_ad(i) = trans_ad(i) + trans_tmp_ad(i) * EXP(-2._jprb * trans(i)) * trans_tmp(i)**2
      DO j = 1, tessem_nin
        new_x_ad(j) = new_x_ad(j) + net%w1(i,j) * trans_ad(i)
      ENDDO
      trans_ad(i) = 0._jprb
    ENDDO
    ! preprocessing
    DO i = 1, tessem_nin
      x_ad(i) = x_ad(i) + 2._jprb * new_x_ad(i) / (net%x_max(i) - net%x_min(i))
    ENDDO

  END SUBROUTINE prop_neuralnet_ad

  SUBROUTINE rttov_tessem_ad(freq, theta, windspeed, tskin, salinity, &
                             windspeed_ad, tskin_ad, salinity_ad, &
                             emis_h_ad, emis_v_ad)

    REAL(jprb), INTENT(IN)    :: freq, theta, windspeed, tskin, salinity
    REAL(jprb), INTENT(INOUT) :: windspeed_ad, tskin_ad, salinity_ad
    REAL(jprb), INTENT(IN)    :: emis_h_ad, emis_v_ad

    REAL(jprb) :: x(tessem_nin), x_ad(tessem_nin), y_ad(tessem_nout)

    x(:)    = (/ freq, theta, windspeed, tskin, salinity /)

    x_ad(:) = 0._jprb
    y_ad(1) = emis_v_ad
    CALL prop_neuralnet_ad(net_v, x, x_ad, y_ad)

    y_ad(1) = emis_h_ad
    CALL prop_neuralnet_ad(net_h, x, x_ad, y_ad)

    windspeed_ad = windspeed_ad + x_ad(3)
    tskin_ad     = tskin_ad     + x_ad(4)
    salinity_ad  = salinity_ad  + x_ad(5)

  END SUBROUTINE rttov_tessem_ad

END MODULE rttov_tessem_mod
