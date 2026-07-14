!
! TELSEM2_Atlas_Module
!
! Module implementing the TELSEM2 microwave land surface emissivity interpolator.
!
! This is a port of the stand-alone RTTOV TELSEM2 module (mod_mwatlas_m2.F90,
! Copyright 2016 EUMETSAT/NWP SAF) to operate on the CRTM TELSEM2Atlas_type. The
! grid geometry (equal-area cell numbering), the multi-linear angular-correction
! regression and the inter-frequency interpolation reproduce the reference code
! so that results match RTTOV. Only the emissivity is computed here; the atlas
! uncertainty (std/covariance) returned by the reference optional arguments is
! not needed by the CRTM surface optics and is omitted.
!
! Reference:
!   D. Wang, C. Prigent, L. Kilic, S. Fox, R. C. Harlow, C. Jimenez, F. Aires,
!   C. Grassotti, F. Karbou, 2017. Surface emissivity at microwaves to
!   millimeter waves over polar regions: parameterization and evaluation with
!   aircraft experiments. J. Atmos. Oceanic Technol.
!

MODULE TELSEM2_Atlas_Module

  ! -----------------
  ! Environment setup
  ! -----------------
  USE Type_Kinds         , ONLY: fp, Long, Double
  USE TELSEM2Atlas_Define, ONLY: TELSEM2Atlas_type, TELSEM2Atlas_Associated
  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  PUBLIC :: TELSEM2_Setup_Grid
  PUBLIC :: TELSEM2_Emissivity
  ! Lower-level routines exposed for unit testing
  PUBLIC :: TELSEM2_CalcCellnum
  PUBLIC :: TELSEM2_GetCoordinates


  ! -----------------
  ! Module parameters
  ! -----------------
  REAL(fp), PARAMETER :: ZERO = 0.0_fp
  REAL(fp), PARAMETER :: ONE  = 1.0_fp
  ! Earth radius used by the equal-area grid generator (km)
  REAL(Double), PARAMETER :: REARTH = 6371.2_Double
  ! Maximum number of cells gathered in the spatial-averaging box
  INTEGER, PARAMETER :: MAX_BOX_CELLS = 400

  ! Frequency-interpolation ratios for water-like classes above 85 GHz
  ! (class2 = 10..13). Indexed by class2-9.
  REAL(fp), PARAMETER :: RAPPORT43_32(4) = [ 0.62_fp, 0.37_fp, 0.46_fp, 0.63_fp ]
  REAL(fp), PARAMETER :: RAPPORT54_43(4) = [ 0.30_fp, 0.60_fp, 0.47_fp, 0.35_fp ]

  ! Angular-correction regression coefficients (3 anchor frequencies x 10 classes).
  ! Transcribed verbatim from the RTTOV TELSEM2 reference (column-major order).
  REAL(fp), PARAMETER :: A0_K0(3,10) = RESHAPE( [ &
     0.11509_fp,0.091535_fp,0.34796_fp,0.10525_fp,0.16627_fp,0.24434_fp, &
     0.29217_fp,0.23809_fp,0.28954_fp,0.17516_fp,0.19459_fp,0.28697_fp, &
     0.10521_fp,0.12126_fp,0.30278_fp,0.18212_fp,0.19625_fp,0.14551_fp, &
     -0.19202_fp,0.5411_fp,0.03739_fp,0.10292_fp,0.5486_fp,-0.058937_fp, &
     -0.022672_fp,0.44492_fp,-0.058448_fp,-0.33894_fp,-0.17621_fp,0.14742_fp ], [3,10] )
  REAL(fp), PARAMETER :: A0_K1(3,10) = RESHAPE( [ &
     0.61168_fp,0.59095_fp,0.7918_fp,0.60271_fp,0.69213_fp,0.62218_fp, &
     0.32728_fp,0.34334_fp,0.37062_fp,0.51217_fp,0.4491_fp,0.50101_fp, &
     0.48913_fp,0.41932_fp,0.29734_fp,0.64474_fp,0.30637_fp,0.031107_fp, &
     1.0405_fp,0.17538_fp,1.3215_fp,0.61819_fp,0.31298_fp,1.7218_fp, &
     0.87761_fp,0.47583_fp,1.2583_fp,1.0959_fp,0.92842_fp,0.51033_fp ], [3,10] )
  REAL(fp), PARAMETER :: A0_K2(3,10) = RESHAPE( [ &
     0.26726_fp,0.32033_fp,-0.14778_fp,0.28547_fp,0.13592_fp,0.13193_fp, &
     0.37178_fp,0.41813_fp,0.33875_fp,0.30203_fp,0.35479_fp,0.20189_fp, &
     0.40663_fp,0.47493_fp,0.40668_fp,0.14811_fp,0.52382_fp,0.86634_fp, &
     0.14286_fp,0.27164_fp,-0.37947_fp,0.2737_fp,0.12001_fp,-0.67315_fp, &
     0.13492_fp,0.065463_fp,-0.19316_fp,0.24905_fp,0.25475_fp,0.34637_fp ], [3,10] )
  REAL(fp), PARAMETER :: A0_EVEH(3,10) = RESHAPE( [ &
     0.9592599869E+00_fp,0.9565299749E+00_fp,0.9511899948E+00_fp, &
     0.9560700059E+00_fp,0.9541199803E+00_fp,0.9483199716E+00_fp, &
     0.9461100101E+00_fp,0.9439799786E+00_fp,0.9387800097E+00_fp, &
     0.9317600131E+00_fp,0.9289000034E+00_fp,0.9236800075E+00_fp, &
     0.9208700061E+00_fp,0.9190599918E+00_fp,0.9105200171E+00_fp, &
     0.9162799716E+00_fp,0.8937299848E+00_fp,0.8014699817E+00_fp, &
     0.9570500255E+00_fp,0.9213600159E+00_fp,0.7893999815E+00_fp, &
     0.9639400244E+00_fp,0.9530599713E+00_fp,0.8850200176E+00_fp, &
     0.9685299993E+00_fp,0.9622600079E+00_fp,0.9118800163E+00_fp, &
     0.8997200131E+00_fp,0.9012699723E+00_fp,0.9107499719E+00_fp ], [3,10] )
  REAL(fp), PARAMETER :: A1_EVEH(3,10) = RESHAPE( [ &
     0.3627802414E-07_fp,-0.7778328204E-08_fp,0.4396108011E-07_fp, &
     0.2503205394E-06_fp,0.1996262995E-06_fp,0.2929977541E-06_fp, &
     0.4190530660E-06_fp,0.3655744649E-06_fp,0.3519195673E-06_fp, &
     0.5574374313E-06_fp,0.5273076340E-06_fp,0.5376484182E-06_fp, &
     0.1026844529E-05_fp,0.9679998811E-06_fp,0.8616486866E-06_fp, &
     0.3180800832E-06_fp,0.2886778532E-06_fp,0.2310362675E-06_fp, &
     -0.1118036366E-06_fp,-0.1502856577E-06_fp,0.4842232926E-07_fp, &
     -0.8410978580E-08_fp,-0.3478669441E-07_fp,0.2209441590E-06_fp, &
     0.2485776633E-06_fp,0.1800235907E-06_fp,0.2510202251E-06_fp, &
     0.2687000915E-06_fp,0.1740325644E-06_fp,0.3562134339E-06_fp ], [3,10] )
  REAL(fp), PARAMETER :: A2_EVEH(3,10) = RESHAPE( [ &
     0.3067140824E-05_fp,0.2520012231E-05_fp,0.4831396382E-05_fp, &
     0.8213598448E-05_fp,0.7378375358E-05_fp,0.1022081960E-04_fp, &
     0.1225889173E-04_fp,0.1165553113E-04_fp,0.1188659007E-04_fp, &
     0.1693615741E-04_fp,0.1648317448E-04_fp,0.1715818144E-04_fp, &
     0.2744720041E-04_fp,0.2642072104E-04_fp,0.2671847506E-04_fp, &
     0.1349592094E-04_fp,0.1261523357E-04_fp,0.5447756394E-05_fp, &
     0.2064244654E-05_fp,0.1919016057E-06_fp,0.5940860319E-06_fp, &
     0.5334760772E-05_fp,0.4130339221E-05_fp,0.4104662821E-05_fp, &
     0.6530796327E-05_fp,0.5727014013E-05_fp,0.7451782039E-05_fp, &
     0.1071246970E-04_fp,0.9539280654E-05_fp,0.1034286015E-04_fp ], [3,10] )
  REAL(fp), PARAMETER :: A3_EVEH(3,10) = RESHAPE( [ &
     -0.2004991551E-07_fp,-0.6895366056E-07_fp,-0.2047409282E-06_fp, &
     -0.7322448425E-07_fp,-0.1273002681E-06_fp,-0.2729916844E-06_fp, &
     -0.9421125213E-07_fp,-0.1683332300E-06_fp,-0.2726891637E-06_fp, &
     -0.1317753799E-06_fp,-0.2107972250E-06_fp,-0.3556060904E-06_fp, &
     -0.1889465580E-06_fp,-0.2757958271E-06_fp,-0.4909850304E-06_fp, &
     0.7339644004E-08_fp,-0.4058669560E-06_fp,-0.4146343997E-06_fp, &
     0.6170279931E-07_fp,-0.1998567996E-06_fp,-0.4713119139E-07_fp, &
     -0.1361754887E-07_fp,-0.1765622955E-06_fp,-0.2348146637E-06_fp, &
     -0.3901189061E-07_fp,-0.1305666189E-06_fp,-0.1533838798E-06_fp, &
     -0.2679148992E-07_fp,-0.4441960044E-07_fp,-0.1815613899E-06_fp ], [3,10] )
  REAL(fp), PARAMETER :: B0_EVEH(3,10) = RESHAPE( [ &
     0.9592599869E+00_fp,0.9565299749E+00_fp,0.9511899948E+00_fp, &
     0.9560700059E+00_fp,0.9541199803E+00_fp,0.9483199716E+00_fp, &
     0.9461100101E+00_fp,0.9439799786E+00_fp,0.9387800097E+00_fp, &
     0.9317600131E+00_fp,0.9289000034E+00_fp,0.9236800075E+00_fp, &
     0.9208700061E+00_fp,0.9190599918E+00_fp,0.9105200171E+00_fp, &
     0.9162799716E+00_fp,0.8937299848E+00_fp,0.8014699817E+00_fp, &
     0.9570500255E+00_fp,0.9213600159E+00_fp,0.7893999815E+00_fp, &
     0.9639400244E+00_fp,0.9530599713E+00_fp,0.8850200176E+00_fp, &
     0.9685299993E+00_fp,0.9622600079E+00_fp,0.9118800163E+00_fp, &
     0.8997200131E+00_fp,0.9012699723E+00_fp,0.9107499719E+00_fp ], [3,10] )
  REAL(fp), PARAMETER :: B1_EVEH(3,10) = RESHAPE( [ &
     0.3626608347E-07_fp,-0.7786279177E-08_fp,0.4393379172E-07_fp, &
     0.2502746099E-06_fp,0.1995944388E-06_fp,0.2929554341E-06_fp, &
     0.4189516289E-06_fp,0.3655020180E-06_fp,0.3518483140E-06_fp, &
     0.5572838404E-06_fp,0.5271903092E-06_fp,0.5375342766E-06_fp, &
     0.1026605219E-05_fp,0.9677979733E-06_fp,0.8614680951E-06_fp, &
     0.3179358714E-06_fp,0.2884899004E-06_fp,0.2308632219E-06_fp, &
     -0.1118781370E-06_fp,-0.1503948681E-06_fp,0.4834672396E-07_fp, &
     -0.8455684153E-08_fp,-0.3485171618E-07_fp,0.2208606134E-06_fp, &
     0.2485595019E-06_fp,0.1799959364E-06_fp,0.2509846695E-06_fp, &
     0.2686167306E-06_fp,0.1739760478E-06_fp,0.3561317214E-06_fp ], [3,10] )
  REAL(fp), PARAMETER :: B2_EVEH(3,10) = RESHAPE( [ &
     0.3065537157E-05_fp,0.2518960400E-05_fp,0.4829731552E-05_fp, &
     0.8209894986E-05_fp,0.7375769655E-05_fp,0.1021809931E-04_fp, &
     0.1225203869E-04_fp,0.1165053800E-04_fp,0.1188218721E-04_fp, &
     0.1692612022E-04_fp,0.1647546378E-04_fp,0.1715117833E-04_fp, &
     0.2743142431E-04_fp,0.2640772436E-04_fp,0.2670711910E-04_fp, &
     0.1348545720E-04_fp,0.1260529825E-04_fp,0.5439695997E-05_fp, &
     0.2058213340E-05_fp,0.1860650656E-06_fp,0.5898303925E-06_fp, &
     0.5330772183E-05_fp,0.4126528893E-05_fp,0.4100859314E-05_fp, &
     0.6528573977E-05_fp,0.5725009032E-05_fp,0.7449450095E-05_fp, &
     0.1070590315E-04_fp,0.9534271157E-05_fp,0.1033751869E-04_fp ], [3,10] )
  REAL(fp), PARAMETER :: B3_EVEH(3,10) = RESHAPE( [ &
     -0.1370247134E-06_fp,-0.1436897747E-06_fp,-0.2954870411E-06_fp, &
     -0.3118435643E-06_fp,-0.2916583242E-06_fp,-0.4311032171E-06_fp, &
     -0.5048401022E-06_fp,-0.4662823869E-06_fp,-0.5206445053E-06_fp, &
     -0.7210980471E-06_fp,-0.6662896794E-06_fp,-0.7548637200E-06_fp, &
     -0.1110204039E-05_fp,-0.1030801400E-05_fp,-0.1140921199E-05_fp, &
     -0.6330818110E-06_fp,-0.9186441048E-06_fp,-0.7947813856E-06_fp, &
     -0.3242539890E-06_fp,-0.5027602583E-06_fp,-0.2777987334E-06_fp, &
     -0.2747250676E-06_fp,-0.3811997260E-06_fp,-0.4102405455E-06_fp, &
     -0.1994112324E-06_fp,-0.2555484855E-06_fp,-0.2842682534E-06_fp, &
     -0.4413041665E-06_fp,-0.3717419474E-06_fp,-0.4975536854E-06_fp ], [3,10] )


CONTAINS


!--------------------------------------------------------------------------------
! TELSEM2_Setup_Grid
!
! Construct the equal-area grid geometry (Cells_Per_Band, First_Cell) and the
! per-month reverse lookup (Correspondence) from the loaded sparse atlas data.
! Must be called once after the atlas data arrays have been populated.
!--------------------------------------------------------------------------------
  SUBROUTINE TELSEM2_Setup_Grid( atlas )
    TYPE(TELSEM2Atlas_type), INTENT(IN OUT) :: atlas
    INTEGER :: m, k, p

    ! Equal-area grid geometry
    CALL equare( atlas )

    ! Reverse lookup: cell number -> stacked data index, per month (0 = empty)
    atlas%Correspondence = 0
    DO m = 1, atlas%n_Months
      DO k = 1, atlas%Month_Data_Count(m)
        p = atlas%Month_Offset(m) + k
        atlas%Correspondence( atlas%Cell_Number(p), m ) = p
      END DO
    END DO
  END SUBROUTINE TELSEM2_Setup_Grid


!--------------------------------------------------------------------------------
! TELSEM2_Emissivity
!
! Return the TELSEM2 V- and H-pol emissivity for a location, month, frequency and
! zenith angle. Valid is .FALSE. when the atlas has no land climatology at the
! requested cell (e.g. open water / permanent ice), in which case the caller
! should fall back to another model.
!--------------------------------------------------------------------------------
  SUBROUTINE TELSEM2_Emissivity( &
    atlas        , &  ! Input,  loaded atlas
    Latitude     , &  ! Input,  degrees, -90..90
    Longitude    , &  ! Input,  degrees (any range; reduced modulo 360)
    Month        , &  ! Input,  1..12
    Frequency    , &  ! Input,  GHz
    Zenith_Angle , &  ! Input,  degrees
    Emissivity_V , &  ! Output, V-pol emissivity
    Emissivity_H , &  ! Output, H-pol emissivity
    Valid        , &  ! Output, .TRUE. if atlas data found
    Resolution     )  ! Optional input, spatial-averaging resolution (deg)
    ! Arguments
    TYPE(TELSEM2Atlas_type), INTENT(IN)  :: atlas
    REAL(fp),                INTENT(IN)  :: Latitude
    REAL(fp),                INTENT(IN)  :: Longitude
    INTEGER,                 INTENT(IN)  :: Month
    REAL(fp),                INTENT(IN)  :: Frequency
    REAL(fp),                INTENT(IN)  :: Zenith_Angle
    REAL(fp),                INTENT(OUT) :: Emissivity_V
    REAL(fp),                INTENT(OUT) :: Emissivity_H
    LOGICAL,                 INTENT(OUT) :: Valid
    REAL(fp),      OPTIONAL, INTENT(IN)  :: Resolution
    ! Local variables
    REAL(fp) :: lon, resol
    REAL(fp) :: ev_a(3), eh_a(3), ev, eh, ev_sum, eh_sum
    INTEGER  :: cell(MAX_BOX_CELLS), nb_cell, ii, ipos, inumb

    Emissivity_V = ZERO
    Emissivity_H = ZERO
    Valid        = .FALSE.

    ! Guard against an unloaded atlas or out-of-range month
    IF ( .NOT. TELSEM2Atlas_Associated(atlas) ) RETURN
    IF ( Month < 1 .OR. Month > atlas%n_Months ) RETURN

    resol = atlas%Resolution
    IF ( PRESENT(Resolution) ) resol = Resolution
    lon = MODULO( Longitude, 360.0_fp )

    ! List of atlas cells contributing at the requested resolution (1 at native)
    CALL calc_cellnum_mult( atlas, Latitude, lon, resol, cell, nb_cell )

    ev_sum = ZERO
    eh_sum = ZERO
    inumb  = 0
    DO ii = 1, nb_cell
      ipos = atlas%Correspondence( cell(ii), Month )
      IF ( ipos > 0 ) THEN
        inumb = inumb + 1
        ev_a(1) = REAL(atlas%Emissivity(ipos,1), fp)
        eh_a(1) = REAL(atlas%Emissivity(ipos,2), fp)
        ev_a(2) = REAL(atlas%Emissivity(ipos,4), fp)
        eh_a(2) = REAL(atlas%Emissivity(ipos,5), fp)
        ev_a(3) = REAL(atlas%Emissivity(ipos,6), fp)
        eh_a(3) = REAL(atlas%Emissivity(ipos,7), fp)
        CALL emis_interp( Zenith_Angle, Frequency, &
                          atlas%Class1(ipos), atlas%Class2(ipos), &
                          ev_a, eh_a, ev, eh )
        ev_sum = ev_sum + ev
        eh_sum = eh_sum + eh
      END IF
    END DO

    IF ( inumb > 0 ) THEN
      Emissivity_V = ev_sum / REAL(inumb, fp)
      Emissivity_H = eh_sum / REAL(inumb, fp)
      Valid        = .TRUE.
    END IF
  END SUBROUTINE TELSEM2_Emissivity


!--------------------------------------------------------------------------------
! TELSEM2_CalcCellnum (public wrapper for testing)
!--------------------------------------------------------------------------------
  FUNCTION TELSEM2_CalcCellnum( atlas, lat, lon ) RESULT( cellnum )
    TYPE(TELSEM2Atlas_type), INTENT(IN) :: atlas
    REAL(fp),                INTENT(IN) :: lat, lon
    INTEGER :: cellnum
    cellnum = calc_cellnum( atlas, lat, lon )
  END FUNCTION TELSEM2_CalcCellnum


!--------------------------------------------------------------------------------
! TELSEM2_GetCoordinates (public wrapper for testing)
!--------------------------------------------------------------------------------
  SUBROUTINE TELSEM2_GetCoordinates( atlas, cellnum, lat, lon )
    TYPE(TELSEM2Atlas_type), INTENT(IN)  :: atlas
    INTEGER,                 INTENT(IN)  :: cellnum
    REAL(fp),                INTENT(OUT) :: lat, lon
    CALL get_coordinates( atlas, cellnum, lat, lon )
  END SUBROUTINE TELSEM2_GetCoordinates


!################################################################################
!##                          ## PRIVATE PROCEDURES ##                          ##
!################################################################################

  ! Equal-area grid: number of cells per latitude band and the first cell number
  ! of each band. Faithful port of the reference equare routine.
  SUBROUTINE equare( atlas )
    TYPE(TELSEM2Atlas_type), INTENT(IN OUT) :: atlas
    INTEGER :: maxlat, maxlon, maxlt2, lat, icellr, lat1, lat2, numcel, numcls, lon, i
    REAL(Double) :: dlat, rcells, pi, aezon, hezon, aecell, xlatb
    REAL(Double) :: rlatb, rlate, xlate, htb, hte, htzone, rcelat, azone
    INTEGER, ALLOCATABLE :: tocell(:,:)

    dlat   = atlas%Resolution
    maxlat = FLOOR(180._Double/dlat)
    maxlon = FLOOR(360._Double/dlat)
    ALLOCATE( tocell(maxlon,maxlat) )

    pi     = 2.0_Double * ASIN(1.0_Double)
    rcelat = (dlat*pi)/180.0_Double
    hezon  = REARTH*SIN(rcelat)
    aezon  = 2.0_Double*pi*REARTH*hezon
    aecell = (aezon*dlat)/360.0_Double
    maxlt2 = maxlat/2
    DO lat = 1, maxlt2
      xlatb = (lat-1)*dlat
      xlate = xlatb+dlat
      rlatb = (2.0_Double*pi*xlatb)/360.0_Double
      rlate = (2.0_Double*pi*xlate)/360.0_Double
      htb    = REARTH*SIN(rlatb)
      hte    = REARTH*SIN(rlate)
      htzone = hte-htb
      azone  = 2.0_Double*pi*REARTH*htzone
      rcells = azone/aecell
      icellr = FLOOR(rcells+0.50_Double)
      lat1 = lat+maxlt2
      lat2 = maxlt2+1-lat
      atlas%Cells_Per_Band(lat1) = icellr
      atlas%Cells_Per_Band(lat2) = icellr
    END DO
    numcel = 0
    DO lat = 1, maxlat
      numcls = atlas%Cells_Per_Band(lat)
      DO lon = 1, numcls
        numcel = numcel + 1
        tocell(lon,lat) = numcel
      END DO
    END DO
    DEALLOCATE( tocell )

    atlas%First_Cell(1) = 1
    DO i = 2, maxlat
      atlas%First_Cell(i) = atlas%First_Cell(i-1) + atlas%Cells_Per_Band(i-1)
    END DO
  END SUBROUTINE equare


  ! Cell number for a given lat/lon (lon already in [0,360)).
  FUNCTION calc_cellnum( atlas, lat, lon ) RESULT( cellnum )
    TYPE(TELSEM2Atlas_type), INTENT(IN) :: atlas
    REAL(fp),                INTENT(IN) :: lat, lon
    INTEGER :: cellnum
    INTEGER :: ilat, ilon
    ilat = MIN( INT((lat+90._fp)/atlas%Resolution)+1, INT(atlas%n_Latitude_Bands) )
    ilon = INT( lon/(360._fp/atlas%Cells_Per_Band(ilat)) ) + 1
    cellnum = atlas%First_Cell(ilat) + ilon - 1
  END FUNCTION calc_cellnum


  ! Cell-center lat/lon for a given cell number (inverse of calc_cellnum).
  SUBROUTINE get_coordinates( atlas, cellnum, lat, lon )
    TYPE(TELSEM2Atlas_type), INTENT(IN)  :: atlas
    INTEGER,                 INTENT(IN)  :: cellnum
    REAL(fp),                INTENT(OUT) :: lat, lon
    REAL(fp) :: res_lat
    INTEGER  :: i, index_lat_max, index_lat, index_lon

    res_lat       = atlas%Resolution
    index_lat_max = INT(180._fp/res_lat)
    lat = ZERO
    lon = ZERO
    IF ( cellnum >= atlas%First_Cell(index_lat_max) ) THEN
      index_lat = index_lat_max
      lat = (index_lat - 0.5_fp)*res_lat - 90._fp
      index_lon = cellnum - atlas%First_Cell(index_lat_max) + 1
      lon = (index_lon - 0.5_fp)*(360._fp/atlas%Cells_Per_Band(index_lat))
    ELSE
      DO i = 1, index_lat_max-1
        IF ( (cellnum >= atlas%First_Cell(i)) .AND. (cellnum < atlas%First_Cell(i+1)) ) THEN
          index_lat = i
          lat = (index_lat - 0.5_fp)*res_lat - 90._fp
          index_lon = cellnum - atlas%First_Cell(i) + 1
          lon = (index_lon - 0.5_fp)*(360._fp/atlas%Cells_Per_Band(index_lat))
        END IF
      END DO
    END IF
  END SUBROUTINE get_coordinates


  ! Cell numbers within a box of size 'resol' centred on (lat,lon).
  SUBROUTINE calc_cellnum_mult( atlas, lat, lon, resol, cell, nb_cell )
    TYPE(TELSEM2Atlas_type), INTENT(IN)  :: atlas
    REAL(fp),                INTENT(IN)  :: lat, lon, resol
    INTEGER,                 INTENT(OUT) :: cell(:)
    INTEGER,                 INTENT(OUT) :: nb_cell
    INTEGER  :: ilat, ilon, nbcel, i2lon, i3lon, i2lat, nbreslat, nbreslon, i4lon
    REAL(fp) :: maxlat

    ilat = MIN( INT((lat+90._fp)/atlas%Resolution)+1, INT(atlas%n_Latitude_Bands) )
    ilon = INT( lon/(360._fp/atlas%Cells_Per_Band(ilat)) ) + 1

    nbcel = 1
    cell(1) = atlas%First_Cell(ilat) + ilon - 1
    maxlat  = 180._fp/atlas%Resolution
    nbreslat = INT( resol/atlas%Resolution/2._fp )
    IF ( nbreslat >= 1 ) THEN
      DO i2lat = ilat-nbreslat, ilat+nbreslat
        IF ( (i2lat < 1) .OR. (i2lat > INT(maxlat)) ) CYCLE
        IF ( ABS((i2lat-0.5_fp)*atlas%Resolution-90._fp-lat) <= resol/2._fp ) THEN
          i2lon = INT( lon/(360._fp/atlas%Cells_Per_Band(i2lat)) ) + 1
          nbreslon = INT( resol/(360._fp/(REAL(atlas%Cells_Per_Band(i2lat),fp)))/2._fp )
          DO i3lon = i2lon-nbreslon, i2lon+nbreslon
            IF ( MOD(ABS((i3lon-0.5_fp)*(360._fp/atlas%Cells_Per_Band(i2lat))-lon),360._fp) <= resol/2._fp ) THEN
              IF ( nbcel >= SIZE(cell) ) CYCLE
              ! Longitude wrap: valid cell indices are 1..Cells_Per_Band, so only
              ! indices strictly beyond the band wrap (i3lon == Cells_Per_Band is
              ! the band's valid last cell, not a wrap case).
              i4lon = i3lon
              IF ( i3lon < 1 )                           i4lon = atlas%Cells_Per_Band(i2lat) + i3lon
              IF ( i3lon > atlas%Cells_Per_Band(i2lat) ) i4lon = i3lon - atlas%Cells_Per_Band(i2lat)
              nbcel = nbcel + 1
              cell(nbcel) = atlas%First_Cell(i2lat) + i4lon - 1
              IF ( cell(nbcel) == cell(1) ) nbcel = nbcel - 1
            END IF
          END DO
        END IF
      END DO
    END IF
    nb_cell = nbcel
  END SUBROUTINE calc_cellnum_mult


  ! Inter-frequency interpolation of emissivity (piecewise linear 19/37/85 GHz,
  ! with frequency dependence above 85 GHz for water-like classes 10..13).
  SUBROUTINE interp_freq2( emiss19, emiss37, emiss85, f, class2, emiss, an, bn, cn )
    REAL(fp),           INTENT(IN)  :: emiss19, emiss37, emiss85, f
    INTEGER,            INTENT(IN)  :: class2
    REAL(fp),           INTENT(OUT) :: emiss
    REAL(fp), OPTIONAL, INTENT(OUT) :: an, bn, cn
    REAL(fp) :: a, b, c

    IF ( f <= 19.35_fp ) THEN
      a = ONE; b = ZERO; c = ZERO
      emiss = emiss19
    ELSE IF ( f <= 37._fp ) THEN
      a = (37._fp-f)/(37._fp-19.35_fp)
      b = (f-19.35_fp)/(37._fp-19.35_fp)
      c = ZERO
      emiss = a*emiss19 + b*emiss37
    ELSE IF ( f < 85.5_fp ) THEN
      a = ZERO
      b = (85.5_fp-f)/(85.5_fp-37._fp)
      c = (f-37._fp)/(85.5_fp-37._fp)
      emiss = b*emiss37 + c*emiss85
    ELSE
      a = ZERO; b = ZERO; c = ONE
      emiss = emiss85
      IF ( (class2 > 9) .AND. (class2 < 14) .AND. (emiss85 > emiss37) ) THEN
        IF ( f <= 150._fp ) THEN
          emiss = emiss85 + (f-85.5_fp)*((emiss85-emiss37)/(85.5_fp-37._fp))*RAPPORT43_32(class2-9)
        ELSE IF ( f <= 190._fp ) THEN
          emiss = emiss85 + (150._fp-85.5_fp)*((emiss85-emiss37)/(85.5_fp-37._fp))*RAPPORT43_32(class2-9)
          emiss = emiss + (f-150._fp)*((emiss-emiss85)/(150._fp-85.5_fp))*RAPPORT54_43(class2-9)
        ELSE
          emiss = emiss85 + (150._fp-85.5_fp)*((emiss85-emiss37)/(85.5_fp-37._fp))*RAPPORT43_32(class2-9)
          emiss = emiss + (190._fp-150._fp)*((emiss-emiss85)/(150._fp-85.5_fp))*RAPPORT54_43(class2-9)
        END IF
        IF ( emiss > ONE ) emiss = ONE
      END IF
    END IF

    IF ( PRESENT(an) ) an = a
    IF ( PRESENT(bn) ) bn = b
    IF ( PRESENT(cn) ) cn = c
  END SUBROUTINE interp_freq2


  ! Angular + frequency interpolation of the 3 anchor emissivities to (theta,freq).
  SUBROUTINE emis_interp( theta, freq, class1, class2, ev, eh, emiss_interp_v, emiss_interp_h )
    REAL(fp), INTENT(IN)  :: theta, freq, ev(3), eh(3)
    INTEGER,  INTENT(IN)  :: class1, class2
    REAL(fp), INTENT(OUT) :: emiss_interp_v, emiss_interp_h
    REAL(fp) :: e0, theta0, theta53, emiss_scal_v(3), emiss_scal_h(3)
    REAL(fp) :: S1_v, S1_h, S2_v, S2_h, S_v, S_h, a0, a1, a2, a3, b0, b1, b2, b3
    REAL(fp) :: em53_v, em53_h, emtheta_v, emtheta_h
    INTEGER  :: j

    DO j = 1, 3
      ! Nadir value from multi-linear regression on the V/H anchor emissivities
      e0 = A0_K0(j,class1) + A0_K1(j,class1)*ev(j) + A0_K2(j,class1)*eh(j)
      a0 = A0_EVEH(j,class1); a1 = A1_EVEH(j,class1)
      a2 = A2_EVEH(j,class1); a3 = A3_EVEH(j,class1)
      b0 = B0_EVEH(j,class1); b1 = B1_EVEH(j,class1)
      b2 = B2_EVEH(j,class1); b3 = B3_EVEH(j,class1)
      theta0  = ZERO
      theta53 = 53._fp
      ! V polarization
      S1_v = ((theta-theta53)/(theta0-theta53)) * ((e0-a0)/a0)
      em53_v = a3*(theta53**3) + a2*(theta53**2) + a1*theta53 + a0
      S2_v = ((theta-theta0)/(theta53-theta0))*((ev(j)-em53_v)/em53_v)
      S_v = ONE + S1_v + S2_v
      emtheta_v = a3*(theta**3) + a2*(theta**2) + a1*theta + a0
      emiss_scal_v(j) = S_v * emtheta_v
      ! H polarization
      S1_h = ((theta-theta53)/(theta0-theta53)) * ((e0-b0)/b0)
      em53_h = b3*(theta53**3) + b2*(theta53**2) + b1*theta53 + b0
      S2_h = ((theta-theta0)/(theta53-theta0))*((eh(j)-em53_h)/em53_h)
      S_h = ONE + S1_h + S2_h
      emtheta_h = b3*(theta**3) + b2*(theta**2) + b1*theta + b0
      emiss_scal_h(j) = S_h * emtheta_h
    END DO

    CALL interp_freq2( emiss_scal_v(1), emiss_scal_v(2), emiss_scal_v(3), freq, class2, emiss_interp_v )
    CALL interp_freq2( emiss_scal_h(1), emiss_scal_h(2), emiss_scal_h(3), freq, class2, emiss_interp_h )

    ! Where V < H, average the two (reference behaviour)
    IF ( emiss_interp_v < emiss_interp_h ) THEN
      emiss_interp_v = (emiss_interp_v + emiss_interp_h)/2._fp
      emiss_interp_h = emiss_interp_v
    END IF
    emiss_interp_h = MIN( ONE, emiss_interp_h )
    emiss_interp_v = MIN( ONE, emiss_interp_v )
  END SUBROUTINE emis_interp

END MODULE TELSEM2_Atlas_Module
