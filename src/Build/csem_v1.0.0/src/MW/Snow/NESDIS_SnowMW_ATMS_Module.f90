!--------------------------------------------------------------------------------
!M+
! NAME:
!       NESDIS_ATMS_SnowEM_Module
!
! PURPOSE:
!       Module containing the snow-typing algorithms. A general interface is used to call one of the 
!       snow-typing algorithms in !terms of the input arguments. This Module is used together with 
!       NESDIS_SnowEM_ATMS_Parameters Module to implement the library-based snow emissivity  model.
!
! REFERENCES:
!       Yan, B., F. Weng and K.Okamoto,2004: "A microwave snow emissivity model, 8th Specialist Meeting on
!
!       Microwave Radiometry and Remote Sension Applications,24-27 February, 2004, Rome, Italy.
!
! CATEGORY:
!       Surface : MW Surface Snow Emissivity of ATMS
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!
!       USE NESDIS_ATMS_SnowEM_Module
!
! MODULES:
!       Type_Kinds:               Module containing definitions for kinds of variable types
!
!       NESDIS_MWEmiss_Snow_Module:     Module containing the microwave land emissivity model
!
!       NESDIS_SnowEM_Parameters: Module containing the predefined microwave snow emissivity spectra
!
! CONTAINS:
!
! PUBLIC SUNPROGRAMS:
!
!       NESDIS_ATMS_SNOWEM:       Subroutine to calculate the microwave snow emissivity from ATMS
!
!
! PRIVATE SUBPROGRAMS:
!       These subroutines are used to determine the snow types from the brightness temperatures(TB) 
!       of five ATMS window channels( 23.8 GHz, 31.4 GHz, 50.3 GHz, 88.2 GHz, 165.5 GHz) and/or
!       surface temperature plus snow depth. The five channels are further divided into two 
!       groups: Group-1 ( 23.8 GHz, 31.4 GHz, 50.3 GHz, 88.2 GHz) and Group-2 (88.2 GHz, 165.5GHz), 
!       corresponding to the window channels of AMSU-A and AMSU-B, respectively.
!       Different combinations of available ATMS window-channel and surface observations result
!       in differenet snow-typing algotrithms:   
!
!       ATMS_SNOW_ByTB_A      : by ATMS TBs of Group-1 channels 
!       ATMS_SNOW_ByTB_B      : by ATMS TBs of Group-2 channels 
!       ATMS_SNOW_ByTBs       : by the TBs  of all the five ATMS channels 
!       ATMS_SNOW_ByTBTs_A    : by ATMS TBs of Group-1 channels and surface temperature
!       ATMS_SNOW_ByTBTs_B    : by ATMS TBs of Group-2 channels and surface temperature
!       ATMS_SNOW_ByTBTs      : by the TBs  of all the five ATMS channels and surface temperature (regression-based)
!       ATMS_SNOW_ByTBTs_D    : by the TBs  of all the five ATMS channels and surface temperature (diagnosis-based)
!       ATMS_SNOW_ByTypes     : bydefault surface type (4)
!       ATMS_ALandEM_Snow     : Subroutine to initilize the vaiables to default values 
!       em_initialization     : Subroutine to initialization snow emissivity
!       em_interpolate        : Subroutine to perform frequency interpolation of snow emissivity
!
! INCLUDE FILES:
!       None.
!
! EXTERNALS:
!       None.
!
! COMMON BLOCKS:
!       None.
!
! FILES ACCESSED:
!       None.
!
! CREATION HISTORY:
!       Written by:     Ming Chen, IMSG Inc., Banghua.Yan@noaa.gov (04-28-2012)
!
!
!       and             Fuzhong Weng, NOAA/NESDIS/ORA, Fuzhong.Weng@noaa.gov
!
!  Copyright (C) 2012 Fuzhong Weng and Ming Chen
!
!  This program is free software; you can redistribute it and/or modify it under the terms of the GNU
!  General Public License as published by the Free Software Foundation; either version 2 of the License,
!  or (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
!  License for more details.
!
!  You should have received a copy of the GNU General Public License along with this program; if not, write
!  to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
!M-
!--------------------------------------------------------------------------------

MODULE NESDIS_ATMS_SnowEM_Module

  USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp, ip => CSEM_ip, &
                            Double => CSEM_Double
  USE NESDIS_MW_SnowEmiss_Util
  USE NESDIS_SnowEM_ATMS_Parameters
  IMPLICIT NONE

! Visibilities
  PRIVATE 
  PUBLIC  :: NESDIS_ATMS_SNOWEM
  REAL(fp), PARAMETER :: ONE = 1.0_fp, ZERO = 0.0_fp
  REAL(fp), PARAMETER :: PI = 3.141592653589793238462643_fp
  

CONTAINS


   !################################################################################
   !################################################################################
   !##                                                                            ##
   !##                         ## PUBLIC MODULE ROUTINES ##                       ##
   !##                                                                            ##
   !################################################################################
   !################################################################################

   !-------------------------------------------------------------------------------------------------------------
   !
   ! NAME:
   !       NESDIS_ATMS_SNOWEM
   !
   ! PURPOSE:
   !       Subroutine to simulate microwave emissivity over snow conditions from ATMS measurements at window
   !       channels.
   !
   !
   ! CATEGORY:
   !       CRTM : Surface : MW SNOWEM
   !
   ! LANGUAGE:
   !       Fortran-95
   !
   ! CALLING SEQUENCE:
   !       CALL NESDIS_ATMS_SNOWEM
   !
   ! INPUT ARGUMENTS:
   !
   !         Frequency                Frequency User defines
   !                                  This is the "I" dimension
   !                                  UNITS:      GHz
   !                                  TYPE:       REAL( fp )
   !                                  DIMENSION:  Scalar
   !
   !
   !         Satellite_Angle          The local zenith angle in degree for ATMS measurements.
   !                                  ** NOTE: THIS IS A MANDATORY MEMBER **
   !                                  **       OF THIS STRUCTURE          **
   !                                  UNITS:      Degrees
   !                                  TYPE:       REAL( fp )
   !                                  DIMENSION:  Rank-1, (I)
   !
   !         User_Angle               The local angle value in degree user defines.
   !                                  ** NOTE: THIS IS A MANDATORY MEMBER **
   !                                  **       OF THIS STRUCTURE          **
   !                                  UNITS:      Degrees
   !                                  TYPE:       REAL( fp )
   !                                  DIMENSION:  Rank-1, (I)
   !
   !
   !         Tbs                      BRIGHTNESS TEMPERATURES AT FIVE ATMS WINDOW CHANNELS
   !                                  UNITS:      Kelvin, K
   !                                  TYPE:       REAL( fp )
   !                                  DIMENSION   4*1 SCALAR
   !
   !                        WHICH ARE
   !                                  tbs[1] = TB at 23.8 GHz
   !                                  tbs[2] = TB at 31.4 GHz
   !                                  tbs[3] = TB at 50.3 GHz
   !                                  tbs[4] = TB at 88.2 GHz
   !                                  tbs[5] = TB at 165.5 GHz
   !
   !
   !         Tss = Land_Temperature:  The land surface temperature.
   !                                  UNITS:      Kelvin, K
   !                                  TYPE:       REAL( fp )
   !                                  DIMENSION:  Scalar
   !
   !
   !         Snow_Depth:              The snow depth.
   !                                  UNITS:      mm
   !                                  TYPE:       REAL( fp )
   !                                  DIMENSION:  Scalar
   !
   ! **** IMPORTANT NOTES ****
   !
   !        When one variable among  Tbs[],  Ts are not available, set -999.0
   !
   !
   !
   !
   ! OUTPUT ARGUMENTS:
   !
   !         Emissivity_H:            The surface emissivity at a horizontal polarization.
   !                                  ** NOTE: THIS IS A MANDATORY MEMBER **
   !                                  **       OF THIS STRUCTURE          **
   !                                  UNITS:      N/A
   !                                  TYPE:       REAL( fp )
   !                                  DIMENSION:  Scalar
   !
   !         Emissivity_V:            The surface emissivity at a vertical polarization.
   !                                  ** NOTE: THIS IS A MANDATORY MEMBER **
   !                                  **       OF THIS STRUCTURE          **
   !                                  UNITS:      N/A
   !                                  TYPE:       REAL( fp )
   !                                  DIMENSION:  Scalar
   !
   !
   !
   !
   ! OPTIONAL OUTPUT ARGUMENTS:
   !
   !       snow_type  -  snow type (not output here)
   !                     1 : Wet Snow
   !                     2 : Grass_after_Snow
   !                     3 : RS_Snow (A)
   !                     4 : Powder Snow
   !                     5 : RS_Snow (B)
   !                     6 : RS_Snow (C)
   !                     7 : RS_Snow (D)
   !                     8 : Thin Crust Snow
   !                     9 : RS_Snow (E)
   !                     10: Bottom Crust Snow (A)
   !                     11: Shallow Snow
   !                     12: Deep Snow
   !                     13: Crust Snow
   !                     14: Medium Snow
   !                     15: Bottom Crust Snow (B)
   !                     16: Thick Crust Snow
   !                   -999: ATMS measurements are not available or over non-snow conditions
   !
   ! CALLS:
   !       ATMS_SNOW_ByTB_A      : by ATMS TBs of Group-1 channels
   !       ATMS_SNOW_ByTB_B      : by ATMS TBs of Group-2 channels
   !       ATMS_SNOW_ByTBs       : by the TBs  of all the five ATMS channels
   !       ATMS_SNOW_ByTBTs_A    : by ATMS TBs of Group-1 channels and surface temperature
   !       ATMS_SNOW_ByTBTs_B    : by ATMS TBs of Group-2 channels and surface temperature
   !       ATMS_SNOW_ByTBTs      : by the TBs  of all the five ATMS channels and surface temperature
   !       ATMS_SNOW_ByTypes     : bydefault surface type (4)
   !       ATMS_ALandEM_Snow     : Subroutine to initilize the vaiables to default values
   !       em_initialization     : Subroutine to initialization snow emissivity
   !       NESDIS_MWEmiss_Snow         : EM physical model for angular interpolations
   !
   ! SIDE EFFECTS:
   !       None.
   !
   ! RESTRICTIONS:
   !       None.
   !
   ! COMMENTS:
   !       Note the INTENT on the output SensorData argument is IN OUT rather than
   !       just OUT. This is necessary because the argument may be defined upon
   !       input. To prevent memory leaks, the IN OUT INTENT is a must.
   !
   ! CREATION HISTORY:
   !       Written by:
   !                       Ming Chen, IMSG Inc., ming.chen@noaa.gov (04-28-2012)
   !                       Fuzhong Weng, NOAA/NESDIS/ORA, Fuzhong.Weng@noaa.gov
   !
   !
   !  Copyright (C) 2005 Fuzhong Weng and Ming Chen
   !
   !  This program is free software; you can redistribute it and/or modify it under the terms of the GNU
   !  General Public License as published by the Free Software Foundation; either version 2 of the License,
   !  or (at your option) any later version.
   !
   !  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
   !  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
   !  License for more details.
   !
   !  You should have received a copy of the GNU General Public License along with this program; if not, write
   !  to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
   !
   !------------------------------------------------------------------------------------------------------------


   SUBROUTINE  NESDIS_ATMS_SNOWEM(Satellite_Angle,                             &  ! INPUT
                                  User_Angle,                                  &  ! INPUT
                                  Frequency,                                   &  ! INPUT
                                  Tbs,                                         &  ! INPUT
                                  Tss,                                         &  ! INPUT
                                  Snow_Depth,                                  &  ! INPUT
                                  Emissivity_H,                                &  ! OUPUT
                                  Emissivity_V)                                   ! OUTPUT



     IMPLICIT NONE

     INTEGER, PARAMETER :: nch = N_FREQ_ATMS, nwch = 5
     REAL(fp),INTENT(IN) :: Satellite_Angle,User_Angle,Frequency
     REAL(fp),INTENT(IN),  OPTIONAL  :: Tbs(:), Tss, Snow_Depth
     REAL(fp),INTENT(OUT)            :: Emissivity_H,Emissivity_V
     REAL(fp) :: em_vector(2),esh1,esv1,esh2,esv2,desh,desv,dem
     REAL(fp) :: Ts = 273.15 ! default skin-surface temperature
     INTEGER :: Snow_Type = 4 ! default snow type
     INTEGER :: i

     LOGICAL  :: VALID_SNOW_DEPTH = .FALSE.
     INTEGER  :: input_type = 0

   ! Analyze the input types and determine which algorithms to be used

     IF (PRESENT(Snow_Depth) .AND. PRESENT(Tss) ) THEN
       VALID_SNOW_DEPTH = .TRUE.
       IF ((Snow_Depth < 0.0_fp) .OR. (Snow_Depth >= 3000.0_fp)) &
       VALID_SNOW_DEPTH = .FALSE.
     ENDIF

     IF (PRESENT(TBs)) THEN
        IF (SIZE(Tbs) >= 4) input_type=IBSET(input_type, 0)
        IF (SIZE(Tbs) >= 5) input_type=IBSET(input_type, 1)
        DO i=1,SIZE(Tbs)
          IF ((Tbs(i) <= 100.0_fp) .OR. (Tbs(i) >= 290.0_fp) ) THEN
            IF (i <= 4) input_type=IBCLR(input_type, 0)
            IF (i >= 4) input_type=IBCLR(input_type, 1)
          ENDIF
        END DO
     ENDIF

     IF (PRESENT(Tss) ) THEN
        input_type=IBSET(input_type, 2)
        Ts=Tss
        IF ((Ts <= 150.0_fp) .OR. (Ts >= 290.0_fp) ) THEN
           input_type=IBCLR(input_type, 2)
           VALID_SNOW_DEPTH = .FALSE.
        ENDIF
     ENDIF


   ! Initialization
     CALL em_initialization(frequency,em_vector)

   ! the above regression-based snow-typing algs are superseded by the diagnosis-based snow-typing
     CALL ATMS_SNOW_ByTBTs_D(Frequency,Tbs,Ts,Snow_Type,em_vector)

   ! Get the emissivity angle dependence
     CALL NESDIS_MWEmiss_Snow(Satellite_Angle,frequency,&
          0.0_fp,0.0_fp,Ts,Ts,0.0_fp,9,13,2.0_fp,esh1,esv1)
     CALL NESDIS_MWEmiss_Snow(User_Angle,frequency,&
          0.0_fp,0.0_fp,Ts,Ts,0.0_fp,9,13,2.0_fp,esh2,esv2)

     desh = esh1 - esh2
     desv = esv1 - esv2
     dem = ( desh + desv ) * 0.5_fp

   ! Emissivity at User's Angle
     Emissivity_H = em_vector(1) - dem;  Emissivity_V = em_vector(2)- dem

     IF (Emissivity_H > one) Emissivity_H = one
     IF (Emissivity_V > one) Emissivity_V = one

     IF (Emissivity_H < 0.3_fp) Emissivity_H = 0.3_fp
     IF (Emissivity_V < 0.3_fp) Emissivity_V = 0.3_fp

     RETURN

   END SUBROUTINE NESDIS_ATMS_SNOWEM




   !##################################################################################
   !##################################################################################
   !##                                                                              ##
   !##                          ## PRIVATE MODULE ROUTINES ##                       ##
   !##                                                                              ##
   !##################################################################################
   !##################################################################################




   SUBROUTINE ATMS_SNOW_ByTBTs_D(frequency,tb,ts,snow_type,em_vector)
   !----------------------------------------------------------------------------------------------------------!
   !$$$  subprogram documentation block
   !
   ! subprogram: Calculate emissivity by diagnosis-based algorithm
   !
   !
   ! abstract:
   !   Diagnose the snow type, and use the emissivity spectrum of the snow type as the first-guess to diagnose
   !   the magnitude of necessary adjustment with respect to window-channel TBs and surface temperature Ts. Perfrom
   !   necessary interpolation/extrapolation a required frequency and user angle.
   !
   !
   ! input argument list:
   !
   !     frequency        -  frequency in GHz
   !     theta            -  local zenith angle (currently, not used here)
   !     tb[1] ~ tb[5]    -  brightness temperature at five ATMS window channels:
   !                              tb[1] : 23.8 GHz
   !                              tb[2] : 31.4 GHz
   !                              tb[3] : 50.3 GHz
   !                              tb[4] : 88.2 GHz
   !                              tb[5] : 165.5 GHz
   !
   ! output argument list:
   !
   !      em_vector[1] and [2]  -  emissivity at two polarizations.
   !                              set esv = esh here and will be updated
   !      snow_type        -  snow type
   !
   !
   ! remarks:
   !
   ! program history log:
   !            Ming Chen, IMSG at NOAA/NESDIS/STAR                 date: 2012-04-28
   !
   !
   !----------------------------------------------------------------------------------------------------------!

     IMPLICIT NONE

     INTEGER , PARAMETER  :: ntype = N_SNOW_TYPES, nch = N_FREQ_ATMS, nwch = 5
     REAL(fp), PARAMETER  :: earthrad = 6374._fp, satheight = 833.4_fp
     INTEGER  :: freq_idx,snow_type
     REAL(fp) :: frequency
     REAL(fp) :: em(nch,ntype), em_vector(:)
     REAL(fp) :: tb(:),freq(nch)
     REAL(fp) :: ts, emissivity
     REAL(fp) :: ediff(ntype), X(nwch),Y(nwch),emw(nwch)
     REAL(fp) :: XX,XY,del,dem,dem2,delta,deltb
     INTEGER  :: minlc(1)
     INTEGER  :: windex(nwch)=(/1,2,3,11,12/)             ! window channel index of the library spectrum
     ! Coefficients of quadratic equations used to estimate the emissivity difference
     ! between Ch-31.4GHz and 88.2GHZ
     REAL(fp) :: coe1(6)=(/ -0.837001E+00_fp, 0.954882E-02_fp, -0.271471E-04_fp, &
                             -0.536112E-02_fp, 0.144279E-04_fp, 0.218317E-02_fp/)
     ! between Ch-31.4GHz and 165.5GHZ
     REAL(fp) :: coe2(6)=(/ -0.854226E+00_fp, 1.203220E-02_fp, -0.216043E-04_fp, &
                             -0.887216E-02_fp, 0.118303E-04_fp, 0.263699E-02_fp/)
     ! Quadratic EQ terms (1.0,tb(4),tb(4)^2,tb(5),tb(5)^2,Ts)
     REAL(fp) :: coe3(6)

   ! Sixteen candidate snow emissivity spectra
     em = SNOW_EMISS_ATMS_LIB
     freq = FREQUENCY_ATMS


     minlc =minloc(ABS(freq-frequency)); freq_idx=minlc(1)

   !*** IDENTIFY SNOW TYPE
     snow_type = 4 !default
     ediff=abs(Tb(1)/em(1,:)-Tb(2)/em(2,:))+abs(Tb(2)/em(2,:)-Tb(4)/em(11,:))
     minlc = minloc(ediff) ; snow_type=minlc(1)

   !*** adjustment from the library values
     emw=em(windex,snow_type)
     X=1.0_fp/emw ; Y=LOG(Tb/(Ts*emw))
     IF(frequency > 100.0_fp) THEN
       XX=DOT_PRODUCT(X((/1,2,4,5/)),X((/1,2,4,5/)))
       XY=DOT_PRODUCT(X((/1,2,4,5/)),Y((/1,2,4,5/)))
       del=XY/XX
       deltb=Tb(3)-Tb(5)
     ELSE
       XX=DOT_PRODUCT(X((/1,2,4/)),X((/1,2,4/)))
       XY=DOT_PRODUCT(X((/1,2,4/)),Y((/1,2,4/)))
       del=XY/XX
       deltb=Tb(3)-Tb(4)
     ENDIF
     dem = 0.0_fp; delta = 0.0_fp
     IF(frequency <= 30.0_fp ) dem = 1.1_fp*del
     IF(frequency > 30._fp .AND. frequency <= 40.0_fp ) dem = 1.05_fp*del
     IF(frequency > 40._fp .AND. frequency <= 50.0_fp ) dem = 1.0_fp*del
     IF(frequency > 50_fp) THEN
        IF(del .LE. 0.0_fp .AND. ABS(deltb) .LT. 30.0_fp) delta=0.5_fp+deltb/50.0_fp
        IF(del .LE. 0.0_fp .AND. ABS(deltb) .GE. 30.0_fp) delta=1.0_fp+deltb/50.0_fp
        IF(del .GT. 0.0_fp .AND. ABS(deltb) .LT. 35.0_fp) delta=1.05_fp-deltb/70.0_fp
        IF(del .GT. 0.0_fp .AND. ABS(deltb) .GE. 35.0_fp) delta=.85_fp-deltb/70.0_fp
        IF(frequency <= 100.0_fp) dem  = del+(delta*del-del)*(frequency-50.0_fp)/(100.0_fp-50.0_fp)
        IF(frequency >  100.0_fp) dem  = 0.65_fp*delta*del
     ENDIF
     dem2=dem
     IF (frequency > 80.0_fp  )THEN
       coe3=(/1.0_fp,tb(4),tb(4)*tb(4),tb(5),tb(5)*tb(5),Ts/)
       IF(del<-0.13_fp)del=-0.13_fp
       IF(frequency <= 100.0_fp )THEN
         dem2=del-SUM(coe1*coe3)
       ELSE
         dem2=del-SUM(coe2*coe3)
       ENDIF
     ENDIF

     emissivity = em(freq_idx,snow_type)+(dem+dem2)/2.0_fp
     IF (emissivity >  1.0_fp )emissivity = 1.0_fp
     IF (emissivity <= 0.3_fp )emissivity = 0.3_fp

     em_vector(1) = emissivity
     em_vector(2) = emissivity

    RETURN
   END SUBROUTINE ATMS_SNOW_ByTBTs_D



   SUBROUTINE em_initialization(frequency,em_vector)

   !----------------------------------------------------------------------------------------------------------!
   !$$$  subprogram documentation block
   !
   ! subprogram: ATMS snow emissivity initialization
   !
   !
   ! abstract:   ATMS snow emissivity initialization
   !
   !
   ! input argument list:
   !
   !      frequency   - frequency in GHz
   !
   ! output argument list:
   !
   !     em_vector[1] and [2]  -  initial emissivity at two polarizations.
   !
   ! important internal variables:
   !
   !      freq[1~10]  - ten frequencies for sixteen snow types of emissivity
   !      em[1~16,*]  - sixteen snow emissivity spectra
   !      snow_type   - snow type
   !                    where it is initialized to as the type 4,i.e, Powder Snow
   !
   ! remarks:
   !
   ! program history log:
   !            Banghua Yan, nesdis                                 date: 2003-08-18
   !            Ming Chen, IMSG at NOAA/NESDIS/STAR                 date: 2012-04-28
   !
   ! attributes:
   !   language: f90
   !   machine:  ibm rs/6000 sp
   !
   !----------------------------------------------------------------------------------------------------------!

     IMPLICIT NONE

     INTEGER,PARAMETER :: nch = N_FREQ_ATMS,ncand=N_SNOW_TYPES
     REAL(fp) :: frequency,em_vector(:),freq(nch)
     REAL(fp) :: em(ncand,nch)
     REAL(fp) :: kratio, bconst,emissivity
     INTEGER :: ich

   ! Sixteen candidate snow emissivity spectra

     em = TRANSPOSE(SNOW_EMISS_ATMS_LIB)
     freq = FREQUENCY_ATMS

   ! Initialization for emissivity at certain frequency
   !    In case of no any inputs available for various options
   !    A constant snow type & snow emissivity spectrum is assumed
   !                    (e.g., powder) snow_type = 4
   emissivity = 1.0_fp
   ! Specify snow emissivity at required frequency
     DO ich = 2, nch
        IF (frequency <  freq(1))   EXIT
        IF (frequency >= freq(nch)) EXIT
        IF (frequency <  freq(ich)) THEN
           emissivity = em(4,ich-1) + (em(4,ich) - em(4,ich-1))     &
                *(frequency - freq(ich-1))/(freq(ich) - freq(ich-1))
           EXIT
        ENDIF
     END DO

   ! Extrapolate to lower frequencies than 4.9GHz
     IF (frequency <= freq(1)) THEN
        kratio = (em(4,2) - em(4,1))/(freq(2) - freq(1))
        bconst = em(4,1) - kratio*freq(1)
        emissivity =  kratio*frequency + bconst
        IF (emissivity >  one)         emissivity = one
        IF (emissivity <= 0.8_fp) emissivity = 0.8_fp
     ENDIF


   ! Assume emissivity = constant at frequencies >= 150 GHz
     IF (frequency >= freq(nch)) emissivity = em(4,nch)
     em_vector(1) = emissivity
     em_vector(2) = emissivity

     RETURN
   END SUBROUTINE em_initialization



END MODULE NESDIS_ATMS_SnowEM_Module
