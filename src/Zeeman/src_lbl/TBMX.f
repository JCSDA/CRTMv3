      SUBROUTINE TBMX(TBRC,TBLC,TBLIN,H,TAV,AL,GP,GM,GPI,CBTH)
C
C  LANGUAGE- FORTRAN 77
C
C  VERSION- 1    DATE- 11/28/87    PROGRAMMER- P. ROSENKRANZ
C
C  FUNCTION-
C  COMPUTES MATRIX RADIATIVE TRANSFER THROUGH A HOMOGENEOUS LAYER
C
C  ARGUMENT SPECIFICATIONS-
      REAL TBRC,TBLC,H,TAV,AL,CBTH
      COMPLEX GP,GM,GPI,TBLIN
C
C  NAME    IN/OUT  UNITS    DESCRIPTON
C
C  TBRC     I,O    KELVIN   RIGHT-CIRCULARLY POLARIZED BRIGHTNESS TEMPERATURE
C  TBLC     I,O    KELVIN   LEFT-CIRCULARLY POLARIZED BRIGHTNESS TEMPERATURE.
C  TBLIN    I,O    KELVIN   COMPLEX OFF-DIAGONAL TERM OF TB COHERENCY MATRIX;
C                           IN TERMS OF STOKES PARAMETERS, Q=2*REAL(TBLIN),
C                           U=2*AIMAG(TBLIN).
C                           ON INPUT, TBRC, TBLC AND TBLIN SHOULD HAVE INITIAL
C                           VALUES AT ONE SIDE OF THE LAYER. ON OUTPUT, 
C                           THEY HAVE BEEN COMPUTED FOR THE OTHER SIDE OF THE LAYER.
C  H         I     KM       PATH LENGTH THROUGH LAYER
C  TAV       I     KELVIN   TEMPERATURE
C  AL        I     1/KM     POWER ABSORPTION COEFF NOT DEPENDING ON POLARIZATION;
C                           I.E., FROM DISTANT LINES.
C  GP        I     1/KM     COMPLEX AMPLITUDE PROPAGATION COEFFICIENT OF
C                           ZEEMAN-SPLIT LINE, FOR DELTA M = +1.
C  GM        I     1/KM     COMPLEX AMPLITUDE PROPAGATION COEFFICIENT OF
C                           ZEEMAN-SPLIT LINE, FOR DELTA M = -1.
C  GPI       I     1/KM     COMPLEX AMPLITUDE PROPAGATION COEFFICIENT OF
C                           ZEEMAN-SPLIT LINE, FOR DELTA M = 0.
C  CBTH      I              COSINE(ANGLE BETWEEN MAGNETIC FIELD AND PROPAGATION
C                           DIRECTION)
C********************************************************************************
C
C  LOCAL VARIABLES-
      COMPLEX C11,C12,C22,CP,CM,CPI,DELTA,
     &  E11,E22,E12,C2,DUM,EL1,EL2,D,C212,AB,BB,CB
      SQMAG(DUM) = REAL(DUM)**2 + AIMAG(DUM)**2
C
      PB = AL*.5*H
      CP = GP*H
      CM = GM*H
      CPI = GPI*H
      SIN2 = .5*(1.-CBTH**2)
      OPC = .5*(1.+CBTH)**2
      OMC = .5*(1.-CBTH)**2
      C11 = OPC*CP + OMC*CM + SIN2*CPI +CMPLX(PB,0.)
      C22 = OMC*CP + OPC*CM + SIN2*CPI + CMPLX(PB,0.)
      C12 = SIN2*(CP+CM-CPI)
      IF(ABS(REAL(C12))+ABS(AIMAG(C12))) 9,9,10
C      MATRIX IS DIAGONAL
9      E11 = CEXP(-C11)
      E22 = CEXP(-C22)
      E12 = CMPLX(0.,0.)
      GOTO 20
10      C2 = C22 - C11
      DELTA = CSQRT(C2*C2 + 4.*C12*C12)
      IF(REAL(C2)*REAL(DELTA).LT.0.) DELTA = -DELTA
      C2 = .5*(C2+DELTA)
      EL2 = CEXP(-C2-C11)
      EL1 = EL2*CEXP(DELTA)
      C212 = C2*C12
      C2 = C2*C2
      C12 = C12*C12
      D = C2 + C12
      E11 = (C2*EL1+C12*EL2)/D
      E22 = (C12*EL1+C2*EL2)/D
      E12 = C212*(EL2-EL1)/D
20      TBR = TBRC - TAV
      TBL = TBLC - TAV
      AA = SQMAG(E11)
      AB = E11*CONJG(E12)
      AC = SQMAG(E12)
      TBRC = TAV + TBR*AA + 2.*REAL(TBLIN*AB) + TBL*AC
      CB = E12*CONJG(E22)
      CC = SQMAG(E22)
      TBLC = TAV + TBR*AC + 2.*REAL(TBLIN*CB) + TBL*CC
      BB = E11*CONJG(E22)
      TBLIN = TBR*AB + TBLIN*BB + CONJG(TBLIN)*AC + TBL*CB
      RETURN
      END
