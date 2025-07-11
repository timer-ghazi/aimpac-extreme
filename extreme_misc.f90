! Generated from extreme.f90
! Contains: ALPHA, AZero, BOND, BLOCK DATA


SUBROUTINE ALPHA
!
   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
   Parameter (MaxOff=200000,MaxAtm=100)
   COMMON CO(MaxOff),IC(MaxOff),NCENT,NMO,NPRIMS
   COMMON /ANG/ ANGLE(3,MaxAtm,MaxAtm)
   COMMON /OFFSET/ ITYPE,ICENT,IMO,IEORB,IE,ICHARG,IXC,IYC,IZC,&
   &IXX, IYY, IZZ,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,IGXX,IGXY,IGXZ,IGYY,&
   &IGYZ,IGZZ,IGXXX,IGXXY,IGXXZ,IGXYY,IGXZZ,IGXYZ,IGYYY,IGYYZ,&
   &IGYZZ,IGZZZ,IGXXXX,IGXXXY,IGXXXZ,IGXXYY,IGXXZZ,IGXYYY,IGXZZZ,&
   &IGXXYZ,IGXYYZ,IGXYZZ,IGYYYY,IGYYYZ,IGYYZZ,IGYZZZ,IGZZZZ,ICOFMx
   COMMON /UNITS/ INPT, IOUT, IWFN, IWLP
   Save One,ThSxty
   Data One/1.d0/,ThSxty/3.6d2/
1000 FORMAT(' HERE IS A LIST OF BOND PATHS THAT HAVE BEEN TRACED ')
1010 FORMAT(/,' INPUT ATOMS FOR WHICH ANGLE IS DESIRED (EG.1,2,3) ',$)
1020 FORMAT(/,' BOND PATH ANGLE = ',F8.4)
1030 FORMAT(/,' GEOMETRIC BOND ANGLE = ',F8.4)
1040 FORMAT(/,' DEVIATION OF FIRST BOND PATH VECTOR FROM FIRST',&
   &' GEOMETRIC BOND VECTOR   ',F8.4)
1050 FORMAT(/,' DEVIATION OF SECOND BOND PATH VECTOR FROM SECOND',&
   &' GEOMETRIC BOND VECTOR ',F8.4)
1060 FORMAT(/,' GEOMETRIC BOND ANGLE MINUS BOND PATH ANGLE ',F8.4)
1070 FORMAT(/,' WOULD YOU LIKE TO DO ANOTHER COMPARISION ? (0/1) ',$)
1080 FORMAT(/,' COMPARISION FOR ATOMS ',3I4)
!
   Pi=DaCos(-One)
   Degree=thsxty/(Two*Pi)
!
100 WRITE (IOUT,1010)
   READ (INPT,*) I,J,K
   WRITE (IOUT,1080) I,J,K
   WRITE (IWLP,1080) I,J,K
!
!    GET BOND PATH ANGLE
!
   S = ANGLE(1,J,K)*ANGLE(1,J,I) +&
   &ANGLE(2,J,K)*ANGLE(2,J,I) +&
   &ANGLE(3,J,K)*ANGLE(3,J,I)
   BANGLE = (DACOS(S))*Degree
   WRITE (IOUT,1020) BANGLE
   WRITE (IWLP,1020) BANGLE
!
!    USE UNIT VECTORS ALONG THE BONDS TO CALCULATE
!    GEOMETRIC BOND ANGLE
!
   C1 = CO(IXC+I) - CO(IXC+J)
   D1 = CO(IXC+K) - CO(IXC+J)
   C2 = CO(IYC+I) - CO(IYC+J)
   D2 = CO(IYC+K) - CO(IYC+J)
   C3 = CO(IZC+I) - CO(IZC+J)
   D3 = CO(IZC+K) - CO(IZC+J)
   S1 = SQRT(C1**2 + C2**2 + C3**2)
   S2 = SQRT(D1**2 + D2**2 + D3**2)
   C1 = C1/S1
   C2 = C2/S1
   C3 = C3/S1
   D1 = D1/S2
   D2 = D2/S2
   D3 = D3/S2
!
!      GET ANGLES BETWEEN BOND PATHS AND BOND DIRECTIONS
!
   S1 = ANGLE(1,J,I)*C1 + ANGLE(2,J,I)*C2 + ANGLE(3,J,I)*C3
   S2 = ANGLE(1,J,K)*D1 + ANGLE(2,J,K)*D2 + ANGLE(3,J,K)*D3
   S3 = C1*D1 + C2*D2 + C3*D3
   ANGLE1 = (DACOS(S1))*Degree
   ANGLE2 = (DACOS(S2))*Degree
   ANGLE3 = (DACOS(S3))*Degree
   WRITE (IOUT,1030) ANGLE3
   WRITE (IWLP,1030) ANGLE3
   WRITE (IOUT,1040) ANGLE1
   WRITE (IWLP,1040) ANGLE1
   WRITE (IOUT,1050) ANGLE2
   WRITE (IWLP,1050) ANGLE2
!
!    DETERMINE THE DIFFERENCE BETWEEN THE BOND PATH ANGLE AND THE
!    GEOMETRIC BOND ANGLE
!
   DANGLE = ANGLE3 - BANGLE
   WRITE (IOUT,1060) DANGLE
   WRITE (IWLP,1060) DANGLE
!
   WRITE (IOUT,1070)
   READ (INPT,*) IGO
   IF (IGO .EQ. 1) GOTO 100
!
   RETURN
!
!    FORMATS
!
END

Subroutine AZero
!
   Implicit Double Precision (A-H,O-Z)
!
   Common /AFG/ A0,AX,AY,AZ,AXX,AXY,AXZ,AYY,AYZ,AZZ,AXXX,&
   &AXXY,AXXZ,AXYY,AXZZ,AXYZ,AYYY,AYYZ,AYZZ,AZZZ,AXXXX,AXXXY,&
   &AXXXZ,AXXYY,AXXZZ,AXYYY,AXZZZ,AXXYZ,AXYYZ,AXYZZ,AYYYY,&
   &AYYYZ,AYYZZ,AYZZZ,AZZZZ,F0,FX,FY,FZ,FXX,FXY,FXZ,FYY,FYZ,&
   &FZZ,FXXX,FXXY,FXXZ,FXYY,FXZZ,FXYZ,FYYY,FYYZ,FYZZ,FZZZ,&
   &FXXXX,FXXXY,FXXXZ,FXXYY,FXXZZ,FXYYY,FXZZZ,FXXYZ,FXYYZ,&
   &FXYZZ,FYYYY,FYYYZ,FYYZZ,FYZZZ,FZZZZ,G0,GX,GY,GZ,GXX,GXY,&
   &GXZ,GYY,GYZ,GZZ,GXXX,GXXY,GXXZ,GXYY,GXZZ,GXYZ,GYYY,GYYZ,&
   &GYZZ,GZZZ,GXXXX,GXXXY,GXXXZ,GXXYY,GXXZZ,GXYYY,GXZZZ,&
   &GXXYZ,GXYYZ,GXYZZ,GYYYY,GYYYZ,GYYZZ,GYZZZ,GZZZZ
   Save Zero
   Data Zero/0.0d0/
!
   A0=Zero
   AX=Zero
   AY=Zero
   AZ=Zero
   AXX=Zero
   AXY=Zero
   AXZ=Zero
   AYY=Zero
   AYZ=Zero
   AZZ=Zero
   AXXX=Zero
   AXXY=Zero
   AXXZ=Zero
   AXYY=Zero
   AXZZ=Zero
   AXYZ=Zero
   AYYY=Zero
   AYYZ=Zero
   AYZZ=Zero
   AZZZ=Zero
   AXXXX=Zero
   AXXXY=Zero
   AXXXZ=Zero
   AXXYY=Zero
   AXXZZ=Zero
   AXYYY=Zero
   AXZZZ=Zero
   AXXYZ=Zero
   AXYYZ=Zero
   AXYZZ=Zero
   Ayyyy=Zero
   Ayyyz=Zero
   Ayyzz=Zero
   Ayzzz=Zero
   Azzzz=Zero
!
   Return
End

SUBROUTINE BOND(X,Y,Z,BL,NCAL,PN,IFunc,Iwhole)
!
   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
   Parameter(MaxOff=200000,MaxAtm=100,MaxStp=141,MinStp=6)
   COMMON CO(MaxOff),IC(MaxOff),NCENT,NMO,NPRIMS
   COMMON /DIST/ BANGLE(6),RMIN,RMAX,MINR
   COMMON /OFFSET/ ITYPE,ICENT,IMO,IEORB,IE,ICHARG,IXC,IYC,IZC,&
   &IXX, IYY, IZZ,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,IGXX,IGXY,IGXZ,IGYY,&
   &IGYZ,IGZZ,IGXXX,IGXXY,IGXXZ,IGXYY,IGXZZ,IGXYZ,IGYYY,IGYYZ,&
   &IGYZZ,IGZZZ,IGXXXX,IGXXXY,IGXXXZ,IGXXYY,IGXXZZ,IGXYYY,IGXZZZ,&
   &IGXXYZ,IGXYYZ,IGXYZZ,IGYYYY,IGYYYZ,IGYYZZ,IGYZZZ,IGZZZZ,ICofMx
   Common /What/ IWhat
   DIMENSION R1(3,MaxStp),PN(4),PR(3,MaxStp),H(3,3),SG(3,3),&
   &PPN(3),GRHO(3),tmpxyz(3)
   Save Zero,dxyz,pt15,fifty,pt2,one,two,three,fvhund,&
   &pt01,pt001,small
   Data Zero/0.0d0/,dxyz/1.d-3/,pt15/0.15d0/,fifty/50.0d0/,&
   &pt2/0.2d0/,One/1.0d0/,Two/2.0d0/,Three/3.0d0/,&
   &fvhund/500.0d0/,pt01/0.01d0/,pt001/0.001d0/,small/1.d-9/,&
   &Hund/1.d2/
!
   Iwhat=0
   Pi=Dacos(-one)
   BL = Zero
   NCAL = 0
   R1(1,1) = X
   R1(2,1) = Y
   R1(3,1) = Z
   tmpxyz(1)=R1(1,1)
   tmpxyz(2)=R1(2,1)
   tmpxyz(3)=R1(3,1)
   CALL GRDRHO(0,tmpxyz,RHO,GRHO,GRAD,H,SG)
   PR(1,1)=GRHO(1)
   PR(2,1)=GRHO(2)
   PR(3,1)=GRHO(3)
   NCAL = NCAL+1
!
110 NN = INT((RMIN-Pt15)*Fifty)
   IF (NN.LT.MinStp) NN = MinStp
   IF (NN.GT.(Maxstp-1)) NN = MaxStp-1
   DS = (RMIN-Pt15)/DFLOAT(NN)
   CALL DES(R1,PR,DS,NN)
   NCAL = NCAL+NN+NN
   BL = BL+DFLOAT(NN)*DS
   DO 100 I = 1,3
      R1(I,1) = R1(I,NN+1)
      PR(I,1) = PR(I,NN+1)
100 CONTINUE
   IF (RMIN.GT.Pt2) GOTO 110
   PPN(1) = CO(IXC+MINR)
   PPN(2) = CO(IYC+MINR)
   PPN(3) = CO(IZC+MINR)
   BL = BL+RMIN+DXYZ
   IFAIL = 1
   CALL NEWTON(PPN,IFAIL,IFUNC,Iwhole,NITER,1)
   PN(1) = PPN(1)
   PN(2) = PPN(2)
   PN(3) = PPN(3)
   PN(4) = RMIN
130 GR = DSQRT((PN(1)-R1(1,1))**2+(PN(2)-R1(2,1))**2+&
   &(PN(3)-R1(3,1))**2)
   IF (GR .GT. pt01) THEN
      NN = INT((GR-dxyz)*fvhund)
      IF (NN.LT.minstp) NN = minstp
      IF (NN.GT.(maxstp-1)) NN = MaxStp-1
      DS = (GR-dxyz)/DFLOAT(NN)
      CALL DES(R1,PR,DS,NN)
      NCAL = NCAL+NN+NN
      DO 120 I = 1,3
         R1(I,1) = R1(I,NN+1)
         PR(I,1) = PR(I,NN+1)
120   CONTINUE
      GOTO 130
   END IF
!
   BANGLE(4) = one
   BANGLE(5) = DACOS(BANGLE(3))
   IF (DABS(BANGLE(1)) .LT. Small) THEN
      BANGLE(6) = Pi/Two
      IF (BANGLE(2).LT.Zero) BANGLE(6) = Three*BANGLE(6)
   ELSE
      BANGLE(6) = DATAN(BANGLE(2)/BANGLE(1))
   ENDIF
   IF (BANGLE(1).LT.Zero) BANGLE(6) = BANGLE(6)+Pi
   IF(BANGLE(1) .GE. Zero .AND. BANGLE(2) .LT. Zero) BANGLE(6)=&
   &BANGLE(6)+Two*Pi
!
   RETURN
END

