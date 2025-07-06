! Generated from extreme.f90
! Contains: GAUS4, GEOM


SUBROUTINE GAUS4
!
   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
   Parameter (MaxOff=200000)
   COMMON /OFFSET/ ITYPE,ICENT,IMO,IEORB,IE,ICHARG,IXC,IYC,IZC,&
   &IXX, IYY, IZZ,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,IGXX,IGXY,IGXZ,IGYY,&
   &IGYZ,IGZZ,IGXXX,IGXXY,IGXXZ,IGXYY,IGXZZ,IGXYZ,IGYYY,IGYYZ,&
   &IGYZZ,IGZZZ,IGXXXX,IGXXXY,IGXXXZ,IGXXYY,IGXXZZ,IGXYYY,IGXZZZ,&
   &IGXXYZ,IGXYYZ,IGXYZZ,IGYYYY,IGYYYZ,IGYYZZ,IGYZZZ,IGZZZZ,ICOFMX
   COMMON CO(MaxOff),IC(MaxOff),NCENT,NMO,NPRIMS
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
   Common /What/ Iwhat
   Common /Options/ Icut, Iprint,Eps,EpsNuc,Dmp,DmpNuc
   Save Zero,Ifill,IEmpty,Ten
   Data Zero/0.0d0/,Ifill/1/,IEmpty/0/,Ten/10.0d0/
!
   IPlac=35
   Cutoff=Ten**(-Icut)
!
   DO 10 J=1,Iplac*NMO
      CO(IPSI+J)=ZERO
10 CONTINUE
!
   DO 20 I=1,NPRIMS
      K=IC(ICENT+I)
      ENT=CO(IE+I)
      X=CO(IXX+K)
      Y=CO(IYY+K)
      Z=CO(IZZ+K)
      EXPON=DEXP(-ENT*CO(IR2+K))
!
      Call PrimeF(Ent,Expon,X,Y,Z)
!
      IWhich=IC(ITYPE+I)
!
      Call PrimeA(IWhich,Ifill,X,Y,Z)
!
      Call PrimeG
!
      Call PrimeA(IWhich,IEmpty,X,Y,Z)
!
      Check=Dmax1(dabs(G0),dabs(Gx),dabs(Gy),dabs(Gz))
      If(check*CO(Icofmx+I).gt.cutoff)Then
         DO 30 J=1,NMO
            CIJ=CO(IMO+NPRIMS*(J-1)+I)
            CO(IPSI+J)=CO(IPSI+J)+CIJ*G0
            CO(IGX+J)=CO(IGX+J)+CIJ*GX
            CO(IGY+J)=CO(IGY+J)+CIJ*GY
30       CO(IGZ+J)=CO(IGZ+J)+CIJ*GZ
      Endif
      If(Iwhat.gt.0)Then
         Check=Dmax1(dabs(GXX),dabs(GXY),dabs(GXZ),dabs(GYY),&
         &dabs(GYZ),dabs(GZZ))
         If(check*CO(Icofmx+I).gt.cutoff)Then
            DO 40 J=1,NMO
               CIJ=CO(IMO+NPRIMS*(J-1)+I)
               CO(IGXX+J)=CO(IGXX+J)+CIJ*GXX
               CO(IGXY+J)=CO(IGXY+J)+CIJ*GXY
               CO(IGXZ+J)=CO(IGXZ+J)+CIJ*GXZ
               CO(IGYY+J)=CO(IGYY+J)+CIJ*GYY
               CO(IGYZ+J)=CO(IGYZ+J)+CIJ*GYZ
40          CO(IGZZ+J)=CO(IGZZ+J)+CIJ*GZZ
         Endif
      Endif
      If(Iwhat.gt.1)Then
         Check=Dmax1(dabs(GXXX),dabs(GXXY),dabs(GXXZ),dabs(GXYY),&
         &dabs(GXZZ),dabs(GXYZ),dabs(GYYY),dabs(GYYZ),dabs(GYZZ),&
         &dabs(GZZZ))
         If(check*CO(Icofmx+I).gt.cutoff)Then
            DO 50 J=1,NMO
               CIJ=CO(IMO+NPRIMS*(J-1)+I)
               CO(IGXXX+J)=CO(IGXXX+J)+CIJ*GXXX
               CO(IGXYY+J)=CO(IGXYY+J)+CIJ*GXYY
               CO(IGXZZ+J)=CO(IGXZZ+J)+CIJ*GXZZ
               CO(IGXXY+J)=CO(IGXXY+J)+CIJ*GXXY
               CO(IGXXZ+J)=CO(IGXXZ+J)+CIJ*GXXZ
               CO(IGXYZ+J)=CO(IGXYZ+J)+CIJ*GXYZ
               CO(IGYYY+J)=CO(IGYYY+J)+CIJ*GYYY
               CO(IGYYZ+J)=CO(IGYYZ+J)+CIJ*GYYZ
               CO(IGYZZ+J)=CO(IGYZZ+J)+CIJ*GYZZ
50          CO(IGZZZ+J)=CO(IGZZZ+J)+CIJ*GZZZ
         Endif
      Endif
      If(Iwhat.gt.2)Then
         Check=Dmax1(dabs(GXXXX),dabs(GXXXY),dabs(GXXXZ),&
         &dabs(GXXYY),dabs(GXXZZ),dabs(GXYYY),dabs(GXZZZ),&
         &dabs(GXYYZ),dabs(GXYZZ),dabs(GYYYY),dabs(GYYYZ),&
         &dabs(GYYZZ),dabs(GYZZZ),dabs(GZZZZ))
         If(check*CO(Icofmx+I).gt.cutoff)Then
            DO 60 J=1,NMO
               CIJ=CO(IMO+NPRIMS*(J-1)+I)
               CO(IGXXXX+J)=CO(IGXXXX+J)+CIJ*GXXXX
               CO(IGXXXY+J)=CO(IGXXXY+J)+CIJ*GXXXY
               CO(IGXXXZ+J)=CO(IGXXXZ+J)+CIJ*GXXXZ
               CO(IGXXYY+J)=CO(IGXXYY+J)+CIJ*GXXYY
               CO(IGXXZZ+J)=CO(IGXXZZ+J)+CIJ*GXXZZ
               CO(IGXYYY+J)=CO(IGXYYY+J)+CIJ*GXYYY
               CO(IGXZZZ+J)=CO(IGXZZZ+J)+CIJ*GXZZZ
               CO(IGXXYZ+J)=CO(IGXXYZ+J)+CIJ*GXXYZ
               CO(IGXYYZ+J)=CO(IGXYYZ+J)+CIJ*GXYYZ
               CO(IGXYZZ+J)=CO(IGXYZZ+J)+CIJ*GXYZZ
               CO(IGYYYY+J)=CO(IGYYYY+J)+CIJ*GYYYY
               CO(IGYYYZ+J)=CO(IGYYYZ+J)+CIJ*GYYYZ
               CO(IGYYZZ+J)=CO(IGYYZZ+J)+CIJ*GYYZZ
               CO(IGYZZZ+J)=CO(IGYZZZ+J)+CIJ*GYZZZ
60          CO(IGZZZZ+J)=CO(IGZZZZ+J)+CIJ*GZZZZ
         Endif
      Endif
20 CONTINUE
!
   RETURN
END

SUBROUTINE GEOM (N,XYZ,RN,AYZ,AXZ,AXY)
!
   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
   Parameter(MaxOff=200000)
   COMMON CO(MaxOff),IC(MaxOff),NCENT,NMO,NPRIMS
   COMMON /OFFSET/ ITYPE,ICENT,IMO,IEORB,IE,ICHARG,IXC,IYC,IZC,&
   &IXX, IYY, IZZ,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,IGXX,IGXY,IGXZ,IGYY,&
   &IGYZ,IGZZ,IGXXX,IGXXY,IGXXZ,IGXYY,IGXZZ,IGXYZ,IGYYY,IGYYZ,&
   &IGYZZ,IGZZZ,IGXXXX,IGXXXY,IGXXXZ,IGXXYY,IGXXZZ,IGXYYY,IGXZZZ,&
   &IGXXYZ,IGXYYZ,IGXYZZ,IGYYYY,IGYYYZ,IGYYZZ,IGYZZZ,IGZZZZ,ICofMx
   Dimension XYZ(3)
   Save Zero,One,Two,thsxty
   Data Zero/0.d0/,One/1.0d0/,Two/2.0d0/,Thsxty/360.0d0/
!
   Pi=Dacos(-one)
   Degree=Thsxty/(Two*pi)
!
   XN2 = (XYZ(1)-CO(IXC+N))*(XYZ(1)-CO(IXC+N))
   YN2 = (XYZ(2)-CO(IYC+N))*(XYZ(2)-CO(IYC+N))
   ZN2 = (XYZ(3)-CO(IZC+N))*(XYZ(3)-CO(IZC+N))
   RN = DSQRT(XN2+YN2+ZN2)
   IF (RN.EQ.Zero) THEN
      AXY = Zero
      AXZ = Zero
      AYZ = Zero
      RETURN
   ELSE
      RYZ = DSQRT(YN2+ZN2)
      RXZ = DSQRT(XN2+ZN2)
      RXY = DSQRT(XN2+YN2)
      AYZ = DACOS(RYZ/RN)*DEGREE
      AXZ = DACOS(RXZ/RN)*DEGREE
      AXY = DACOS(RXY/RN)*DEGREE
      RETURN
   ENDIF
END
