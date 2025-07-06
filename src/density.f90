module density
  implicit none
contains
SUBROUTINE DIVSTR(dsig,dsigm)
!
   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
   Parameter(MaxOff=200000)
   COMMON CO(MaxOff),IC(MaxOff),NCENT,NMO,NPRIMS
   COMMON /OFFSET/ ITYPE,ICENT,IMO,IEORB,IE,ICHARG,IXC,IYC,IZC,&
   &IXX, IYY, IZZ,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,IGXX,IGXY,IGXZ,IGYY,&
   &IGYZ,IGZZ,IGXXX,IGXXY,IGXXZ,IGXYY,IGXZZ,IGXYZ,IGYYY,IGYYZ,&
   &IGYZZ,IGZZZ,IGXXXX,IGXXXY,IGXXXZ,IGXXYY,IGXXZZ,IGXYYY,IGXZZZ,&
   &IGXXYZ,IGXYYZ,IGXYZZ,IGYYYY,IGYYYZ,IGYYZZ,IGYZZZ,IGZZZZ,ICOFMx
   DIMENSION dsig(3)
   Save Zero,Pt5
   Data Zero/0.d0/,Pt5/0.5d0/
!
   DO 100 I = 1,3
100 Dsig(I) = Dsig(I)
!
   DO 110 J = 1,NMO
      Temp = -(CO(IGXX+J)+CO(IGYY+J)+CO(IGZZ+J))
      DSig(1) = DSig(1)+CO(IP+J)*(CO(IPSI+J)*(CO(IGXXX+J)+&
      &CO(IGXYY+J)+CO(IGXZZ+J)) + temp*CO(IGX+J))
      DSig(2) = DSig(2)+CO(IP+J)*(CO(IPSI+J)*(CO(IGXXY+J)+&
      &CO(IGYYY+J)+CO(IGYZZ+J)) + temp*CO(IGY+J))
110 DSig(3) = DSig(3)+CO(IP+J)*(CO(IPSI+J)*(CO(IGXXZ+J)+&
   &CO(IGYYZ+J)+CO(IGZZZ+J)) + temp*CO(IGZ+J))
!
   DO 120 I = 1,3
120 Dsig(I) = pt5*Dsig(I)
!
   Dsigm = dsqrt(dsig(1)**2+dsig(2)**2+dsig(3)**2)

   RETURN
END
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
SUBROUTINE GRDD2R(Iopt,R,VALUE,GRADD2,GRADD,HDEL2)
!
   IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!
   Parameter(MaxOff=200000)
   COMMON CO(MaxOff),IC(MaxOff),NCENT,NMO,NPRIMS
   COMMON /DIST/ BANGLE(6),RMIN,RMAX,MINR
   COMMON /OFFSET/ ITYPE,ICENT,IMO,IEORB,IE,ICHARG,IXC,IYC,IZC,&
   &IXX, IYY, IZZ,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,IGXX,IGXY,IGXZ,IGYY,&
   &IGYZ,IGZZ,IGXXX,IGXXY,IGXXZ,IGXYY,IGXZZ,IGXYZ,IGYYY,IGYYZ,&
   &IGYZZ,IGZZZ,IGXXXX,IGXXXY,IGXXXZ,IGXXYY,IGXXZZ,IGXYYY,IGXZZZ,&
   &IGXXYZ,IGXYYZ,IGXYZZ,IGYYYY,IGYYYZ,IGYYZZ,IGYZZZ,IGZZZZ,ICofMx
   DIMENSION GRADD2(3),R(3),HDEL2(3,3)
   Save Zero,Small,Two,Thou
   Data Zero/0.d0/,Small/1.d-13/,Two/2.d0/,Thou/1.d3/
!
   If(Iopt.eq.1)Goto 30
!
   RMIN = Thou
   RMAX = Zero
!
   DO 10 I = 1,NCENT
      CO(IXX+I) = R(1)-CO(IXC+I)
      CO(IYY+I) = R(2)-CO(IYC+I)
      CO(IZZ+I) = R(3)-CO(IZC+I)
      CO(IR2+I) = CO(IXX+I)*CO(IXX+I) +&
      &CO(IYY+I)*CO(IYY+I) +&
      &CO(IZZ+I)*CO(IZZ+I)
      CO(IRR+I) = DSQRT(CO(IR2+I))
      IF(CO(IRR+I) .LT. RMIN) THEN
         RMIN = CO(IRR+I)
         MINR = I
         BANGLE(1) = CO(IXX+I)/(CO(IRR+I)+Small)
         BANGLE(2) = CO(IYY+I)/(CO(IRR+I)+Small)
         BANGLE(3) = CO(IZZ+I)/(CO(IRR+I)+Small)
         IF (CO(IRR+I) .GT. RMAX) THEN
            RMAX = CO(IRR+I)
         END IF
      END IF
10 CONTINUE
!
   CALL GAUS4
!
   DO 20 I = 1,3
      GRADD2(I) = Zero
      DO 20 J = 1,3
         HDEL2(I,J) = Zero
20 CONTINUE
!
30 VALUE=Zero
   DO 40 I = 1,NMO
40 VALUE=VALUE+Two*CO(IP+I)*(CO(IPSI+I)*&
   &(CO(IGXX+I)+CO(IGYY+I)+CO(IGZZ+I))+&
   &CO(IGX+I)**2+CO(IGY+I)**2+CO(IGZ+I)**2)
!
   If(Iopt.eq.1)Goto 70
!
   DO 50 I = 1,NMO
      PSI=CO(IPSI+I)
      COMM1=Two*CO(IP+I)
      COMM2=CO(IGXX+I)+CO(IGYY+I)+CO(IGZZ+I)
      COMM3=CO(IGXXX+I)+CO(IGXYY+I)+CO(IGXZZ+I)
      COMM4=CO(IGXXY+I)+CO(IGYYY+I)+CO(IGYZZ+I)
      COMM5=CO(IGXXZ+I)+CO(IGYYZ+I)+CO(IGZZZ+I)
!
      GRADD2(1) = GRADD2(1) + COMM1*(Two*(CO(IGX+I)*CO(IGXX+I)+&
      &CO(IGY+I)*CO(IGXY+I)+CO(IGZ+I)*CO(IGXZ+I))+&
      &CO(IGX+I)*COMM2+PSI*COMM3)
      GRADD2(2)=GRADD2(2)+COMM1*(Two*(CO(IGX+I)*CO(IGXY+I)+&
      &CO(IGY+I)*CO(IGYY+I)+CO(IGZ+I)*CO(IGYZ+I))+&
      &CO(IGY+I)*COMM2+PSI*COMM4)
      GRADD2(3)=GRADD2(3)+COMM1*(Two*(CO(IGX+I)*CO(IGXZ+I)+&
      &CO(IGY+I)*CO(IGYZ+I)+CO(IGZ+I)*CO(IGZZ+I))+&
      &CO(IGZ+I)*COMM2+PSI*COMM5)
!
      HDEL2(1,1)=HDEL2(1,1)+COMM1*(Two*(CO(IGXX+I)*CO(IGXX+I)+&
      &CO(IGXY+I)*CO(IGXY+I)+CO(IGXZ+I)*CO(IGXZ+I)+&
      &CO(IGX+I)*CO(IGXXX+I)+CO(IGY+I)*CO(IGXXY+I)+&
      &CO(IGZ+I)*CO(IGXXZ+I)+CO(IGX+I)*COMM3)+CO(IGXX+I)*COMM2+&
      &PSI*(CO(IGXXXX+I)+CO(IGXXYY+I)+CO(IGXXZZ+I)))
      HDEL2(1,2)=HDEL2(1,2)+COMM1*(Two*(CO(IGXY+I)*&
      &(CO(IGXX+I)+CO(IGYY+I))+CO(IGXZ+I)*CO(IGYZ+I)+&
      &CO(IGX+I)*CO(IGXXY+I)+CO(IGY+I)*CO(IGXYY+I)+&
      &CO(IGZ+I)*CO(IGXYZ+I))+CO(IGXY+I)*COMM2+&
      &CO(IGX+I)*COMM4+CO(IGY+I)*COMM3+&
      &PSI*(CO(IGXXXY+I)+CO(IGXYYY+I)+CO(IGXYZZ+I)))
      HDEL2(1,3)=HDEL2(1,3)+COMM1*(Two*(CO(IGXZ+I)*&
      &(CO(IGXX+I)+CO(IGZZ+I))+CO(IGXY+I)*CO(IGYZ+I)+&
      &CO(IGX+I)*CO(IGXXZ+I)+CO(IGY+I)*CO(IGXYZ+I)+&
      &CO(IGZ+I)*CO(IGXZZ+I))+CO(IGXZ+I)*COMM2+&
      &CO(IGX+I)*COMM5+CO(IGZ+I)*COMM3+&
      &PSI*(CO(IGXXXZ+I)+CO(IGXYYZ+I)+CO(IGXZZZ+I)))
      HDEL2(2,2)=HDEL2(2,2)+COMM1*(Two*(CO(IGYY+I)*CO(IGYY+I)+&
      &CO(IGXY+I)*CO(IGXY+I)+CO(IGYZ+I)*CO(IGYZ+I)+&
      &CO(IGX+I)*CO(IGXYY+I)+CO(IGY+I)*CO(IGYYY+I)+&
      &CO(IGZ+I)*CO(IGYYZ+I)+CO(IGY+I)*COMM4)+CO(IGYY+I)*COMM2+&
      &PSI*(CO(IGYYYY+I)+CO(IGXXYY+I)+CO(IGYYZZ+I)))
      HDEL2(2,3)=HDEL2(2,3)+COMM1*(Two*(CO(IGYZ+I)*&
      &(CO(IGYY+I)+CO(IGZZ+I))+CO(IGXY+I)*CO(IGXZ+I)+&
      &CO(IGX+I)*CO(IGXYZ+I)+CO(IGY+I)*CO(IGYYZ+I)+&
      &CO(IGZ+I)*CO(IGYZZ+I))+CO(IGYZ+I)*COMM2+&
      &CO(IGY+I)*COMM5+CO(IGZ+I)*COMM4+&
      &PSI*(CO(IGXXYZ+I)+CO(IGYYYZ+I)+CO(IGYZZZ+I)))
      HDEL2(3,3)=HDEL2(3,3)+COMM1*(Two*(CO(IGZZ+I)*CO(IGZZ+I)+&
      &CO(IGXZ+I)*CO(IGXZ+I)+CO(IGYZ+I)*CO(IGYZ+I)+&
      &CO(IGX+I)*CO(IGXZZ+I)+CO(IGY+I)*CO(IGYZZ+I)+&
      &CO(IGZ+I)*CO(IGZZZ+I)+CO(IGZ+I)*COMM5)+CO(IGZZ+I)*COMM2+&
      &PSI*(CO(IGZZZZ+I)+CO(IGXXZZ+I)+CO(IGYYZZ+I)))
!
50 CONTINUE
!
   HDEL2(2,1)=HDEL2(1,2)
   HDEL2(3,1)=HDEL2(1,3)
   HDEL2(3,2)=HDEL2(2,3)
   GRADD=Zero
   DO 60 I=1,3
60 GRADD=GRADD+GRADD2(I)*GRADD2(I)
   GRADD=DSQRT(GRADD)
!
70 Continue
!
   RETURN
END
SUBROUTINE GRDKEG(IOpt,R,VALUE,W,GRAD,H)
!
   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
   Parameter(MaxOff=200000)
   COMMON CO(MaxOff),IC(MaxOff),NCENT,NMO,NPRIMS
   COMMON /DIST/ BANGLE(6),RMIN,RMAX,MINR
   COMMON /OFFSET/ ITYPE,ICENT,IMO,IEORB,IE,ICHARG,IXC,IYC,IZC,&
   &IXX, IYY, IZZ,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,IGXX,IGXY,IGXZ,IGYY,&
   &IGYZ,IGZZ,IGXXX,IGXXY,IGXXZ,IGXYY,IGXZZ,IGXYZ,IGYYY,IGYYZ,&
   &IGYZZ,IGZZZ,IGXXXX,IGXXXY,IGXXXZ,IGXXYY,IGXXZZ,IGXYYY,IGXZZZ,&
   &IGXXYZ,IGXYYZ,IGXYZZ,IGYYYY,IGYYYZ,IGYYZZ,IGYZZZ,IGZZZZ,ICofMx
   DIMENSION W(3),R(3),H(3,3)
   Save Zero,Small,Pt5,Thou
   Data Zero/0.d0/,Small/1.d-13/,Pt5/0.5d0/,Thou/1.d3/
!
   If(Iopt.eq.1)Goto 30
!
   RMIN = Thou
   RMAX = Zero
!
   DO 10 I=1,NCENT
      CO(IXX+I) = R(1)-CO(IXC+I)
      CO(IYY+I) = R(2)-CO(IYC+I)
      CO(IZZ+I) = R(3)-CO(IZC+I)
      CO(IR2+I) = CO(IXX+I)*CO(IXX+I) +&
      &CO(IYY+I)*CO(IYY+I) +&
      &CO(IZZ+I)*CO(IZZ+I)
      CO(IRR+I) = DSQRT(CO(IR2+I))
      IF(CO(IRR+I) .LT. RMIN) THEN
         RMIN = CO(IRR+I)
         MINR = I
         BANGLE(1) = CO(IXX+I)/(CO(IRR+I)+Small)
         BANGLE(2) = CO(IYY+I)/(CO(IRR+I)+Small)
         BANGLE(3) = CO(IZZ+I)/(CO(IRR+I)+Small)
      ELSE IF (CO(IRR+I) .GT. RMAX) THEN
         RMAX = CO(IRR+I)
      END IF
10 CONTINUE
!
   CALL GAUS4
!
   DO 20 I = 1,3
      W(I) = Zero
      DO 20 J = I,3
         H(I,J) = Zero
20 CONTINUE
!
30 Value=Zero
   DO 40 I = 1,NMO
40 Value=Value+CO(IP+I)*&
   &(CO(IGX+I)**2+CO(IGY+I)**2+CO(IGZ+I)**2)
   VALUE=Pt5*VALUE
!
   If(Iopt.eq.1)Goto 70
!
   DO 50 I = 1,NMO
      P0 = CO(IP+I)
      W(1) = W(1) + P0*(CO(IGX+I)*CO(IGXX+I)+&
      &CO(IGY+I)*CO(IGXY+I)+&
      &CO(IGZ+I)*CO(IGXZ+I))
      W(2) = W(2) + P0*(CO(IGY+I)*CO(IGYY+I)+&
      &CO(IGX+I)*CO(IGXY+I)+&
      &CO(IGZ+I)*CO(IGYZ+I))
      W(3) = W(3) + P0*(CO(IGZ+I)*CO(IGZZ+I)+&
      &CO(IGY+I)*CO(IGYZ+I)+&
      &CO(IGX+I)*CO(IGXZ+I))
      H(1,1) = H(1,1) + P0*(CO(IGXX+I)**2+&
      &CO(IGX+I)*CO(IGXXX+I)+CO(IGXY+I)**2+&
      &CO(IGY+I)*CO(IGXXY+I)+&
      &CO(IGXZ+I)**2+CO(IGZ+I)*CO(IGXXZ+I))
      H(1,2) = H(1,2) + P0*(CO(IGXX+I)*CO(IGXY+I)+&
      &CO(IGX+I)*CO(IGXXY+I)+CO(IGXY+I)*CO(IGYY+I)+&
      &CO(IGY+I)*CO(IGXYY+I)+&
      &CO(IGXZ+I)*CO(IGYZ+I)+CO(IGZ+I)*CO(IGXYZ+I))
      H(1,3) = H(1,3) + P0*(CO(IGXX+I)*CO(IGXZ+I)+&
      &CO(IGX+I)*CO(IGXXZ+I)+CO(IGXY+I)*CO(IGYZ+I)+&
      &CO(IGY+I)*CO(IGXYZ+I)+&
      &CO(IGXZ+I)*CO(IGZZ+I)+CO(IGZ+I)*CO(IGXZZ+I))
      H(2,2) = H(2,2) + P0*(CO(IGYY+I)**2+&
      &CO(IGY+I)*CO(IGYYY+I)+CO(IGXY+I)**2+&
      &CO(IGX+I)*CO(IGXYY+I)+&
      &CO(IGYZ+I)**2+CO(IGZ+I)*CO(IGYYZ+I))
      H(2,3) = H(2,3) + P0*(CO(IGXY+I)*CO(IGXZ+I)+&
      &CO(IGX+I)*CO(IGXYZ+I)+CO(IGYY+I)*CO(IGYZ+I)+&
      &CO(IGY+I)*CO(IGYYZ+I)+&
      &CO(IGYZ+I)*CO(IGZZ+I)+CO(IGZ+I)*CO(IGYZZ+I))
      H(3,3) = H(3,3) + P0*(CO(IGZZ+I)**2+&
      &CO(IGZ+I)*CO(IGZZZ+I)+CO(IGYZ+I)**2+&
      &CO(IGY+I)*CO(IGYZZ+I)+&
      &CO(IGXZ+I)**2+CO(IGX+I)*CO(IGXZZ+I))
50 CONTINUE
!
   GRAD=Zero
!
   DO 60 I = 1,3
      GRAD = GRAD+W(I)*W(I)
      DO 60 J = I,3
         H(J,I) = H(I,J)
60 CONTINUE
!
   GRAD = DSQRT(GRAD)
!
70 Continue
!
   RETURN
END
SUBROUTINE GRDKEK(Iopt,R,VALUE,W,GRADK,H)
!
   IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!
   Parameter(MaxOff=200000)
   COMMON CO(MaxOff),IC(MaxOff),NCENT,NMO,NPRIMS
   COMMON /DIST/ BANGLE(6),RMIN,RMAX,MINR
   COMMON /OFFSET/ ITYPE,ICENT,IMO,IEORB,IE,ICHARG,IXC,IYC,IZC,&
   &IXX, IYY, IZZ,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,IGXX,IGXY,IGXZ,IGYY,&
   &IGYZ,IGZZ,IGXXX,IGXXY,IGXXZ,IGXYY,IGXZZ,IGXYZ,IGYYY,IGYYZ,&
   &IGYZZ,IGZZZ,IGXXXX,IGXXXY,IGXXXZ,IGXXYY,IGXXZZ,IGXYYY,IGXZZZ,&
   &IGXXYZ,IGXYYZ,IGXYZZ,IGYYYY,IGYYYZ,IGYYZZ,IGYZZZ,IGZZZZ,ICofMx
   DIMENSION RS(3),W(3),R(3),H(3,3)
   Save Zero,Small,Pt5,Two,Thou
   Data Zero/0.d0/,Small/1.d-13/,Pt5/0.5d0/,Two/2.d0/,Thou/1.d3/
!
   If(Iopt.eq.1)Goto 30
!
   RMIN = Thou
   RMAX = Zero
!
   DO 10 I = 1,NCENT
      CO(IXX+I) = R(1)-CO(IXC+I)
      CO(IYY+I) = R(2)-CO(IYC+I)
      CO(IZZ+I) = R(3)-CO(IZC+I)
      CO(IR2+I) = CO(IXX+I)*CO(IXX+I) +&
      &CO(IYY+I)*CO(IYY+I) +&
      &CO(IZZ+I)*CO(IZZ+I)
      CO(IRR+I) = DSQRT(CO(IR2+I))
      IF(CO(IRR+I) .LT. RMIN) THEN
         RMIN = CO(IRR+I)
         MINR = I
         BANGLE(1) = CO(IXX+I)/(CO(IRR+I)+Small)
         BANGLE(2) = CO(IYY+I)/(CO(IRR+I)+Small)
         BANGLE(3) = CO(IZZ+I)/(CO(IRR+I)+Small)
      ELSE IF (CO(IRR+I) .GT. RMAX) THEN
         RMAX = CO(IRR+I)
      END IF
10 CONTINUE
!
   CALL GAUS4
!
   DO 20 I = 1,3
      W(I) = Zero
      DO 20 J = 1,3
         H(I,J) = Zero
20 CONTINUE
!
30 Value=Zero
   DO 40 I = 1,NMO
40 Value=Value+CO(IP+I)*CO(IPSI+I)*(CO(IGXX+I)+&
   &CO(IGYY+I)+CO(IGZZ+I))
   Value=-Pt5*Value
!
   If(Iopt.eq.1)Goto 70
!
   DO 50 I = 1,NMO
      PSI=CO(IPSI+I)
      COMM1=-CO(IP+I)/Two
      COMM2=CO(IGXX+I)+CO(IGYY+I)+CO(IGZZ+I)
      COMM3=CO(IGXXX+I)+CO(IGXYY+I)+CO(IGXZZ+I)
      COMM4=CO(IGXXY+I)+CO(IGYYY+I)+CO(IGYZZ+I)
      COMM5=CO(IGXXZ+I)+CO(IGYYZ+I)+CO(IGZZZ+I)
!
      W(1) = W(1) + COMM1*(CO(IGX+I)*COMM2+PSI*COMM3)
      W(2) = W(2) + COMM1*(CO(IGY+I)*COMM2+PSI*COMM4)
      W(3) = W(3) + COMM1*(CO(IGZ+I)*COMM2+PSI*COMM5)
!
      H(1,1)=H(1,1)+COMM1*(TWO*CO(IGX+I)*COMM3+CO(IGXX+I)*COMM2+&
      &PSI*(CO(IGXXXX+I)+CO(IGXXYY+I)+CO(IGXXZZ+I)))
!
      H(1,2)=H(1,2)+COMM1*(CO(IGXY+I)*COMM2+CO(IGX+I)*COMM4+&
      &CO(IGY+I)*COMM3+PSI*(CO(IGXXXY+I)+CO(IGXYYY+I)+CO(IGXYZZ+I)))
!
      H(1,3)=H(1,3)+COMM1*(CO(IGXZ+I)*COMM2+CO(IGX+I)*COMM5+&
      &CO(IGZ+I)*COMM3+PSI*(CO(IGXXXZ+I)+CO(IGXYYZ+I)+CO(IGXZZZ+I)))
!
      H(2,2)=H(2,2)+COMM1*(TWO*CO(IGY+I)*COMM4+CO(IGYY+I)*COMM2+&
      &PSI*(CO(IGXXYY+I)+CO(IGYYYY+I)+CO(IGYYZZ+I)))
!
      H(2,3)=H(2,3)+COMM1*(CO(IGYZ+I)*COMM2+CO(IGY+I)*COMM5+&
      &CO(IGZ+I)*COMM4+PSI*(CO(IGXXYZ+I)+CO(IGYYYZ+I)+CO(IGYZZZ+I)))
!
      H(3,3)=H(3,3)+COMM1*(TWO*CO(IGZ+I)*COMM5+CO(IGZZ+I)*COMM2+&
      &PSI*(CO(IGXXZZ+I)+CO(IGYYZZ+I)+CO(IGZZZZ+I)))
!
50 CONTINUE
!
   H(2,1)=H(1,2)
   H(3,1)=H(1,3)
   H(3,2)=H(2,3)
   GRADK=Zero
   DO 60 I=1,3
60 GRADK=GRADK+W(I)*W(I)
   GRADK=DSQRT(GRADK)
!
70 Continue
!
   RETURN
END
SUBROUTINE GRDRHO(Iopt,R,VALUE,W,GRAD,H,SG)
!
   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
   Parameter(MaxOff=200000)
   COMMON CO(MaxOff),IC(MaxOff),NCENT,NMO,NPRIMS
   COMMON /DIST/ BANGLE(6),RMIN,RMAX,MINR
   COMMON /OFFSET/ ITYPE,ICENT,IMO,IEORB,IE,ICHARG,IXC,IYC,IZC,&
   &IXX, IYY, IZZ,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,IGXX,IGXY,IGXZ,IGYY,&
   &IGYZ,IGZZ,IGXXX,IGXXY,IGXXZ,IGXYY,IGXZZ,IGXYZ,IGYYY,IGYYZ,&
   &IGYZZ,IGZZZ,IGXXXX,IGXXXY,IGXXXZ,IGXXYY,IGXXZZ,IGXYYY,IGXZZZ,&
   &IGXXYZ,IGXYYZ,IGXYZZ,IGYYYY,IGYYYZ,IGYYZZ,IGYZZZ,IGZZZZ,ICofMx
   DIMENSION W(3),R(3),H(3,3),SG(3,3)
   Save Zero,Small,pt5,Two,Thou
   Data Zero/0.d0/,Small/1.d-13/,pt5/0.5d0/,Two/2.d0/,Thou/1.d3/
!
   If(Iopt.eq.2)Goto 30
!
   RMIN = Thou
   RMAX = Zero
!
   DO 10 I=1,NCENT
      CO(IXX+I) = R(1)-CO(IXC+I)
      CO(IYY+I) = R(2)-CO(IYC+I)
      CO(IZZ+I) = R(3)-CO(IZC+I)
      CO(IR2+I) = CO(IXX+I)*CO(IXX+I) +&
      &CO(IYY+I)*CO(IYY+I) +&
      &CO(IZZ+I)*CO(IZZ+I)
      CO(IRR+I) = DSQRT(CO(IR2+I))
      IF(CO(IRR+I) .LT. RMIN) THEN
         RMIN = CO(IRR+I)
         MINR = I
         BANGLE(1) = CO(IXX+I)/(CO(IRR+I)+Small)
         BANGLE(2) = CO(IYY+I)/(CO(IRR+I)+Small)
         BANGLE(3) = CO(IZZ+I)/(CO(IRR+I)+Small)
      ELSE IF (CO(IRR+I) .GT. RMAX) THEN
         RMAX = CO(IRR+I)
      END IF
10 CONTINUE
!
   CALL GAUS4
!
   DO 20 I = 1,3
      W(I) = Zero
      DO 20 J = I,3
         H(I,J) = Zero
         SG(I,J) = Zero
20 CONTINUE
!
30 Value=Zero
   DO 40 I = 1,NMO
      P0 = CO(IP+I)
      P1 = P0*CO(IPSI+I)
      Value=Value+P1*CO(IPSI+I)
      W(1) = W(1) + P1*CO(IGX+I)
      W(2) = W(2) + P1*CO(IGY+I)
40 W(3) = W(3) + P1*CO(IGZ+I)
!
   If(Iopt.eq.2)Goto 70
!
   If(IOpt.eq.1)Then
      DO 50 I = 1,NMO
         P0 = CO(IP+I)
         P1 = P0*CO(IPSI+I)
         H(1,1) = H(1,1) + P0*CO(IGX+I)**2
         H(1,2) = H(1,2) + P0*CO(IGX+I)*CO(IGY+I)
         H(1,3) = H(1,3) + P0*CO(IGX+I)*CO(IGZ+I)
         H(2,2) = H(2,2) + P0*CO(IGY+I)**2
         H(2,3) = H(2,3) + P0*CO(IGY+I)*CO(IGZ+I)
         H(3,3) = H(3,3) + P0*CO(IGZ+I)**2
         SG(1,1) = SG(1,1) + P1*CO(IGXX+I)
         SG(1,2) = SG(1,2) + P1*CO(IGXY+I)
         SG(1,3) = SG(1,3) + P1*CO(IGXZ+I)
         SG(2,2) = SG(2,2) + P1*CO(IGYY+I)
         SG(2,3) = SG(2,3) + P1*CO(IGYZ+I)
50    SG(3,3) = SG(3,3) + P1*CO(IGZZ+I)
   Endif
!
   Do 60 I=1,3
      DO 60 J = I,3
         DM = H(I,J)
         H(I,J) = (DM+SG(I,J))*Two
         SG(I,J) = (-DM+SG(I,J))*pt5
         H(J,I) = H(I,J)
         SG(J,I) = SG(I,J)
60 CONTINUE
!
70 GRAD=Zero
   DO 80  I = 1,3
      W(I) = W(I)*Two
80 GRAD = GRAD+W(I)*W(I)
   GRAD = DSQRT(GRAD)
!
   RETURN
END
SUBROUTINE GRDV(Iopt,R,VALUE,GRADD2,GRADD,HDEL2)
!
   IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!
   Parameter(MaxOff=200000)
   COMMON CO(MaxOff),IC(MaxOff),NCENT,NMO,NPRIMS
   COMMON /DIST/ BANGLE(6),RMIN,RMAX,MINR
   COMMON /OFFSET/ ITYPE,ICENT,IMO,IEORB,IE,ICHARG,IXC,IYC,IZC,&
   &IXX, IYY, IZZ,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,IGXX,IGXY,IGXZ,IGYY,&
   &IGYZ,IGZZ,IGXXX,IGXXY,IGXXZ,IGXYY,IGXZZ,IGXYZ,IGYYY,IGYYZ,&
   &IGYZZ,IGZZZ,IGXXXX,IGXXXY,IGXXXZ,IGXXYY,IGXXZZ,IGXYYY,IGXZZZ,&
   &IGXXYZ,IGXYYZ,IGXYZZ,IGYYYY,IGYYYZ,IGYYZZ,IGYZZZ,IGZZZZ,ICofMx
   DIMENSION GRADD2(3),R(3),HDEL2(3,3)
   Save Zero,Small,Pt5,Two,Thou
   Data Zero/0.d0/,Small/1.d-13/,Pt5/0.5d0/,Two/2.d0/,Thou/1.d3/
!
   If(Iopt.eq.1)Goto 30
!
   RMIN = Thou
   RMAX = Zero
!
   DO 10 I = 1,NCENT
      CO(IXX+I) = R(1)-CO(IXC+I)
      CO(IYY+I) = R(2)-CO(IYC+I)
      CO(IZZ+I) = R(3)-CO(IZC+I)
      CO(IR2+I) = CO(IXX+I)*CO(IXX+I) +&
      &CO(IYY+I)*CO(IYY+I) +&
      &CO(IZZ+I)*CO(IZZ+I)
      CO(IRR+I) = DSQRT(CO(IR2+I))
      IF(CO(IRR+I) .LT. RMIN) THEN
         RMIN = CO(IRR+I)
         MINR = I
         BANGLE(1) = CO(IXX+I)/(CO(IRR+I)+Small)
         BANGLE(2) = CO(IYY+I)/(CO(IRR+I)+Small)
         BANGLE(3) = CO(IZZ+I)/(CO(IRR+I)+Small)
         IF (CO(IRR+I) .GT. RMAX) THEN
            RMAX = CO(IRR+I)
         END IF
      END IF
10 CONTINUE
!
   CALL GAUS4
!
   DO 20 I = 1,3
      GRADD2(I) = Zero
      DO 20 J = 1,3
         HDEL2(I,J) = Zero
20 CONTINUE
!
30 VALUE=Zero
   DO 40 I = 1,NMO
40 VALUE=VALUE+Pt5*CO(IP+I)*(CO(IPSI+I)*&
   &(CO(IGXX+I)+CO(IGYY+I)+CO(IGZZ+I))-&
   &(CO(IGX+I)**2+CO(IGY+I)**2+CO(IGZ+I)**2))
!
   If(Iopt.eq.1)Goto 70
!
   DO 50 I = 1,NMO
      PSI=CO(IPSI+I)
      COMM1=Pt5*CO(IP+I)
      COMM2=CO(IGXX+I)+CO(IGYY+I)+CO(IGZZ+I)
      COMM3=CO(IGXXX+I)+CO(IGXYY+I)+CO(IGXZZ+I)
      COMM4=CO(IGXXY+I)+CO(IGYYY+I)+CO(IGYZZ+I)
      COMM5=CO(IGXXZ+I)+CO(IGYYZ+I)+CO(IGZZZ+I)
!
      GRADD2(1) = GRADD2(1) + COMM1*(-Two*(CO(IGX+I)*CO(IGXX+I)+&
      &CO(IGY+I)*CO(IGXY+I)+CO(IGZ+I)*CO(IGXZ+I))+&
      &CO(IGX+I)*COMM2+PSI*COMM3)
      GRADD2(2)=GRADD2(2)+COMM1*(-Two*(CO(IGX+I)*CO(IGXY+I)+&
      &CO(IGY+I)*CO(IGYY+I)+CO(IGZ+I)*CO(IGYZ+I))+&
      &CO(IGY+I)*COMM2+PSI*COMM4)
      GRADD2(3)=GRADD2(3)+COMM1*(-Two*(CO(IGX+I)*CO(IGXZ+I)+&
      &CO(IGY+I)*CO(IGYZ+I)+CO(IGZ+I)*CO(IGZZ+I))+&
      &CO(IGZ+I)*COMM2+PSI*COMM5)
!
      HDEL2(1,1)=HDEL2(1,1)+COMM1*(-Two*(CO(IGXX+I)*CO(IGXX+I)+&
      &CO(IGXY+I)*CO(IGXY+I)+CO(IGXZ+I)*CO(IGXZ+I)+&
      &CO(IGX+I)*CO(IGXXX+I)+CO(IGY+I)*CO(IGXXY+I)+&
      &CO(IGZ+I)*CO(IGXXZ+I)-CO(IGX+I)*COMM3)+CO(IGXX+I)*COMM2+&
      &PSI*(CO(IGXXXX+I)+CO(IGXXYY+I)+CO(IGXXZZ+I)))
      HDEL2(1,2)=HDEL2(1,2)+COMM1*(-Two*(CO(IGXY+I)*&
      &(CO(IGXX+I)+CO(IGYY+I))+CO(IGXZ+I)*CO(IGYZ+I)+&
      &CO(IGX+I)*CO(IGXXY+I)+CO(IGY+I)*CO(IGXYY+I)+&
      &CO(IGZ+I)*CO(IGXYZ+I))+CO(IGXY+I)*COMM2+&
      &CO(IGX+I)*COMM4+CO(IGY+I)*COMM3+&
      &PSI*(CO(IGXXXY+I)+CO(IGXYYY+I)+CO(IGXYZZ+I)))
      HDEL2(1,3)=HDEL2(1,3)+COMM1*(-Two*(CO(IGXZ+I)*&
      &(CO(IGXX+I)+CO(IGZZ+I))+CO(IGXY+I)*CO(IGYZ+I)+&
      &CO(IGX+I)*CO(IGXXZ+I)+CO(IGY+I)*CO(IGXYZ+I)+&
      &CO(IGZ+I)*CO(IGXZZ+I))+CO(IGXZ+I)*COMM2+&
      &CO(IGX+I)*COMM5+CO(IGZ+I)*COMM3+&
      &PSI*(CO(IGXXXZ+I)+CO(IGXYYZ+I)+CO(IGXZZZ+I)))
      HDEL2(2,2)=HDEL2(2,2)+COMM1*(-Two*(CO(IGYY+I)*CO(IGYY+I)+&
      &CO(IGXY+I)*CO(IGXY+I)+CO(IGYZ+I)*CO(IGYZ+I)+&
      &CO(IGX+I)*CO(IGXYY+I)+CO(IGY+I)*CO(IGYYY+I)+&
      &CO(IGZ+I)*CO(IGYYZ+I)-CO(IGY+I)*COMM4)+CO(IGYY+I)*COMM2+&
      &PSI*(CO(IGYYYY+I)+CO(IGXXYY+I)+CO(IGYYZZ+I)))
      HDEL2(2,3)=HDEL2(2,3)+COMM1*(-Two*(CO(IGYZ+I)*&
      &(CO(IGYY+I)+CO(IGZZ+I))+CO(IGXY+I)*CO(IGXZ+I)+&
      &CO(IGX+I)*CO(IGXYZ+I)+CO(IGY+I)*CO(IGYYZ+I)+&
      &CO(IGZ+I)*CO(IGYZZ+I))+CO(IGYZ+I)*COMM2+&
      &CO(IGY+I)*COMM5+CO(IGZ+I)*COMM4+&
      &PSI*(CO(IGXXYZ+I)+CO(IGYYYZ+I)+CO(IGYZZZ+I)))
      HDEL2(3,3)=HDEL2(3,3)+COMM1*(-Two*(CO(IGZZ+I)*CO(IGZZ+I)+&
      &CO(IGXZ+I)*CO(IGXZ+I)+CO(IGYZ+I)*CO(IGYZ+I)+&
      &CO(IGX+I)*CO(IGXZZ+I)+CO(IGY+I)*CO(IGYZZ+I)+&
      &CO(IGZ+I)*CO(IGZZZ+I)-CO(IGZ+I)*COMM5)+CO(IGZZ+I)*COMM2+&
      &PSI*(CO(IGZZZZ+I)+CO(IGXXZZ+I)+CO(IGYYZZ+I)))
!
50 CONTINUE
!
   HDEL2(2,1)=HDEL2(1,2)
   HDEL2(3,1)=HDEL2(1,3)
   HDEL2(3,2)=HDEL2(2,3)
   GRADD=Zero
   DO 60 I=1,3
60 GRADD=GRADD+GRADD2(I)*GRADD2(I)
   GRADD=DSQRT(GRADD)
!
70 Continue
!
   RETURN
END
SUBROUTINE GRDVNE(Iopt,R,VALUE,W,GRAD,H)
!
   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
   Parameter(MaxOff=200000)
   COMMON CO(MaxOff),IC(MaxOff),NCENT,NMO,NPRIMS
   COMMON /DIST/ BANGLE(6),RMIN,RMAX,MINR
   COMMON /OFFSET/ ITYPE,ICENT,IMO,IEORB,IE,ICHARG,IXC,IYC,IZC,&
   &IXX, IYY, IZZ,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,IGXX,IGXY,IGXZ,IGYY,&
   &IGYZ,IGZZ,IGXXX,IGXXY,IGXXZ,IGXYY,IGXZZ,IGXYZ,IGYYY,IGYYZ,&
   &IGYZZ,IGZZZ,IGXXXX,IGXXXY,IGXXXZ,IGXXYY,IGXXZZ,IGXYYY,IGXZZZ,&
   &IGXXYZ,IGXYYZ,IGXYZZ,IGYYYY,IGYYYZ,IGYYZZ,IGYZZZ,IGZZZZ,ICofMx
   DIMENSION W(3),R(3),H(3,3)
   Save Zero,Small,One,Three,Thou
   Data Zero/0.d0/,Small/1.d-13/,One/1.d0/,Three/3.d0/,Thou/1.d3/
!
   If(Iopt.eq.1)Goto 30
!
   RMIN = Thou
   RMAX = Zero
!
   DO 10 I=1,NCENT
      CO(IXX+I) = R(1)-CO(IXC+I)
      CO(IYY+I) = R(2)-CO(IYC+I)
      CO(IZZ+I) = R(3)-CO(IZC+I)
      CO(IR2+I) = CO(IXX+I)*CO(IXX+I) +&
      &CO(IYY+I)*CO(IYY+I) +&
      &CO(IZZ+I)*CO(IZZ+I)
      CO(IRR+I) = DSQRT(CO(IR2+I))
      IF(CO(IRR+I) .LT. RMIN) THEN
         RMIN = CO(IRR+I)
         MINR = I
         BANGLE(1) = CO(IXX+I)/(CO(IRR+I)+Small)
         BANGLE(2) = CO(IYY+I)/(CO(IRR+I)+Small)
         BANGLE(3) = CO(IZZ+I)/(CO(IRR+I)+Small)
      ELSE IF (CO(IRR+I) .GT. RMAX) THEN
         RMAX = CO(IRR+I)
      END IF
10 CONTINUE
!
   DO 20 I = 1,3
      W(I) = Zero
      DO 20 J = I,3
         H(I,J) = Zero
20 CONTINUE
!
30 Value=Zero
   DO 40 I=1,NCENT
      OVERR=One/CO(IRR+I)
40 Value=Value+CO(ICHARG+I)*OVERR
!
   If(Iopt.eq.1)Goto 70
!
   DO 50 I=1,NCENT
      OVERR=One/CO(IRR+I)
      OVERR3=OVERR**3
      OVERR5=OVERR**5
      W(1) = W(1) - CO(ICHARG+I)*CO(IXX+I)*OVERR3
      W(2) = W(2) - CO(ICHARG+I)*CO(IYY+I)*OVERR3
      W(3) = W(3) - CO(ICHARG+I)*CO(IZZ+I)*OVERR3
      H(1,1) = H(1,1) - CO(ICHARG+I)*(OVERR3-&
      &Three*CO(IXX+I)**2*OVERR5)
      H(2,2) = H(2,2) - CO(ICHARG+I)*(OVERR3-&
      &Three*CO(IYY+I)**2*OVERR5)
      H(3,3) = H(3,3) - CO(ICHARG+I)*(OVERR3-&
      &Three*CO(IZZ+I)**2*OVERR5)
      H(1,2) = H(1,2) + Three*CO(ICHARG+I)*CO(IXX+I)*CO(IYY+I)*&
      &OVERR5
      H(1,3) = H(1,3) + Three*CO(ICHARG+I)*CO(IXX+I)*CO(IZZ+I)*&
      &OVERR5
50 H(2,3) = H(2,3) + Three*CO(ICHARG+I)*CO(IYY+I)*CO(IZZ+I)*&
   &OVERR5
!
   GRAD=Zero
!
   DO 60 I = 1,3
      GRAD = GRAD+W(I)*W(I)
      DO 60 J = I,3
         H(J,I) = H(I,J)
60 CONTINUE
!
   GRAD = DSQRT(GRAD)
!
70 Continue
!
   RETURN
END
INTEGER FUNCTION IDAMAX(N,DX,INCX)
   DOUBLE PRECISION DX(1),DMAX
   INTEGER I,INCX,IX,N
   IDAMAX = 0
   IF( N .LT. 1 ) RETURN
   IDAMAX = 1
   IF(N.EQ.1)RETURN
   IF(INCX.EQ.1)GO TO 20
   IX = 1
   DMAX = DABS(DX(1))
   IX = IX + INCX
   DO 10 I = 2,N
      IF(DABS(DX(IX)).LE.DMAX) GO TO 5
      IDAMAX = I
      DMAX = DABS(DX(IX))
5     IX = IX + INCX
10 CONTINUE
   RETURN
20 DMAX = DABS(DX(1))
   DO 30 I = 2,N
      IF(DABS(DX(I)).LE.DMAX) GO TO 30
      IDAMAX = I
      DMAX = DABS(DX(I))
30 CONTINUE
   RETURN
END
Subroutine PrimeA(IWhich,IMode,X,Y,Z)
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
   Save Zero,One,Two,Three,Six
   Data Zero/0.0d0/,One/1.0d0/,Two/2.0d0/,Three/3.0d0/,Six/6.0d0/
!
   If(IWhich.eq.1)Then
      If(IMode.eq.1)Then
         A0=One
      ElseIf(Imode.eq.0)Then
         A0=Zero
      Endif
!
   ElseIf(Iwhich.eq.2)Then
      If(IMode.eq.1)Then
         A0=X
         Ax=One
      ElseIf(Imode.eq.0)Then
         A0=Zero
         Ax=Zero
      Endif
!
   ElseIf(Iwhich.eq.3)Then
      If(IMode.eq.1)Then
         A0=Y
         Ay=One
      ElseIf(Imode.eq.0)Then
         A0=Zero
         Ay=Zero
      Endif
!
   ElseIf(Iwhich.eq.4)Then
      If(IMode.eq.1)Then
         A0=Z
         AZ=One
      ElseIf(Imode.eq.0)Then
         A0=Zero
         AZ=Zero
      Endif
!
   ElseIf(Iwhich.eq.5)Then
      If(IMode.eq.1)Then
         A0=X**2
         AX=Two*X
         AXX=Two
      ElseIf(Imode.eq.0)Then
         A0=Zero
         AX=Zero
         AXX=Zero
      Endif
!
   ElseIf(Iwhich.eq.6)Then
      If(IMode.eq.1)Then
         A0=Y**2
         AY=Two*Y
         AYY=Two
      ElseIf(Imode.eq.0)Then
         A0=Zero
         AY=Zero
         AYY=Zero
      Endif
!
   ElseIf(Iwhich.eq.7)Then
      If(IMode.eq.1)Then
         A0=Z**2
         AZ=Two*Z
         AZZ=Two
      ElseIf(Imode.eq.0)Then
         A0=Zero
         AZ=Zero
         AZZ=Zero
      Endif
!
   ElseIf(Iwhich.eq.8)Then
      If(IMode.eq.1)Then
         A0=X*Y
         AX=Y
         AY=X
         AXY=One
      ElseIf(Imode.eq.0)Then
         A0=Zero
         AX=Zero
         AY=Zero
         AXY=Zero
      Endif
!
   ElseIf(Iwhich.eq.9)Then
      If(IMode.eq.1)Then
         A0=X*Z
         AX=Z
         AZ=X
         AXZ=One
      ElseIf(Imode.eq.0)Then
         A0=Zero
         AX=Zero
         AZ=Zero
         AXZ=Zero
      Endif
!
   ElseIf(Iwhich.eq.10)Then
      If(IMode.eq.1)Then
         A0=Y*Z
         AY=Z
         AZ=Y
         AYZ=One
      ElseIf(Imode.eq.0)Then
         A0=Zero
         AY=Zero
         AZ=Zero
         AYZ=Zero
      Endif
!
   ElseIf(Iwhich.eq.11)Then
      If(IMode.eq.1)Then
         A0=X**3
         AX=Three*X**2
         AXX=Six*X
         AXXX=Six
      ElseIf(Imode.eq.0)Then
         A0=Zero
         AX=Zero
         AXX=Zero
         AXXX=Zero
      Endif
!
   ElseIf(Iwhich.eq.12)Then
      If(IMode.eq.1)Then
         A0=Y**3
         AY=Three*Y**2
         AYY=Six*Y
         AYYY=Six
      ElseIf(Imode.eq.0)Then
         A0=Zero
         AY=Zero
         AYY=Zero
         AYYY=Zero
      Endif
!
   ElseIf(Iwhich.eq.13)Then
      If(IMode.eq.1)Then
         A0=Z**3
         AZ=Three*Z**2
         AZZ=Six*Z
         AZZZ=Six
      ElseIf(Imode.eq.0)Then
         A0=Zero
         AZ=Zero
         AZZ=Zero
         AZZZ=Zero
      Endif
!
   ElseIf(Iwhich.eq.14)Then
      If(IMode.eq.1)Then
         A0=X**2*Y
         AX=Two*X*Y
         AY=X**2
         AXX=Two*Y
         AXY=Two*X
         AXXY=Two
      ElseIf(Imode.eq.0)Then
         A0=Zero
         AX=Zero
         AY=Zero
         AXX=Zero
         AXY=Zero
         AXXY=Zero
      Endif
!
   ElseIf(Iwhich.eq.15)Then
      If(IMode.eq.1)Then
         A0=X**2*Z
         AX=Two*X*Z
         AZ=X**2
         AXX=Two*Z
         AXZ=Two*X
         AXXZ=Two
      ElseIf(Imode.eq.0)Then
         A0=Zero
         AX=Zero
         AZ=Zero
         AXX=Zero
         AXZ=Zero
         AXXZ=Zero
      Endif
!
   ElseIf(Iwhich.eq.16)Then
      If(IMode.eq.1)Then
         A0=Y**2*Z
         AY=Two*Y*Z
         AZ=Y**2
         AYY=Two*Z
         AYZ=Two*Y
         AYYZ=Two
      ElseIf(Imode.eq.0)Then
         A0=Zero
         AY=Zero
         AZ=Zero
         AYY=Zero
         AYZ=Zero
         AYYZ=Zero
      Endif
!
   ElseIf(Iwhich.eq.17)Then
      If(IMode.eq.1)Then
         A0=X*Y**2
         AX=Y**2
         AY=Two*Y*X
         AYY=Two*X
         AXY=Two*Y
         AXYY=Two
      ElseIf(Imode.eq.0)Then
         A0=Zero
         AX=Zero
         AY=Zero
         AYY=Zero
         AXY=Zero
         AXYY=Zero
      Endif
!
   ElseIf(Iwhich.eq.18)Then
      If(IMode.eq.1)Then
         A0=X*Z**2
         AX=Z**2
         AZ=Two*Z*X
         AZZ=Two*X
         AXZ=Two*Z
         AXZZ=Two
      ElseIf(Imode.eq.0)Then
         A0=Zero
         AX=Zero
         AZ=Zero
         AZZ=Zero
         AXZ=Zero
         AXZZ=Zero
      Endif
!
   ElseIf(Iwhich.eq.19)Then
      If(IMode.eq.1)Then
         A0=Y*Z**2
         AY=Z**2
         AZ=Two*Z*Y
         AZZ=Two*Y
         AYZ=Two*Z
         AYZZ=Two
      ElseIf(Imode.eq.0)Then
         A0=Zero
         AY=Zero
         AZ=Zero
         AZZ=Zero
         AYZ=Zero
         AYZZ=Zero
      Endif
!
   ElseIf(Iwhich.eq.20)Then
      If(IMode.eq.1)Then
         A0=X*Y*Z
         AX=Y*Z
         AY=X*Z
         AZ=X*Y
         AXY=Z
         AXZ=Y
         AYZ=X
         AXYZ=One
      ElseIf(Imode.eq.0)Then
         A0=Zero
         AX=Zero
         AY=Zero
         AZ=Zero
         AXY=Zero
         AXZ=Zero
         AYZ=Zero
         AXYZ=Zero
      Endif
   Endif
!
   Return
End
Subroutine PrimeF(Ent,EXPON,X,Y,Z)
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
   Common /What/ IWhat
!
   Save Two,Three
   Data Two/2.0d0/,Three/3.d0/
!
   toalp=-Two*Ent
   toalpe=toalp*expon
!
   F0=EXPON
   FX=toalpe*X
   FY=toalpe*Y
   FZ=toalpe*Z
   If(Iwhat.gt.0)Then
      FXX=toalp*X*FX+toalpe
      FXY=toalp*Y*FX
      FXZ=Toalp*Z*FX
      FYY=toalp*Y*FY+toalpe
      FYZ=Toalp*Z*FY
      FZZ=toalp*Z*FZ+toalpe
   Endif
   If(Iwhat.gt.1)Then
      FXXX=Toalp*(Two*FX+X*FXX)
      FYYY=Toalp*(Two*FY+Y*FYY)
      FZZZ=Toalp*(Two*FZ+Z*FZZ)
      FXXY=toalp*Y*FXX
      FXXZ=toalp*Z*FXX
      FXYY=toalp*X*FYY
      FXZZ=toalp*X*FZZ
      FXYZ=toalp*Z*FXY
      FYYZ=toalp*Z*FYY
      FYZZ=toalp*Y*FZZ
   Endif
   If(Iwhat.gt.2)Then
      FXXXX=Toalp*(X*FXXX+Three*FXX)
      FXXXY=FXXX*Toalp*Y
      FXXXZ=FXXX*Toalp*Z
      FXXYY=Toalp*(FXX+Y*FXXY)
      FXXZZ=Toalp*(FXX+Z*FXXZ)
      FXYYY=FYYY*ToAlp*X
      FXZZZ=FZZZ*ToAlp*X
      FXXYZ=FXXY*ToAlp*Z
      FXYYZ=FXYY*ToAlp*Z
      FXYZZ=FXZZ*ToAlp*Y
      FYYYY=Toalp*(Y*FYYY+Three*FYY)
      FYYYZ=FYYY*ToAlp*Z
      FYYZZ=ToAlp*(FYY+Z*FYYZ)
      FYZZZ=FZZZ*ToAlp*Y
      FZZZZ=Toalp*(Z*FZZZ+Three*FZZ)
   Endif
!
   Return
End
Subroutine PrimeG
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
   Common /What/ IWhat
   Save Two,Three,Four,Six
   Data Two/2.0d0/,Three/3.0d0/,Four/4.0d0/,Six/6.0d0/
!
   G0=A0*F0
   GX=AX*F0+A0*FX
   GY=AY*F0+A0*FY
   GZ=AZ*F0+A0*FZ
   If(IWhat.gt.0)Then
      GXX=AXX*F0+TWO*AX*FX+A0*FXX
      GXY=AXY*F0+AX*FY+AY*FX+A0*FXY
      GXZ=AXZ*F0+AX*FZ+AZ*FX+A0*FXZ
      GYY=AYY*F0+TWO*AY*FY+A0*FYY
      GYZ=AYZ*F0+AY*FZ+AZ*FY+A0*FYZ
      GZZ=AZZ*F0+TWO*AZ*FZ+A0*FZZ
   Endif
   If(IWhat.gt.1)Then
      GXXX=AXXX*F0+THREE*AXX*FX+THREE*AX*FXX+A0*FXXX
      GXXY=AXXY*F0+AXX*FY+TWO*AXY*FX+TWO*AX*FXY+AY*FXX+A0*FXXY
      GXXZ=AXXZ*F0+AXX*FZ+TWO*AXZ*FX+TWO*AX*FXZ+AZ*FXX+A0*FXXZ
      GYYY=AYYY*F0+THREE*AYY*FY+THREE*AY*FYY+A0*FYYY
      GXYY=AXYY*F0+AYY*FX+TWO*AXY*FY+TWO*AY*FXY+AX*FYY+A0*FXYY
      GYYZ=AYYZ*F0+AYY*FZ+TWO*AYZ*FY+TWO*AY*FYZ+AZ*FYY+A0*FYYZ
      GZZZ=AZZZ*F0+THREE*AZZ*FZ+THREE*AZ*FZZ+A0*FZZZ
      GXZZ=AXZZ*F0+AZZ*FX+TWO*AXZ*FZ+TWO*AZ*FXZ+AX*FZZ+A0*FXZZ
      GYZZ=AYZZ*F0+AZZ*FY+TWO*AYZ*FZ+TWO*AZ*FYZ+AY*FZZ+A0*FYZZ
      GXYZ=AXYZ*F0+AXY*FZ+AYZ*FX+AY*FXZ+AXZ*FY+AX*FYZ+AZ*FXY+&
      &A0*FXYZ
   Endif
   If(IWhat.gt.2)Then
      GXXXX=AXXXX*F0+FOUR*AXXX*FX+SIX*AXX*FXX+FOUR*AX*FXXX+A0*FXXXX
      GXXXY=AXXXY*F0+AXXX*FY+THREE*AXXY*FX+THREE*AXX*FXY+THREE*AXY*FXX+&
      &THREE*AX*FXXY+AY*FXXX+A0*FXXXY
      GXXXZ=AXXXZ*F0+AXXX*FZ+THREE*AXXZ*FX+THREE*AXX*FXZ+THREE*AXZ*FXX+&
      &THREE*AX*FXXZ+AZ*FXXX+A0*FXXXZ
      GXXYY=AXXYY*F0+TWO*AXXY*FY+AXX*FYY+TWO*AXYY*FX+FOUR*AXY*FXY+&
      &TWO*AX*FXYY+AYY*FXX+TWO*AY*FXXY+A0*FXXYY
      GXXZZ=AXXZZ*F0+TWO*AXXZ*FZ+AXX*FZZ+TWO*AXZZ*FX+FOUR*AXZ*FXZ+&
      &TWO*AX*FXZZ+AZZ*FXX+TWO*AZ*FXXZ+A0*FXXZZ
      GXYYY=AXYYY*F0+AYYY*FX+THREE*AXYY*FY+THREE*AYY*FXY+THREE*AXY*FYY+&
      &THREE*AY*FXYY+AX*FYYY+A0*FXYYY
      GXZZZ=AXZZZ*F0+AZZZ*FX+THREE*AXZZ*FZ+THREE*AZZ*FXZ+THREE*AXZ*FZZ+&
      &THREE*AZ*FXZZ+AX*FZZZ+A0*FXZZZ
      GXXYZ=AXXYZ*F0+AXXY*FZ+AXXZ*FY+AXX*FYZ+TWO*AXYZ*FX+TWO*AXY*FXZ+&
      &TWO*AXZ*FXY+TWO*AX*FXYZ+AYZ*FXX+AY*FXXZ+AZ*FXXY+A0*FXXYZ
      GXYYZ=AXYYZ*F0+AYYZ*FX+AXYY*FZ+AYY*FXZ+TWO*AXYZ*FY+TWO*AYZ*FXY+&
      &TWO*AXY*FYZ+TWO*AY*FXYZ+AXZ*FYY+AZ*FXYY+AX*FYYZ+A0*FXYYZ
      GXYZZ=AXYZZ*F0+AYZZ*FX+AXZZ*FY+AZZ*FXY+TWO*AXYZ*FZ+TWO*AYZ*FXZ+&
      &TWO*AXZ*FYZ+TWO*AZ*FXYZ+AXY*FZZ+AY*FXZZ+AX*FYZZ+A0*FXYZZ
      GYYYY=AYYYY*F0+FOUR*AYYY*FY+SIX*AYY*FYY+FOUR*AY*FYYY+A0*FYYYY
      GYYYZ=AYYYZ*F0+AYYY*FZ+THREE*AYYZ*FY+THREE*AYY*FYZ+THREE*AYZ*FYY+&
      &THREE*AY*FYYZ+AZ*FYYY+A0*FYYYZ
      GYYZZ=AYYZZ*F0+TWO*AYYZ*FZ+AYY*FZZ+TWO*AYZZ*FY+FOUR*AYZ*FYZ+&
      &TWO*AY*FYZZ+AZZ*FYY+TWO*AZ*FYYZ+A0*FYYZZ
      GYZZZ=AYZZZ*F0+AZZZ*FY+THREE*AYZZ*FZ+THREE*AZZ*FYZ+THREE*AYZ*FZZ+&
      &THREE*AZ*FYZZ+AY*FZZZ+A0*FYZZZ
      GZZZZ=AZZZZ*F0+FOUR*AZZZ*FZ+SIX*AZZ*FZZ+FOUR*AZ*FZZZ+A0*FZZZZ
   Endif
!
   Return
End
Subroutine Props(XYZ,IFUNC,ICrit,Iwhole,J1,J2,INP,Evsave)
!
   Implicit Double Precision (A-H,O-Z)
!
   Character*8 Atnam
   Character*80 WFNTTL,JOBTTL
   Parameter(MaxAtm=100,MaxOff=200000,MaxCrt=500)
   COMMON CO(MaxOff),IC(MaxOff),NCENT,NMO,NPRIMS
   COMMON /ANG/ ANGLE(3,MaxAtm,MaxAtm)
   COMMON /DIST/ BANGLE(6),RMIN,RMAX,MINR
   COMMON /OFFSET/ ITYPE,ICENT,IMO,IEORB,IE,ICHARG,IXC,IYC,IZC,&
   &IXX, IYY, IZZ,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,IGXX,IGXY,IGXZ,IGYY,&
   &IGYZ,IGZZ,IGXXX,IGXXY,IGXXZ,IGXYY,IGXZZ,IGXYZ,IGYYY,IGYYZ,&
   &IGYZZ,IGZZZ,IGXXXX,IGXXXY,IGXXXZ,IGXXYY,IGXXZZ,IGXYYY,IGXZZZ,&
   &IGXXYZ,IGXYYZ,IGXYZZ,IGYYYY,IGYYYZ,IGYYZZ,IGYZZZ,IGZZZZ,ICofMx
   COMMON /STRING/ WFNTTL,JOBTTL,ATNAM(MaxAtm)
   COMMON /UNITS/  INPT,IOUT,IWFN,IWLP
   COMMON /WHAT/ Iwhat
   Common /Options/ Icut,Iprint,Eps,Epsnuc,Dmp,DmpNuc
   Dimension XYZ(3),HRHO(3,3),HD2RHO(3,3),HGEE(3,3),HQUAY(3,3),&
   &HVNE(3,3),WRHO(3),WD2RHO(3),WGEE(3),WQUAY(3),WVNE(3),&
   &SG(3,3),EV(3),EU(3),DSIG(3),Evsave(maxCrt,3),SV(3),PN(4),&
   &WORK(3,3),HVEE(3,3),WVEE(3)
   Save Zero,dxyz,One,Four
   Data Zero/0.d0/,dxyz/1.d-3/,One/1.d0/,Four/4.0d0/
1000 FORMAT(/,' COORDINATES OF CRITICAL POINT AND',&
   &' DISTANCE FROM MOLECULAR ORIGIN')
1010 FORMAT(/,' COORDINATES OF POINT AND DISTANCE',&
   &' FROM MOLECULAR ORIGIN')
1020 FORMAT(10X,'X = ',1PE16.8,/,10X,'Y = ',1PE16.8,/,&
   &10X,'Z = ',1PE16.8,/,10X,'R = ',1PE16.8)
1030 FORMAT(/,' VECTORS AND DISTANCES FROM NUCLEI TO CRITICAL POINT')
1040 FORMAT(/,' VECTORS AND DISTANCES FROM NUCLEI TO POINT')
1050 FORMAT(/,' NUCLEUS',8X,'   X   ',8X,'   Y   ',8X,'   Z   ',&
   &8X,'  Dist  ')
1060 FORMAT(1X,A8,1X,1P4E16.8)
1070 FORMAT(/,' ANGLES OF VECTORS WITH RESPECT TO YZ/XZ/XY PLANES',&
   &' OF MCS')
1080 FORMAT(/,' NUCLEUS',8X,'YZ ANGLE',8X,'XZ ANGLE',8X,'XY ANGLE')
1090 FORMAT(1X,A8,1X,1P3E16.8)
1100 FORMAT(/,' EIGENVALUES OF THE HESSIAN ')
1110 FORMAT(1X,1P3E18.8)
1120 FORMAT(/,' THE ELLIPTICITY IS ',1PE16.8)
1130 FORMAT(/,' EIGENVECTORS OF THE HESSIAN ')
1140 FORMAT(/,' EIGENVALUES OF THE STRESS TENSOR ')
1150 FORMAT(/,' THE TRACE OF THE STRESS TENSOR IS ',1PE16.8)
1160 FORMAT(/,' EIGENVECTORS OF THE STRESS TENSOR ')
1170 FORMAT(/,' COMPONENTS OF THE DIVERGENCE OF THE STRESS TENSOR ')
1180 FORMAT(/,' MAGNITUDE OF THE DIVERGENCE OF THE STRESS TENSOR ',&
   &1PE16.8)
1190 FORMAT(/,' VALUES ',/,' Rho(r)',18X,1PE17.10,/,&
   &' |GRAD(Rho(r))|',10X,1PE17.10,/,' GRAD(Rho(r))x',&
   &11X,1PE17.10,/,&
   &' GRAD(Rho(r))y',11X,1PE17.10,/,' GRAD(Rho(r))z',11X,1PE17.10,/,&
   &' DEL**2(Rho(r))',10X,1PE17.10,/,' G(r)',20X,1PE17.10,/,&
   &' K(r)',20X,1PE17.10,/,' L(r)',20X,1PE17.10,/,&
   &' Vnuc(r)',17X,1PE17.10,/,' V(r)',20X,1PE17.10)
1200 FORMAT(/,' DO YOU WANT TO TRACE THE BOND PATH ? (0=no/1=yes) ',$)
1210 FORMAT(/,' BOND PATH LINKED TO ',A8,' (',I3,' CALLS)',/' EPS=',&
   &1PE16.8,'   LENGTH=',1PE16.8)
1220 FORMAT(/,' CHARGE DENSITY MAXIMUM OCCURS AT ',/,&
   &12X,' X = ',1PE16.8,/,12X,' Y = ',1PE16.8,/,12X,' Z = ',1PE16.8,&
   &//,'  DISPLACEMENT = ',1PE16.8)
1230 FORMAT(/,' UNIT VECTOR FROM NUCLEUS IN BOND PATH DIRECTION ',&
   &/,'    CARTESIAN COORD.    SPHERICAL COORD.')
1240 FORMAT(1H ,3(2X,F16.8,4X,F16.8,/,1X))
1250 FORMAT(/,' TOTAL BOND PATH LENGTH = ',1PE16.8)
1260 FORMAT(' GEOMETRIC BOND LENGTH = ',1PE16.8)
1270 FORMAT(' BOND PATH LENGTH MINUS GEOMETRICAL BOND LENGTH = ',&
   &1PE16.8)
1280 FORMAT(/,' VALUES ',/,' Del**2*(Rho(r))',14X,1PE17.10,/,&
   &' |GRAD(Del**2(Rho(r)))|',7X,1PE17.10,/,&
   &' GRAD(Del**2(Rho(r)))x',8X,1PE17.10,/,&
   &' GRAD(Del**2(Rho(r)))y',8X,1PE17.10,/,&
   &' GRAD(Del**2(Rho(r)))z',8X,1PE17.10,/,&
   &' DEL**2(DEL**2(Rho(r)))',7X,1PE17.10,/,&
   &' Rho(r)',23X,1PE17.10,/,' G(r)',25X,1PE17.10,/,&
   &' K(r)',25X,1PE17.10,/,' L(r)',25X,1PE17.10,/,&
   &' Vnuc(r)',22X,1PE17.10,/,' V(r)',25X,1PE17.10)
1290 FORMAT(/,' VALUES ',/,' G(r)',20X,1PE17.10,/,&
   &' |GRAD(G(r)|',13X,1PE17.10,/,' GRAD(G(r))x',&
   &13X,1PE17.10,/,' GRAD(G(r))y',13X,1PE17.10,/,&
   &' GRAD(G(r))z',13X,1PE17.10,/,' DEL**2(G(r))',12X,1PE17.10,/,&
   &' Rho(r)',18X,1PE17.10,/,' DEL**2(Rho(r))',10X,1PE17.10,/,&
   &' K(r)',20X,1PE17.10,/,' L(r)',20X,1PE17.10,/,' Vnuc(r)',&
   &17X,1PE17.10,/,' V(r)',20X,1PE17.10)
1300 FORMAT(/,' VALUES ',/,' K(r)',20X,1PE17.10,/,&
   &' |GRAD(K(r))|',12X,1PE17.10,/,' GRAD(K(r))x',&
   &13X,1PE17.10,/,' GRAD(K(r))y',13X,1PE17.10,/,&
   &' GRAD(K(r))z',13X,1PE17.10,/,' Del**2(G(r))',12X,1PE17.10,/,&
   &' Rho(r)',18X,1PE17.10,/,' DEL**2(Rho(r))',10X,1PE17.10,/,&
   &' G(r)',20X,1PE17.10,/,' L(r)',20X,1PE17.10,/,' Vnuc(r)',&
   &17X,1PE17.10,/,' V(r)',20X,1PE17.10)
1310 FORMAT(/,' VALUES ',/,' Vnuc(r)',17X,1PE17.10,/,&
   &' |GRAD(Vnuc(r))|',9X,1PE17.10,/,' GRAD(Vnuc(r))x',&
   &10X,1PE17.10,/,' GRAD(Vnuc(r))y',10X,1PE17.10,/,&
   &' GRAD(Vnuc(r))z',10X,1PE17.10,/,' Del**2(Vnuc(r))',&
   &9X,1PE17.10,/,' Rho(r)',18X,1PE17.10,/,' DEL**2(Rho(r))',&
   &10X,1PE17.10,/,' G(r)',20X,1PE17.10,/,' K(r)',20X,&
   &1PE17.10,/,' L(r)',20X,1PE17.10,/,' V(r)',20X,1PE17.10)
1320 FORMAT(/,' VALUES ',/,' V(r)',20X,1PE17.10,/,&
   &' |GRAD(V(r))|',12X,1PE17.10,/,' GRAD(V(r))x',&
   &13X,1PE17.10,/,' GRAD(V(r))y',13X,1PE17.10,/,&
   &' GRAD(V(r))z',13X,1PE17.10,/,' Del**2(V(r))',&
   &12X,1PE17.10,/,' Rho(r)',18X,1PE17.10,/,' DEL**2(Rho(r))',&
   &10X,1PE17.10,/,' G(r)',20X,1PE17.10,/,' K(r)',20X,&
   &1PE17.10,/,' L(r)',20X,1PE17.10,/,' Vnuc(r)',17X,1PE17.10)
!
   R = DSQRT(XYZ(1)**2+XYZ(2)**2+XYZ(3)**2)
   If(Icrit.eq.1)Then
      If(Iprint.eq.1)WRITE (IOUT,1000)
      WRITE (IWLP,1000)
   Else
      If(Iprint.eq.1)WRITE (IOUT,1010)
      WRITE (IWLP,1010)
   Endif
   If(Iprint.eq.1)WRITE (IOUT,1020) XYZ(1),XYZ(2),XYZ(3),R
   WRITE (IWLP,1020) XYZ(1),XYZ(2),XYZ(3),R
!
   If(Icrit.eq.1)Then
      If(Iprint.eq.1)WRITE (IOUT,1030)
      WRITE (IWLP,1030)
   Else
      If(Iprint.eq.1)WRITE (IOUT,1040)
      WRITE (IWLP,1040)
   Endif
!
   If(Iprint.eq.1)WRITE (IOUT,1050)
   WRITE (IWLP,1050)
   DO 10 I = 1,NCENT
      XREL=XYZ(1)-CO(IXC+I)
      YREL=XYZ(2)-CO(IYC+I)
      ZREL=XYZ(3)-CO(IZC+I)
      Dist=dsqrt(xrel**2+yrel**2+zrel**2)
      If(Iprint.eq.1)WRITE(IOUT,1060)ATnam(i),Xrel,Yrel,Zrel,Dist
      WRITE(IWLP,1060)ATnam(i),Xrel,Yrel,Zrel,Dist
10 Continue
!
   If(Iprint.eq.1)WRITE (IOUT,1070)
   WRITE (IWLP,1070)
   If(Iprint.eq.1)WRITE (IOUT,1080)
   WRITE (IWLP,1080)
   DO 20 I = 1,NCENT
      CALL GEOM (I,XYZ,RN,AYZ,AXZ,AXY)
      If(Iprint.eq.1)WRITE(IOUT,1090) ATNAM(I),AYZ,AXZ,AXY
      WRITE(IWLP,1090) ATNAM(I),AYZ,AXZ,AXY
20 CONTINUE
!
   If(IFUNC.eq.1)Then
      IWhat=2
      CALL GRDRHO(1,XYZ,RHO,WRHO,GRHO,HRHO,SG)
      CALL GRDD2R(1,XYZ,D2Rho,WD2Rho,GD2Rho,HD2Rho)
      CALL GRDKEG(1,XYZ,GEE,WGEE,GGEE,HGEE)
      CALL GRDKEK(1,XYZ,QUAY,WQUAY,GQUAY,HQUAY)
      CALL GRDVNE(1,XYZ,VNE,WVNE,GVNE,HVNE)
      CALL GRDV(1,XYZ,VEE,WVEE,GVEE,HVEE)
      CALL TRACE(HRHO,EV,WORK,3,IFAIL)
      If(ICRIT.EQ.1)Then
         EVSAVE(INP,1)=EV(1)
         EVSAVE(INP,2)=EV(2)
         EVSAVE(INP,3)=EV(3)
      Endif
      CALL TRACE(SG,EU,WORK,3,IFAIL)
      SV(1) = HRHO(1,3)
      SV(2) = HRHO(2,3)
      SV(3) = HRHO(3,3)
      XLag = -D2Rho/Four
      If(Iprint.eq.1)WRITE (IOUT,1100)
      WRITE (IWLP,1100)
      If(Iprint.eq.1)WRITE (IOUT,1110) (EV(I),I=1,3)
      WRITE (IWLP,1110) (EV(I),I=1,3)
      IF ((EV(1).LT.Zero).AND.(EV(2).LT.Zero)) THEN
         ELLIPT = (EV(1)/EV(2))-One
         If(Iprint.eq.1)WRITE (IOUT,1120) ELLIPT
         WRITE (IWLP,1120) ELLIPT
      ENDIF
      If(Iprint.eq.1)WRITE (IOUT,1130)
      WRITE (IWLP,1130)
      DO 150 I = 1,3
         If(Iprint.eq.1)WRITE(IOUT,1110)(HRHO(I,J),J=1,3)
         WRITE (IWLP,1110) (HRHO(I,J),J=1,3)
150   CONTINUE
!
      If(Iprint.eq.1)WRITE (IOUT,1140)
      WRITE (IWLP,1140)
      If(Iprint.eq.1)WRITE (IOUT,1110) (EU(I),I=1,3)
      WRITE (IWLP,1110) (EU(I),I=1,3)
      TRSG = EU(1)+EU(2)+EU(3)
      If(Iprint.eq.1)WRITE (IOUT,1150) TRSG
      WRITE (IWLP,1150) TRSG
      If(Iprint.eq.1)WRITE (IOUT,1160)
      WRITE (IWLP,1160)
      DO 160 I = 1,3
         If(Iprint.eq.1)WRITE (IOUT,1110) (SG(I,J),J=1,3)
         WRITE (IWLP,1110) (SG(I,J),J=1,3)
160   CONTINUE
      CALL DIVSTR(DSIG,DSIGM)
      If(Iprint.eq.1)WRITE (IOUT,1170)
      WRITE (IWLP,1170)
      If(Iprint.eq.1)WRITE (IOUT,1110) (DSIG(I),I=1,3)
      WRITE (IWLP,1110) (DSIG(I),I=1,3)
      If(Iprint.eq.1)WRITE (IOUT,1180) DSIGM
      WRITE (IWLP,1180) DSIGM
      If(Iprint.eq.1)Then
         WRITE (IOUT,1190) RHO,GRHO,WRHO(1),WRHO(2),WRHO(3),D2RHO,&
         &GEE,QUAY,XLag,VNE,VEE
      Endif
      WRITE (IWLP,1190) RHO,GRHO,WRHO(1),WRHO(2),WRHO(3),D2RHO,&
      &GEE,QUAY,XLag,VNE,VEE
!
      IF ((EV(1).LT.Zero) .AND. (EV(2).LT.Zero)&
      &.AND. (EV(3).GT.Zero).and.iwhole.eq.0.and.icrit.eq.1)THEN
         WRITE (IOUT,1200)
         READ (INPT,*) IWK
         IF (IWK .EQ. 1) THEN
            X1 = XYZ(1) + DXYZ*SV(1)
            Y1 = XYZ(2) + DXYZ*SV(2)
            Z1 = XYZ(3) + DXYZ*SV(3)
            CALL BOND(X1,Y1,Z1,BL1,NCAL,PN,IFunc,IWhole)
            If(Iprint.eq.1)WRITE (IOUT,1210) ATNAM(MINR),NCAL,EPSD,BL1
            WRITE (IWLP,1210) ATNAM(MINR),NCAL,EPSD,BL1
            If(Iprint.eq.1)WRITE (IOUT,1220) (PN(I),I=1,4)
            WRITE (IWLP,1220) (PN(I),I=1,4)
end module density
