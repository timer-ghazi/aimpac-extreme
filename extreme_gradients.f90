! Generated from extreme.f90
! Contains: GRDRHO, GRDD2R, GRDKEG, GRDKEK, GRDV, GRDVNE, MULTI1, DES, NEWTON


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

SUBROUTINE MULTI1(W,R,W1,I,DS,S)
!
   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
   COMMON /C2/ EPSD,AB(9,9),AM(9,10)
   DIMENSION W1(3,I),W(3),R(3),R1(3),HESS(3,3),SG(3,3)
!
!     PERFORM A PREDICTOR STEP OF ORDER I.
!
   DO 100 L = 1,3
      H = 0.0D0
      DO 110 J = 1,I
         H = H+AB(I,J)*W1(L,J)
110   CONTINUE
      R(L) = W(L)+H*DS
100 CONTINUE
!
!     PERFORM A CORRECTOR STEP OF ORDER I+1.
!
   CALL GRDRHO(0,R,RHO,R1,Gnorm,HESS,SG)
   DO 120 L = 1,3
      H = AM(I,1)*R1(L)/Gnorm
      DO 130 J = 1,I
         H = H+AM(I,J+1)*W1(L,J)
130   CONTINUE
      R(L) = W(L)+H*DS
120 CONTINUE
   RETURN
END

SUBROUTINE DES(R1,PR,DS,NN)
!
   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
   Parameter (MaxStp=141)
   COMMON /C2/ EPSD,AB(9,9),AM(9,10)
   DIMENSION R1(3,MaxStp),W1(3,9),R(3),W(3),PR(3,MaxStp),&
   &H(3,3), SG(3,3),tmpxyz(3),grho(3)
   Save Zero,One
   Data Zero/0.0d0/,One/1.0d0/
!
!     IF IT IS THE FIRST STEP, START WITH A SINGLE STEP METHOD.
!
   DO 100 I = 1,5
      DO 110 L = 1,3
         W(L) = R1(L,I)
         DO 110 J = 1,I
            G = DSQRT(PR(1,I-J+1)**2+PR(2,I-J+1)**2+PR(3,I-J+1)**2)
            W1(L,J) = PR(L,I-J+1)/G
110   CONTINUE
      CALL MULTI1(W,R1(1,I+1),W1,I,DS,DFLOAT(I)*DS)
      tmpxyz(1)=R1(1,I+1)
      tmpxyz(2)=R1(2,I+1)
      tmpxyz(3)=R1(3,I+1)
      CALL GRDRHO(0,tmpxyz,RHO,GRHO,GRAD,H,SG)
      PR(1,I+1)=GRHO(1)
      PR(2,I+1)=GRHO(2)
      PR(3,I+1)=GRHO(3)
100 CONTINUE
!
!     PERFORM THE REMAINING STEPS WITH A 6 STEP METHOD OF ORDER 7.
!
   DO 130 I = 6,NN
      DO 140 L = 1,3
         W(L) = R1(L,I)
         DO 140 J = 1,6
            G = DSQRT(PR(1,I-J+1)**2+PR(2,I-J+1)**2+PR(3,I-J+1)**2)
            W1(L,J) = PR(L,I-J+1)/G
140   CONTINUE
      CALL MULTI1(W,R1(1,I+1),W1,6,DS,DFLOAT(I)*DS)
      tmpxyz(1)=R1(1,I+1)
      tmpxyz(2)=R1(2,I+1)
      tmpxyz(3)=R1(3,I+1)
      CALL GRDRHO(0,tmpxyz,RHO,GRHO,GRAD,H,SG)
      PR(1,I+1)=GRHO(1)
      PR(2,I+1)=GRHO(2)
      PR(3,I+1)=GRHO(3)
      IF (I .NE. (I/20)*20) GOTO 130
      DO 150 L = 1,3
         W(L) = R1(L,I)
150   CONTINUE
      CALL MULTI1(W,R,W1,5,DS,DFLOAT(I)*DS)
!
!     ESTIMATE THE ERROR.
!
      EPS1 = Zero
      DO 160 L = 1,3
         EPS1 = EPS1+(R(L)-R1(L,I+1))**2
160   CONTINUE
      EPS1 = DSQRT(EPS1)
      RH = DMAX1(One,DSQRT(R1(1,I+1)**2+R1(2,I+1)**2+R1(3,I+1)**2))
      EPSD = DMAX1(EPSD,EPS1/RH)
130 CONTINUE
   RETURN
END

SUBROUTINE NEWTON(R,IFAIL,IFUNC,IWHOLE,NITER,INuc)
!
   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
   Parameter(MaxSch=20)
   COMMON /UNITS/ INPT, IOUT, IWFN, IWLP
   Common /What/ IWhat
   Common /Options/ ICut,IPrint,Eps,EpsNuc,Dmp,Dmpnuc
   INTEGER INFO
!     Use LAPACK routine DGESV for solving the Hessian system
   DIMENSION W(3),R1(3),IW(3),XX(2),R(3),SG(3,3),H(3,3)
   Save Zero,pt1,pt05,pt025,One
   Data Zero/0.0d0/,pt1/0.1d0/,pt05/0.05d0/,pt025/0.025d0/,&
   &One/1.0d0/,pt25/0.25d0/
1000 FORMAT(' STEP',3X,'CURRENT COORDINATES ',40X,'|GRAD(Rho)|',/)
1001 FORMAT(' STEP',3X,'CURRENT COORDINATES ',40X,&
   &'|GRAD(DEL**2(Rho))|',/)
1002 FORMAT(' STEP',3X,'CURRENT COORDINATES ',40X,&
   &'|GRAD(G)|',/)
1003 FORMAT(' STEP',3X,'CURRENT COORDINATES ',40X,&
   &'|GRAD(K)|',/)
1004 FORMAT(' STEP',3X,'CURRENT COORDINATES ',40X,&
   &'|GRAD(Vnuc)|',/)
1005 FORMAT(' STEP',3X,'CURRENT COORDINATES ',40X,&
   &'|GRAD(V)|',/)
1010 FORMAT(1X,I4,3X,1P4E18.10)
1020 FORMAT(' NUMBER OF NEWTON-RAPHSON ITERATIONS : ',I4)
1030 FORMAT(' CRITICAL POINT NOT YET FOUND.  CONTINUE SEARCH',&
   &' ? (0=no/1=yes) ',$)
1040 FORMAT(' TRY AGAIN ')
!
   Cutoff=eps
   If(Inuc.eq.1)cutoff=epsnuc
   Damp=Dmp
   If(Inuc.eq.1)Damp=Dmpnuc
!
   IWhat=1
!
   CALL GRDRHO(1,R,RHO,W,GRAD,H,SG)
!
   IS = 0
   If(Iprint.eq.1)Then
      IF (IFAIL .EQ. 0) Then
         If(Iwhole.eq.0)Then
            WRITE (IOUT,1000)
            WRITE (IOUT,1010) IS,R(1),R(2),R(3),GRAD
         Endif
      Endif
   Endif
!
!    BEGIN NEWTON RAPHSON SEARCH
!
99 DO 100 I = 1,MaxSch
      IS = IS + 1
      DO 105 J = 1,3
         R1(J) = W(J)
105   CONTINUE
      CALL DGESV(3,1,H,3,IW,R1,3,INFO)
      IF(INFO.NE.0) THEN
         IFAIL = 1
         GOTO 150
      ENDIF
      DO 110 J = 1,3
         DJ = DABS(R1(J))
         IF (DJ.GT.Damp) R1(J) = Damp*R1(J)/DJ
         R(J) = R(J) - R1(J)
110   CONTINUE
!
      CALL GRDRHO(1,R,RHO,W,GRAD,H,SG)
!
      IF (IFAIL .EQ. 0.and.Iwhole.eq.0.and.iprint.eq.1) Then
         WRITE (IOUT,1010) IS,R(1),R(2),R(3),GRAD
      Endif
!
      IF (GRAD .LE. Cutoff) THEN
         NITER=I
         Goto 150
      END IF
100 CONTINUE
!
!    SHALL WE CONTINUE ?
!
   IF (IFAIL .EQ. 1) Goto 150
!
130 CONTINUE
   IYN=0
   IF(IWHOLE.EQ.0)Then
      WRITE (IOUT,1030)
      READ (INPT,*) IYN
   Endif
   IF (IYN .EQ. 0) THEN
      IFAIL = 1
      Goto 150
   ELSE IF (IYN .EQ. 1) THEN
      GOTO 99
   ELSE
      WRITE (IOUT,1040)
      GOTO 130
   END IF
!
150 Continue
!
   RETURN
END
