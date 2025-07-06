module io
  implicit none
contains
   SUBROUTINE RDWFN
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      CHARACTER*80 WFNTTL,JOBTTL
      CHARACTER*8 ATNAM
      CHARACTER*4 MODE
      PARAMETER (NTYPE=20,MaxOff=2000000,MaxAtm=500)
      COMMON CO(MaxOff),IC(MaxOff),NCENT,NMO,NPRIMS
      COMMON /OFFSET/ ITYPE,ICENT,IMO,IEORB,IE,ICHARG,IXC,IYC,IZC,&
      &IXX, IYY, IZZ,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,IGXX,IGXY,IGXZ,IGYY,&
      &IGYZ,IGZZ,IGXXX,IGXXY,IGXXZ,IGXYY,IGXZZ,IGXYZ,IGYYY,IGYYZ,&
      &IGYZZ,IGZZZ,IGXXXX,IGXXXY,IGXXXZ,IGXXYY,IGXXZZ,IGXYYY,IGXZZZ,&
      &IGXXYZ,IGXYYZ,IGXYZZ,IGYYYY,IGYYYZ,IGYYZZ,IGYZZZ,IGZZZZ,ICOFMX
      COMMON /STRING/ WFNTTL,JOBTTL,ATNAM(MaxAtm)
      COMMON /UNITS/  INPT,IOUT,IWFN,IWLP
      DATA ENDATA /8HEND DATA/,ZERO/0.d0/
!
      READ(IWFN,101) WFNTTL
!
      READ(IWFN,102) MODE,NMO,NPRIMS,NCENT
!
!     Set OffSets For Storage in IC:
!
!                             FUNCTION TYPES
      ITYPE = 0
!                             FUNCTION CENTRES
      ICENT = ITYPE+NPRIMS
!
!     Set OffSets For Storage in CO - MO Coefficients, Then others:
!
      IMO=0
!                             ORBITAL ENERGIES OF EACH M.O.
      IEORB = IMO+NMO*NPRIMS
!                             EXPONENTS OF THE FUNCTIONS
      IE = IEORB+NMO
!                             NUCLEAR CHARGE OF EACH CENTRE
      ICHARG = IE+NPRIMS
!
!                             XYZ COORDS OF THE CENTRES IN
!                             ORIGINAL CARTESIAN SYSTEM
      IXC = ICHARG+NCENT
      IYC = IXC+NCENT
      IZC = IYC+NCENT
!
!                             XYZ COORDS OF CENTRES RELATIVE TO A
!                             CURRENT TEST POINT IN THE INTEGRATION
      IXX = IZC+NCENT
      IYY = IXX+NCENT
      IZZ = IYY+NCENT
!
!                             SQUARE OF DISTANCE FROM CENTRES TO
!                             CURRENT POINT
      IRR = IZZ+NCENT
!                             DISTANCE FROM CENTRES TO CURRENT POINT
      IR2 = IRR+NCENT
!                             OCCUPATION NUMBER OF EACH M.O
      IP = IR2+NCENT
!                             PSI VALUES FOR EACH M.O.
      IPSI = IP+NMO
!
      IGX = IPSI+NMO
      IGY = IGX+NMO
      IGZ = IGY+NMO
!
      IGXX = IGZ+NMO
      IGXY = IGXX+NMO
      IGXZ = IGXY+NMO
      IGYY = IGXZ+NMO
      IGYZ = IGYY+NMO
      IGZZ = IGYZ+NMO
!
      IGXXX=IGZZ+NMO
      IGXXY=IGXXX+NMO
      IGXXZ=IGXXY+NMO
      IGXYY=IGXXZ+NMO
      IGXZZ=IGXYY+NMO
      IGXYZ=IGXZZ+NMO
      IGYYY=IGXYZ+NMO
      IGYYZ=IGYYY+NMO
      IGYZZ=IGYYZ+NMO
      IGZZZ=IGYZZ+NMO
!
      IGXXXX=IGZZZ+NMO
      IGXXXY=IGXXXX+NMO
      IGXXXZ=IGXXXY+NMO
      IGXXYY=IGXXXZ+NMO
      IGXXZZ=IGXXYY+NMO
      IGXYYY=IGXXZZ+NMO
      IGXZZZ=IGXYYY+NMO
      IGXXYZ=IGXZZZ+NMO
      IGXYYZ=IGXXYZ+NMO
      IGXYZZ=IGXYYZ+NMO
      IGYYYY=IGXYZZ+NMO
      IGYYYZ=IGYYYY+NMO
      IGYYZZ=IGYYYZ+NMO
      IGYZZZ=IGYYZZ+NMO
      IGZZZZ=IGYZZZ+NMO
      ICOFMX=IGZZZZ+NMO
!
      DO 100 I = 1,NCENT
         READ (IWFN,103) ATNAM(I),J,CO(IXC+J),CO(IYC+J),CO(IZC+J),&
         &CO(ICHARG+J)
100   CONTINUE
      READ (IWFN,104) (IC(ICENT+I),I=1,NPRIMS)
      READ (IWFN,104) (IC(ITYPE+I),I=1,NPRIMS)
      READ (IWFN,105) (CO(IE+I),I=1,NPRIMS)
!
      DO 130 I = 1,NPRIMS
         IF(IC(ITYPE+I).GT.NTYPE) GOTO 999
130   CONTINUE
!
      DO 120 I = 1,NMO
         READ (IWFN,106) CO(IP+I),CO(IEORB+I)
         K = NPRIMS*(I-1)+IMO
         READ (IWFN,107) (CO(K+J),J=1,NPRIMS)
120   Continue
!
      Do 125 I=1,NPrims
         check=zero
         Do 126 J=1,NMO
            temp=dabs(CO(IMO+NPRIMS*(J-1)+I))
            If(temp.gt.check)check=temp
126      Continue
         CO(ICOFMX+I)=check
125   Continue

      READ (IWFN,108) CHECK
      IF (CHECK .NE. ENDATA) STOP ' RDWFN : END CARD NOT FOUND '
!
!    READ IN TOTAL SCF ENERGY AND -V/T
!
      READ (IWFN,109) TOTE,GAMMA
      GOTO 9999
999   STOP 'EXTREME CANNOT WORK WITH g-, h- OR HIGHER PRIMITIVES'
9999  CONTINUE
      RETURN
!
!
101   FORMAT (A80)
102   FORMAT (4X,A4,12X,3(I3,17X))
103   FORMAT(A8,11X,I3,2X,3F12.8,10X,F5.1)
104   FORMAT (20X,20I3)
105   FORMAT (10X,5E14.7)
106   FORMAT(35X,F12.8,15X,F12.8)
107   FORMAT (5E16.8)
108   FORMAT(A8)
109   FORMAT(17X,F20.12,18X,F13.8)
   END
   Subroutine Store(XYZ,X,Y,Z,J1,J2,Icancl,INP,JI1,JI2,NITER)
!
      Implicit Double Precision (A-H,O-Z)
!
      Parameter(MaxCrt=500)
      COMMON /UNITS/  INPT,IOUT,IWFN,IWLP
      Common /options/ Icut,Iprint,Eps,Epsnuc,Dmp,DmpNuc
      Dimension XYZ(3),X(MaxCrt),Y(MaxCrt),Z(MaxCrt),JI1(MaxCrt),&
      &JI2(MaxCrt)
      Save Close
      Data Close/1.d-6/
1000  FORMAT(' NEW CRITICAL POINT FOUND:  NUMBER ',I4)
1010  FORMAT(' REDUNDANT CRITICAL POINT FOUND:  SAME AS NUMBER ',I4)
1020  FORMAT(' NUMBER OF NEWTON-RAPHSON ITERATIONS : ',I4)
!
      ICancl=0
      If(INP.gt.0)Then
         DO 10 I=1,INP
            Dist=dsqrt((xyz(1)-X(I))**2+(xyz(2)-Y(I))**2+(XYZ(3)-Z(I))**2)
            If(Dist.lt.close)Icancl=1
            IF(Dist.lt.Close)ISame=I
10       Continue
      Endif
!
      If(Icancl.eq.0)Then
         INP = INP + 1
         X(INP) = XYZ(1)
         Y(INP) = XYZ(2)
         Z(INP) = XYZ(3)
         JI1(INP) = J1
         JI2(INP) = J2
         Write(iout,*)
         Write(IOUT,1000)INP
         If(Iprint.eq.1)Write(iout,*)
         If(Iprint.eq.1)Write(IOUT,1020)NITER
         Write(iwlp,*)
         Write(IWLP,1000)INP
         Write(iwlp,*)
         Write(IWLP,1020)NITER
      Else
         If(Iprint.eq.1)Write(iout,*)
         If(Iprint.eq.1)Write(iwlp,*)
         If(Iprint.eq.1)Write(IOUT,1010)Isame
         If(Iprint.eq.1)Write(IWLP,1010)Isame
      Endif
!
      Return
   End
end module io
