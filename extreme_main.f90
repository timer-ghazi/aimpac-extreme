! Generated from extreme.f90
! Main program

PROGRAM EXTREME
!
!    VERSION 94 - Revision B
!
!    EXTREME SEARCHES 3-D MOLECULAR SPACE FOR CRITICAL POINTS OF THE
!    FOLLOWING SCALAR FIELDS:
!
!    1) THE ELECTRON DENSITY DISTRIBUTION (RHO).
!    2) THE LAPLACIAN OF THE ELECTRON DENSITY DISTRIBUTION (DEL**2(RHO)).
!    3) THE LAGRANGIAN KINETIC ENERGY DENSITY DISTRIBUTION (KEG).
!    4) THE HAMILTONIAN KINETIC ENERGY DENSITY DISTRIBUTION (KEK).
!    5) THE NUCLEAR POTENTIAL DISTRIBUTION (Vnuc).
!    6) THE POTENTIAL ENERGY DENSITY - TRACE OF THE STRESS TENSOR. (V).
!
!    AS INPUT, AN AIMPAC WAVEFUNCTION FILE IS REQUIRED.  THE REST OF THE
!    INPUT IS SUPPLIED BY THE USER INTERACTIVELY.
!
!    BASICALLY, THE CRITICAL POINTS ARE FOUND USING A NEWTON-RAPHSON
!    SEARCH, WITH INITIAL GUESSES SUPPLIED BY THE USER, OR AUTOMATICALLY.
!
!    EXTREME CAN HANDLE S, P(3), D(6) and F(10) TYPE GAUSSIAN BASIS
!    FUNCTIONS.
!
!    ORIGINAL EXTREME (SADDLE) DEVELOPED BY VARIOUS MEMBERS OF
!    RICHARD BADER'S RESEARCH GROUP - MOST EXTENSIVELY BY KEITH E
!    LAIDIG APRIL 1989.
!
!    HEAVILY MODIFIED, EXTENDED AND CLEANED UP:  ADDITION OF THE FIELDS
!    DEL**2(RHO), G, K, NUCLEAR POTENTIAL AND POTENTIAL ENERGY DENSITY.
!    INCORPORATION OF AUTOMATED SEARCHING OPTIONS, AND MANY OTHER
!    NICETIES.
!    TAK 8/94
!
!    Questions and Suggestions Should be Directed to:
!    R.F.W. Bader:  bader@mcmail.cis.mcmaster.ca
!    or
!    T.A. Keith:    keith@babbage.chemistry.mcmaster.ca
!
   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
   CHARACTER*8  ATNAM
   CHARACTER*80 JOBTTL,WFNTTL
   CHARACTER*4 FWFN /'.wfn'/, FWLP /'.crt'/
   CHARACTER*8 BLANK /'        '/
   CHARACTER*40 WFN,WLP
   Parameter(MaxAtm=100,MaxOff=200000,MaxCrt=500)
   COMMON CO(MaxOff),IC(MaxOff),NCENT,NMO,NPRIMS
   COMMON /ANG/ ANGLE(3,MaxAtm,MaxAtm)
   COMMON /OFFSET/ ITYPE,ICENT,IMO,IEORB,IE,ICHARG,IXC,IYC,IZC,&
   &IXX, IYY, IZZ,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,IGXX,IGXY,IGXZ,IGYY,&
   &IGYZ,IGZZ,IGXXX,IGXXY,IGXXZ,IGXYY,IGXZZ,IGXYZ,IGYYY,IGYYZ,&
   &IGYZZ,IGZZZ,IGXXXX,IGXXXY,IGXXXZ,IGXXYY,IGXXZZ,IGXYYY,IGXZZZ,&
   &IGXXYZ,IGXYYZ,IGXYZZ,IGYYYY,IGYYYZ,IGYYZZ,IGYZZZ,IGZZZZ,ICofMx
   COMMON /STRING/ WFNTTL,JOBTTL,ATNAM(MaxAtm)
   COMMON /UNITS/  INPT,IOUT,IWFN,IWLP
   COMMON /OPTIONS/ ICut,IPrint,Eps,EpsNuc,Dmp,DmpNuc
   DIMENSION JI1(MaxCrt),JI2(MaxCrt),X(MaxCrt),Y(MaxCrt),Z(MaxCrt),&
   &XYZ(3),EVSAVE(MaxCrt,3)
   DATA Zero/0.d0/,One/1.d0/,Two/2.d0/,Three/3.d0/,Ten/10.0d0/&
   &Four/4.d0/,pt5/0.5d0/,degen/1.d-8/
1000 FORMAT(' ENTER THE TITLE OF THIS RUN ')
1010 FORMAT(A80)
1020 FORMAT(' EXTREME ')
1030 Format(' NUCLEAR COORDINATES')
1040 FORMAT(6X,A8,4X,3F15.8)
1050 FORMAT(' WHICH SCALAR FIELD TO ANALYZE:  ',/,' ELECTRON',&
   &' DENSITY, RHO (1)',/,' LAPLACIAN OF RHO, DEL**2(RHO) (2)',/,&
   &' LAGRANGIAN KINETIC ENERGY DENSITY, G (3)',/,' HAMILTONIAN ',&
   &'KINETIC ENERGY DENSITY, K (4)',/,&
   &' NUCLEAR POTENTIAL, Vnuc (5)',/,&
   &' POTENTIAL ENERGY DENSITY, V (6)',/,' Enter (1-6):  ', $)
1060 FORMAT(/,' TRY AGAIN ')
1070 FORMAT(' ANALYZING RHO, THE ELECTRON DENSITY')
1080 FORMAT(' ANALYZING DEL**2(RHO), THE LAPLACIAN OF THE ELECTRON',&
   &' DENSITY')
1090 FORMAT(' ANALYZING G, THE LAGRANGIAN KINETIC ENERGY DENSITY')
1100 FORMAT(' ANALYZING K, THE HAMILTONIAN KINETIC ENERGY DENSITY')
1110 FORMAT(' ANALYZING Vnuc, THE NUCLEAR POTENTIAL FIELD')
1115 FORMAT(' ANALYZING V, THE POTENTIAL ENERGY DENSITY ')
1120 FORMAT(/,' COORDS or NUCLEAR (0), BOND (1), RING (2),',&
   &' CAGE (3), ANGLE (4), POINT (5),',/,&
   &' MEGA (6), OPTIONS (7), LIST CURRENT CPs (8), STOP (9)',/,&
   &' Enter (1-9)  ',$)
1121 Format(/,' Nuclear Coordinates (1) or Arbitrary',&
   &' Coordinates (0) ? ',$)
1122 Format(/' Enter Nucleus Number for Nuclear Critical Point',&
   &' Search  ',$)
1130 FORMAT(/,' INPUT STARTING COORDINATES FOR NEWTON-RAPHSON',&
   &' SEARCH ',/,' ',$)
1140 FORMAT(/,' STARTING COORDINATES FOR NEWTON-RAPHSON SEARCH ',&
   &/,1X,3F15.8)
1141 FORMAT(/,' NUCLEUS ',I4,' STARTING COORDINATES FOR',&
   &' NEWTON-RAPHSON SEARCH ',&
   &/,1X,3F15.8)
1150 FORMAT(/,' INPUT NUMBERS OF TWO ATOMS TO SEARCH BETWEEN: ',$)
1160 FORMAT(/,' ONE OR MORE OF THE ATOMS DOES NOT EXIST:',&
   &' INPUT NUMBERS AGAIN ',$)
1170 FORMAT(/,' START HALFWAY BETWEEN ATOMS (1) OR ANOTHER FRACTION',&
   &' OF DISTANCE (0) ? ',$)
1180 FORMAT(/,' FRACTION (eg. 0.5) OF DISTANCE FROM ATOM ',I4,&
   &' TO ATOM ', I4,' TO START SEARCH ? ',/,$)
1190 FORMAT(/,' SEARCHING BETWEEN ATOMS ',2I4)
1200 FORMAT(/,' INPUT NUMBERS OF THREE ATOMS TO SEARCH BETWEEN: ',$)
1210 FORMAT(/,' SEARCHING BETWEEN ATOMS ',3I4)
1220 FORMAT(/,' INPUT NUMBERS OF FOUR ATOMS TO SEARCH BETWEEN: ',$)
1230 FORMAT(/,' SEARCHING BETWEEN ATOMS ',4I4)
1240 FORMAT(/,' INPUT COORDINATES FOR PROPERTY EVALUATION ',/,' ')
1250 FORMAT(/,' COORDINATES FOR PROPERTY EVALUATION ',/,' ',3F9.6)
1260 FORMAT(/,' Using Avg of Corresponding Nuclear Coordinates as',&
   &' Initial Guess,',&
   &/,' Search:  For Nuclear Critical Points (1); Between',&
   &' Atom Pairs (2);',/,' Between Atom Triads (3);',&
   &' Between Atom Quads (4); All of the Above (5);',/,&
   &' or None of the Above (0) ? ',$)
1265 FORMAT(/,' Search Along Whole Line Between Two Atoms (1)',/,&
   &' Search Another Line (2)',/,&
   &' None (0)',/,' Enter (0-2)  ',$)
1266 Format(/,' Enter Two Atom Numbers and the Number of Starting',&
   &' Points  ',$)
1267 Format(/,' Enter Coordinates of One End of Desired Line',&
   &/,$)
1268 Format(/,' Enter Coordinates of Other End of ',&
   &' Line and Number of Starting Points',/,$)
1270 FORMAT(/,' Search Over one or all atom-centered spheres ?',&
   &' (1=yes/0=no) ',$)
1280 FORMAT(/,' What Atom Number ?  Type 0 for all atoms ',$)
1290 FORMAT(/,' Input Spherical Grid Search Parameters: ',/,&
   &' NPhi,NTheta,RadMin,RadMax and NradPt (eg. 6 6 0.5 2.0 10) ',&
   &/,$)
1295 FORMAT(' Change Accuracy to be Used in Evaluating',&
   &' Functions (1)',/,&
   &' Change Print Options (2)',/,' Change Gradient Threshold',&
   &' For Calling',&
   &' a Point a Critical Point (3) ',/,&
   &' Change Damping Factor for Newton-Raphson Search (4)',/,&
   &' Finished With Options (0)',/,&
   &' Enter (0-4)  ',$)
1296 FORMAT(' Enter Integer N (Default is 20) to Change the Function',&
   &/,' Accuracy to 10**(-N)  ',$)
1297 FORMAT(' Print Output to Screen ? (1=yes/0=no(default))  ',$)
1298 FORMAT(' Function Accuracy Set to ',1PE16.8)
1299 FORMAT(' Print Option Set to ',I3)
1300 FORMAT(/,' CRITICAL POINTS: ')
1301 FORMAT(' Enter Integer N (Default is 10) to Change the',/,&
   &' Gradient Threshhold Value to 10**(-N) for Non-Nuclear',/,&
   &' Critical Points and 10**(-(N-2)) for Nuclear Critical',&
   &' Points  ',$)
1302 FORMAT(' Gradient Threshold Set to ',1PE16.8,&
   &' for Non-Nuclear',/, ' Critical Points and ',1PE16.8,&
   &' for Nuclear Critical Points')
1303 FORMAT(' Enter Damping Factor F (Default is 0.25) for',&
   &' Newton-Raphson Search.',/, ' The Value of F for',&
   &' Nuclear Searches Will be (F/10.0)  ',$)
1304 FORMAT(' Newton-Raphson Damping Factor Set to ',1PE16.8,&
   &' for Non-Nuclear',/, ' Critical Points and ',1PE16.8,&
   &' for Nuclear Critical Points')
1310 FORMAT(I4,1X,3(1PE16.8),2A8,2X,'(',I1,',',I2,')')
!
   CALL MAKNAME(1,WFN,ILEN,FWFN)
   IF (ILEN .EQ. 0) STOP 'usage:  extreme wfnfile '
!     Generate output filename from input filename
   CALL GETARG(1,WLP)
   CALL MAKNAME_DIRECT(WLP,ILEN,FWLP)
   IF (ILEN .EQ. 0) STOP 'usage:  extreme wfnfile '
!
   OPEN (IWFN,FILE=WFN,status='unknown')
   OPEN (IWLP,FILE=WLP,status='unknown')
!
   INP=0
!
   CALL RDWFN
!
   WRITE (IOUT,1010) WFNTTL
   WRITE (IOUT,1000)
   READ (INPT,1010) JOBTTL
   WRITE (IWLP,1020)
   WRITE (IWLP,1010) WFNTTL
   WRITE (IWLP,1010) JOBTTL
!
   Write(IWLP,1030)
   DO 10 I = 1,NCENT
      WRITE(IWLP,1040) ATNAM(I),CO(IXC+I),CO(IYC+I),CO(IZC+I)
10 CONTINUE
!
   IFUNC=1
   WRITE(IWLP,1070)
   WRITE(IOUT,1070)
   WRITE(IOUT,*)
!
30 WRITE (IOUT,1120)
   READ (INPT,*) IST
!
   IF (IST .EQ. 0) THEN
!
      write(iout,1121)
      Read(Inpt,*) nucarb
      If(nucarb.eq.0)Then
         WRITE (IOUT,1130)
         READ (INPT,*) XYZ(1),XYZ(2),XYZ(3)
         If(Iprint.eq.1)WRITE (IOUT,1140) XYZ(1),XYZ(2),XYZ(3)
         WRITE (IOUT,*)
         WRITE (IWLP,*)
         WRITE (IWLP,1140) XYZ(1),XYZ(2),XYZ(3)
         J1 = 0
         J2 = 0
         Inuc=0
      ElseIf(nucarb.eq.1)Then
         J2=0
         Inuc=1
         write(iout,1122)
         Read(Inpt,*)J1
         XYZ(1)=CO(IXC+J1)
         XYZ(2)=CO(IYC+J1)
         XYZ(3)=CO(IZC+J1)
         If(Iprint.eq.1)WRITE (IOUT,1141) XYZ(1),XYZ(2),XYZ(3)
         WRITE (IOUT,*)
         WRITE (IWLP,*)
         WRITE (IWLP,1141) XYZ(1),XYZ(2),XYZ(3)
      Endif
      IWHOLE=0
      Ifail=0
      CALL NEWTON(XYZ,IFAIL,IFUNC,IWHOLE,NITER,Inuc)
      If(Ifail.eq.0)Then
         Icrit=1
         CALL STORE(XYZ,X,Y,Z,J1,J2,ICancl,INP,Ji1,Ji2,NITER)
         If(Icancl.eq.0)Then
            Call Props(XYZ,IFUNC,ICrit,Iwhole,J1,J2,INP,EVSAVE)
         Endif
      Endif
      Goto 30
!
   ELSE IF (IST .EQ. 1) THEN
!
      WRITE (IOUT,1150)
40    READ (INPT,*) J1,J2
      IF(J1.GT.NCENT.OR.J2.GT.NCENT) THEN
         WRITE (IOUT, 1160)
         GOTO 40
      ENDIF
      Write(iout,1170)
      Read(inpt,*) IHALF
      Fract=pt5
      If(IHALF.eq.0)Then
         Write(Iout,1180) J1, J2
         Read(Inpt,*)Fract
      Endif
      WRITE (IOUT,1190) J1,J2
      WRITE (IWLP,1190) J1,J2
      XYZ(1) = CO(IXC+J1)*fract+CO(IXC+J2)*(one-fract)
      XYZ(2) = CO(IYC+J1)*fract+CO(IYC+J2)*(one-fract)
      XYZ(3) = CO(IZC+J1)*fract+CO(IZC+J2)*(one-fract)
      IWHOLE=0
      Ifail=0
      If(Iprint.eq.1)WRITE (IOUT,1140) XYZ(1),XYZ(2),XYZ(3)
      WRITE (IWLP,1140) XYZ(1),XYZ(2),XYZ(3)
      CALL NEWTON(XYZ,IFAIL,IFUNC,IWHOLE,NITER,0)
      If(Ifail.eq.0)Then
         Icrit=1
         CALL STORE(XYZ,X,Y,Z,J1,J2,ICancl,INP,Ji1,Ji2,NITER)
         If(Icancl.eq.0)Then
            Call Props(XYZ,IFUNC,ICrit,Iwhole,J1,J2,INP,EVSAVE)
         Endif
      Endif
      Goto 30
!
   ELSE IF (IST .EQ. 2) THEN
!
      WRITE (IOUT,1200)
50    READ (INPT,*) J1,J2,J3
      IF(J1.GT.NCENT.OR.J2.GT.NCENT.OR.J3.GT.NCENT) THEN
         WRITE (IOUT, 1160)
         GOTO 50
      ENDIF
      WRITE (IOUT,1210) J1,J2,J3
      WRITE (IWLP,1210) J1,J2,J3
      XYZ(1) = (CO(IXC+J1)+CO(IXC+J2)+CO(IXC+J3))/THREE
      XYZ(2) = (CO(IYC+J1)+CO(IYC+J2)+CO(IYC+J3))/THREE
      XYZ(3) = (CO(IZC+J1)+CO(IZC+J2)+CO(IZC+J3))/THREE
      J1=0
      J2=0
      IWHOLE=0
      IFAIL=0
      If(Iprint.eq.1)WRITE (IOUT,1140) XYZ(1),XYZ(2),XYZ(3)
      WRITE (IWLP,1140) XYZ(1),XYZ(2),XYZ(3)
      CALL NEWTON(XYZ,IFAIL,IFUNC,IWHOLE,NITER,0)
      If(IFAIL.eq.0)Then
         Icrit=1
         CALL STORE(XYZ,X,Y,Z,J1,J2,ICancl,INP,Ji1,Ji2,NITER)
         If(Icancl.eq.0)Then
            Call Props(XYZ,IFUNC,ICrit,Iwhole,J1,J2,INP,EVSAVE)
         Endif
      Endif
      Goto 30
!
   ELSE IF (IST .EQ. 3) THEN
!
      WRITE (IOUT,1220)
60    READ (INPT,*) J1,J2,J3,J4
      IF(J1.GT.NCENT.OR.J2.GT.NCENT.OR.J3.GT.NCENT&
      &.OR.J4.GT.NCENT) THEN
         WRITE (IOUT, 1160)
         GOTO 60
      ENDIF
      WRITE (IOUT,1230) J1,J2,J3,J4
      WRITE (IWLP,1230) J1,J2,J3,J4
      XYZ(1) = (CO(IXC+J1)+CO(IXC+J2)+CO(IXC+J3)+CO(IXC+J4))/FOUR
      XYZ(2) = (CO(IYC+J1)+CO(IYC+J2)+CO(IYC+J3)+CO(IYC+J4))/FOUR
      XYZ(3) = (CO(IZC+J1)+CO(IZC+J2)+CO(IZC+J3)+CO(IZC+J4))/FOUR
      J1=0
      J2=0
      IWHOLE=0
      IFAIL=0
      If(Iprint.eq.1)WRITE (IOUT,1140) XYZ(1),XYZ(2),XYZ(3)
      WRITE (IWLP,1140) XYZ(1),XYZ(2),XYZ(3)
      CALL NEWTON(XYZ,IFAIL,IFUNC,IWHOLE,NITER,0)
      If(IFAIL.eq.0)Then
         Icrit=1
         CALL STORE(XYZ,X,Y,Z,J1,J2,ICancl,INP,Ji1,Ji2,NITER)
         If(Icancl.eq.0)Then
            Call Props(XYZ,IFUNC,ICrit,Iwhole,J1,J2,INP,EVSAVE)
         Endif
      Endif
      Goto 30
!
   ELSE IF (IST .EQ. 4) THEN
!
      CALL ALPHA
      GOTO 30
!
   ELSE IF (IST .EQ. 5) THEN
!
      WRITE (IOUT,1240)
      READ (INPT,*) (XYZ(II),II=1,3)
      IF(IPRINT.eq.1)WRITE (IOUT,1250) (XYZ(II),II=1,3)
      WRITE (IWLP,1250) (XYZ(II),II=1,3)
      IWHOLE=0
      Icrit=0
      J1=0
      J2=0
      Call Props(XYZ,IFUNC,ICrit,Iwhole,J1,J2,INP,EVSAVE)
      GOTO 30
!
   ELSEIF (IST.EQ.6)Then
!
71    WRITE(IOUT,1260)
      READ(INPT,*) IWHOLE
!
      If(Iwhole.eq.1.or.Iwhole.eq.5)Then
         Do 75 IAtom=1,Ncent
            XYZ(1) = CO(IXC+Iatom)
            XYZ(2) = CO(IYC+Iatom)
            XYZ(3) = CO(IZC+Iatom)
            Ifail=0
            J1=Iatom
            J2=0
            If(Iprint.eq.1)WRITE (IOUT,1140) XYZ(1),XYZ(2),XYZ(3)
            If(Iprint.eq.1)WRITE (IWLP,1140) XYZ(1),XYZ(2),XYZ(3)
            CALL NEWTON(XYZ,IFAIL,IFUNC,IWHOLE,NITER,1)
            If(IFAIL.eq.0)Then
               Icrit=1
               CALL STORE(XYZ,X,Y,Z,J1,J2,ICancl,INP,Ji1,Ji2,NITER)
               If(Icancl.eq.0)Then
                  Call Props(XYZ,IFUNC,ICrit,Iwhole,J1,J2,INP,EVSAVE)
               Endif
            Endif
75       Continue
         If(Iwhole.eq.1)Goto 71
      Endif
!
      If(Iwhole.eq.2.or.Iwhole.eq.5)Then
         Do 80 IAtom=1,Ncent
            Do 90 Katom=1,Iatom-1
               J1 = IATOM
               J2 = KATOM
               WRITE (IWLP,1190) J1,J2
               XYZ(1) = (CO(IXC+J1)+CO(IXC+J2))/TWO
               XYZ(2) = (CO(IYC+J1)+CO(IYC+J2))/TWO
               XYZ(3) = (CO(IZC+J1)+CO(IZC+J2))/TWO
               Ifail=0
               If(Iprint.eq.1)WRITE (IOUT,1140) XYZ(1),XYZ(2),XYZ(3)
               If(Iprint.eq.1)WRITE (IWLP,1140) XYZ(1),XYZ(2),XYZ(3)
               CALL NEWTON(XYZ,IFAIL,IFUNC,IWHOLE,NITER,0)
               If(IFAIL.eq.0)Then
                  Icrit=1
                  CALL STORE(XYZ,X,Y,Z,J1,J2,ICancl,INP,Ji1,Ji2,NITER)
                  If(Icancl.eq.0)Then
                     Call Props(XYZ,IFUNC,ICrit,Iwhole,J1,J2,INP,EVSAVE)
                  Endif
               Endif
90          Continue
80       Continue
         If(Iwhole.eq.2)Goto 71
      Endif
!
      If(IWhole.eq.3.or.Iwhole.eq.5)Then
         Do 100 IAtom=1,Ncent
            Do 110 Jatom=1,Iatom-1
               Do 120 Katom=1,Jatom-1
                  J1 = Iatom
                  J2 = Jatom
                  J3 = Katom
                  WRITE (IWLP,1210) J1,J2,J3
                  XYZ(1) = (CO(IXC+J1)+CO(IXC+J2)+CO(IXC+J3))/THREE
                  XYZ(2) = (CO(IYC+J1)+CO(IYC+J2)+CO(IYC+J3))/THREE
                  XYZ(3) = (CO(IZC+J1)+CO(IZC+J2)+CO(IZC+J3))/THREE
                  J1=0
                  J2=0
                  Ifail=0
                  If(Iprint.eq.1)WRITE (IOUT,1140) XYZ(1),XYZ(2),XYZ(3)
                  If(Iprint.eq.1)WRITE (IWLP,1140) XYZ(1),XYZ(2),XYZ(3)
                  CALL NEWTON(XYZ,IFAIL,IFUNC,IWHOLE,NITER,0)
                  If(IFAIL.eq.0)Then
                     Icrit=1
                     CALL STORE(XYZ,X,Y,Z,J1,J2,ICancl,INP,Ji1,Ji2,NITER)
                     If(Icancl.eq.0)Then
                        Call Props(XYZ,IFUNC,ICrit,Iwhole,J1,J2,INP,EVSAVE)
                     Endif
                  Endif
120            Continue
110         Continue
100      Continue
         If(Iwhole.eq.3)Goto 71
      Endif
!
      If(Iwhole.eq.4.or.iwhole.eq.5)Then
         Do 130 IAtom=1,Ncent
            Do 140 Jatom=1,Iatom-1
               Do 150 Katom=1,Jatom-1
                  Do 160 Latom=1,Katom-1
                     J1 = Iatom
                     J2 = Jatom
                     J3 = Katom
                     J4 = Latom
                     WRITE (IWLP,1230) J1,J2,J3,J4
                     XYZ(1) = (CO(IXC+J1)+CO(IXC+J2)+CO(IXC+J3)+CO(IXC+J4))/FOUR
                     XYZ(2) = (CO(IYC+J1)+CO(IYC+J2)+CO(IYC+J3)+CO(IYC+J4))/FOUR
                     XYZ(3) = (CO(IZC+J1)+CO(IZC+J2)+CO(IZC+J3)+CO(IYC+J4))/FOUR
                     J1=0
                     J2=0
                     Ifail=0
                     If(Iprint.eq.1)WRITE (IOUT,1140) XYZ(1),XYZ(2),XYZ(3)
                     If(Iprint.eq.1)WRITE (IWLP,1140) XYZ(1),XYZ(2),XYZ(3)
                     CALL NEWTON(XYZ,IFAIL,IFUNC,IWHOLE,NITER,0)
                     If(IFAIL.eq.0)Then
                        Icrit=1
                        CALL STORE(XYZ,X,Y,Z,J1,J2,ICancl,INP,Ji1,Ji2,NITER)
                        If(Icancl.eq.0)Then
                           Call Props(XYZ,IFUNC,ICrit,Iwhole,J1,J2,INP,EVSAVE)
                        Endif
                     Endif
160               Continue
150            Continue
140         Continue
130      Continue
         If(Iwhole.eq.4)Goto 71
      Endif
!
72    Write(Iout,1265)
      Read(Inpt,*) ILine
      Iwhole=6
      If(Iline.eq.1.or.iline.eq.2)Then
         If(ILine.eq.1)Then
            Write(iout,1266)
            Read(Inpt,*) J1,J2,Npts
            Dist=dsqrt((CO(IXC+J1)-CO(IXC+J2))**2+&
            &(CO(IYC+J1)-CO(IYC+J2))**2+&
            &(CO(IZC+J1)-CO(IZC+J2))**2)
            Unitx = (CO(IXC+J2)-CO(IXC+J1))/Dist
            Unity = (CO(IYC+J2)-CO(IYC+J1))/Dist
            Unitz = (CO(IZC+J2)-CO(IZC+J1))/Dist
            X0=CO(IXC+J1)
            Y0=CO(IYC+J1)
            Z0=CO(IZC+J1)
         ElseIf(Iline.eq.2)Then
            Write(Iout,1267)
            Read(Inpt,*) X1,Y1,Z1
            Write(Iout,1268)
            Read(Inpt,*) X2,Y2,Z2,NPts
            Dist=Dsqrt((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2)
            Unitx=(X2-X1)/Dist
            Unity=(Y2-Y1)/Dist
            Unitz=(Z2-Z1)/Dist
            X0=X1
            Y0=Y1
            Z0=Z1
         Endif
         Dinc=Dist/NPts
         Do 162 I=1,Npts
            XYZ(1)= X0 + I*Dinc*UNITX
            XYZ(2)= Y0 + I*Dinc*UNITY
            XYZ(3)= Z0 + I*Dinc*UNITZ
            J1=0
            J2=0
            Ifail=0
            If(Iprint.eq.1)WRITE (IOUT,1140) XYZ(1),XYZ(2),XYZ(3)
            If(Iprint.eq.1)WRITE (IWLP,1140) XYZ(1),XYZ(2),XYZ(3)
            CALL NEWTON(XYZ,IFAIL,IFUNC,IWHOLE,NITER,0)
            If(IFAIL.eq.0)Then
               Icrit=1
               CALL STORE(XYZ,X,Y,Z,J1,J2,ICancl,INP,Ji1,Ji2,NITER)
               If(Icancl.eq.0)Then
                  Call Props(XYZ,IFUNC,ICrit,Iwhole,J1,J2,INP,EVSAVE)
               Endif
            Endif
162      Continue
         Goto 72
      Endif
!
73    Write(IOUT,1270)
      READ(INPT,*) Isphere
      If(ISphere.ne.0)Then
         Iwhole=7
         Pi=Dacos(-one)
         Write(Iout,1280)
         Read(inpt,*)Iall
         Write(Iout,1290)
         Read(inpt,*)Iphi,ITheta,RadMin,Radmax,Nrad
         If(iall.ne.0)Then
            iatom=iall
            natoms=1
         Else
            iatom=1
            natoms=ncent
         Endif
         Phinc=Two*pi/Iphi
         HPhinc=Pt5*Phinc
         Thinc=Two*pi/ITheta
         HThinc=Pt5*Thinc
         Radinc=Radmax/Nrad
         Do 170 jatom=1,natoms
            Do 180 I=1,Iphi
               Phi=Phinc*I-HPhinc
               Do 190 J=1,Itheta
                  Theta=THinc*J-HThinc
                  unitx=dsin(theta)*dcos(phi)
                  unity=dsin(theta)*dsin(phi)
                  unitz=dcos(theta)
                  Do 200 K=1,Nrad
                     XYZ(1)=CO(IXC+jatom+iatom-1)+(RadMin+RADINC*K)*UNITX
                     XYZ(2)=CO(IYC+jatom+iatom-1)+(RadMin+RADINC*K)*UNITY
                     XYZ(3)=CO(IZC+jatom+iatom-1)+(RadMin+RADINC*K)*UNITZ
                     J1=0
                     J2=0
                     Ifail=0
                     CALL NEWTON(XYZ,IFAIL,IFUNC,IWHOLE,NITER,0)
                     If(IFAIL.eq.0)Then
                        Icrit=1
                        CALL STORE(XYZ,X,Y,Z,J1,J2,ICancl,INP,Ji1,Ji2,NITER)
                        If(Icancl.eq.0)Then
                           Call Props(XYZ,IFUNC,ICrit,Iwhole,J1,J2,INP,EVSAVE)
                        Endif
                     Endif
200               Continue
190            Continue
180         Continue
170      Continue
         Goto 73
      Endif
      Goto 30
!
   ELSEIF(IST.EQ.7)Then
!
220   Write(iout,*)
      Write(iout,1295)
      Read(inpt,*) Iopt
      If(Iopt.eq.1)Then
         Write(Iwlp,*)
         Write(iout,*)
         Write(iout,1296)
         Read(inpt,*)ICut
         Cut=Ten**(-Icut)
         Write(iwlp,1298)cut
         Goto 220
      ElseIf(Iopt.eq.2)Then
         Write(Iwlp,*)
         Write(iout,*)
         Write(iout,1297)
         Read(Inpt,*) IPrint
         Write(iwlp,1299)IPrint
         Goto 220
      ElseIf(Iopt.eq.3)Then
         Write(Iwlp,*)
         Write(iout,*)
         Write(iout,1301)
         Read(Inpt,*) Ieps
         eps=Ten**(-ieps)
         iepnuc=ieps-2
         epsnuc=Ten**(-iepnuc)
         Write(iwlp,1302)eps,epsnuc
         Goto 220
      ElseIf(Iopt.eq.4)Then
         Write(Iwlp,*)
         Write(iout,*)
         Write(iout,1303)
         Read(Inpt,*) Dmp
         DmpNuc=Dmp/Ten
         Write(iwlp,1304)Dmp,DmpNuc
         Goto 220
      Endif
      Write(iout,*)
      Goto 30
!
   ELSE IF ((IST .EQ. 8).or.(Ist.eq.9)) THEN
!
      WRITE (IOUT,1300)
      If(IST.eq.9)WRITE (IWLP,1300)
      DO 210 I = 1,INP
         ev1=evsave(i,1)
         ev2=evsave(i,2)
         ev3=evsave(i,3)
         aev1=dabs(ev1)
         aev2=dabs(ev2)
         aev3=dabs(ev3)
         ineg=0
         ipos=0
         If(ev1.lt.zero.and.aev1.gt.degen)ineg=ineg+1
         If(ev2.lt.zero.and.aev2.gt.degen)ineg=ineg+1
         If(ev3.lt.zero.and.aev3.gt.degen)ineg=ineg+1
         If(ev1.gt.zero.and.aev1.gt.degen)ipos=ipos+1
         If(ev2.gt.zero.and.aev2.gt.degen)ipos=ipos+1
         If(ev3.gt.zero.and.aev3.gt.degen)ipos=ipos+1
         irank=ineg+ipos
         isig=ipos-ineg
         IF(JI1(I).NE.0.AND.JI2(I).NE.0) THEN
            WRITE (IOUT,1310) I,X(I),Y(I),Z(I),ATNAM(JI1(I)),&
            &ATNAM(JI2(I)),irank,isig
            If(Ist.eq.9)WRITE (IWLP,1310) I,X(I),Y(I),Z(I),ATNAM(JI1(I)),&
            &ATNAM(JI2(I)),irank,isig
         ELSEIF(JI1(I).NE.0.AND.JI2(I).EQ.0) THEN
            WRITE (IOUT,1310) I,X(I),Y(I),Z(I),ATNAM(JI1(I)),&
            &BLANK,irank,isig
            If(Ist.eq.9)WRITE (IWLP,1310) I,X(I),Y(I),Z(I),ATNAM(JI1(I)),&
            &BLANK,irank,isig
         ELSE
            WRITE (IOUT,1310) I,X(I),Y(I),Z(I),BLANK,BLANK,irank,isig
            If(Ist.eq.9)WRITE (IWLP,1310) I,X(I),Y(I),Z(I),BLANK,BLANK,&
            &irank,isig
         ENDIF
210   CONTINUE
      If(Ist.eq.8)GOTO 30
!
      STOP ' EXTREME IS DONE '
!
   ELSE
      WRITE (IOUT,1100)
      GOTO 30
   END IF
!
END
