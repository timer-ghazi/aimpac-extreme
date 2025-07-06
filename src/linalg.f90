module linalg
  implicit none
contains
   SUBROUTINE      TQLGRM	(N, D, E, Z, IERR)
      IMPLICIT        DOUBLE PRECISION (A-H, O-Z)
      DIMENSION       D(*), E(*), Z(N,N)
      PARAMETER (AMACH = 16.0D-13)
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
!
      IERR    = 0
      IF (N .EQ. 1) RETURN
!
      DO 30 I = 2,N
         E(I-1)  = E(I)
30    CONTINUE

      F       = ZERO
      B       = ZERO
      E(N)    = ZERO

      DO 31 L = 1,N
         J       = 0
         H       = AMACH*(DABS(D(L)) + DABS(E(L)))
         IF (B .LT. H) B = H

105      DO 32 M = L,N
            IF (DABS(E(M)) .LE. B) GOTO 120
32       CONTINUE

120      IF (M .EQ. L) GOTO 220

130      IF (J .EQ. 30) THEN
            IERR    = L
            RETURN
         END IF

         J       = J + 1
         L1      = L + 1
         G       = D(L)
         P       = (D(L1) - G)/(2*E(L))
         IF (DABS(P*AMACH) .GT. ONE) THEN
            R       = P
         ELSE
            R       = DSQRT(P*P + 1)
         END IF
         D(L)    = E(L)/(P + DSIGN(R,P))
         H       = G - D(L)

         DO 33 I = L1,N
            D(I)    = D(I) - H
33       CONTINUE

         F       = F + H
         P       = D(M)
         C       = ONE
         S       = ZERO
         MML     = M - L

         DO 34 II = 1,MML
            I       = M - II
            G       = C*E(I)
            H       = C*P
            IF (DABS(P) .GE. DABS(E(I))) THEN
               C       = E(I)/P
               R       = DSQRT(C*C + 1)
               E(I+1)  = S*P*R
               S       = C/R
               C       = ONE/R
            ELSE
               C       = P/E(I)
               R       = DSQRT(C*C + 1)
               E(I+1)  = S*E(I)*R
               S       = 1.D0/R
               C       = C*S
            END IF
            P       = C*D(I) - S*G
            D(I+1)  = H + S*(C*G + S*D(I))

            DO 35 K = 1,N
               H       = Z(K,I+1)
               Z(K,I+1)= S*Z(K,I) + C*H
               Z(K,I)  = C*Z(K,I) - S*H
35          CONTINUE

34       CONTINUE

         E(L)    = S*P
         D(L)    = C*P
         IF (DABS(E(L)) .GT. B) GOTO 130

220      D(L)    = D(L) + F
31    CONTINUE

      DO 300 II = 2,N
         I       = II - 1
         K       = I
         P       = D(I)

         DO 260 J = II,N
            IF (D(J) .GE. P) GOTO 260
            K       = J
            P       = D(J)
260      CONTINUE

         IF (K .EQ. I) GOTO 300
         D(K)    = D(I)
         D(I)    = P

         DO 37 J = 1,N
            P       = Z(J,I)
            Z(J,I)  = Z(J,K)
            Z(J,K)  = P
37       CONTINUE

300   CONTINUE
      RETURN
   END
   SUBROUTINE	TRACE	(H, E, W, N, IERR)
!
! TRACE CALLS TREDIG AND TLQGRM TO DIAGONALIZE A SYMMETRIC REAL MATRIX.
! THE MATRIX IS PASSED DOWN IN H AND IS REPLACED BY THE EIGENVECTORS.
! THE EIGENVALUES IN E ARE STORED SMALLEST FIRST.
! THE WORK STORE W SHOULD BE AT LEAST OF DIMENSION N.
! SKK ==================================================================

      IMPLICIT        DOUBLE PRECISION (A-H, O-Z)
      DIMENSION       H(N,N), E(*), W(*)
!
      CALL TREDIG	(N, E, W, H)
      CALL TQLGRM	(N, E, W, H, IERR)
!
      RETURN
   END
   SUBROUTINE      TREDIG	(N, D, E, Z)
      IMPLICIT        DOUBLE PRECISION (A-H, O-Z)
      DIMENSION       D(*), E(*), Z(N,N)
      PARAMETER	(ZERO = 0.0D0, ONE = 1.0D0)

      IF (N .EQ. 1) GOTO 320

      DO 30 II = 2,N
         I       = N + 2 - II
         L       = I - 1
         H       = ZERO
         SCALE   = ZERO

         IF (L .LT. 2) GOTO 130

         DO 31 K = 1,L
            SCALE   = SCALE + DABS(Z(I,K))
31       CONTINUE

         IF (SCALE .NE. ZERO) GOTO 140
130      E(I)    = Z(I,L)
         GOTO 290

140      RSCALE	= ONE/SCALE
         DO 32 K = 1,L
            Z(I,K)  = Z(I,K)*RSCALE
            H       = H + Z(I,K)*Z(I,K)
32       CONTINUE
         F       = Z(I,L)
         G       = -DSIGN(DSQRT(H),F)
         E(I)    = SCALE*G
         H       = H - F*G
         Z(I,L)  = F - G
         F       = ZERO
         RH	= ONE/H
         RHSCALE	= RH*RSCALE

         DO 33 J = 1,L
            Z(J,I)  = Z(I,J)*RHSCALE
            G       = ZERO

            DO 34 K = 1,J
               G       = G + Z(J,K)*Z(I,K)
34          CONTINUE

            JP1     = J + 1
            IF (L .LT. JP1) GOTO 220

            DO 35 K = JP1,L
               G       = G + Z(K,J)*Z(I,K)
35          CONTINUE

220         E(J)    = G*RH
            F       = F + E(J)*Z(I,J)
33       CONTINUE

         HH      = F/(H + H)

         DO 36 J = 1,L
            F       = Z(I,J)
            G       = E(J) - HH*F
            E(J)    = G
            DO 37 K = 1,J
               Z(J,K)  = Z(J,K) - F*E(K) - G*Z(I,K)
37          CONTINUE
36       CONTINUE

         DO 38 K	= 1,L
            Z(I,K)  =  SCALE*Z(I,K)
38       CONTINUE

290      D(I)    = H
30    CONTINUE

320   D(1)    = ZERO
      E(1)    = ZERO

      DO 500 I = 1,N
         L       = I - 1
         IF (D(I) .EQ. ZERO) GOTO 380

         DO 40 J	= 1,L
            G       = ZERO

            DO 41 K	= 1,L
               G       = G + Z(I,K)*Z(K,J)
41          CONTINUE

            DO 42 K = 1,L
               Z(K,J)  = Z(K,J) - G*Z(K,I)
42          CONTINUE

40       CONTINUE

380      D(I)    = Z(I,I)
         Z(I,I)  = ONE
         IF(L .LT. 1) GOTO 500

         DO 43 J	= 1,L
            Z(J,I)  = ZERO
            Z(I,J)  = ZERO
43       CONTINUE

500   CONTINUE
      RETURN
   END
end module linalg
