module utils
  implicit none
contains
SUBROUTINE MAKNAME(I,STRING,L,EXT)
   CHARACTER*(*) STRING,EXT
   INTEGER I,J,L,N,EXTLEN,DOT_POS

   CALL GETARG(I,STRING)
   J = LEN(STRING)
   EXTLEN = LEN(EXT)

!     Find the length of the actual argument (before spaces)
   L = 0
   DO 10 N = 1,J
      IF(STRING(N:N) .EQ. ' ') THEN
         L = N - 1
         GOTO 20
      ENDIF
10 CONTINUE
   L = J

20 CONTINUE
   IF (L .EQ. 0) RETURN

!     Check if filename already has the desired extension
   IF (L .GE. EXTLEN) THEN
      IF (STRING(L-EXTLEN+1:L) .EQ. EXT) THEN
!         Already has correct extension, just trim spaces
         STRING = STRING(1:L)
         RETURN
      ENDIF
   ENDIF

!     Check if filename has any extension (contains a dot)
   DOT_POS = 0
   DO 30 N = L,1,-1
      IF (STRING(N:N) .EQ. '.') THEN
         DOT_POS = N
         GOTO 40
      ENDIF
30 CONTINUE

40 CONTINUE
   IF (DOT_POS .GT. 0) THEN
!       Has an extension, keep the full filename as-is
      STRING = STRING(1:L)
   ELSE
!       No extension, add the default one
      STRING = STRING(1:L)//EXT
      L = L + EXTLEN
   ENDIF

   RETURN
END
   SUBROUTINE MAKNAME_DIRECT(STRING,L,EXT)
      CHARACTER*(*) STRING,EXT

      INTEGER L,N,J,EXTLEN,DOT_POS

      J = LEN(STRING)
      EXTLEN = LEN(EXT)

!     Find the length of the actual argument (before spaces)
      L = 0
      DO 10 N = 1,J
         IF(STRING(N:N) .EQ. ' ') THEN
            L = N - 1
            GOTO 20
         ENDIF
10    CONTINUE
      L = J

20    CONTINUE
      IF (L .EQ. 0) RETURN

!     Check if filename already has the desired extension
      IF (L .GE. EXTLEN) THEN
         IF (STRING(L-EXTLEN+1:L) .EQ. EXT) THEN
!         Already has correct extension, just trim spaces
            STRING = STRING(1:L)
            RETURN
         ENDIF
      ENDIF

!     Check if filename has any extension (contains a dot)
      DOT_POS = 0
      DO 30 N = L,1,-1
         IF (STRING(N:N) .EQ. '.') THEN
            DOT_POS = N
            GOTO 40
         ENDIF
30    CONTINUE

40    CONTINUE
      IF (DOT_POS .GT. 0) THEN
!       Replace existing extension with new one
         STRING = STRING(1:DOT_POS-1)//EXT
         L = DOT_POS - 1 + EXTLEN
      ELSE
!       No extension, add the default one
         STRING = STRING(1:L)//EXT
         L = L + EXTLEN
      ENDIF

      RETURN
   END
end module utils
