      SUBROUTINE VARFMT(Vx)
C***********************************************************************
C SET DECIMAL PLACES ACCORDING TO NUMBER SIZE FOR F-FORMAT IN
C VARIABLE FORMAT FMT.
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C DUMMY ARGUMENTS
      REAL*8 Vx(NCOL)
C LOCAL VARIABLES
      INTEGER i,k
      REAL*8 vi
      REAL*8 DABS
      SAVE i,k,vi
C
      DO i = 1,Npt
        vi = DABS(Vx(i))
        k = 2*i + 3
        Fmt(k) = '5,'
        IF ( vi.GE.1. ) Fmt(k) = '4,'
        IF ( vi.GE.10. ) Fmt(k) = '3,'
        IF ( vi.GE.100. ) Fmt(k) = '2,'
        IF ( vi.GE.10000. ) Fmt(k) = '1,'
        IF ( vi.GE.1000000. ) Fmt(k) = '0,'
      ENDDO
      Fmt(29)(2:) = ' '
      END
