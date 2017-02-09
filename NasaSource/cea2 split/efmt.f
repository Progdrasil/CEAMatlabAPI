      SUBROUTINE EFMT(Fone,Aa,Vx)
C***********************************************************************
C WRITE OUTPUT RECORD WITH NUMERICAL VALUES IN SPECIAL EXPONENT FORM.
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C DUMMY ARGUMENTS
      CHARACTER*15 Aa
      CHARACTER*4 Fone
      REAL*8 Vx(MAXMAT)
C LOCAL VARIABLES
      CHARACTER*4 fmix(5),frmt(8)
      INTEGER i,j,j1,ne(NCOL)
      INTEGER IABS
      REAL*8 ee,fe,w(NCOL)
      REAL*8 DABS,DLOG10
      SAVE ee,fe,i,j,j1,ne,w
C
      DATA frmt/'(1H ',',A15',',','9X,','13(F','6.4,','I2,','1X))'/
      DATA fmix/'I3,','6.4,','I2,','9X,','5.3,'/
      frmt(6) = fmix(2)
      frmt(7) = fmix(3)
      j1 = 1
      frmt(4) = '1x,'
      IF ( Fone.EQ.'9X,' ) THEN
        j1 = 2
        frmt(4) = fmix(4)
      ENDIF
      DO i = j1,Npt
        IF ( Vx(i).NE.0. ) THEN
          ee = DLOG10(DABS(Vx(i)))
          ne(i) = ee
          fe = ne(i)
          IF ( ee.LT.-.2181E-05.AND.fe.NE.ee ) ne(i) = ne(i) - 1
          IF ( IABS(ne(i)).GE.10 ) THEN
            frmt(6) = fmix(5)
            frmt(7) = fmix(1)
          ENDIF
          w(i) = Vx(i)/10.**ne(i)
        ELSE
          w(i) = 0.
          ne(i) = 0.
        ENDIF
      ENDDO
      WRITE (IOOUT,frmt) Aa,(w(j),ne(j),j=j1,Npt)
      END
