      SUBROUTINE UTRAN(Readok)
C***********************************************************************
C READ TRANSPORT PROPERTIES FORM I/O UNIT 7 IN RECORD FORMAT AND WRITE
C UNFORMATTED ON I/O UNIT IOTRN.  USES SCRATCH I/O UNIT IOSCH.
C
C UTRAN IS CALLED FROM SUBROUTINE INPUT AFTER A RECORD WITH 'tran'
C IN COLUMNS 1-4 HAS BEEN READ.
C
C NOTE:  THIS ROUTINE MAY BE CALLED DIRECTLY  AND USED BY ITSELF TO
C PROCESS THE TRANSPORT PROPERTY DATA.
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C DUMMY ARGUMENTS
      LOGICAL Readok
C LOCAL VARIABLES
      CHARACTER*16 tname(2)
      CHARACTER*1 vorc
      INTEGER i,ic,in,iv,j,k,ncc,nn,ns,nv
      REAL*8 cc,tc(36),tcin(6),trcoef(6,3,2),vvl
C
      EQUIVALENCE (tc(1),trcoef(1,1,1))
      ns = 0
      REWIND IOSCH
 100  DO i = 1,36
        tc(i) = 0.
      ENDDO
      READ (IOINP,99001) tname,vvl,nv,cc,ncc
      IF ( tname(1).EQ.'end'.OR.tname(1).EQ.'LAST' ) THEN
        WRITE (IOTRN) ns
        REWIND IOSCH
        DO i = 1,ns
          READ (IOSCH,ERR=200) tname,trcoef
          WRITE (IOTRN) tname,trcoef
        ENDDO
        GOTO 300
      ELSE
        ic = 0
        iv = 0
        nn = nv + ncc
        IF ( nv.LE.3.AND.ncc.LE.3 ) THEN
          DO in = 1,nn
            READ (IOINP,99002) vorc,tcin
            IF ( vorc.EQ.'C' ) THEN
              k = 2
              ic = ic + 1
              j = ic
            ELSE
              k = 1
              iv = iv + 1
              j = iv
            ENDIF
            IF ( j.GT.3 ) GOTO 200
            DO i = 1,6
              trcoef(i,j,k) = tcin(i)
            ENDDO
          ENDDO
          ns = ns + 1
          WRITE (IOSCH) tname,trcoef
          GOTO 100
        ENDIF
      ENDIF
 200  WRITE (IOOUT,99003) tname
      Readok = .FALSE.
 300  RETURN
99001 FORMAT (2A16,2X,A1,I1,A1,I1)
99002 FORMAT (1X,A1,2F9.2,4E15.8)
99003 FORMAT (/' ERROR IN PROCESSING trans.inp AT OR NEAR (UTRAN)',/1X,
     &        2A16)
      END
