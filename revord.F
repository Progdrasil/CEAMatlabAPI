#include "fintrf.h"
C
#if 0
C     
C     revord.F
C     .F file needs to be preprocessed to generate .for equivalent
C     
#endif
C     
C     revord.f
C     Example for illustrating how to copy the string data from 
C     MATLAB to a FORTRAN-style string and back again 
C
C     Takes a row vector string and returns a string in reverse order.
C
C	Limit of 100 characters. 
C
C     This is a MEX-file for MATLAB.
C     Copyright 1984-2011 The MathWorks, Inc.
C     
C======================================================================
C     Computational routine

      subroutine revord(input_buf, strlen, output_buf)
      character  input_buf(*), output_buf(*)
      mwSize    i, strlen     
      do 10 i=1,strlen
         output_buf(i) = input_buf(strlen-i+1)
 10   continue
      return
      end

C-----------------------------------------------------------------------
C     Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

C     Declarations
      implicit none

C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs
 
C     Function declarations:
      mwPointer mxCreateString
      mwPointer mxGetString       
      mwPointer mxGetM, mxGetN
      integer mxIsChar
 
	  integer status
 
C     Arguments for computational routine:
C     Maximum string size = maxbuf
      mwSize  maxbuf
      parameter(maxbuf = 100)
      character*100 input_buf, output_buf 
      mwPointer  strlen

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if (nrhs .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:revord:nInput',
     +                           'One input required.')
      elseif (nlhs .gt. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:revord:nOutput',
     +                           'Too many output arguments.')
C     The input must be a string.
      elseif(mxIsChar(prhs(1)) .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:revord:NonString',
     +                           'Input must be a string.')
C     The input must be a row vector.
      elseif (mxGetM(prhs(1)) .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:revord:NonRowVector',
     +                           'Input must be a row vector.')
      endif
    
C     Get the length of the input string.
      strlen = mxGetM(prhs(1))*mxGetN(prhs(1))
      if (strlen .gt. maxbuf) then
         call mexErrMsgIdAndTxt ('MATLAB:revord:maxbuf',
     +                           'Max string length 100.')
      endif

C     Get the string contents (dereference the input pointer).
      status = mxGetString(prhs(1), input_buf, maxbuf)

C     Check if mxGetString is successful.
      if (status .ne. 0) then 
         call mexErrMsgIdAndTxt ('MATLAB:revord:readError',
     +                           'Error reading string.')
      endif

C     Initialize outbuf_buf to blanks. This is necessary on some compilers.
      output_buf = ' '

C     Call the computational subroutine.
      call revord(input_buf, strlen, output_buf)
    
C     Set output_buf to MATLAB mexFunction output.
      plhs(1) = mxCreateString(output_buf)

      return
      end
