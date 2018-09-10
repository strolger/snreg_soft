      SUBROUTINE QSORT8(N,X)
* Sorting program that uses a quicksort algorithm
      parameter (maxstack=256)
      REAL*8 X(N)
      REAL*8 KEY, KL, KR, KM, TEMP
      INTEGER L, R, M, LSTACK(maxstack), RSTACK(maxstack), SP
      INTEGER NSTOP
      LOGICAL MGTL, LGTR, RGTM
      DATA NSTOP /15/

      IF(N.LE.NSTOP) GOTO 100
      SP = 0
      SP = SP + 1
      LSTACK(SP) = 1
      RSTACK(SP) = N

* Sort a subrecord off the stack
* Set KEY = median of X(L), X(M), X(R)
1     L = LSTACK(SP)
      R = RSTACK(SP)
      SP = SP - 1
      M = (L + R) / 2
      KL = X(L)
      KM = X(M)
      KR = X(R)
      MGTL = KM .GT. KL
      RGTM = KR .GT. KM
      LGTR = KL .GT. KR
      IF(MGTL .EQV. RGTM) THEN
          IF(MGTL .EQV. LGTR) THEN
              KEY = KR
          ELSE
              KEY = KL
          ENDIF
      ELSE
          KEY = KM
      ENDIF

      I = L
      J = R

* Find a big record on the left
10    IF(X(I).GE.KEY) GOTO 11
      I = I + 1
      GOTO 10
11    CONTINUE
* Find a small record on the right
20    IF(X(J).LE.KEY) GOTO 21
      J = J - 1
      GOTO 20
21    CONTINUE
      IF(I.GE.J) GOTO 2
* Exchange records
      TEMP = X(I)
      X(I) = X(J)
      X(J) = TEMP
      I = I + 1
      J = J - 1
      GOTO 10

* Subfile is partitioned into two halves, left .le. right
* Push the two halves on the stack
2     IF(J-L+1 .GT. NSTOP) THEN
          SP = SP + 1
          LSTACK(SP) = L
          RSTACK(SP) = J
      ENDIF
      IF(R-J .GT. NSTOP) THEN
          SP = SP + 1
          LSTACK(SP) = J+1
          RSTACK(SP) = R
      ENDIF
      IF(SP.GT.MAXSTACK) THEN
         WRITE(6,*) 'QSORT4: Fatal error from stack overflow'
         WRITE(6,*) 'Fall back on sort by insertion'
         GOTO 100
      END IF

* Anything left to process?
      IF(SP.GT.0) GOTO 1

* Sorting routine that sorts the N elements of double precision
* array X by straight insertion between previously sorted numbers
100   DO 110 J = N-1,1,-1
      K = J
      DO 120 I = J+1,N
      IF(X(J).LE.X(I)) GOTO 121
120   K = I
121   CONTINUE
      IF(K.EQ.J) GOTO 110
      TEMP = X(J)
      DO 130 I = J+1,K
130   X(I-1) = X(I)
      X(K) = TEMP
110   CONTINUE
      RETURN
      END

      SUBROUTINE QSORT4(N,X)
* Sorting program that uses a quicksort algorithm
      parameter (maxstack=256)
      REAL X(N)
      REAL KEY, KL, KR, KM, TEMP
      INTEGER L, R, M, LSTACK(maxstack), RSTACK(maxstack), SP
      INTEGER NSTOP
      LOGICAL MGTL, LGTR, RGTM
      DATA NSTOP /15/

      IF(N.LE.NSTOP) GOTO 100
      SP = 0
      SP = SP + 1
      LSTACK(SP) = 1
      RSTACK(SP) = N

* Sort a subrecord off the stack
* Set KEY = median of X(L), X(M), X(R)
1     L = LSTACK(SP)
      R = RSTACK(SP)
      SP = SP - 1
      M = (L + R) / 2
      KL = X(L)
      KM = X(M)
      KR = X(R)
      MGTL = KM .GT. KL
      RGTM = KR .GT. KM
      LGTR = KL .GT. KR
      IF(MGTL .EQV. RGTM) THEN
          IF(MGTL .EQV. LGTR) THEN
              KEY = KR
          ELSE
              KEY = KL
          ENDIF
      ELSE
          KEY = KM
      ENDIF

      I = L
      J = R

* Find a big record on the left
10    IF(X(I).GE.KEY) GOTO 11
      I = I + 1
      GOTO 10
11    CONTINUE
* Find a small record on the right
20    IF(X(J).LE.KEY) GOTO 21
      J = J - 1
      GOTO 20
21    CONTINUE
      IF(I.GE.J) GOTO 2
* Exchange records
      TEMP = X(I)
      X(I) = X(J)
      X(J) = TEMP
      I = I + 1
      J = J - 1
      GOTO 10

* Subfile is partitioned into two halves, left .le. right
* Push the two halves on the stack
2     IF(J-L+1 .GT. NSTOP) THEN
          SP = SP + 1
          LSTACK(SP) = L
          RSTACK(SP) = J
      ENDIF
      IF(R-J .GT. NSTOP) THEN
          SP = SP + 1
          LSTACK(SP) = J+1
          RSTACK(SP) = R
      ENDIF
      IF(SP.GT.MAXSTACK) THEN
         WRITE(6,*) 'QSORT4: Fatal error from stack overflow'
         WRITE(6,*) 'Fall back on sort by insertion'
         GOTO 100
      END IF

* Anything left to process?
      IF(SP.GT.0) GOTO 1

* Sorting routine that sorts the N elements of single precision
* array X by straight insertion between previously sorted numbers
100   DO 110 J = N-1,1,-1
      K = J
      DO 120 I = J+1,N
      IF(X(J).LE.X(I)) GOTO 121
120   K = I
121   CONTINUE
      IF(K.EQ.J) GOTO 110
      TEMP = X(J)
      DO 130 I = J+1,K
130   X(I-1) = X(I)
      X(K) = TEMP
110   CONTINUE
      RETURN
      END

      SUBROUTINE QSORT2(N,X,Y)
* Sorting program that uses a quicksort algorithm, sorts both x and y by x
      parameter (maxstack=256)
      REAL X(N), Y(N)
      REAL KEY, KL, KR, KM, TEMP
      INTEGER L, R, M, LSTACK(maxstack), RSTACK(maxstack), SP
      INTEGER NSTOP
      LOGICAL MGTL, LGTR, RGTM
      DATA NSTOP /15/

      IF(N.LE.NSTOP) GOTO 100
      SP = 0
      SP = SP + 1
      LSTACK(SP) = 1
      RSTACK(SP) = N

* Sort a subrecord off the stack
* Set KEY = median of X(L), X(M), X(R)
1     L = LSTACK(SP)
      R = RSTACK(SP)
      SP = SP - 1
      M = (L + R) / 2
      KL = X(L)
      KM = X(M)
      KR = X(R)
      MGTL = KM .GT. KL
      RGTM = KR .GT. KM
      LGTR = KL .GT. KR
      IF(MGTL .EQV. RGTM) THEN
          IF(MGTL .EQV. LGTR) THEN
              KEY = KR
          ELSE
              KEY = KL
          ENDIF
      ELSE
          KEY = KM
      ENDIF

      I = L
      J = R

* Find a big record on the left
10    IF(X(I).GE.KEY) GOTO 11
      I = I + 1
      GOTO 10
11    CONTINUE
* Find a small record on the right
20    IF(X(J).LE.KEY) GOTO 21
      J = J - 1
      GOTO 20
21    CONTINUE
      IF(I.GE.J) GOTO 2
* Exchange records
      TEMP = X(I)
      X(I) = X(J)
      X(J) = TEMP
      TEMP = Y(I)
      Y(I) = Y(J)
      Y(J) = TEMP
      I = I + 1
      J = J - 1
      GOTO 10

* Subfile is partitioned into two halves, left .le. right
* Push the two halves on the stack
2     IF(J-L+1 .GT. NSTOP) THEN
          SP = SP + 1
          LSTACK(SP) = L
          RSTACK(SP) = J
      ENDIF
      IF(R-J .GT. NSTOP) THEN
          SP = SP + 1
          LSTACK(SP) = J+1
          RSTACK(SP) = R
      ENDIF
      IF(SP.GT.MAXSTACK) THEN
         WRITE(6,*) 'QSORT4: Fatal error from stack overflow'
         WRITE(6,*) 'Fall back on sort by insertion'
         GOTO 100
      END IF

* Anything left to process?
      IF(SP.GT.0) GOTO 1

* Sorting routine that sorts the N elements of single precision
* array X by straight insertion between previously sorted numbers
100   DO 110 J = N-1,1,-1
      K = J
      DO 120 I = J+1,N
      IF(X(J).LE.X(I)) GOTO 121
120   K = I
121   CONTINUE
      IF(K.EQ.J) GOTO 110
      TEMPX = X(J)
      TEMPY = Y(J)
      DO 130 I = J+1,K
      X(I-1) = X(I)
130   Y(I-1) = Y(I)
      X(K) = TEMPX
      Y(K) = TEMPY
110   CONTINUE
      RETURN
      END
