CCC routine fills 2 GB of RAM, may be virtual
CCC 2147483647= max length of 32-bit integer
CCC 1073741824= 1 GB of RAM...
CCC ignore-- doesnt seem to matter if real/integer*8 or not
CCC complie using gfortran -fno-range-check -m64 -o memflush memflush.f
CCC 34359738367= max length of 64-bit integer
      real*16 x(1073741824)
      integer*16 i
      do i=1,1073741824
         x(i)=i*0.928372451
      end do
      end
