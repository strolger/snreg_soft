CCC     It appears there is a 26-bit limit 
      real*8 x(671088)
      integer i
      do i=1,671088
         x(i)=i*0.928372451
         write (6,*) x(i)
      end do
      end
