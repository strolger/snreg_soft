C     Routines to read and write various format disk files
C     

C Read a Vista format floating point file
      subroutine rvfile(nhead,header,nx,ny,data,fname)
      real data(1)
      character*5760 header
      character*(*) fname
      
      open(unit=19,file=fname,status='old',form='unformatted')
      read(19) header
      nhead = 72
      i = index(header,'NAXIS1')
      read(header(i+10:),'(i20)') nx
      i = index(header,'NAXIS2')
      read(header(i+10:),'(i20)') ny
      dx = 0
      dy = 0
      i = index(header,'CRVAL1')
      if(i.gt.0) read(header(i+10:),'(g20.0)') dx
      i = index(header,'CRVAL2')
      if(i.gt.0) read(header(i+10:),'(g20.0)') dy
      read(19) (data(i),i=1,nx*ny)
      close(19)
      
      return
      end
      
C Write a Vista format floating point file
      subroutine wvfile(nhead,header,nx,ny,data,fname)
      real data(1)
      character*(*) fname
      character*80 header(nhead)
      character*5760 vistaheader
      
      open(unit=19,file=fname,status='unknown',form='unformatted')
      nhline = 36*((nhead+35)/36)
      do 21 i = 1,5760
         j = (i-1)/80 + 1
         k = mod(i-1,80) + 1
         if(j.le.nhead) then
            vistaheader(i:i) = header(j)(k:k)
         else
            vistaheader(i:i) = ' '
         end if
 21   continue
      write(19) vistaheader
      write(19) (data(i),i=1,nx*ny)
      close(19)
      return
      end
      
C Program to read IRAF "floating" (i.e. scaled 32 bit integer) FITS files
      subroutine rirafile(nhead,header,nx,ny,data,fname)
      real data(1)
      character*(*) header(nhead), fname
      equivalence (r4, i4)

      nhead = 500
      call rffint(nhead,header,nx,ny,data,fname)
      bscale = 1
      bzero = 0
      do 5 i = 1,nhead
         if(header(i)(:6) .eq. 'BSCALE') read(header(i),1000) bscale
         if(header(i)(:5) .eq. 'BZERO') read(header(i),1000) bzero
 1000    format(10x,g20.0)
 5    continue

C      write(6,*) nx, ny, bscale, bzero
      do 10 i = 1,nx*ny
         r4 = data(i)
         data(i) = bscale*i4+bzero
 10   continue

      return
      end

C Read a bitmap file into an I*2 array
      subroutine rbitfile(nhead,header,nx,ny,data,fname)
      integer*2 data(1)
      character*(*) fname
      character*(*) header
      
      call rfbitfile(nhead,header,nx,ny,data,fname)
      call b1toi2(data,nx*ny,data)

      return
      end

C Write an I*2 array into a bitmap file
      subroutine wbitfile(nhead,header,nx,ny,data,fname)
      integer*2 data(1)
      character*(*) fname
      character*(*) header
      
      call i2tob1(data,nx*ny,data)
      call wfbitfile(nhead,header,nx,ny,data,fname)

      return
      end
      
C Convert an I2 -> I4 array
      subroutine i2toi4(i2,n,i4)
      integer*2 i2(n)
      integer*4 i4(n)
      do 20 j = n,1,-1
         i4(j) = i2(j)
 20   continue
      return
      end
      
C Convert an I4 -> I2 array
      subroutine i4toi2(i4,n,i2)
      integer*2 i2(n)
      integer*4 i4(n)
      do 20 j = 1,n
         i2(j) = i4(j)
 20   continue
      return
      end
      
C Convert an I4 -> R4 array
      subroutine i4tor4(i4,n,r4)
      integer*4 i4(n)
      real*4 r4(n)
      do 20 j = 1,n
         r4(j) = i4(j)
 20   continue
      return
      end
      
C Convert an R4 -> I4 array
      subroutine r4toi4(r4,n,i4)
      integer*4 i4(n)
      real*4 r4(n)
      do 20 j = 1,n
         i4(j) = nint(r4(j))
 20   continue
      return
      end

C Convert a B1 -> I2 array
      subroutine b1toi2(b1,n,i2)
      integer*2 b1(1), i2(n), k2
      integer k4, bits(16)
      equivalence (k2,k4)
      data bits /1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,
     $     16384,32768/
      
      k = (n+15) / 16
      k2 = b1(k)
      do 20 j = n,1,-1
         i = mod(j-1,16)
         if(and(k4,bits(i+1)) .eq. 0) then
            i2(j) = 0
         else
            i2(j) = 1
         end if
         if(i.eq.0) then
            k = k - 1
            k2 = b1(k)
         end if
 20   continue
      return
      end

C Convert an I2 -> B1 array
      subroutine i2tob1(i2,n,b1)
      integer*2 b1(1), i2(n), k2
      integer k4, bits(16)
      equivalence (k2,k4)
      data bits /1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,
     $     16384,32768/
      
      k = 1
      k4 = 0
      do 20 j = 1,n
         i = mod(j-1,16)
         if(i2(j).ne.0) k4 = k4 + bits(i+1)
         if(i.eq.15) then
            b1(k) = k2
            k = k + 1
            k4 = 0
         end if
 20   continue
      if(i.ne.15) b1(k) = k2
      return
      end

