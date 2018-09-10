* Programs to read and write things to a header
      subroutine newfhead(n,header,npix,nline)
* Creates a new minimal FITS header
      character*80 header(6)
      header(1) = 'SIMPLE  =                    T'
      header(2) = 'BITPIX  =                   16'
      header(3) = 'NAXIS   =                    2'
      write(header(4),'(a,i20)') 'NAXIS1  = ', npix
      write(header(5),'(a,i20)') 'NAXIS2  = ', nline
      header(6) = 'END'
      n = 6
      return
      end

      subroutine afhead(n,header,line)
* Adds line to the FITS header just before END
      character*80 header(n+1), line
      n = n + 1
      header(n) = header(n-1)
      header(n-1) = line
      return
      end

      subroutine cfhead(n,header,line)
* Replaces the line in the FITS header with the same keyword as line
      character*80 header(n), line
      do 900 i = 1,n
          if(header(i)(1:8).eq.line(1:8)) then
              header(i) = line
              return
          end if
900   continue
      write(6,*) 'CFHEAD: keyword',line(1:8),'not found'
      return
      end

      subroutine rfitem(n,header,keyword,nbyte,nr,item)
* Reads the item with key keyword from the FITS header. nbyte is
* 1 for a character string, 2 for short integer, 4 for long integer,
* -4 for s.p. floating, and -8 for d.p. floating.
* Returns the number of bytes returned, nr
      parameter (ifitspar1=11)
      character*80 header(n)
      character*(*) keyword
      character*1 item(8), b(8)
      integer*2 i2
      integer*4 i4
      real r4
      real*8 r8
      equivalence (b,i2), (b,i4), (b,r4), (b,r8)
C      WRITE(6,*) 'RFITEM: KEYWORD ',KEYWORD,' SOUGHT ',LEN(KEYWORD)
* Locate the item
      m = 0
      do 900 i = 1,n
          if(header(i)(1:8).eq.keyword) then
              m = i
              goto 10
          end if
900   continue
      write(6,*) 'RFITEM: keyword ',keyword,' not found'
      nr = 0
      return
10    continue
* Locate the beginnings and ends of the item
      i1 = ifitspar1
* Is it a character string?
      if(header(m)(i1:i1).eq.''''.or.nbyte.eq.1) then
          i1 = i1 + 1
          i2 = i1
          do 901 i = i1+1,80
              if(header(m)(i:i).eq.'''') goto 11
              i2 = i
901       continue
11        continue
* Strip trailing blanks
          nr = lenc(header(m)(i1:i2))
          do 903 i = 0,nr-1
              item(i+1) = header(m)(i1+i:i1+i)
903       continue
          return
      end if
* Translate a number
      i1 = 0
      i2 = 0
      do 904 i = ifitspar1,80
          if(i1.eq.0.and.header(m)(i:i).ne.' ') i1 = i
          if(i1.ne.0 .and.
     1		(header(m)(i:i).eq.' '.or.header(m)(i:i).eq.'/')) then
              i2 = i - 1
              goto 905
          end if
904   continue
905   continue
      if(nbyte.eq.2) then
          read(header(m)(i1:i2),'(i70)') i2
      else if(nbyte.eq.4) then
          read(header(m)(i1:i2),'(i70)') i4
      else if(nbyte.eq.-4) then
          read(header(m)(i1:i2),'(g70.0)') r4
      else if(nbyte.eq.-8) then
          read(header(m)(i1:i2),'(g70.0)') r8
      end if
      nr = iabs(nbyte)
      do 906 i = 1,nr
          item(i) = b(i)
906   continue
      return
      end

      function lenc(s)
      character*(*) s
      lenc = len(s)
      do 900 i = len(s),1,-1
          if(s(i:i).ne.' ') then
              lenc = i
              return
          end if
900   continue
      lenc = 0
      return
      end

      subroutine rh2dfit(n,head,key,nx,ny,scale,coeff)
* Reads header quantities xxxx1, xxxx2 = nx, ny (xxxx = KEY)
* xxxxS10, S11, S20, S21 = scale x, scale y: (z = Sx0 * (x - Sx1))
* xxxxCmn, coefficients of polynomial
      character*80 head(n)
      character*(*) key
      real scale(4)
      real*8 coeff(1)
      character index*2, ky*8
      ky = key
      k = lenc(key)
      call rfitem(n,head,ky(:k)//'D1',4,nr,nx)
      if(nr.ne.4) write(6,*) 'RH2DFIT: NX not found'
      call rfitem(n,head,ky(:k)//'D2',4,nr,ny)
      if(nr.ne.4) write(6,*) 'RH2DFIT: NY not found'
      call rfitem(n,head,ky(:k)//'S10',-4,nr1,scale(1))
      call rfitem(n,head,ky(:k)//'S11',-4,nr2,scale(2))
      call rfitem(n,head,ky(:k)//'S20',-4,nr3,scale(3))
      call rfitem(n,head,ky(:k)//'S21',-4,nr4,scale(4))
      if(nr1+nr2+nr3+nr4.ne.16) write(6,*) 'RH2DFIT: SCALE not found'
      do 900 j = 0,ny-1
          do 901 i = 0,nx-1
              write(index,'(2i1)') i,j
              call rfitem(n,head,ky(:k)//'C'//index,
     1				-8,nr,coeff(ny*i+j+1))
              if(nr.ne.8) write(6,*) 'RH2DFIT: coefficient not found ',
     1                ky(:k)//'C'//index
901       continue
900   continue
      return
      end

      subroutine wh2dfit(n,head,key,nx,ny,scale,coeff)
* Writes header quantities xxxx1, xxxx2 = nx, ny (xxxx = KEY)
* xxxxS10, S11, S20, S21 = scale x, scale y: (z = Sx0 * (x - Sx1))
* xxxxCmn, coefficients of polynomial
      character*80 head(n), line
      character*(*) key
      real scale(4)
      real*8 coeff(1)
      character index*2, ky*8
      ky = key
      k = lenc(key)
      write(line,'(a8,''= '',i20)') ky(:k)//'D1    ', nx
      call afhead(n,head,line)
      write(line,'(a8,''= '',i20)') ky(:k)//'D2    ', ny
      call afhead(n,head,line)
      write(line,'(a8,''= '',1pe20.7)') ky(:k)//'S10   ', scale(1)
      call afhead(n,head,line)
      write(line,'(a8,''= '',1pe20.7)') ky(:k)//'S11   ', scale(2)
      call afhead(n,head,line)
      write(line,'(a8,''= '',1pe20.7)') ky(:k)//'S20   ', scale(3)
      call afhead(n,head,line)
      write(line,'(a8,''= '',1pe20.7)') ky(:k)//'S21   ', scale(4)
      call afhead(n,head,line)
      do 900 j = 0,ny-1
          do 901 i = 0,nx-1
              write(index,'(2i1)') i,j
              write(line,'(a8,''= '',1pe30.15)')
     1            ky(:k)//'C'//index//'   ', coeff(ny*i+j+1)
              call afhead(n,head,line)
901       continue
900   continue
      return
      end

      subroutine rhreal(n,head,key,var)
* Reads real header quantity var
      character*80 head(n)
      character*(*) key
      call rfitem(n,head,key,-4,nr,var)
      if(nr.ne.4) write(6,*) 'RHREAL: key not found', key
      return
      end

      subroutine whreal(n,head,key,var)
* Writes real header quantity var
      character*80 head(n), line
      character*(*) key
      character ky*8
      ky = key
      k = lenc(key)
      write(line,'(a8,''= '',1pe20.7)') ky(:k)//'      ', var
      call afhead(n,head,line)
      return
      end

      subroutine rhint(n,head,key,ival)
* Reads integer header quantity ival
      character*80 head(n)
      character*(*) key
      call rfitem(n,head,key,4,nr,ival)
      if(nr.ne.4) write(6,*) 'RHINT: key not found', key
      return
      end

      subroutine whint(n,head,key,ival)
* Writes integer header quantity ival
      character*80 head(n), line
      character*(*) key
      character ky*8
      ky = key
      k = lenc(key)
      write(line,'(a8,''= '',i20)') ky(:k)//'      ', ival
      call afhead(n,head,line)
      return
      end

      subroutine chint(n,head,key,int)
* Changes an integer header quantity
      character*80 head(n), line
      character*(*) key
      character ky*8
      ky = key
      k = lenc(key)
      write(line,'(a8,''= '',i20)') ky(:k)//'      ', int
      call cfhead(n,head,line)
      return
      end

      subroutine chreal(n,head,key,val)
* Changes an real header quantity
      character*80 head(n), line
      character*(*) key
      character ky*8
      ky = key
      k = lenc(key)
      write(line,'(a8,''= '',1pg20.4)') ky(:k)//'      ', val
      call cfhead(n,head,line)
      return
      end

      subroutine chchar(n,head,key,string)
* Changes an character header quantity
      character*80 head(n), line, temp
      character*(*) key, string
      character ky*8
      ky = key
      k = lenc(key)
      temp = string
      write(line,'(a8,''= '',a)') ky(:k)//'      ',
     1				''''//temp(:lenc(temp))//''''
      call cfhead(n,head,line)
      return
      end

      subroutine rhchar(n,head,key,nr,string)
      character*(*) key
      character*(*) string
      character*80 head(n)
      call rfitem(n,head,key,1,nr,string)
      if(nr.eq.0) write(6,*) 'RHCHAR: key not found', key
      return
      end

      subroutine hcvcorr(n,head,ra,dec,day,ut,hcv)
* Program to return the velocity of the sun
      parameter (pi=3.14159265, dr=pi/180.)
      integer mday(12)
      character*10 date
      character*80 head(n), string
      data mday /31,28,31,30,31,30,31,31,30,31,30,31/

* Assume first off coords in 'DD:MM:SS' format; 
* failing that look for floating numbers.
      call rfitem(n,head,'RA',1,nr,string)
      if(index(string(:nr),':').gt.0) then
         irh = nextract(string(:nr),i)
         irm = nextract(string(i+1:nr),j)
         irs = nextract(string(i+j+2:nr),k)
         irt = nextract(string(i+j+k+3:nr),k)
         ra = 15*(irh+(irm+(irs+irt/100.)/60.)/60.)
C         WRITE(6,*) IRH, IRM, IRS, IRT, RA
      else
         call rfitem(n,head,'RA',-4,nr,ra)
         if(nr.ne.4) write(6,*) 'HCVCORR: RA unknown format'
      end if

      call rfitem(n,head,'DEC',1,nr,string)
      if(index(string(:nr),':').gt.0) then
         idd = nextract(string(:nr),i)
         idm = nextract(string(i+1:nr),j)
         ids = nextract(string(i+j+2:nr),k)
         idt = nextract(string(i+j+k+3:nr),k)
         dec = idd+(idm+(ids+idt/100.)/60.)/60.
C         WRITE(6,*) IDD, IDM, IDS, IDT, DEC
         if(index(string(:nr),'-').gt.0) dec = -dec
      else
         call rfitem(n,head,'DEC',-4,nr,dec)
         if(nr.ne.4) write(6,*) 'HCVCORR: DEC unknown format'
      end if

      call rfitem(n,head,'UT',1,nr,string)
      if(index(string(:nr),':').gt.0) then
         iuh = nextract(string(:nr),i)
         ium = nextract(string(i+1:nr),j)
         ius = nextract(string(i+j+2:nr),k)
         iut = nextract(string(i+j+k+3:nr),k)
         ut = iuh+(ium+(ius+iut/100.)/60.)/60.
C         WRITE(6,*) IUH, IUM, IUS, IUT, UT
      else
         call rfitem(n,head,'UT',-4,nr,ut)
         if(nr.ne.4) write(6,*) 'HCVCORR: UT unknown format'
      end if

* Assume dd/mm/yy, and date is local time of start of night
      call rfitem(n,head,'DATE-OBS',1,nr,date)
      if(nr.gt.0) then
         i1 = nextract(date(1:nr),i)
         i2 = nextract(date(i+1:nr),j)
         i3 = nextract(date(i+j+2:nr),iend)
      else

C* Assume mm/dd/yy, and date is local time of start of night
* Assume mm/dd/yy, and date is local time of start of night
         call rfitem(n,head,'DATE    ',1,nr,date)
C         i2 = nextract(date(1:nr),iend)
C         i1 = nextract(date(iend+1:nr),inew)
         i1 = nextract(date(1:nr),iend)
         i2 = nextract(date(iend+1:nr),inew)
         i3 = nextract(date(iend+inew+2:nr),iend)
      end if

      if(i1.lt.1.or.i1.gt.mday(i2).or.i2.lt.1.or.i2.gt.12)
     1    write(6,*) 'Bad date, day/month/year =',i1,i2,i3,'    ',date
      day = i1
      do 900 i = 1,i2-1
          day = day + mday(i)
900   continue
      if(mod(i3,4).eq.0.and.mod(i3,100).ne.0.and.i2.gt.2) day = day + 1
      day = day + ut/24 - .5
      x = cos(dr*(190.0042 + .9863*day - ra)) * cos(dr*dec)
      y = -23.4542 * cos(dr*((day-80)*.9863))
      hcv = 29.79 * (x * cos(dr*y) + sin(dr*dec) * sin(dr*y))
      return
      end

      function nextract(s,i2)
      character*(*) s
      i1 = 0
      do 10 i = 1,len(s)
          if(i1.eq.0.and.s(i:i).ge.'0'.and.s(i:i).le.'9') i1 = i
          if(i1.ne.0.and.(s(i:i).lt.'0'.or.s(i:i).gt.'9')) then
              i2 = i - 1
              goto 11
          end if
10    continue
11    continue
      read(s(i1:i2),'(i20)',err=666) nextract
      return
 666  nextract = -1
      return
      end

