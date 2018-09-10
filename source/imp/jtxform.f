* Q+D to fit a coordinate transform x'(x,y) and y'(x,y)
* 
* Syntax: jtxform coord_file map_file [NFIT = 3 | 6 | 10]'
*
* coord_file has one header line, then
*               four columns with src_x  src_y  dest_x' dest_y'
*
* rev 1.2 (JT 12 Jul 2000): Added fitinfo and maxres option
* rev 1.1 (JT 22 Sep 1999): Added -halfpix option
*
* John Tonry - 18 June 1999
*
      parameter (maxpt=1000, maxfit=10)
      real x1(maxpt), y1(maxpt), x2(maxpt), y2(maxpt)
      real xfit(maxpt), yfit(maxpt)
      real xpar(maxfit), ypar(maxfit)
      character*80 line, mapfile
      character*10 terms(10)
      data terms /' +','*x +','*y +','*x*x +','*x*y +','*y*y +',
     $     '*x*x*x +','*x*x*y +','*x*y*y +','*y*y*y'/

      nfit = maxfit

      ihalfpix = 0
      ifitinfo = 0
      resmax = -1

      do 77 i = 1,maxfit
         xpar(i) = 0
         ypar(i) = 0
 77   continue

      if(iargc().lt.1) then
         write(6,*) 'Syntax: jtxform coord_file map_file'
         write(6,*) '             [-halfpix] [NFIT: 3|6|10]'
         write(6,*) '             [fitinfo] [maxres=x]'
         write(6,*)
         write(6,*) '        map_file is header then four columns '
         write(6,*) '        with src_x  src_y  dest_x dest_y'
         write(6,*)
         write(6,*) '        Coordinates are assumed to have 0.5,0.5'
         write(6,*) '        at the center of the lower left pixel'
         stop
      end if

      if(iargc().ge.2) then
         call getarg(2,mapfile)
      else
         mapfile = 'jtxform.map'
      end if

      do 7 i = 3,iargc()
         call getarg(i,line)
         if(line(:8).eq.'-halfpix') then
            ihalfpix = 1
         else if(line(:7).eq.'fitinfo') then
            ifitinfo = 1
         else if(line(:7).eq.'maxres=') then
            read(line(8:),*) resmax
         else
            read(line,*) nfit
         end if
 7    continue

      call getarg(1,line)
      open(1,file=line,status='old')

* Is the first line a header or data?
      n = 0
      do 1 i = 1,maxpt
         read(1,*,end=2,err=3) x1(n+1), y1(n+1), x2(n+1), y2(n+1)
         n = n + 1
         if(ihalfpix.eq.1) then
            x1(n) = x1(n) + 0.5
            y1(n) = y1(n) + 0.5
            x2(n) = x2(n) + 0.5
            y2(n) = y2(n) + 0.5
         end if
         if(n.eq.1) then
            xm = x1(n)
            xp = x1(n)
            ym = y1(n)
            yp = y1(n)
         else
            xm = amin1(xm,x1(n))
            xp = amax1(xp,x1(n))
            ym = amin1(ym,y1(n))
            yp = amax1(yp,y1(n))
         end if
         goto 1
 3       if(i.ne.1) then
            write(6,*) 'Error reading line ',i,' from ',line(:30)
            stop
         end if
 1    continue
 2    close(1)

      write(6,*) n, ' coordinates read'

* Chop out bad points by maximum residual
      if(resmax.gt.0) then
         ntry = min(5, n/5)
         do 101 i = 1,ntry
            call cfit(nfit,n,x1,y1,x2,y2,xm,xp,ym,yp,xfit,yfit,
     $           xpar,ypar,xrms,yrms)
            norig = n
            write(6,6101) n, xrms, yrms, resmax
 6101       format('N =',i3,'  Xrms: ',f7.2,'  Yrms: ',f7.2,
     $          '  chopping if >',f5.2 )

            call chopres(n,x1,y1,x2,y2,xfit,yfit,sig,resmax)
            if(n.eq.norig) goto 102
 101     continue
 102     continue
      end if

* Chop out bad points by rms
      do 10 ntry = 0,4
         call cfit(nfit,n,x1,y1,x2,y2,xm,xp,ym,yp,xfit,yfit,
     $        xpar,ypar,xrms,yrms)

         nbad = 0
         sig = 5 - 0.5*ntry
         call prune(n,x1,y1,x2,y2,xfit,yfit,nbad,sig,xrms,yrms)

         write(6,6100) n, xrms, yrms, sig, nbad
 6100    format('N =',i3,'  Xrms: ',f7.2,'  Yrms: ',f7.2,
     $        '  Sig =', f4.1,'  #bad =',i3)
 10   continue

* Get rough scale, rotation, and parity
      pi = 4*atan(1.0)

* This is the current fit of p1 as a function of p2
      xrng = (xp-xm) / 2
      yrng = (yp-ym) / 2
      det = xpar(2)/xrng*ypar(3)/yrng - xpar(3)/yrng*ypar(2)/xrng
      iparity = det / abs(det)
      scale = sqrt(abs(det))
      theta = 90/pi * (atan2(-xpar(3),ypar(3))+atan2(ypar(2),xpar(2)))
      if(theta.lt.0) theta = theta + 360

      write(6,*)
      write(6,6163) n, iparity, scale, theta, sqrt(xrms**2+yrms**2)
 6163 format('Astrometry fit: N =',i3,'  Parity =',i3,
     $     '   sc/th =',f7.4,f7.1,'  rms =',f8.3)

      write(6,*)
      write(6,*) 'C code for transformation: '

      write(6,6103) 'x', '(2*(i+0.5)-(', xp,xm,xp,xm
      write(6,6103) 'y', '(2*(j+0.5)-(', yp,ym,yp,ym
 6103 format(6x,a,' = ',a,f7.1,'+',f7.1,') / (',f7.1,'-',f7.1,');')
      if(nfit.eq.3) then
         write(6,6203) 'xp', (xpar(i),i=1,nfit)
         write(6,6203) 'yp', (ypar(i),i=1,nfit)
 6203    format(6x,a,' = ',f9.3,' +',f8.3,'*x +',f8.3,'*y;')
      else if(nfit.eq.6) then
         write(6,6206) 'xp', (xpar(i),i=1,nfit)
         write(6,6206) 'yp', (ypar(i),i=1,nfit)
 6206    format(6x,a,' = ',f9.3,' +',f8.3,'*x +',f8.3,'*y +',/,
     $        7x,f7.3,'*x*x +',f7.3,'*x*y +',f7.3,'*y*y;')
      else if(nfit.eq.10) then
         write(6,6210) 'xp', (xpar(i),i=1,nfit)
         write(6,6210) 'yp', (ypar(i),i=1,nfit)
 6210    format(6x,a,' = ',f9.3,' +',f8.3,'*x +',f8.3,'*y +',/,
     $        7x,f7.3,'*x*x +',f7.3,'*x*y +',f7.3,'*y*y +',/,
     $        7x,f7.3,'*x*x*x +',f7.3,'*x*x*y +',
     $        f7.3,'*x*y*y +',f7.3,'*y*y*y;')
      end if

      open(unit=2,file=mapfile,status='unknown')
      write(2,2001) xm, xp, ym, yp
      do 45 i = 1,maxfit
         write(2,2002) xpar(i), ypar(i), terms(i)
 45   continue
 2001 format(4f15.5,4x,'Range: x,y = (2*(i,j+0.5)-(p-m))/(p+m)')
 2002 format(2f15.5,4x,a)
      close(2)

      if(ifitinfo.eq.1) then
         l = index(mapfile, '.map')
         if(l.eq.0) l = index(mapfile, ' ')
         line = mapfile(:l-1) // '.mapinfo'
         open(unit=2,file=line,status='unknown')
         write(2,2163) nfit, n, iparity*scale, theta, 
     $        xrms, yrms, sqrt(xrms**2+yrms**2)
 2163    format('N =',i3,i4,' Map =',f7.4,f7.1,' rms =',3f8.3)

      end if

      open(unit=2,file='jtxform.out',status='unknown')
      write(2,2007)
 2007 format('     x1      y1       x2       y2      ',
     $     'xfit     yfit     x2-fit   y2-fit')
      do 50 i = 1,n
         write(2,2000) x1(i), y1(i), x2(i), y2(i), xfit(i), yfit(i),
     $        x2(i)-xfit(i), y2(i)-yfit(i)
 2000    format(8f9.3)
 50   continue

      stop
      end
