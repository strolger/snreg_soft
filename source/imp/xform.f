* Q+D to fit a coordinate transform x'(x,y) and y'(x,y)
* 
* Syntax: xform coord_file map_file [NFIT = 3 | 6 | 10]'
*
* coord_file has one header line, then
*               four columns with src_x  src_y  dest_x' dest_y'
*
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

      do 77 i = 1,maxfit
         xpar(i) = 0
         ypar(i) = 0
 77   continue

      if(iargc().lt.1) then
         write(6,*) 'Syntax: xform coord_file map_file'
         write(6,*) '             [-halfpix] [NFIT: 3|6|10]'
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
         mapfile = 'xform.map'
      end if

      if(iargc().ge.3) then
         call getarg(3,line)
         if(line(:8).eq.'-halfpix') then
            ihalfpix = 1
            if(iargc().ge.4) then
               call getarg(3,line)
               read(line,*) nfit
            end if
         else
            read(line,*) nfit
         end if
      end if

      call getarg(1,line)
      open(1,file=line,status='old')

      n = 0
      read(1,*)
      xmid = 0
      ymid = 0
 1    read(1,*,end=2) x1(n+1), y1(n+1), x2(n+1), y2(n+1)
      if(ihalfpix.eq.1) then
         x1(n+1) = x1(n+1) + 0.5
         y1(n+1) = y1(n+1) + 0.5
         x2(n+1) = x2(n+1) + 0.5
         y2(n+1) = y2(n+1) + 0.5
      end if
      xmid = xmid + x1(n+1)
      ymid = ymid + y1(n+1)
      n = n + 1
      goto 1
 2    close(1)
      write(6,*) n, ' coordinates read'
      xmid = nint(xmid / n)
      ymid = nint(ymid / n)

      do 10 ntry = 0,4
         call cfit(nfit,n,x1,y1,x2,y2,xmid,ymid,xfit,yfit,
     $        xpar,ypar,xrms,yrms)

         nbad = 0
         sig = 5 - 0.5*ntry
         call prune(n,x1,y1,x2,y2,xfit,yfit,nbad,sig,xrms,yrms)

         write(6,6100) n, xmean, xrms, ymean, yrms, sig, nbad
 6100    format('N =',i3,'   X: ',2f7.2,'   Y: ',2f7.2,'  Sig =', f4.1,
     $        '  #bad =',i3)
 10   continue

* Get rough scale, rotation, and parity
      pi = 4*atan(1.0)

* This is the current fit of p1 as a function of p2
      det = xpar(2)/xmid*ypar(3)/ymid - xpar(3)/xmid*ypar(2)/ymid
      iparity = det / abs(det)
      scale = sqrt(abs(det))
      theta = atan2(-xpar(3)/xmid,ypar(3)/ymid) / pi * 180
      if(theta.lt.0) theta = theta + 360

      write(6,*)
      write(6,6163) n, iparity, scale, theta, sqrt(xrms**2+yrms**2)
 6163 format('Astrometry fit: N =',i3,'  Parity =',i3,
     $     '   sc/th =',f7.4,f7.1,'  rms =',f8.3)

      write(6,*)
      write(6,*) 'C code for transformation: '

      write(6,6103) 'x', '(i+0.5)', nint(xmid)
      write(6,6103) 'y', '(j+0.5)', nint(ymid)
 6103 format(6x,a,' = ',a,' /',i4,'.0 - 2;')
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
      write(2,2001) xmid, ymid, 'Mid: x,y = (i,j+0.5)/mid-2'
      do 45 i = 1,maxfit
         write(2,2001) xpar(i), ypar(i), terms(i)
 45   continue
 2001 format(2f15.5,4x,a)
      close(2)

      open(unit=2,file='xform.out',status='unknown')
      do 50 i = 1,n
         write(2,2000) x1(i), y1(i), x2(i), y2(i), xfit(i), yfit(i),
     $        x2(i)-xfit(i), y2(i)-yfit(i)
 2000    format(8f9.3)
 50   continue

      stop
      end

      subroutine prune(n,x1,y1,x2,y2,xfit,yfit,nbad,sig,xrms,yrms)
      real x1(n), y1(n), x2(n), y2(n), xfit(n), yfit(n)
      nbad = 0
      do 10 i = n,1,-1
         if(abs(x2(i)-xfit(i)).gt.sig*xrms .or.
     $        abs(y2(i)-yfit(i)).gt.sig*yrms) then
            nbad = nbad + 1
            do 11 k = i,n-1
               x1(k) = x1(k+1)
               y1(k) = y1(k+1)
               x2(k) = x2(k+1)
               y2(k) = y2(k+1)
               xfit(k) = xfit(k+1)
               yfit(k) = yfit(k+1)
 11         continue
         end if
 10   continue
      n = n - nbad
      return
      end


      subroutine cfit(nfit,n,x1,y1,x2,y2,xmid,ymid,xfit,yfit,
     $     xpar,ypar,xrms,yrms)
      real x1(n), y1(n), x2(n), y2(n), xfit(n), yfit(n)
      real xpar(nfit), ypar(nfit)
      parameter (maxpt=1000, maxfit=10)
      real s(maxfit*maxpt)
      
      do 10 i = 1,n
         x = x1(i) / xmid - 2
         y = y1(i) / ymid - 2
         s( 1 + nfit*(i-1)) = 1
         s( 2 + nfit*(i-1)) = x
         s( 3 + nfit*(i-1)) = y
         if(nfit.ge.6) then
            s( 4 + nfit*(i-1)) = x*x
            s( 5 + nfit*(i-1)) = x*y
            s( 6 + nfit*(i-1)) = y*y
         end if
         if(nfit.ge.10) then
            s( 7 + nfit*(i-1)) = x*x*x
            s( 8 + nfit*(i-1)) = x*x*y
            s( 9 + nfit*(i-1)) = x*y*y
            s(10 + nfit*(i-1)) = y*y*y
         end if
 10   continue

      call linearfit(n, x2, nfit, s, xpar)
      call linearfit(n, y2, nfit, s, ypar)

      do 20 i = 1,n
         xfit(i) = 0
         yfit(i) = 0
         do 21 k = 1,nfit
            xfit(i) = xfit(i) + s(k+nfit*(i-1)) * xpar(k)
            yfit(i) = yfit(i) + s(k+nfit*(i-1)) * ypar(k)
 21      continue
         xmean = xmean + x2(i) - xfit(i)
         xrms = xrms + (x2(i) - xfit(i))**2
         ymean = ymean + y2(i) - yfit(i)
         yrms = yrms + (y2(i) - yfit(i))**2
 20   continue

      xmean = xmean / n
      xrms = sqrt(xrms / n - xmean*xmean)
      ymean = ymean / n
      yrms = sqrt(yrms / n - ymean*ymean)

      return
      end

* Routine to fit a linear function
      subroutine linearfit(npt, y, nx, x, param)
      parameter (maxx=20)
      real x(nx,npt), y(npt), param(nx)
C Y = data being fitted as a linear function of X(I,*)
C PARAM = returned parameters of fit:
C	Y(*) = PARAM(1)*X(1,*)
C	     + PARAM(2)*X(2,*)
C	     + PARAM(3)*X(3,*)   etc.
C
      real*8 v(maxx), am(maxx*maxx), par8(maxx)
      if(nx.gt.maxx) then
         write(6,*) 'Too many parameters requested:', nx, maxx
         return
      end if

      if(npt.lt.nx) then
         write(6,*) 'Too few points provided:', npt, nx
         return
      end if

      do 10 i = 1,nx*nx
         if(i.le.nx) v(i) = 0
         am(i) = 0
 10   continue

C Accumulate sums
      do 20 i = 1,npt
         do 22 j = 1,nx
            v(j) = v(j) + dble(y(i))*dble(x(j,i))
            do 23 k = 1,j
               am(k+(j-1)*nx) = am(k+(j-1)*nx) + 
     $              dble(x(j,i))*dble(x(k,i))
 23         continue
 22      continue
 20   continue

      call solve(nx, v, par8, am)

      do 40 i = 1,nx
         param(i) = sngl(par8(i))
 40   continue

      return
      end

* Solve a system of linear equations
      subroutine solve(n, y, x, a)
      implicit real*8 (a-h,o-z)
      real*8 x(n), y(n), a(n,n)

* Fill in symmetrical matrix (assume bottom half full)
      do 30 l = 1,n
         do 31 k = l+1,n
            a(k,l) = a(l,k)
 31      continue
 30   continue

C Solve matrix by Gaussian elimination
      do 40 j = 1,n-1
         if(a(j,j).eq.0) then
            write(6,*) 'WHOA: singular matrix'
            return
         end if
         do 41 i = j+1,n
            r = a(i,j) / a(j,j)
            y(i) = y(i) - r*y(j)
            do 42 k = j+1,n
               a(i,k) = a(i,k) - r*a(j,k)
 42         continue
 41      continue
 40   continue

C Back substitute to solve for parameters
      x(n) = y(n) / a(n,n)
      do 50 i = 1,n-1
         j = n - i
         x(j) = y(j)
         do 51 k = j+1,n
            x(j) = x(j) - a(j,k)*x(k)
 51      continue
         x(j) = x(j) / a(j,j)
 50   continue

      return
      end
