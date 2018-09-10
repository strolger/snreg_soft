* Routines fit a coordinate transform x'(x,y) and y'(x,y)
*
* John Tonry - 18 June 1999
*
      subroutine cfit(nfit,n,x1,y1,x2,y2,xm,xp,ym,yp,xfit,yfit,
     $     xpar,ypar,xrms,yrms)
      real x1(n), y1(n), x2(n), y2(n), xfit(n), yfit(n)
      real xpar(nfit), ypar(nfit)
      parameter (maxpt=1000, maxfit=10)
      real s(maxfit*maxpt)
      
      do 10 i = 1,n
         x = (2*x1(i) - (xm+xp)) / (xp-xm)
         y = (2*y1(i) - (ym+yp)) / (yp-ym)
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

      subroutine chopres(n,x1,y1,x2,y2,xfit,yfit,sig,resmax)
      real x1(n), y1(n), x2(n), y2(n), xfit(n), yfit(n)
      worst = abs(x2(1)-xfit(1))
      iworst = 1
      do 10 i = 1,n
         if(abs(x2(i)-xfit(i)).gt.worst .or.
     $        abs(y2(i)-yfit(i)).gt.worst) then
            worst = amax1(abs(x2(i)-xfit(i)),abs(y2(i)-yfit(i)))
            iworst = i
         end if
 10   continue

      if(worst.gt.resmax) then
         do 11 k = iworst,n-1
            x1(k) = x1(k+1)
            y1(k) = y1(k+1)
            x2(k) = x2(k+1)
            y2(k) = y2(k+1)
            xfit(k) = xfit(k+1)
            yfit(k) = yfit(k+1)
 11      continue
      end if
      n = n - 1
      return
      end
