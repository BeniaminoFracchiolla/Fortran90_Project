      program nozzle
c
c  This program solves the 1D Euler equations for a gamma-law gas
c  global dt evaluation per time-accurate solution and wall condition at left and right boundaries
c
      implicit double precision (a-h,o-z)
      parameter (maxmx = 1000)
      parameter (meqn = 3)
      parameter (mwaves = 3)
      parameter (mbc = 2)
c
      dimension q(1-mbc:maxmx+mbc, meqn)
      dimension qold(1-mbc:maxmx+mbc, meqn)
      dimension qtmp(1-mbc:maxmx+mbc, meqn)
      dimension rhs(1-mbc:maxmx+mbc, meqn)
      dimension dt(1-mbc:maxmx+mbc)
      dimension dt1(1-mbc:maxmx+mbc)
      dimension dt2(1-mbc:maxmx+mbc)
      dimension qpr(1-mbc:maxmx+mbc, meqn)
      dimension qm(1-mbc:maxmx+mbc, meqn),qp(1-mbc:maxmx+mbc, meqn)
      dimension fm(1-mbc:maxmx+mbc, meqn),fp(1-mbc:maxmx+mbc, meqn)
      dimension fluxt(1-mbc:maxmx+mbc, meqn)
      dimension dw(1-mbc:maxmx+mbc, meqn,2)
      dimension x(1-mbc:maxmx+mbc)
      dimension area(1-mbc:maxmx+mbc)
      dimension alph(4)
      common /param/  gamma,gamma1
c
c
c     # Runge-Kutta parameters
      irkmax = 2
      alph(1) = 0.42
      alph(2) = 1.  
c     irkmax = 1
c     alph(1) = 1.
c
      write(6,*) 'input gamma'
      read(5,*) gamma
      gamma1 = gamma - 1.d0
      write(6,*) 'input dx'
      read(5,*) dx
      write(6,*) 'input cfl, itmax, itstmp'
      read(5,*) cfl,itmax,itstmp
      write(6,*) 'input tmax'
      read(5,*) tmax
      write(6,*) 'MUSCL coefficient'
      read(5,*) rkap
      write(6,*) 'input epsilon'
      read(5,*) epsilon
c
c
c     # set grid and initial conditions
c
c     # domain is  0 <= x <= 1:
      ax = 0.d0
!     ax = -1.d0
!     bx = 10.d0
      bx = 1.d0
      mx= (bx-ax + 1d-12)/dx 
c
c     # check for errors in data:
c
      if (mx .gt. maxmx) then
         write(6,*) "attento al dimensionamento"
         call exit
      endif
c
      do 10 i=1-mbc,mx+mbc
         x(i) = ax + (i-0.5d0) * dx
         write(20,1020) x(i)
 1020    format(e16.6)
   10    continue
c
c define the area function
      do i=1-mbc,mx+mbc
         area(i) = 1.0d0
      end do

c
c     # initial condition
      call ic(maxmx,meqn,mbc,mx,x,dx,q)
c
      call out1eu(maxmx,meqn,mbc,mx,x,q,area)

      t = 0.d0
c 
c     # main loop
      do 100 n=1,itmax

         dtgl = 1.d10
c
c        # extend data from grid to bordering boundary cells:
         call bcw(maxmx,meqn,mbc,mx,q)

c
c        # store the solution
         do m=1,meqn
            do i=1-mbc,mx+mbc 
              qold(i,m) = q(i,m)
            end do
         end do
c
c        # Runge-Kutta
         do irk=1,irkmax    !RK
c
c        # compute primitive variables (rho, u, p)
         do i=1-mbc,mx+mbc
           qpr(i,1) = q(i,1)
           qpr(i,2) = q(i,2)/q(i,1)
           qpr(i,3) = gamma1*(q(i,3)-0.5*q(i,2)*q(i,2)/q(i,1))
           if(qpr(i,3) .le. 0.d0) then
             write(6,*) "Attento: pressione negativa in i = ",i,n
             write(6,*) "pressione = ", qpr(i,3)
             call exit
           end if
           if(irk .eq. 1) then
             u = qpr(i,2)
             a = sqrt(gamma*qpr(i,3)/qpr(i,1))
             dt1(i) = cfl*dx/dabs(u+a)
             dt2(i) = cfl*dx/dabs(u-a)
             dtgl = min(dt1(i),dt2(i),dtgl)
           end if
         end do
c
c        # compute left (m)  and right (p) state
c        1st order scheme
         do m=1,meqn
           do i=1,mx+1
             qp(i,m) = qpr(i,m)
             qm(i,m) = qpr(i-1,m)
           end do
         end do

c        higher-order scheme
         do m=1,meqn
           do i=0,mx+1
             dw(i,m,1) = qpr(i,m) - qpr(i-1,m)
             dw(i,m,2) = qpr(i+1,m) - qpr(i,m)
           end do
         end do
         rkp=1.+rkap
         rkm=1.-rkap
         do m=1,meqn
           do i=1,mx+1
c limiter
            erre = dw(i,m,2)/(dw(i,m,1)+1.e-12)
c           phi1=dmax1(0.d0, dmin1(1.d0,erre))        !minmod
            phi1=(erre*erre+erre)/(1.d0+erre*erre)    !van Albada

            erre = dw(i,m,1)/(dw(i,m,2)+1.e-12)
c           phi2=dmax1(0.d0, dmin1(1.d0,erre))        !minmod
            phi2=(erre*erre+erre)/(1.d0+erre*erre)    !van Albada

            erre = dw(i,m,1)/(dw(i-1,m,1)+1.e-12)
c           phi3=dmax1(0.d0, dmin1(1.d0,erre))        !minmod
            phi3=(erre*erre+erre)/(1.d0+erre*erre)    !van Albada

            erre = dw(i-1,m,1)/(dw(i,m,1)+1.e-12)
c           phi4=dmax1(0.d0, dmin1(1.d0,erre))        !minmod
            phi4=(erre*erre+erre)/(1.d0+erre*erre)    !van Albada
c
!           qp(i,m)=qpr(i,m)-0.25*(rkp*phi1*dw(i,m,1)
!    .                            +rkm*phi2*dw(i,m,2))
!           qm(i,m)=qpr(i-1,m)+0.25*(rkm*phi3*dw(i-1,m,1)
!    .                              +rkp*phi4*dw(i-1,m,2))
           end do
         end do
c
c       #compute flux at interface
         do i=1,mx+1
c          compute right state using p
           rhor=qp(i,1)
           if(rhor.lt.0) then
             write(6,*) "Negative density right in i = ",i,rhor
             call exit
           end if
           ur=qp(i,2)
           ppr=qp(i,3)
           q2r=ur*ur
           aar = sqrt(gamma*ppr/rhor)
           e0r=ppr/(gamma1*rhor)+.5*q2r
           h0r=e0r+ppr/rhor
           fp(i,1)=rhor*ur
           fp(i,2)=(rhor*ur*ur+ppr)
           fp(i,3)=(gamma/gamma1*ppr+0.5*rhor*q2r)*ur
c          compute left state using m
           rhol=qm(i,1)
           if(rhol.lt.0) then
             write(6,*) "Negative density left in i = ",i
             call exit
           end if
           ul=qm(i,2)
           ppl=qm(i,3)
           q2l=ul*ul
           aal = sqrt(gamma*ppl/rhol)
           e0l=ppl/(gamma1*rhol)+.5*q2l
           h0l=e0l+ppl/rhol
           fm(i,1)=rhol*ul
           fm(i,2)=(rhol*ul*ul+ppl)
           fm(i,3)=(gamma/gamma1*ppl+0.5*rhol*q2l)*ul
c          calculate delta vector
           drho=qp(i,1)-qm(i,1)
           dp=ppr-ppl
           du=ur-ul
c          calculate Roe averaged values
           rhols=sqrt(rhol)
           rhors=sqrt(rhor)
           rho=rhols*rhors
           denom=rhols+rhors
           u=(rhols*ul+rhors*ur)/denom
           h0=(rhols*h0l+rhors*h0r)/denom
           q2=u*u
           as=gamma1*(h0-q2/2.)
           if (as .le. 0.d0) then
             write(6,*) "Negative square speed of sound in i = ", i
             call exit
           end if
           a=sqrt(as)
c          calculate delta flux terms
           rm1=abs(u)
           rm2=abs(u+a)
           rm3=abs(u-a)
c  entropy fix
           epp1 = dmax1(1.d-20, (u-ul), (ur-u))
           epp2 = dmax1(1.d-20, ((u+a)-(ul+aal)), ((ur+aar)-(u+a)))
           epp3 = dmax1(1.d-20, ((u-a)-(ul-aal)), ((ur-aar)-(u-a)))
           if ( rm1 .lt. epp1 ) then
             rm1 = 0.5 * (rm1*rm1/epp1 + epp1)
           end if
           if ( rm2 .lt. epp2 ) then
             rm2 = 0.5 * (rm2*rm2/epp2 + epp2)
           end if
           if ( rm3 .lt. epp3 ) then
             rm3 = 0.5 * (rm3*rm3/epp3 + epp3)
           end if
c
           alp1=(drho-dp/as)*rm1
           alp2=(dp/2./as+rho*du/2./a)*rm2
           alp3=(dp/2./as-rho*du/2./a)*rm3
           df1=alp1+alp2+alp3
           df2=(u*alp1)+(u+a)*alp2+(u-a)*alp3
           df3=(.5*q2*alp1)+(h0+a*u)*alp2+(h0-a*u)*alp3
c     calculate total flux at interface
           fluxt(i,1)=.5*(fp(i,1)+fm(i,1)-df1)
           fluxt(i,2)=.5*(fp(i,2)+fm(i,2)-df2)
           fluxt(i,3)=.5*(fp(i,3)+fm(i,3)-df3)
         end do
c
c       #compute residuals and update the solution
         do m=1,meqn
           do i=1,mx
             rhs(i,m) = (fluxt(i+1,m)-fluxt(i,m))/dx
             qtmp(i,m) = qold(i,m) - alph(irk)*
c    .                dt(i) * (fluxt(i+1,m)-fluxt(i,m)) / dx
     .                dtgl * (fluxt(i+1,m)-fluxt(i,m)) / dx
           end do
         end do        
c     source terms for quasi-1d flow
         do i = 1, mx
           rho = q(i,1)
           rhou = q(i,2)
           u = rhou/rho
           e = q(i,3)
           h = gamma*e - gamma1*0.5d0*rhou*u
           so1 = - rhou * (area(i+1)-area(i-1))*0.5d0/dx/area(i)
           so2 = - rhou*u * (area(i+1)-area(i-1))*0.5d0/dx/area(i)
           so3 = - h*u * (area(i+1)-area(i-1))*0.5d0/dx/area(i)
           rhs(i,1) = rhs(i,1) - so1
           rhs(i,2) = rhs(i,2) - so2
           rhs(i,3) = rhs(i,3) - so3
c          q(i,1) = qtmp(i,1) + alph(irk)* dt(i) * so1
c          q(i,2) = qtmp(i,2) + alph(irk)* dt(i) * so2
c          q(i,3) = qtmp(i,3) + alph(irk)* dt(i) * so3
           q(i,1) = qtmp(i,1) + alph(irk)* dtgl * so1
           q(i,2) = qtmp(i,2) + alph(irk)* dtgl * so2
           q(i,3) = qtmp(i,3) + alph(irk)* dtgl * so3
         end do

       end do     !RK
c
c       # compute residual norms
         rhs1_l1=0.
         rhs2_l1=0.
         rhs3_l1=0.
         do i=1,mx
           rhs1_l1=rhs1_l1+dabs(rhs(i,1))
           rhs2_l1=rhs2_l1+dabs(rhs(i,2))
           rhs3_l1=rhs3_l1+dabs(rhs(i,3))
         end do
         rhs1_l1=rhs1_l1/mx
         rhs2_l1=rhs2_l1/mx
         rhs3_l1=rhs3_l1/mx
c    
         rhs_inf=-1.
         do i=1,mx
           do m=1,meqn
             rhs_inf=max(rhs_inf,dabs(rhs(i,m)))
             i_max = i
           end do
         end do

         t = t + dtgl

         if(n/itstmp*itstmp.eq.n) then
            write(6,*) n,t,rhs1_l1,rhs2_l1,rhs3_l1,i_max,rhs_inf
            call out1eu(maxmx,meqn,mbc,mx,x,q,area)
         end if

         if (rhs_inf .lt. epsilon) go to 101
         if(t.ge.tmax) go to 101
c
  100    continue
c
  101    continue

         write(6,*) n,t,rhs1_l1,rhs2_l1,rhs3_l1,i_max,rhs_inf
         call out1eu(maxmx,meqn,mbc,mx,x,q,area)
c
      stop 
      end
c
c =========================================================
       subroutine ic(maxmx,meqn,mbc,mx,x,dx,q)
c =========================================================
c
c     # Set initial conditions for q.
c
c     # quasi 1d nozzle
c
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, meqn)
      dimension x(1-mbc:maxmx+mbc)
      common /param/  gamma,gamma1
c
c
       write(6,*) 'input left state:p,rho,u'
       read(5,*) pL,rhoL,uL
       write(6,*) "valori L",pL,rhoL,uL
       write(6,*) 'input right state:p,rho,u '
       read(5,*) pR,rhoR,uR
       write(6,*) "Valori R",pR,rhoR,uR

c      write(6,*) pL,rhoL,uL
c      write(6,*) pR,rhoR,uR

c  
       do 150 i=1-mbc,mx+mbc
!        if(x(i).le.5.0d0) then
!         if(x(i).le.0.8d0) then
         if(x(i).le.0.3d0) then
           q(i,1) = rhoL
           q(i,2) = rhoL*uL
           q(i,3) = pL/gamma1 + 0.5d0*rhoL*uL**2
         else
           q(i,1) = rhoR
           q(i,2) = rhoR*uR
           q(i,3) = pR/gamma1 + 0.5d0*rhoR*uR**2
         end if
  150  continue
c
      do 20 i=1,mx
          do m=1,meqn
c            # exponents with more than 2 digits cause problems reading
c            # into matlab... reset tiny values to zero:
             if (dabs(q(i,m)) .lt. 1d-99) q(i,m) = 0.d0
             enddo
          rho = q(i,1)
          pr = gamma1*(q(i,3) - 0.5d0*q(i,2)**2 / rho)
          u = q(i,2)/rho
          c = dsqrt(gamma*pr/rho)
          ama = u/c
c ORIG    write(10,1002) rho,u,pr
          write(10,*) x(i),rho,u,pr
 1002     format(3e14.6)
 20       continue

      return
      end
c
c =========================================================
      subroutine out1eu(maxmx,meqn,mbc,mx,x,q,area)
c =========================================================
c
c     # Output the results in primitive variables
c     # (density, velocity, pressure) for the Euler equations.
c
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, meqn)
      dimension x(1-mbc:maxmx+mbc)
      dimension area(1-mbc:maxmx+mbc)
      character*11 fname
      common /param/  gamma,gamma1
c
c
         fname = "outfile"
c        open(unit=100,file=fname,status="unknown",
c    .       form="formatted")
         open(10,file=fname,status="unknown",
     .       form="formatted")
         fname = "out"
         open(11,file=fname,status="unknown",
     .       form="formatted")


 1001 format(e16.8, 2i5, /)
      write(10,*) 'VARIABLES = "x", "A", "rho", "u", "p", "M"'
      write(10,*) 'ZONE F=POINT, I=',mx
      do 20 i=1,mx
          do m=1,meqn
c            # exponents with more than 2 digits cause problems reading
c            # into matlab... reset tiny values to zero:
             if (dabs(q(i,m)) .lt. 1d-99) q(i,m) = 0.d0
             enddo
          rho = q(i,1)
          pr = gamma1*(q(i,3) - 0.5d0*q(i,2)**2 / rho)
          u = q(i,2)/rho
          c = dsqrt(gamma*pr/rho)
          ama = u/c
c ORIG    write(10,1002) rho,u,pr
          write(10,*) x(i),area(i),rho,u,pr,ama
          write(11,*) x(i),rho,u,pr
 1002     format(3e14.6)
 20       continue
c
      close(10)
      close(11)
      return
      end
c
c =========================================================
      subroutine bcw(maxmx,meqn,mbc,mx,q)
c =========================================================
c
c     # Extend the data from the computational region 
c     #      i = 1, 2, ..., mx2
c     # to the virtual cells outside the region, with
c     #      i = 1-ibc  and   i = mx+ibc   for ibc=1,...,mbc
c
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, meqn)
      common /param/  gamma,gamma1
c
c     # left boundary wall
c     # right boundary wall
c     ----------------------------------------------
c

c Left Boundary  
      rhol = q(1,1)
      ul = q(1,2)/q(1,1)
      pl = (q(1,3)-0.5d0*rhol*ul**2)*gamma1
      al = dsqrt(gamma*pl/rhol)
      riel = ul - 2.d0*al/gamma1
      aml = ul/al
      pt = pl * (1.d0 + gamma1*0.5d0 * aml**2)**(gamma/gamma1)
      rhot = rhol * (1.d0 + gamma1*0.5d0 * aml**2)**(1.d0/gamma1)

      u = 0.
      a = - riel * gamma1 * 0.5
      p = pt 
      rho = rhot 

C closed wall condition at left boundary/extrapolation condition      
c      do 10 ibc=1,mbc
c            q(1-ibc,1) = rho
c            q(1-ibc,2) = rho*u
c            q(1-ibc,3) = p/gamma1 + 0.5d0*rho*u**2

C supersonic inlet condition at left boundary
      do 10 ibc=1,mbc
          q(1-ibc,1)=q(1,1)
          q(1-ibc,2)=q(1,2)
          q(1-ibc,3)=q(1,3)
   10       continue

c Right Boundary  

      rhor = q(mx,1)
      ur = q(mx,2)/q(mx,1)
      pr = (q(mx,3)-0.5d0*rhor*ur**2)*gamma1
      ar = dsqrt(gamma*pr/rhor)
      rier = ur + 2.d0*ar/gamma1
      amr = ur/ar
      pt = pr * (1.d0 + gamma1*0.5d0 * amr**2)**(gamma/gamma1)
      rhot = rhor * (1.d0 + gamma1*0.5d0 * amr**2)**(1.d0/gamma1)

      u = 0.
      a = rier * gamma1 * 0.5
      p = pt 
      rho = rhot 

C closed wall condition at right boundary
c      do 20 ibc=1,mbc
c            q(mx+ibc,1) = rho
c            q(mx+ibc,2) = rho*u
c            q(mx+ibc,3) = p/gamma1 + 0.5d0*rho*u**2

C supersonic outlet condition at right boundary/extrapolation condition
       do 20 ibc=1,mbc
          q(mx+ibc,1)=q(mx,1)
          q(mx+ibc,2)=q(mx,2)
          q(mx+ibc,3)=q(mx,3)
   20       continue

C subsonic outlet at right boundary

c     rhor = q(mx,1)
c     ur = q(mx,2)/q(mx,1)
c     pr = (q(mx,3)-0.5d0*rhor*ur**2)*gamma1
c     ar = dsqrt(gamma*pr/rhor)
c     rier = ur + 2.d0*ar/gamma1
c     porgr = pr/rhor**gamma

c     pL = 1.d5
c     pR = pL/10.d0
c     pexit = pR
c     p = pexit
c     rho = (p/porgr)**(1.d0/gamma)
c     a = dsqrt(gamma*p/rho)
c     u = rier - 2.d0*a/gamma1
c

c     do 20 ibc=1,mbc
c           q(mx+ibc,1) = rho
c           q(mx+ibc,2) = rho*u
c           q(mx+ibc,3) = p/gamma1 + 0.5d0*rho*u**2
c  20       continue

      return
      end
c
