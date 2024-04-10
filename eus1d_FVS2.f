      program eus1d_FVS2
c
c  This program solves the 1D Euler equations for a gamma-law gas
c  global dt evaluation per time-accurate solution and wall condition at left and right boundaries
c  This is done through Flux-Vector-Splitting method, in the specific
c  case with the Liou-Steffen scheme.
c
c     # Variables declaration
      implicit double precision (a-h,o-z)
      parameter (maxmx = 1000)
      parameter (meqn = 3)
      parameter (mbc = 2)
c
      dimension q(1-mbc:maxmx+mbc, meqn)          !conservative vectors
      dimension qold(1-mbc:maxmx+mbc, meqn)       !old conservative vec.
      dimension rho(1-mbc:maxmx+mbc)              !density vector
      dimension u(1-mbc:maxmx+mbc)                !velocity vector
      dimension pr(1-mbc:maxmx+mbc)               !pressure vector
      dimension a(1-mbc:maxmx+mbc)                !sound speed vector
      dimension dt(1-mbc:maxmx+mbc)
      dimension dt1(1-mbc:maxmx+mbc)
      dimension dt2(1-mbc:maxmx+mbc)
      dimension qpr(1-mbc:maxmx+mbc, meqn)        !primitive vectors
      dimension amai(1-mbc:maxmx+mbc)             !cell-interface MACH
      dimension pri(1-mbc:maxmx+mbc)              !cell-interface Pr.
      dimension fluxt(1-mbc:maxmx+mbc, meqn)      !cell-interface Flux
      dimension x(1-mbc:maxmx+mbc)                !cell-centered x-axis
      dimension area(1-mbc:maxmx+mbc)
      dimension alph(4)                           !RK coefficients vector
      common /param/  gamma,gamma1
c
c
c     # Runge-Kutta parameters
      irkmax = 2
      alph(1) = 0.42
      alph(2) = 1.  
c      irkmax = 1
c      alph(1) = 1.
c
c     # Input section
      write(6,*) 'input position of diaphram'
      read(5,*) diaph
      write(6,*) 'input gamma'
      read(5,*) gamma
      gamma1 = gamma - 1.d0
      write(6,*) 'input dx'
      read(5,*) dx
      write(6,*) 'input cfl, itmax, nstmp'       !itmax: max iteration n
      read(5,*) cfl,itmax,nstmp                  !nstmp: out stamp ratio
      write(6,*) 'input tmax'                    !tmax: max time allowed
      read(5,*) tout
c
c     # Set grid and initial conditions
c
c     # Domain is  0 <= x <= 1:
      ax = 0.d0
      bx = 1.d0
      mx = (bx-ax + 1d-12)/dx 
c
c     # Check for errors in data:
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
c     # Define the area function
      do i=1-mbc,mx+mbc
         area(i) = 1.0d0
      end do
c
c     # Initial condition
      call ic(maxmx,meqn,mbc,mx,x,dx,q,diaph)
c
      call out1eu(maxmx,meqn,mbc,mx,x,q,area)
c
      t = 0.d0
c 
c     # Main loop
      do 100 n=1,itmax
c
         dtgl = 1.d10
c
c        # Extend data from grid to wall boundary cells
         call bcw(maxmx,meqn,mbc,mx,q)
c
c        # Store the solutions in terms of conservative variables
         do m=1,meqn
            do i=1-mbc,mx+mbc 
              qold(i,m) = q(i,m)
            end do
         end do
c
c        # Runge-Kutta
         do irk=1,irkmax    !RK
c
c        # Compute primitive variables (rho, u, p) and CFL condition
           do i=1-mbc,mx+mbc
             qpr(i,1) = q(i,1)
             qpr(i,2) = q(i,2)/q(i,1)
             qpr(i,3) = gamma1*(q(i,3)-0.5d0*q(i,2)*q(i,2)/q(i,1))
c            # Rename of primitive variables in different vectors              
             rho(i)=qpr(i,1)
             u(i)=qpr(i,2)
             pr(i)=qpr(i,3)
             a(i)=dsqrt(gamma*qpr(i,3)/qpr(i,1))
c            # Check of non-physical solutions         
             if(qpr(i,3) .le. 0.d0) then
               write(6,*) "Attento: pressione negativa in i = ",i,n
               write(6,*) "pressione = ", qpr(i,3)
               call exit
             end if
c            # Evaluation of global dt
             if(irk.eq.1) then
               dt1(i) = cfl*dx/dabs(u(i)+a(i))
               dt2(i) = cfl*dx/dabs(u(i)-a(i))
               dtgl = min(dt1(i),dt2(i),dtgl)
             end if
           end do
c  
c        # Lius Steffen scheme
           do i=0,mx+1
c
             ama = u(i)/a(i)                     ! local Mach number
c
c        # RECOSTRUCTION OF LEFT CELL STATE AND RIGHT CELL STATE
c             
c            # Splitting of the interface Mach number and interface Pr.
             if(dabs(ama).le.1.0d0)then
c                   SUBSONIC CASE
c                   right/left Mach component
                    amar = +0.25d0*(ama+1.d0)**2
                    amal = -0.25d0*(ama-1.d0)**2
c                   right/left Pressure component 
                    prr = 0.5d0*pr(i)*(1.d0+ama)
                    prl = 0.5d0*pr(i)*(1.d0-ama)
             else
c                   SUPERSONIC CASE
c                   right/left Mach component
                    amar = 0.5d0*(ama+dabs(ama))
                    amal = 0.5d0*(ama-dabs(ama))
c                   right/left Pressure component
                    prr = 0.5d0*pr(i)*(ama+dabs(ama))/ama
                    prl = 0.5d0*pr(i)*(ama-dabs(ama))/ama
             end if
c             
c            # Compute interface Mach number and interface pressure
c             M(i+1/2)=M(i)+ + M(i+1)-
c             p(i+1/2)=p(i)+ + p(i+1)-
             if(i.ge.0.and.i.le.mx)then
                    amai(i) = amar
                    pri(i) = prr
             end if
c
             if(i.ge.1)then
                    amai(i-1) = amai(i-1) + amal
                    pri(i-1) = pri(i-1) + prl
             end if
           end do  
c
c        # COMPUTE INTERFACE NUMERICAL FLUXES
c
           do i=0,mx
             if(amai(i).ge.0.0)then
c                advection from left to right (evalueted in i)
                 ene = 0.5d0*u(i)*u(i)+a(i)*a(i)/(gamma1)
c             
                 fluxt(i,1) = amai(i)*rho(i)*a(i)
                 fluxt(i,2) = amai(i)*rho(i)*a(i)*u(i)+pri(i)
                 fluxt(i,3) = amai(i)*rho(i)*a(i)*ene
c            
             else
c                advection from right to left (evaluated in i+1)
                 ene = 0.5d0*u(i+1)*u(i+1)+a(i+1)*a(i+1)/(gamma1)
c
                 fluxt(i,1) = amai(i)*rho(i+1)*a(i+1)
                 fluxt(i,2) = amai(i)*rho(i+1)*a(i+1)*u(i+1)+pri(i)
                 fluxt(i,3) = amai(i)*rho(i+1)*a(i+1)*ene
c
             end if
           end do
c
c       # Update the solution
           do m=1,meqn
             do i=1,mx
               q(i,m) = qold(i,m) - alph(irk) *
     .                 dtgl * (fluxt(i,m)-fluxt(i-1,m)) / dx
             end do
           end do        
c
         end do     !RK
c
c       # Update the time
         t = t + dtgl
c        
         if(n/nstmp*nstmp.eq.n) then
            write(6,*) n,t,tout
         end if
c
         if(t.ge.tout) go to 101
c
  100    continue
c
  101    continue

         write(6,*) n,t,tout
         call out1eu(maxmx,meqn,mbc,mx,x,q,area)
c
      stop 
      end
c
c =========================================================
       subroutine ic(maxmx,meqn,mbc,mx,x,dx,q,diaph)
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

c      set the conservative value from wall boundary cells
c      rho, rho*u, p/rho+0.5*rho*u^2
       do 150 i=1-mbc,mx+mbc
         if(x(i).le.diaph) then
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
c     # left boundary open wall
c     # right boundary open wall
c     ----------------------------------------------
c
C supersonic inlet condition at left boundary
      do 10 ibc=1,mbc
          q(1-ibc,1)=q(1,1)
          q(1-ibc,2)=q(1,2)
          q(1-ibc,3)=q(1,3)
   10       continue
c
C supersonic outlet condition at right boundary/extrapolation condition
       do 20 ibc=1,mbc
          q(mx+ibc,1)=q(mx,1)
          q(mx+ibc,2)=q(mx,2)
          q(mx+ibc,3)=q(mx,3)
   20       continue
c
      return
      end
c
