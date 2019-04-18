module Precision
       implicit none
       integer, parameter :: dl = KIND(1.d0)
       integer, parameter :: sp = KIND(1.0)
end module Precision
module Pub
    use Precision
    implicit none
    real(dl), parameter :: tol=1.d-6
    real(dl), parameter :: PI=3.1415926d0
    real(dl) :: rs_sat,c_sat,rand2
    integer, parameter :: nhalo=176330227
    
    contains
       subroutine init_random_seed()

       INTEGER :: i, n, clock
       INTEGER, DIMENSION(:), ALLOCATABLE :: seed

       CALL RANDOM_SEED(size = n)
       ALLOCATE(seed(n))

       CALL SYSTEM_CLOCK(COUNT=clock)

       seed = clock + 37 * (/ (i - 1, i = 1, n) /)
       CALL RANDOM_SEED(PUT = seed)

       DEALLOCATE(seed)
       end


end module Pub
Subroutine bisection(f,x1,x2,eps,Root,flag)
!============================================================
! Solutions of equation f(x)=0 on [x1,x2] interval
! Method: Bisectional (closed domain) (a single root)
! Alex G. January 2010
!------------------------------------------------------------
! input ...
! f   - function - evaluates f(x) for any x in [x1,x2]
! x1  - left endpoint of initial interval
! x2  - right endpoint of initial interval
! eps - desired uncertainity of the root as |b-a|<eps
! output ...
! Root  - root of the equation f(x)=0
! flag  - indicator of success
!         >0 - a single root found, flag=number of iterations
!          0 - no solutions for the bisectional method
!         <0 - not a root but singularity, flag=number of iterations
!
! Comments: Function f(x) has to change sign between x1 and x2
!           Max number of iterations - 200 (accuracy (b-a)/2**200)
!====================================================================
use Precision
implicit none
real(dl):: f, x1, x2, eps, Root
real(dl):: a, b, c
integer i, flag
integer, parameter:: iter=200

!* check the bisection condition

if(f(x1)*f(x2)>0.0) then
  flag = 0
  return
end if

!* initialize calculations
a=x1
b=x2

!* Iterative refining the solution 
do i=1,iter
  c=(b+a)/2.0
  if(f(a)*f(c).le.0.0) then
      b = c
    else
      a=c
  end if
! condition(s) to stop iterations)
  if(abs(b-a)<= eps) exit  
end do
Root=(b+a)/2.0

!* check if it is a root or singularity
if (abs(f(Root)) < 1.0) then
  flag=i
  else
  flag = -i
end if
end subroutine bisection


function nfw(r)
use Precision
use Pub
implicit none
real(dl) nfw,r
nfw=(dlog((rs_sat+r)/rs_sat)-r/(r+rs_sat))/(dlog(1.d0+c_sat)-c_sat/(1+c_sat))-rand2
end function nfw


       
    
program hod_selection
      use Precision
      use Pub
      implicit none
      real(dl), dimension(:),allocatable :: mass, x, y, z, rvir, rs 
      real(dl) :: m_cut, sigma, kappa, m1, alpha !baseline hod parameters
      real(dl) :: N_cent, N_sat
      real(dl) :: rand, rand_mu, rand_phi
      real(dl) :: r,dx,dy,dz
      real(dl) :: costheta,sintheta,cosphi,sinphi
      real(dl) :: g_x, g_y, g_z
      real(dl) :: starttime, stoptime
      real(dl), external :: nfw 
      
      integer :: i,j,count,flag,count_sat
      
      m_cut=10.d0**13.08d0
      m1=10.d0**14.06d0
      sigma=0.98d0
      kappa=1.13d0
      alpha=0.9d0
      
      count=0
      count_sat=0
      call init_random_seed() 
      call cpu_time(starttime) 
      allocate(mass(nhalo))
      allocate(x(nhalo))
      allocate(y(nhalo))
      allocate(z(nhalo))
      allocate(rvir(nhalo))
      allocate(rs(nhalo))
      open(unit=100, file='/home/hyzhang/Documents/data/hod/MDhalos.txt')
     
      do i=1,nhalo

          read(100,*) mass(i),rvir(i),rs(i),x(i),y(i),z(i)
          rvir(i)=rvir(i)/1000.d0
          rs(i)=rs(i)/1000.d0
      end do
      close(100)
      call cpu_time(stoptime)
      write(*,*) 'Reading Halo File Complete',' reading time =',stoptime-starttime
      open(unit=200, file='/home/hyzhang/Documents/data/hod/central_output_2.txt')
      open(unit=300, file='/home/hyzhang/Documents/data/hod/satellite_output_2.txt')
      open(unit=400, file='/home/hyzhang/Documents/data/hod/all_output_2.txt')
        do i=1,nhalo
          N_cent = 1.d0/2.d0*erfc(dlog10(m_cut/mass(i))/(dsqrt(2.d0)*sigma))
          N_sat = ((mass(i)-kappa*m_cut)/m1)**alpha
          call RANDOM_NUMBER(rand)
          if(rand.lt.N_cent) then
              if(x(i).le.2500.d0.and.y(i).le.2500.d0.and.z(i).le.2500.d0) then
              count=count+1
              g_x=x(i)
              g_y=y(i)
              g_z=z(i)
              end if
              write(200,"(3F9.3)") g_x,g_y,g_z
              write(400,"(3F9.3)") g_x,g_y,g_z
!          if(g_x.ge.2500.d0.or.g_y.ge.2500.d0.or.g_z.ge.2500.d0) write(*,*) i,count

          end if
          if(N_sat.ge.1.d0) then
               !write(*,*) N_sat
               count=count+int(N_sat)
               count_sat=count_sat+int(N_sat)
               do j=1,int(N_sat)
               call RANDOM_NUMBER(rand2)
               rs_sat=rs(i)
               c_sat=rvir(i)/rs(i)
               call bisection(nfw,0.d0,rvir(i),1.d-4,r,flag)
!               write(*,*) r
               call RANDOM_NUMBER(rand_mu)
               costheta=rand_mu*2.d0-1.d0
               sintheta=dsqrt(1-rand_mu**2)
               call RANDOM_NUMBER(rand_phi)
               cosphi=dcos(rand_phi*2*PI)
               sinphi=dsin(rand_phi*2*PI)
               dx=r*sintheta*cosphi
               dy=r*sintheta*sinphi
               dz=r*costheta
               g_x=x(i)+dx               
               g_y=y(i)+dy               
               g_z=z(i)+dz
               if(g_x.gt.2500.d0) g_x=g_x-2500.d0
               if(g_y.gt.2500.d0) g_y=g_y-2500.d0
               if(g_z.gt.2500.d0) g_z=g_z-2500.d0
               if(g_x.lt.0.d0) g_x=g_x+2500.d0
               if(g_y.lt.0.d0) g_y=g_y+2500.d0
               if(g_z.lt.0.d0) g_z=g_z+2500.d0
               end do               
               write(300,"(3F9.3)") g_x,g_y,g_z
               write(400,"(3F9.3)") g_x,g_y,g_z
           end if
              
      end do
      write(*,*) 'satellite number = ',count_sat
      write(*,*) 'total number = ',count
      close(300)
      close(400)
      deallocate(mass)
      deallocate(x)
      deallocate(y)
      deallocate(z)
      deallocate(rs)
      deallocate(rvir)
      call cpu_time(stoptime)
      write(*,*) 'Writing Galaxy File Complete',' total time =',stoptime-starttime
    end program hod_selection
    
              

