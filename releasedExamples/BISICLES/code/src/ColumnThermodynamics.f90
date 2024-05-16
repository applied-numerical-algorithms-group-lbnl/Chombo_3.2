! column_themodynamics.f90 : f90 subroutines used to carry 
! out thermodynamics calculations that do not depend on the horizontal mesh
! Copyright (C) 2012-2015 BISICLES contributors
!
! Please refer to Copyright.txt, in Chombo's root directory.
!
module column_thermodynamics
  implicit none

  !physical constants, all SI units
  real(kind=8) sput,  rhoi , rhoo,  grav , shci, lhci, coni , conm,  pmlt, trpt
  ! seconds per time unit, density of ice, density of water, acceleration due to gravity 
  ! specific heat capacity of ice, latent heat of fusion of water, thermal conductivity of cold ice,  
  ! regularizing conductivity for temperate ice, 
  ! factor for dependence of melting point on pressure, triple point of water (K)
  
  !real(kind = 8)  water_fraction_drain = 0.01d0, water_fraction_max = 0.05d0, drain_factor = 0.5d0,  ground_water_drain_factor = 0.01d0
  real(kind = 8)  water_fraction_drain, water_fraction_max, water_drain_factor
  real(kind = 8)  deprecated_till_water_drain_factor, till_water_max
  real(kind = 8)  floating_base_max_heat_flux
  integer, parameter :: groundedmaskval = 1, floatingmaskval = 2, openseamaskval = 4, openlandmaskval = 8
  real(kind=8), parameter :: temp_eps = 1.0e-3, max_delta_energy = 200.0
  real(kind=8), parameter :: zero_debug = 0.0d0
  real(kind=8), parameter :: temperature_min = 200.0d0
contains

  real(kind=8) function fo_upwind(u,sminus,splus)
    !first-order upwind scheme.
    real(kind=8), intent(in) :: u,sminus,splus
    if (u .ge. 0.0) then
       fo_upwind = u * sminus
    else
       fo_upwind = u * splus 
    end if
    return
  end function fo_upwind




  subroutine tdmasolve(n,x,a,b,c,d)
    ! solve the n x n tridiagonal system
    ! overwriting the coefficients a,b,c and d
    ! diagonal elements must not be zero
    ! |b1 c1 0 0 0  ...    0| |x1|   |d1|
    ! |a2 b2 c2 0 0 0  ... 0| |x2| = |d2|
    ! |0  a3 b3 c3 0 ......0| |x3|   |d3|
    ! etc      
    integer, intent(in) :: n
    real(kind=8), dimension(1:n), intent(inout) :: x,a,b,c,d
    !locals
    integer :: i
    real(kind=8) :: t

    c(1) = c(1)/b(1)
    d(1) = d(1)/b(1)

    do i = 2,n
       t = b(i) - a(i)*c(i-1)
       c(i) = c(i) / t
       d(i) = (d(i)-a(i)*d(i-1))/t
    end do

    x(n) = d(n)

    do i = n-1, 1,-1
       x(i) = (d(i)-c(i)*x(i+1))
    end do


  end subroutine tdmasolve

  subroutine sigma_advect(rhs, energy, senergy, benergy, husig ,  & 
       fsig, dt, mask, n)
    !increment rhs with interlayer advection
    !TODO : replace first order upwind fluxes with something nice
    integer :: n
    real(kind=8), dimension(1:n), intent(inout) :: rhs
    real(kind=8), dimension(1:n), intent(in) :: energy ! internal energy density
    real(kind=8), dimension(1:n+1), intent(in) :: husig ! thickness * u^{\sigma}
    real(kind=8), dimension(1:n+1), intent(in) :: fsig ! sigma at layer faces
    real(kind=8), intent(in) :: senergy, benergy, dt ! surface & basal temperature, time step
    integer, intent(in) :: mask
    !locals
    integer l
    real(kind=8) flux

    !flux across the surface : Dirichlett condition
    flux = dt * (fo_upwind(husig(1) ,senergy,energy(1)) )
    rhs(1) = rhs(1) + flux / (fsig(2)-fsig(1))

    do l = 2, n
       flux = dt * (fo_upwind(husig(l) ,energy(l-1),energy(l)) )
       rhs(l-1) = rhs(l-1) - flux / (fsig(l)-fsig(l-1))
       rhs(l) = rhs(l) + flux / (fsig(l+1)-fsig(l))
    end do

    !flux across the base
    flux =  dt * (fo_upwind(husig(n+1) ,energy(n),  benergy) )
    rhs(n) = rhs(n) - flux / (fsig(n+1)-fsig(n))

  end subroutine sigma_advect
  


  subroutine moisture_transport(energy, epmp, cdsig, n)
    !!!  redistribute moisture such that w = water_fration_drain up to the water level
    implicit none
    integer :: n
    real(kind=8), dimension(1:n), intent(inout) :: energy, epmp, cdsig
    !locals
    real(kind=8), dimension(1:n) :: energy2,latent,emax
    real(kind=8) :: excess,emin
    integer :: i
    
    !record for a conservation check
    energy2 = energy
    ! decompose into sensible and latent heat
    where (energy .gt. epmp)
       latent = energy - epmp
       energy = epmp
    elsewhere
       latent = 0.0d0
    end where
    ! total moisture (per unit thickness)
    excess = sum(latent*cdsig)
    ! distribute latent heat (allowing refreeze where C*Tpmp - C*T <  Lw)
    do i = n, 1, -1
       if (excess .gt. 0.0d0) then
          energy(i) = energy(i) + min(excess/cdsig(i),lhci*water_fraction_drain)
          excess = excess - lhci*water_fraction_drain*cdsig(i)
       end if
    end do
   
    !put any excess in the bottom layer 
    excess = sum( (energy2 - energy) * cdsig)
    energy(n) = energy(n) + excess/cdsig(n)
    
    !bodge upper and lower limits. ugly, but hopefully only important in stupid zones
    emax = epmp + lhci * water_fraction_max
    where (energy .gt. emax)
       energy = emax
    end where
    emin = shci * temperature_min
    where (energy .lt. emin)
       energy = emin
    end where

   
  end subroutine moisture_transport
  
  subroutine fo_diffusive_advance(energy,tillwaterdepth,senergy,sflux,sdiric, & 
       benergy,bflux,rhs,thckold,thcknew,tillwaterdrainfactor,fsig,dt,mask,n)
    ! solve the equation 
    ! H*E - dt * 1/H * d/dsigma (q) = Ho*Eo + rhs
    !
    !Non-advective fluxes follow Aschwanden 2012, with
    !     { - (Ki + Kr) dE/dz               if E < Epmp
    ! q = {   
    !     { - Ki d/dz(Epmp) - Kr dE/dz  if E > Epmp
    !   
    ! Ki = ki/ci and Kr = km/ci
    ! ki is the ice conductivity, km is the (small) mositure conductivity ci is the specific heat capacity of ice, 
    ! Tpmp is the pressure melting point, Epmp = ci * Tpmp, kr is a regularization parameter 
    ! with either a Dirichlett boundary condition E = sE  
    ! or a flux boundary condition - chi(z) dE/dsigma = sflux at sigma = 0
    ! and a flux boundary condition - chi(sigma) dE/dsigma = bflux at sigma = 1
    !
    !
    ! (first order - backward Euler integration)
    implicit none
    integer :: n
    real(kind=8), intent(inout) :: tillwaterdepth,senergy,benergy,sflux,bflux
    real(kind=8), intent(in) :: thckold,thcknew,dt,tillwaterdrainfactor
    real(kind=8), dimension(1:n), intent(in) :: rhs
    real(kind=8), dimension(1:n), intent(inout) :: energy
    real(kind=8), dimension(1:n+1), intent(in) :: fsig ! sigma, conductivity at layer faces
    integer, intent(in) :: mask
    logical, intent(in) :: sdiric
    !locals
    real(kind=8), dimension(1:n) :: a,c,b,r,a2,c2,b2,r2 ! TDMA coefficients
    real(kind=8), dimension(1:n) :: csig,cdsig,epmp,drain,energy2,energy0,emax,latent
    real(kind=8), dimension(1:n+1) :: fdsig,  kcc, ktt, kct, ktc
    real(kind=8) :: bmb, eps,k, kr,bepmp,sepmp, kb,  emin, sigmaw, excess, lhf, tmp
    
    integer :: i
    logical :: btemperate
   
    eps = 1.0e-3 
    energy0 = energy

    do i = 1,n
       cdsig(i) = fsig(i+1)-fsig(i)
       csig(i) = 0.5d0*(fsig(i+1)+fsig(i))
    end do
    do i = 2,n
       fdsig(i) = 0.5d0 * cdsig(i) + 0.5d0 * cdsig(i-1)
    end do
    fdsig(1) = cdsig(1)
    fdsig(n+1) = cdsig(n)

    !pressure melting point
    epmp = shci * (trpt - pmlt * 0.5d0 * (thckold+thcknew) * csig * rhoi * grav) ! at layer midpoints
    sepmp = shci * trpt !at surface
    bepmp = shci * (trpt - pmlt * 0.5d0 * (thckold+thcknew) * rhoi * grav) !at base

    ! conduction coeffs
    kr =  2.0 * dt / (thcknew+thckold) * conm/(shci * rhoi) *  sput
    k = 2.0 * dt / (thcknew+thckold) * coni/(shci * rhoi) *  sput - kr
    !face centered conduction coeffcients
    kcc = 0.0d0
    ktt = 0.0d0 
    kct = 0.0d0
    ktc = 0.0d0

    !regularization
    kcc(1:n+1) =  kcc(1:n+1) + kr /  fdsig(1:n+1)

    !cold-cold interfaces
    where ( (energy(1:n-1).lt.epmp(1:n-1)) .and. (energy(2:n).lt.epmp(2:n)) )
       kcc(2:n) = kcc(2:n) + k /  fdsig(2:n)
    end where

    !cold-temperate interfaces
    where ( (energy(1:n-1).lt.epmp(1:n-1)) .and. (energy(2:n).ge.epmp(2:n)) )
       kct(2:n) =  kct(2:n) +  k /  fdsig(2:n)
    end where
    
    !temperate-cold interfaces
    where ( (energy(1:n-1).ge.epmp(1:n-1)) .and. (energy(2:n).lt.epmp(2:n)) )
       ktc(2:n) =  ktc(2:n) +  k /  fdsig(2:n)
    end where

    !temperate-temperate interfaces
    where ( (energy(1:n-1).ge.epmp(1:n-1)) .and. (energy(2:n).ge.epmp(2:n)) )
       ktt(2:n) =  ktt(2:n) +  k /  fdsig(2:n)
    end where

    !TDMA coefficients
    !interior ice layers 2 to n-1
    do i = 2, n-1
       a(i) =  (- kcc(i)   - kct(i)  )/cdsig(i)
       c(i) =  (- kcc(i+1) - ktc(i+1))/cdsig(i)
       b(i) = thcknew + (kcc(i) + ktc(i) + kcc(i+1) + kct(i+1) )/cdsig(i)
       r(i) = thckold * energy(i) + rhs(i) &
            - epmp(i-1)/cdsig(i)* (-ktt(i) - ktc(i)) &
            - epmp(i+1)/cdsig(i)* (-ktt(i+1) - kct(i+1)) &
!            - epmp(i)/cdsig(i) * (ktt(i) + ktc(i) + ktt(i+1) + kct(i+1))
            - epmp(i)/cdsig(i) * (ktt(i) + kct(i) + ktt(i+1) + ktc(i+1))
    end do

    ! top layer
    i = 1
    if (sdiric) then
       !Dirichlett boundary condition, E(0) + E(1) = 2 Es. 
       if (senergy.lt.sepmp) then
          !cold surface
          if ( energy(i).lt.epmp(i) ) then
             !cold-cold
             kcc(i) = kcc(i) + k /  fdsig(i)
          else
             !cold-temperate
             kct(i) =  kct(i) +  k /  fdsig(i)
          end if
       else
          !temperate surface
           if ( energy(i).lt.epmp(i) ) then
             !temperate-cold
             ktc(i) = ktc(i) + k /  fdsig(i)
          else
             !temperate-temperate
             ktt(i) =  ktt(i) +  k /  fdsig(i)
          end if   
       end if
      
       a(i) =  0.0d0
       c(i) =  (- kcc(i+1) - ktc(i+1))/cdsig(i)
       b(i) = thcknew + (2.0d0 * kcc(i) + ktc(i) + kct(i) & 
            + kcc(i+1) + kct(i+1) )/cdsig(i)
       r(i) = thckold * energy(i) + rhs(i) &
            + 2.0d0*senergy/cdsig(i)  * (kcc(i) + kct(i)) & 
            - sepmp/cdsig(i)* (-ktt(i) - ktc(i)) &
            - epmp(i+1)/cdsig(i)* (-ktt(i+1) - kct(i+1)) ! &
            !- epmp(i)/cdsig(i) * (ktt(i) + kct(i) + ktt(i+1) + ktc(i+1)) &
                
    else
       !Flux boundary condition
       a(i) =  0.0d0
       c(i) =  (- kcc(i+1) - ktc(i+1))/cdsig(i)
       b(i) = thcknew + (kcc(i+1) + kct(i+1) )/cdsig(i)
       r(i) = thckold * energy(i) + rhs(i) &
            - epmp(i+1)/cdsig(i)* (-ktt(i+1) - kct(i+1)) &
            - epmp(i)/cdsig(i) * ( ktt(i+1) + ktc(i+1)) & 
            +  dt * sflux/cdsig(1)
    end if

    ! bottom layer, flux boundary condition
    i = n 
   
    btemperate = (mask.eq.floatingmaskval)
    if (btemperate) then
       !basal flux will be kb/dt * (bepmp - energy(i))
       kb = 2.0d0 * (k + kr) / ( fdsig(i+1))
       kb = min(floating_base_max_heat_flux,kb)       
       a(i) = ( -kcc(i) - kct(i) )/cdsig(i)
       c(i) = 0.0d0
       b(i) = thcknew + (kcc(i) + ktc(i) + kb)/ cdsig(i) 
       r(i) = thckold * energy(i) + rhs(i) &
            - epmp(i-1)/cdsig(i)* (-ktt(i) - ktc(i)) &
            - epmp(i)/cdsig(i) * (ktt(i) + kct(i)) &
            + bepmp/cdsig(i) * kb
    else
       a(i) = ( -kcc(i) - kct(i) )/cdsig(i)
       c(i) = 0.0d0
       b(i) = thcknew + (kcc(i) + ktc(i))/cdsig(i) 
       r(i) = thckold * energy(i) + rhs(i) &
            - epmp(i-1)/cdsig(i)* (-ktt(i) - ktc(i)) &
            - epmp(i)/cdsig(i) * (ktt(i) + kct(i)) &
            +  dt * bflux/cdsig(i)
       
    end if
    
    a2 = a
    b2 = b
    c2 = c
    r2 = r
    call tdmasolve(n,energy2,a2,b2,c2,r2)
    a2 = a
    b2 = b
    c2 = c
    r2 = r


    !drainage.  \fixme
    where (energy2 .gt. epmp + water_fraction_drain * lhci)
       drain = water_drain_factor * dt
    elsewhere
       drain = 0.0d0
    end where
   
    !re-solve with drainage sink.
    b = b + thcknew*drain
    r = r + thcknew*drain*epmp
    call tdmasolve(n,energy,a,b,c,r)

    !bodge upper and lower limits. ugly, but hopefully only important in stupid zones
    emax = epmp + lhci * water_fraction_max
    where (energy .gt. emax)
       energy = emax
    end where
    emin = shci * temperature_min
    where (energy .lt. emin)
       energy = emin
    end where

    drain = (energy2 - energy)/lhci 

    
    !compute surface energy density or heat flux
    if (sdiric) then
       sflux = 0.0 ! FIX ME
    else
       senergy = 1.5d0 * energy(1) - 0.5d0 * energy(2)
    end if

    !compute basal energy density and modify heat flux
    if (btemperate) then
       bflux = kb/dt * (bepmp - energy(i))
       benergy = bepmp
    else
       kb = coni/(shci * rhoi) *  sput
       benergy = energy(n) + bflux / kb * thcknew * 0.5d0 * (fsig(n+1)-fsig(n))
    end if

    !basal melt rate (negative outward)
    bmb = - sum(drain) / dt * cdsig(n) * thcknew
 
    !update till water fraction (Crank-Nicolson)
    tmp = 0.5d0 * tillwaterdrainfactor * dt 
    tillwaterdepth = (tillwaterdepth * (1.0d0 - tmp) - bmb*dt)/( 1.0d0 + tmp)

    !ice shelf / open sea regions - set till water to max (seems more sensible than zero)
    if ( (mask.eq.floatingmaskval).or.(mask.eq.openseamaskval) )  then
       tillwaterdepth = till_water_max
    end if
    !limit
    tillwaterdepth = max(0.0d0,min(tillwaterdepth, till_water_max))
    return

  end subroutine fo_diffusive_advance

end module column_thermodynamics


subroutine column_thermodynamics_set_constants( asput , arhoi , arhoo, agrav , ashci, alhci, aconi , aconm, apmlt, atrpt)
  use column_thermodynamics
  implicit none
  real(kind=8) asput,  arhoi , arhoo,  agrav , ashci, alhci, &
       aconi , aconm, apmlt, atrpt
  sput = asput
  rhoi = arhoi 
  rhoo = arhoo
  grav = agrav 
  shci = ashci
  lhci = alhci
  coni = aconi 
  conm = aconm
  pmlt = apmlt
  trpt = atrpt
  
end subroutine column_thermodynamics_set_constants

subroutine column_thermodynamics_set_water_constants ( &
     a_water_fraction_drain, a_water_fraction_max, &
     a_water_drain_factor, a_till_water_drain_factor, &
     a_till_water_max, a_floating_base_max_heat_flux  )
  use column_thermodynamics
  implicit none
  
  real(kind=8) a_water_fraction_drain, a_water_fraction_max, &
       a_water_drain_factor, a_till_water_drain_factor, &
       a_till_water_max, a_floating_base_max_heat_flux
  
  water_fraction_drain = a_water_fraction_drain
  water_drain_factor = a_water_drain_factor
  water_fraction_max = a_water_fraction_max
  deprecated_till_water_drain_factor = a_till_water_drain_factor
  till_water_max = a_till_water_max
  floating_base_max_heat_flux = a_floating_base_max_heat_flux
end subroutine column_thermodynamics_set_water_constants
  
subroutine column_thermodynamics_update_internal_energy(energy, tillwaterdepth, &
     senergy, sflux, sdiric, benergy, mask, &
     bflux, rhs, thckold, thcknew,tillwaterdrainfactor, usig, fsig, dsig, time, dt, n)

  !update the mid-layer internal energy density (energy), 
  !the surface  energy density (senergy), and the basal energy density (benergy)
  !
  !Advective fluxes are computed from the cross-layer contravariant component of the velocity, usig
 
  !
  !There are n layers, with sigma = 0 at the surface, 
  !sigma = sigma(i) at the middle of layer i and sigma = 1
  !at the base of the ice
  ! 
  use column_thermodynamics

  implicit none
  integer :: n
  real(kind=8), dimension(1:n), intent(inout) :: energy,rhs,dsig
  real(kind=8), dimension(1:n+1), intent(inout) :: usig,fsig
  real(kind=8), intent(inout) :: tillwaterdepth,senergy,benergy,bflux,sflux
  real(kind=8), intent(in) :: time,dt,thckold,thcknew,tillwaterdrainfactor
  integer, intent(in) :: mask
  logical, intent(in) :: sdiric
  !locals

  real(kind=8) :: dtcfl,tthcknew,tthckold,tt,melt,bepmp
  real(kind=8), dimension(1:n) :: rhsl,drain,moist,csig
  ! face-centered conductivity, internal energy, pressure-melting-point, additional heat flux
  integer l,nt,it,i,npicard,ipicard
  

  if (thcknew.gt.1.0d0 .and. thckold.gt.1.0d0) then
     
     do i = 1,n
        csig(i) = 0.5*(fsig(i+1)+fsig(i))
     end do
    
     !work out a stable time step

     !advection cfl
     dtcfl = dt
     do i = 1,n
        dtcfl = min(dtcfl,(thckold+thcknew)*(fsig(i+1)-fsig(i))/(abs(usig(i+1)) + abs(usig(i))))
     end do

!!$     !a simplified criterion for temperate ice conduction, which is essentially explicit diffusion
!!$     !hopefully the scaling will not be a major issue. TODO revisit and limit to cases where
!!$     !there is temperate ice above the bottom layer if it proves a problem
!!$     !using real(kind=8) tt as a workspace
     tt = (fsig(2)-fsig(1))**2
     do i = 2,n
        tt = min(tt, (fsig(i+1)-fsig(i))**2)
     end do
     tt = tt * min(thckold,thcknew)**2 / (2.0* coni)

     !maybe have a factor < 1 here? So far has not been needed
     dtcfl = min(dtcfl, tt)

     nt = ceiling(dt/dtcfl)
     dtcfl = dt/dble(nt)

     tt = 0.0d0
     do it = 1,nt

        rhsl = rhs/dble(nt)

        !modify rhs to include interlayer advection across constant sigma faces
        !
        
        bepmp = shci * (trpt - pmlt * 0.5 * (thckold+thcknew) * rhoi * grav)
        benergy = energy(n) + (energy(n)-energy(n-1))/(csig(n)-csig(n-1)) &
             * (fsig(n+1) - csig(n))
        call sigma_advect(rhsl, energy, senergy, benergy, usig , &
             fsig, dtcfl, mask,  n)

        !compute vertical diffusion & update the internal energy
        tthckold = tt/dt * thcknew + (1.0d0-tt/dt) *thckold
        tt = tt + dtcfl
        tthcknew = tt/dt * thcknew + (1.0d0-tt/dt) *thckold

        
        call fo_diffusive_advance(energy, tillwaterdepth, senergy, sflux, sdiric, benergy, bflux, rhsl, &
             tthckold, tthcknew,tillwaterdrainfactor,fsig,dtcfl,mask,n)
        
     end do
  else
     !no ice
     senergy = max(senergy, shci*(trpt - temp_eps - 0.0))
     if (sdiric) then
        sflux = bflux
     end if
     energy = senergy 
     benergy = senergy
     if ( (mask.eq.floatingmaskval).or.(mask.eq.openseamaskval) ) then
        tillwaterdepth = till_water_max
     end if
  end if

  return

end subroutine column_thermodynamics_update_internal_energy





subroutine column_compute_sigma_vel( husig,ux,uy, divhu, &
     dsig, n, dht, smb, bmb)

  !compute the contravariant \sigma component of the velocity u^{\sigma} = u.e^{\sigma}
  !This appears in cross-layer advection

  !This means solving a first order ODE w' = f with
  !known w(s) and w(b), which is overdetermined unless 
  !w(s) - w(b) = integral(f dz,s,b) - which should 
  !be the case.
  
  implicit none
  integer, intent(in) :: n
  real(kind=8), dimension(1:n+1), intent(in) :: ux,uy ! horizontal velocity at layer faces
  real(kind=8), dimension(1:n+1), intent(out) :: husig ! h*sigma velocity at   layer faces
  real(kind=8), dimension(1:n), intent(in) :: divhu ! horizontal part of div(Hu)
  real(kind=8), dimension(1:n), intent(in) :: dsig ! dsigma at layer midoints
  real(kind=8), intent(in) :: dht ! dh/dt
  real(kind=8), intent(in) :: smb, bmb ! surface & basal mass balance rates
  !locals
  integer :: i

  !velocity at the base
  husig(n+1) =  - bmb! + gft(n+1)

  !integration
  do i = n , 1, -1
     husig(i) =  husig(i+1) + (divhu(i) + dht)*dsig(i) 
  end do

  return

end subroutine column_compute_sigma_vel

  
subroutine column_compute_z_vel(uz,uzs,ux,uy, divhu, &
     fsig, dsig, n, dsx, dhx, dsy, dhy, dst, dht, & 
     smb, bmb)
  !compute the vertical velocity from the divergence
  !of horizontal velocity and the kinematic boundary
  !conditions. Also compute the contravariant \sigma 
  !component of the velocity u^{\sigma} = u.e^{\sigma}
  !which appears in cross-layer advection
  !
  !This means solving a first order ODE w' = f with
  !known w(s) and w(b), which is ill-posed unless 
  !w(s) - w(b) = integral(f dz,s,b) - which should 
  !be the case. To check this, we compute both 
  !w(1) - the surface velocity computed by integration, and
  !ws - the surface velocity boundary condition
  implicit none
  integer, intent(in) :: n
  real(kind=8), dimension(1:n+1), intent(in) :: ux,uy ! horizontal velocity at layer faces
  real(kind=8), dimension(1:n+1), intent(out) :: uz ! vertical velocity at layer faces
  real(kind=8), intent(out) :: uzs ! uz at surface by kinematic bc
  real(kind=8), dimension(1:n+1), intent(in) :: fsig ! sigma at layer faces
  real(kind=8), dimension(1:n), intent(in) :: divhu ! horizontal part of div(Hu)
  real(kind=8), dimension(1:n), intent(in) :: dsig ! sigma, dsigma at layer midoints
  real(kind=8), intent(in) :: dsx,dsy,dhx,dhy,dst,dht !H, ds/dx, ds/dy, dH/dx, dH/dy. dH/dt
  real(kind=8), intent(in) :: smb, bmb ! surface & basal mass balance rates
  !locals
  integer :: i
  real(kind=8), dimension(1:n+1) :: gfx,gfy,gft
  real(kind=8) :: uzref
  !geoemtry factors
  gfx(1:n+1) = (dsx-fsig(1:n+1)*dhx)
  gfy(1:n+1) = (dsy-fsig(1:n+1)*dhy)
  gft(1:n+1) = (dst-fsig(1:n+1)*dht)

  !velocity at the base
  uz(n+1) =  gft(n+1) + ux(n+1)*gfx(n+1) + uy(n+1)*gfy(n+1) + bmb

  !velocity at the surface
  uzs = dst + ux(1)*dsx + uy(1)*dsy - smb

  !integration
  do i = n , 1, -1
     uz(i) = uz(i+1) - divhu(i)*dsig(i) & 
          + (gfx(i)*ux(i) - gfx(i+1)*ux(i+1)) &
          + (gfy(i)*uy(i) - gfy(i+1)*uy(i+1)) 
  end do

  return

end subroutine column_compute_z_vel
