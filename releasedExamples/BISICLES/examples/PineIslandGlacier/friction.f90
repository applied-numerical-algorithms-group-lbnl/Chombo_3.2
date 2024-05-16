module ncdump
 
 use netcdf
 use mgrelax
 implicit none
  
 contains
 subroutine nccheck(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       stop 2
    end if
  end subroutine nccheck

 subroutine ncsavebike(x,y,thck,topg,beta, uvel, vvel, velc, &
      divuh, divuhc, temp,bheatflux, bdiss, n,m,upn,file)
    ! create a netcdf file on a bisicles grid
    ! given bisicles grid cell-centered data
    implicit none
    character(len=*),intent(in) :: file
    integer,intent(in) :: n,m,upn
    real(kind=8), dimension(1:n) ,intent(in):: x
    real(kind=8), dimension(1:m) ,intent(in):: y
    real(kind=8), dimension(1:n,1:m) ,intent(in):: thck, topg, beta, uvel, vvel,&
         velc, divuh, divuhc, bheatflux, bdiss
    real(kind=8), dimension(1:n,1:m,1:upn) ,intent(in):: temp
    character(len=10) :: tempname
    real(kind=8) :: time;
    integer i,j,nc_id, var_id, cell_dim_id(2), cell_3dim_id(3)
    
    
    call nccheck( nf90_create(file,  NF90_CLOBBER, nc_id) )

    !cell centred dimensions
    call nccheck( nf90_def_dim(nc_id, "x", n , cell_dim_id(1)) )
    call nccheck( nf90_def_dim(nc_id, "y", m , cell_dim_id(2)) )
    call nccheck( nf90_def_dim(nc_id, "sigma", upn , cell_3dim_id(3)) )

    !dimension definitions
    call nccheck( nf90_def_var(nc_id, "x", nf90_real8, cell_dim_id(1), var_id))
    call nccheck( nf90_put_att(nc_id, var_id, "units", "m") )
    call nccheck( nf90_def_var(nc_id, "y", nf90_real8, cell_dim_id(2), var_id))
    call nccheck( nf90_put_att(nc_id, var_id, "units", "m") )
    
    cell_3dim_id(1:2) = cell_dim_id
    call nccheck( nf90_def_var(nc_id, "sigma", nf90_real8, cell_3dim_id(3), var_id))
    
    !cell centered array defintions
    call nccheck( nf90_def_var(nc_id, "thk", nf90_real8, cell_dim_id, var_id) )
    call nccheck( nf90_def_var(nc_id, "topg", nf90_real8, cell_dim_id, var_id) )
    call nccheck( nf90_def_var(nc_id, "beta", nf90_real8, cell_dim_id, var_id) )
    call nccheck( nf90_def_var(nc_id, "xvel", nf90_real8, cell_dim_id, var_id) )
    call nccheck( nf90_def_var(nc_id, "yvel", nf90_real8, cell_dim_id, var_id) )
    call nccheck( nf90_def_var(nc_id, "velc", nf90_real8, cell_dim_id, var_id) )
    call nccheck( nf90_def_var(nc_id, "divuh", nf90_real8, cell_dim_id, var_id) )
    call nccheck( nf90_def_var(nc_id, "divuhc", nf90_real8, cell_dim_id, var_id) )

    do i = 1,upn
       write(tempname,'("temp",i6.6)') i-1
       write(*,*) tempname
       call nccheck( nf90_def_var(nc_id, tempname, nf90_real8, cell_dim_id, var_id) )
    end do

    call nccheck( nf90_def_var(nc_id, "bheatflux", nf90_real8, cell_dim_id, var_id) )
    call nccheck( nf90_def_var(nc_id, "bdiss", nf90_real8, cell_dim_id, var_id) )
    !end of definition section
    call nccheck( nf90_enddef(nc_id) )

    !cell centred x and y
    call nccheck( nf90_inq_varid(nc_id, "x", var_id) )
    call nccheck( nf90_put_var(nc_id, var_id , x) )
    call nccheck( nf90_inq_varid(nc_id, "y", var_id) )
    call nccheck( nf90_put_var(nc_id, var_id , y) )

    !cell centred arrays
    call nccheck( nf90_inq_varid(nc_id, "thk", var_id) )
    call nccheck( nf90_put_var(nc_id, var_id , thck) )
    call nccheck( nf90_inq_varid(nc_id, "topg", var_id) )
    call nccheck( nf90_put_var(nc_id, var_id , topg) )
    call nccheck( nf90_inq_varid(nc_id, "beta", var_id) )
    call nccheck( nf90_put_var(nc_id, var_id , beta) )
    call nccheck( nf90_inq_varid(nc_id, "xvel", var_id) )
    call nccheck( nf90_put_var(nc_id, var_id , uvel) )
    call nccheck( nf90_inq_varid(nc_id, "yvel", var_id) )
    call nccheck( nf90_put_var(nc_id, var_id , vvel) )
    call nccheck( nf90_inq_varid(nc_id, "velc", var_id) )
    call nccheck( nf90_put_var(nc_id, var_id , velc) )
    call nccheck( nf90_inq_varid(nc_id, "divuh", var_id) )
    call nccheck( nf90_put_var(nc_id, var_id , divuh) )
    call nccheck( nf90_inq_varid(nc_id, "divuhc", var_id) )
    call nccheck( nf90_put_var(nc_id, var_id , divuhc) )
    
    do i = 1,upn
       write(tempname,'("temp",i6.6)') i-1
       call nccheck( nf90_inq_varid(nc_id, tempname, var_id) )
       call nccheck( nf90_put_var(nc_id, var_id , temp(:,:,i) ) )
    end do

    call nccheck( nf90_inq_varid(nc_id, "bheatflux", var_id) )
    call nccheck( nf90_put_var(nc_id, var_id , bheatflux) )

    call nccheck( nf90_inq_varid(nc_id, "bdiss", var_id) )
    call nccheck( nf90_put_var(nc_id, var_id , bdiss) )

    call nccheck( nf90_close(nc_id) )
  end subroutine ncsavebike


end module ncdump




subroutine applydrop(drop,usrf,lsrf,topg,thck,ewn,nsn)
  implicit none
  integer, intent(in) :: ewn,nsn
  real(kind=8),parameter :: rhoi = 918.0d0, rhoo = 1028.0d0, grav = 9.81d0, cm = 1.0e-5
  real (kind = 8), dimension(1:ewn,1:nsn), intent(inout) :: topg, lsrf, usrf,  thck
  real (kind = 8), intent(in) :: drop
  
  usrf = usrf - drop
  topg = topg - drop
  lsrf = lsrf - drop
  ! some fiddling to ensure the division into sheet/shelf is
  ! as the data says, and usrf is correct also. topography may change
  ! near the grounding line.
  usrf = max(usrf,0.0)
  lsrf = max(lsrf,topg)
  thck = max(usrf - lsrf,0.0)
  where (lsrf .gt. topg + 1)
     thck = usrf / (1.0 - rhoi/rhoo) 
     lsrf = usrf - thck
     topg = min(topg,lsrf - cm)
  elsewhere
     thck = usrf - lsrf
     where (topg+thck.le.(1.0 - rhoi/rhoo)*thck)
        thck = - rhoo * topg / rhoi + cm
     end where
  end where
  !now use the bisicles expressions, 
  usrf = max(topg + thck, thck * (1.0 - rhoi/rhoo))
  lsrf = usrf - thck
 
end subroutine applydrop

subroutine yflip(a,n,m)
  integer :: i,n,m
  real(kind=8), dimension(1:n,1:m) :: a,b
  do i = 1,m
     b(1:n,m-i+1) = a(1:n,i)
  end do
  a = b
end subroutine yflip

subroutine readone(a,file,ewn,nsn,ewnin,nsnin,ewoff,nsoff)
  integer ewn,nsn,ewnin,nsnin,ewoff,nsoff,ew,ns,unin
  real(kind=8), dimension(1:ewn,1:nsn) :: a
  real(kind=8), dimension(1:ewnin,1:nsnin) :: in
  character(len=*),intent(in) :: file
  integer w,e,s,n
  
  w = 1 + ewoff
  e = ewn + w - 1
  s = 1 + nsoff
  n = nsn + s - 1
  
  unin = 9
  open(unin,file=file,status='old')
  read(unin,*) ((in(ew,ns), ew=1,ewnin), ns=1,nsnin)
  close(unin)
  a(1:ewn,1:nsn) = in(w:e,s:n)
  call yflip(a,ewn,nsn)
end subroutine readone

subroutine readdata(usrf,lsrf,topg,uvel,vvel,ewn,nsn)
 integer, parameter :: ewnin = 292, nsnin = 465, ewoff = 34, nsoff = 48
 integer :: ewn,nsn
 real(kind=8), dimension(1:ewn,1:nsn) :: usrf,lsrf,topg,uvel,vvel
 character(len=12) :: file

 !file = 'topg_5pg.dat'
 call readone(topg,'topg_5pg.dat',ewn,nsn,ewnin,nsnin,ewoff,nsoff)
 !file = 'usrf_5pg.dat'
 call readone(usrf,'usrf_5pg.dat',ewn,nsn,ewnin,nsnin,ewoff,nsoff)
 !file = 'lsrf_5pg.dat'
 call readone(lsrf,'lsrf_5pg.dat',ewn,nsn,ewnin,nsnin,ewoff,nsoff)
 !file = 'uvel_5pg.dat'
 call readone(uvel,'uvel_5pg.dat',ewn,nsn,ewnin,nsnin,ewoff,nsoff)
 !file = 'vvel_5pg.dat'
 call readone(vvel,'vvel_5pg.dat',ewn,nsn,ewnin,nsnin,ewoff,nsoff)
 vvel = -vvel

 

end subroutine readdata

subroutine computeA(flwa,temp,thck,n,ewn,nsn,upn)
  implicit none
  integer :: ewn,nsn,upn
!  real(kind=8), parameter :: Cgas = 0.16612d0, Q = 7.88d4, R=8.31d0, K=1.17d0, B0=2.207, &
!       tR = 273.39d0, enhance=1.0d0
  real(kind=8) :: Cgas, Q, R, K , B0, tR, enhance, rhog, icepmeltf
  real(kind=8), dimension(1:ewn,1:nsn) :: thck
  real(kind=8), dimension(1:upn) :: sigma
  real(kind=8), dimension(1:ewn,1:nsn,1:upn) :: temp,flwa
  real(kind=8) :: n, factor, theta, tmp, pmp, pfactor
  integer i, j, up
  
  Cgas = 0.16612d0
  Q = 7.88d4
  R=8.31d0
  K=1.17d0
  B0=2.207
  tR = 273.39d0
  enhance=1.0d0
  factor = enhance* ( (1.0d0/B0)**n )
  pfactor = 9.81d0 * 918d0 * 9.7456d-8 
  do i = 1,ewn
     do j = 1,nsn
        do up = 1,upn
        theta = (temp(i,j,up)) * (1.0d0 + sigma(up) * thck(i,j) * pfactor)
        pmp = 273.0
        theta = (min(pmp,theta))

        tmp = 3.0*Cgas/(tR - theta)**K - (Q/R) / theta
        flwa(i,j,up) = factor * exp(tmp)
     end do
     end do
  end do
  return
end subroutine computeA

subroutine growbeta(beta,typ,ewn,nsn,niter,maxseabeta)
  implicit none
  integer, parameter :: typ_grounded = 0, typ_iceshelf = 1, typ_opensea = 2
  integer :: ewn,nsn,niter
  real(kind=8), dimension(1:ewn,1:nsn) :: beta
  real(kind=8) :: maxseabeta
  integer, dimension(1:ewn,1:nsn) :: typ, typ2
  integer :: ew,ns,iter,k,bw,be,bn,bs

  do iter = 1,niter
     do ew = 2,ewn-1
        do ns = 2,nsn-1
           if (typ(ew,ns).ne.typ_grounded) then
              beta(ew,ns) = 0.25 * &
                   (beta(ew+1,ns) + beta(ew-1,ns) &
                   + beta(ew,ns+1) + beta(ew,ns-1))
              beta(ew,ns) = max(beta(ew,ns),maxseabeta)

              if (typ(ew-1,ns).eq.typ_grounded) then
                 beta(ew,ns) = min(beta(ew,ns),beta(ew-1,ns))
              end if
              if (typ(ew+1,ns).eq.typ_grounded) then
                 beta(ew,ns) = min(beta(ew,ns),beta(ew+1,ns))
              end if
              if (typ(ew,ns-1).eq.typ_grounded) then
                 beta(ew,ns) = min(beta(ew,ns),beta(ew,ns-1))
              end if
              if (typ(ew,ns+1).eq.typ_grounded) then
                 beta(ew,ns) = min(beta(ew,ns),beta(ew,ns+1))
              end if

           end if
        end do
     end do
  end do
  return
end subroutine growbeta



subroutine smooth(v,a,dx,ewn,nsn)
  ! smooth a field by solving v(out) - a^2 v(out)" = v(in)
  use mgrelax
  implicit none
  integer :: ewn,nsn
  integer :: iter,niter
  real (kind = 8), dimension(1:ewn,1:nsn) :: v
  real (kind = 8), dimension(0:ewn+1,0:nsn+1) :: umg,vmg,dumg,rmg
  real (kind=8) :: mu, dx, a, resNorm
  
  mu = (a**2) / (dx**2)
  write(*,*)  mu
  rmg = 0.0d0
  vmg = 0.0d0
  vmg(1:ewn,1:nsn) = v
  umg = vmg
  call resid(umg,vmg,rmg,mu,ewn,nsn)
  resNorm = sum(abs(rmg) )
  write(*,*) 'initial res= ', resNorm
  niter = 10
  do iter = 1, niter
     dumg = 0.0
     call vcycle(dumg,rmg,mu,ewn,nsn,10,4)
     umg = umg + dumg
     call resid(umg,vmg,rmg,mu,ewn,nsn)
     resNorm = sum(abs(rmg) )
     write (*,*) 'MG iter', iter, " norm = ", resNorm
     
  end do
  v = umg(1:ewn,1:nsn) 

end subroutine smooth

subroutine inventtemp(temp,thck,sigma,ewn,nsn,upn)
  
  !make up a 3d temperature field
  implicit none
  integer :: ewn,nsn,upn,k,i,j
  real(kind=8), dimension(1:upn) :: sigma
  real(kind=8), dimension(1:ewn,1:nsn) :: thck
  real(kind=8), dimension(1:ewn,1:nsn,1:upn) :: temp

  real(kind=8), dimension(1:ewn,1:nsn) :: stemp
  real(kind=8) :: ts,tb, dt

  stemp = 268.0d0
  dt = 0.0
  temp = 268.0d0
  !do i = 1,ewn
  !!   do j = 1,nsn
  !      ts = stemp(i,j)
  !      tb = min(272.0, ts + dt)
  !      do k = 1,upn     
  !         temp(i,j,k) = (tb-ts)*sigma(k)**2 + ts
  !      end do
  !   end do
 ! end do

end subroutine inventtemp


program t
   use ncdump
  implicit none
  
  integer, parameter :: ewn = 256, nsn = 384, upn = 10
  
  real(kind=8),parameter :: scyr = 31556926.0d0, rhoi = 918.0d0, rhoo = 1028.0d0,  &
       grav = 9.81d0, dx = 1.0d+3, dy = dx,  tol = 1.0d-3, minbeta = 20.0, &
       maxbeta = 1.0e+4, seabeta = 100.0, shelfbeta = 1.0e-5,cm = 1.0e-2, &
       glen_n = 3.0, maxseabeta = 100.0, lambda = 4.0e+3
  
  real (kind = 8), dimension(1:ewn,1:nsn) :: topg, lsrf, usrf, thck, uvel, vvel, velc,divuh, &
       divuhc, beta, betar, dsx, dsy, umod, umodsia, bheatflux, bdiss, tmp
  
  real (kind = 8), dimension(1:ewn,1:nsn,1:upn) :: temp, flwa

  character(len=64) :: filename
  
  integer , dimension(1:ewn,1:nsn) ::  typ
 
  integer, parameter :: typ_grounded = 0, typ_iceshelf = 1, typ_opensea = 2
   
  real (kind = 8), dimension(1:ewn) :: x
  real (kind = 8), dimension(1:nsn) :: y
  real(kind=8), dimension(1:upn) :: sigma
  real(kind=8) :: ds,drop,vlambda,vdx
  integer iter, radius, i, j, k, ew, ns

  typ = typ_grounded ! 0 for grounded ice

  x(1) = -1707000.0d0 + 0.5d0*dx
  do ew = 2,ewn
     x(ew) = x(ew-1) + dx
  end do
  y(1) =  -384000.0d0 + 0.5d0*dx
  do ns = 2,nsn
     y(ns) = y(ns-1) + dy
  end do

  call readdata(usrf,lsrf,topg,uvel,vvel,ewn,nsn)

  ds = upn
  ds = 1.0/ds
  sigma(1) = ds/2.0
  do i = 2,upn
     sigma(i) = sigma(i-1) + ds
  end do

 
  

  !raw data - need to set up typ,mask,thck
  typ = typ_grounded ! grounded ice
  thck = usrf - lsrf

  drop =  15.2! Value from Anne Le Brocq
  !call applydrop(drop,usrf,lsrf,topg,thck,ewn,nsn)

  thck = thck - drop
  where (thck.lt.cm)
     thck = cm
  end where
  !now use the bisicles expressions, 
  usrf = max(topg + thck, thck * (1.0 - rhoi/rhoo))
  lsrf = usrf - thck


  where (lsrf - topg .gt. cm)
     typ = typ_iceshelf ! ice shelf
  end where

  where ((lsrf - topg .gt. 120.0) .and. (thck .lt. 120.0))
     typ = typ_opensea ! open sea
     thck = 0.0
     usrf = 0.0
     lsrf = 0.0
  end where


  


  call inventtemp(temp,thck,sigma,ewn,nsn,upn)
 
  

  !assume same 1/ (s * sigma) everywhere
  velc = 1.0d0;
  divuhc = 1.0d0;

  ! at the moment, we don't know what div(uh) should be, but
  ! we don't mind if it is large (and positive) in shelf regions 
  divuh = 0.0d0; 
  where (typ .ne. typ_grounded)
      divuhc = 0.0d0
  end where

  where (uvel .lt. -9000.0)
     velc = 0.0d0
     uvel = 0
     vvel = 0
  end where

  where (typ .eq. typ_opensea)
     uvel = 0.0
     vvel = 0.0
     velc = 0.0d0
  end where
  radius = 2
  !set velc = 0 close to open sea
  do ew = 1+radius,ewn-radius
     do  ns = 1+radius,nsn-radius
        if (typ(ew,ns).eq.typ_opensea) then
           do i = -radius,radius
              do j = -radius,radius
                 velc(ew+i,ns+j) = 0.0;
              end do
           end do
        end if
     end do
  end do


  !expand the region where velc = 0 by a few cells
  tmp = velc
  do iter = 1,1
     do ew = 2,ewn-1
        do ns = 2,nsn-1
           if (velc(ew,ns).gt.0.1) then
           tmp(ew,ns) = 0.25 * ( &
                velc(ew-1,ns) + velc(ew+1,ns) &
                +velc(ew,ns-1) + velc(ew,ns+1))
           end if
        end do
     end do
  end do
  velc = tmp

  where (velc.lt.0.95)
     velc = 0.0d0
  end where

  umod = sqrt(uvel**2 + vvel**2)

  !compute surface gradients
  dsx = 0
  dsy = 0
  
  where ((typ (2:ewn-1,1:nsn) .eq. typ_grounded) &
       .and. (typ (3:ewn,1:nsn) .eq. typ_grounded) &
       .and. (typ (1:ewn-2,1:nsn) .eq. typ_grounded))
     dsx(2:ewn-1,1:nsn) = (usrf(3:ewn, 1:nsn) - usrf(1:ewn-2, 1:nsn)) &
          /(2.0d0 * dx)
  end where


  where ((typ (1:ewn,2:nsn-1) .eq. typ_grounded) &
       .and. (typ (1:ewn,3:nsn) .eq. typ_grounded) &
       .and. (typ (1:ewn,1:nsn-2) .eq. typ_grounded))
     dsy(1:ewn,2:nsn-1) = (usrf(1:ewn,3:nsn) - usrf(1:ewn,1:nsn-2)) &
          /(2.0d0 * dx)
  end where

  
  !A in Glen's law
  call computeA(flwa,temp,thck,glen_n,ewn,nsn,upn)
  flwa = 4.0e-17
  !basal traction coefficient
  where (typ.eq.typ_grounded)
     where (velc.gt.0.0)
        umodsia = (2.0d0*flwa(1:ewn,1:nsn,5) *thck**(glen_n+1)) / (glen_n+1.0d0) * (rhoi* grav)**glen_n  & 
             * sqrt(dsx**2 + dsy**2)**glen_n
        beta = (1 + rhoi* grav * thck * sqrt(dsx**2 + dsy**2)) / ( max(1.0d-6,umod - umodsia))

         !beta is the effective drag f(u) (Tb = f(u) * u), and wewant to run with Tb = beta * |u|^1/3
        !beta = beta * umod**(2.0/3.0)
        !beta
     elsewhere
        beta = 100.0
     end where
  end where

  beta = min(beta,maxbeta * 1.0e+4**(2.0/3.0))
  beta = max(beta,minbeta)

  call growbeta(beta,typ,ewn,nsn,20,maxseabeta)
  vlambda = 4.0e+3
  vdx = dx
  beta = log(beta)
  call smooth(beta,vlambda,vdx,ewn,nsn)
  beta = exp(beta)

!!$  where (typ.eq.typ_grounded)
!!$     beta = 1000
!!$  end where
!!$
  where (typ.eq.typ_opensea)
     beta = 100
  end where

  where (typ.eq.typ_iceshelf)
     beta = minbeta
  end where

  call growbeta(beta,typ,ewn,nsn,1,maxseabeta)

 
  
  !set basal heat flux to 100 mW m^2. Units requires are J a^-1 m^-2
  bheatflux = 100*1e-3 * 365 * 24 * 3600

  !set basal dissipation to beta * (umod - 100.0)**2
  where (typ.eq.typ_grounded)
     bdiss = beta * (max(1.0d-6,umod - 100.0))**2
  elsewhere
     bdiss = 0.0d0
  end where

  filename =  'pig-bisicles-1km.nc'
  call ncsavebike(x,y,thck,topg,beta,uvel,vvel,velc,divuh,divuhc,temp,bheatflux,bdiss,ewn,nsn,upn,filename)

  
end program t 
