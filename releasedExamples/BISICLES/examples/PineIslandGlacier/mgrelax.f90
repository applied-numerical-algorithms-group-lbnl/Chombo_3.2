! multigrid for constant ceoefficient u - l^2 u'' = v with Neumann BCs
! on uniform, square nodal grids
module mgrelax
  implicit none
contains

  subroutine warning()
    integer code
    code = 0
    return
  end subroutine warning


 subroutine restrictnc(fine,coarse,nfx,nfy,ncx,ncy)
    ! restriction : full weighting
    integer nfx,nfy,ncx,ncy
    real(kind=8) :: fine(0:nfx+1,0:nfy+1)
    real(kind=8) :: coarse(0:ncx+1,0:ncy+1), t

    integer i,j, ii, jj
    !ensure that nf = 2*(nc)
    if (nfx .ne. 2*(ncx)) then
       write(*,*) ' nfx =  != 2*(ncx)'
       stop
    end if
    if (nfy .ne. 2*(ncy)) then
       write(*,*) ' nfx =  != 2*(ncx) '
       stop
    end if

!    write(*,*) 'in restrict'

    coarse(0,:) = 0.0
    coarse(ncx+1,:) = 0.0
    coarse(:,0) = 0.0
    coarse(:,ncy+1) = 0.0
     do j = 1, ncy
        do i = 1, ncx
          ii = 2*i
          jj = 2*j
          t = fine(ii,jj) + fine(ii,jj+1) + fine(ii+1,jj) + fine(ii+1,jj+1)
          coarse(i,j) = t / 4.0

       end do
    end do
    return
  end subroutine restrictnc


  subroutine prolongnc(fine,coarse,nfx,nfy,ncx,ncy)
    !prolongation : linear interpolation from coarse to fine, nodal
    integer nfx,nfy,ncx,ncy
    real(kind=8) :: fine(0:nfx+1,0:nfy+1)
    real(kind=8) :: coarse(0:ncx+1,0:ncy+1)

    integer i,j, ii, jj
    !ensure that nf = 2*(nc)
    if (nfx .ne. 2*(ncx)) then
       write(*,*) ' nfx =  != 2*(ncx) +1 '
       stop
    end if

!    write(*,*) 'in prolong' 
    if (nfy .ne. 2*(ncy)) then
       write(*,*) ' nfx =  != 2*(ncx) +1 '
       stop
    end if


    fine(nfx,nfy) = coarse(ncx,ncy)
    do j = 1, ncy
       do i = 1, ncx
          ii = 2*i
          jj = 2*j
          fine(ii,jj) = coarse(i,j)
          fine(ii+1,jj) = coarse(i,j)
          fine(ii,jj+1) = coarse(i,j) 
          fine(ii+1,jj+1) = coarse(i,j)
       end do
    end do
    return
  end subroutine prolongnc

  subroutine apply(u,lu,a,ni,nj)
    integer :: ni,nj,i,j
    real(kind=8) :: a
    real(kind=8) :: u(0:ni+1,0:nj+1)
    real(kind=8) :: lu(0:ni+1,0:nj+1)
    real(kind=8) :: b

!    write(*,*) 'in apply'
    !Neumann BC
    u(0,:) =  u(2,:)
    u(ni+1,:) =  u(ni-1,:)
    u(:,0) =  u(:,2)
    u(:,nj+1) =  u(:,nj-1)
    lu(0,:) = 0.0
    lu(ni+1,:) = 0.0
    lu(:,nj+1) = 0.0
    lu(:,0) = 0.0

     do j = 1, nj
        do i = 1, ni
           lu(i,j) = u(i,j)*(1+4.0*a) &
                - a * (u(i+1,j) +  u(i-1,j) + u(i,j-1) +  u(i,j+1))

        end do
     end do
     i = 0
     return
  end subroutine apply

  subroutine resid(u,v,r,a,ni,nj)
    integer :: ni,nj,i,j
    real(kind=8) :: a
    real(kind=8) :: u(0:ni+1,0:nj+1)
    real(kind=8) :: v(0:ni+1,0:nj+1)
    real(kind=8) :: r(0:ni+1,0:nj+1)
    real(kind=8) :: b
    
!    write(*,*) 'in resid'

    !Neumann BC
    u(0,:) =  u(2,:)
    u(ni+1,:) =  u(ni-1,:)
    u(:,0) =  u(:,2)
    u(:,nj+1) =  u(:,nj-1)
    r(0,:) = 0.0
    r(ni+1,:) = 0.0
    r(:,nj+1) = 0.0
    r(:,0) = 0.0

     do j = 1, nj
        do i = 1, ni
           r(i,j) = v(i,j) - (u(i,j)*(1+4.0*a) &
                - a * (u(i+1,j) +  u(i-1,j) + u(i,j-1) +  u(i,j+1)))

        end do
     end do
     return
  end subroutine resid

  subroutine gsrb(u,v,a,ni,nj,n)

    ! n gsrb passes. u is assumed to have one layer of ghost cells,
    ! a = lambda**2 / dx **2 
    integer :: n,  ni, nj
    real(kind=8) :: u(0:ni+1,0:nj+1)
    real(kind=8) :: v(0:ni+1,0:nj+1)
    real(kind=8) :: a
    
    integer :: iter, i, j , color
    
!    write (*,*) 'in gsrb, n =', n

    do iter = 1,n

       !Neumann BC
       u(0,:) =  u(2,:)
       u(ni+1,:) =  u(ni-1,:)
       u(:,0) =  u(:,2)
       u(:,nj+1) =  u(:,nj-1)

       ! red
       do j = 1, nj
          i = 1 +  mod(j,2)
          do while (i .le. ni)
             u(i,j) = (v(i,j) + a * (u(i+1,j) +  u(i-1,j) &
                  +  u(i,j + 1) + u(i,j -1)))/ (1.0 + 4.0*a)
             

             i = i + 2
          end do
       end do

      

       ! black
       do j = 1, nj
          i = 1 +  mod(j + 1,2)
          
          do while (i .le. ni)
             u(i,j) = (v(i,j) + a * (u(i+1,j) +  u(i-1,j) &
                  +  u(i,j + 1) + u(i,j -1)))/ (1.0 + 4.0*a)
             
             i = i + 2
             
          end do
       end do

       
    end do
    return
  end subroutine gsrb

  recursive subroutine vcycle(u,v,a,nx,ny,depth,nsmooth)

    integer :: nx, ny, nsmooth, depth
    real(kind=8) :: u(0:nx+1,0:ny+1)
    real(kind=8) :: v(0:nx+1,0:ny+1)
    real(kind=8) :: a, resnorm


    integer ncx,ncy,i,j
    real(kind=8) :: ac ! a on the coarse mesh
    real(kind=8), allocatable :: du(:,:)
    real(kind=8), allocatable :: res(:,:)
    real(kind=8), allocatable :: lu(:,:)
    real(kind=8), allocatable :: duc(:,:)
    real(kind=8), allocatable :: resc(:,:)

!    write (*,*) 'vcycle -- depth = ', depth
!    write(*,*) '    nsmooth = ', nsmooth
    
    if (( mod(nx,2) .ne. 0) .or. ( mod(ny,2) .ne. 0 )) then
       ! no more coarsenings possible
!       write(*,*) 'at bottom: nx,ny = ', nx, ny
!       write(*,*) 'mod(nx,2) = ', mod(nx,2)
!       write(*,*) 'mod(ny,2) = ', mod(ny,2)
       call gsrb(u,v,a,nx,ny,nsmooth*nx*ny)
       return
    end if

    ncx = (nx + 1)/2
    ncy = (ny + 1)/2

    if (depth > 1) then

       allocate(lu(0:nx+1,0:ny+1))
       allocate(du(0:nx+1,0:ny+1))
       allocate(res(0:nx+1,0:ny+1))
       allocate(duc(0:ncx+1,0:ncy+1))
       allocate(resc(0:ncx+1,0:ncy+1))
       
       call gsrb(u,v,a,nx,ny,nsmooth)
       call resid(u,v,res,a,nx,ny)
       call restrictnc(res,resc,nx,ny,ncx,ncy)
       ac = a / 4.0
       duc = 0.0
       call vcycle(duc , resc, ac, ncx, ncy,  depth -1, nsmooth)
       call prolongnc(du, duc , nx, ny,  ncx, ncy ) 
       u = u + du
       call gsrb(u,v,a,nx,ny,nsmooth)
    else
       ! at the bottom of the cycle
       call gsrb(u,v,a,nx,ny,nsmooth*nx*ny)
    end if
    return
  end subroutine vcycle


end module mgrelax
