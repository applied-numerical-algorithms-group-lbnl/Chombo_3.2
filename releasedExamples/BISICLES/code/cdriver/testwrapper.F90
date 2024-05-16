
#include "cdriverconstants.h"

program fwrapper
  ! Want to know how to run BISICLES from FORTRAN 90? Avoid learning
  ! C++ ? Follow the example below, but hang your head in shame
#ifdef CH_MPI
  use mpi
#endif	
  implicit none
  character(len=25) :: file

  integer instance_id,  it, nt, max_step
  integer, dimension(1:2) :: dims, boxlo, boxhi
  real(kind=8) :: dx, max_time,start_time, snow2ice
  integer, parameter :: nx = 64, ny = 96
  real(kind=8), dimension(:,:), allocatable :: smb, bmbf, bmbg, usrf, seb, melange, snow
  integer ixlo,ixhi,iylo,iyhi
  
  integer :: rank,nrank,ierr,comm

#ifdef CH_MPI
  call MPI_Init ( ierr )
  call MPI_Comm_rank ( mpi_comm_world, rank, ierr )
  call MPI_Comm_size ( mpi_comm_world, nrank, ierr )
  comm = mpi_comm_world
#else
  rank = 0
  nrank = 1
  comm = 0
#endif

  !example domain decomposition scheme : in serial, one block of data and in parallel two blocks
  !owned by ranks 0,1, with all the other processors contributing nothing 

  if (nrank.eq.1) then
     ixlo = 0
     ixhi = nx - 1
     iylo = 0
     iyhi = ny - 1
  else 
     if (rank.eq.0) then
         ixlo = 0
         ixhi = nx/2 - 1
         iylo = 0
         iyhi = ny - 1
      else if (rank.eq.1) then
         ixlo = nx/2
         ixhi = nx - 1
         iylo = 0
         iyhi = ny - 1
     else
        !no data : specify a block outside the ice sheet domain, which bisicles then ignores
        ixlo = nx + 1
        iylo = ny + 1
        ixhi = nx + 1
        iyhi = ny + 1
     end if

  end if	

  !memory allocation
  allocate (smb(ixlo:ixhi,iylo:iyhi))
  allocate (seb(ixlo:ixhi,iylo:iyhi))
  allocate (usrf(ixlo:ixhi,iylo:iyhi))
  allocate (bmbf(ixlo:ixhi,iylo:iyhi))
  allocate (bmbg(ixlo:ixhi,iylo:iyhi))
  allocate (snow(ixlo:ixhi,iylo:iyhi))
  !name of the configuration file
  file = "inputs.pigv5.1km.l1l2.l1" 
  
  !create an instance: instance_id will be used to identify it from now on
  call f_bisicles_new_instance(instance_id, file, len(file), comm)
 
  
  !now set up some uniform mesh data to read from / write to the interface
  !mesh spacing
  dx = 4.0e+3
  !data dimensions
  dims(1) = nx
  dims(2) = ny
  !lower left corner on the BISICLES level with resolution dx
  boxlo(1) = ixlo
  boxlo(2) = iylo
  !top right corner
  boxhi(1) = ixhi
  boxhi(2) = iyhi
     
     
  !tell BISICLES to read various data from the arrays smb,bmbg,bmbf,seb
  !surface flux  (m/a)
  ! Automatically added to whatever is specified in the configuration file
  smb = 1.0d0 +  dble(mod(rank,2)) ! a strange smb that makes it clear which rank is being read from
  call f_bisicles_set_2d_data(instance_id, smb, BISICLES_FIELD_SURFACE_FLUX, dx, dims, boxlo, boxhi)

  !grounded ice basal flux, m/a
  ! Automatically added to whatever is specified in the configuration file
  bmbg = -0.0d0
  call f_bisicles_set_2d_data(instance_id, bmbg, BISICLES_FIELD_GROUNDED_ICE_BASAL_FLUX, dx, dims, boxlo, boxhi)

  !floating ice basal flux m/a
  ! Automatically added to whatever is specified in the configuration file
  bmbf = -10.0d0
  call f_bisicles_set_2d_data(instance_id, bmbf, BISICLES_FIELD_FLOATING_ICE_BASAL_FLUX, dx, dims, boxlo, boxhi)

  !energy flux, J/a
  !surfaceHeatBoundaryData.type = MemoryLevelData must be set in the configuration file
  seb = 1.0d+7
  call f_bisicles_set_2d_data(instance_id, seb, BISICLES_FIELD_SURFACE_HEAT_FLUX, dx, dims, boxlo, boxhi)
  

  !At this point, we have given BISICLES as much data as it needs to run.
  !So, initialize it, at which point the AMR velocity problem gets solved (or a checkpoint loaded)
  !After this step we can still change the data in (say) smb but we cannot change its location
  call f_bisicles_init_instance(instance_id)


  
  !now do some time stepping
  nt = 3
  max_step = 0
  max_time = 0.0d0
  start_time = 0.0d0
  do it = 1, nt

     
     !snow -> ice conversion (testing a nasty UKESM feature)
     snow = 55.0 - dble(it-1)*5.0
     snow2ice = 50.0
     call f_bisicles_push_thin_ice(instance_id, snow, snow2ice, dx, dims, boxlo, boxhi)
     write(*,*) 'snow to ice sum,max ', sum(snow), maxval(snow)
     
     !read the surface elevation into a regular array (usrf)
     call f_bisicles_get_2d_data(instance_id, usrf, BISICLES_FIELD_SURFACE_ELEVATION, dx, dims, boxlo, boxhi)
     
     !re-compute the smb - just change the data. this is obviously a weird field that is supposed 
     !to show the domain decomposed reads/writes are working (look at surfaceThicknessBalance in the output)
     smb = usrf / 1000.0 +  dble(mod(rank,2))

    
     
     !advance in time
     start_time = max_time
     max_time = start_time + 1.0d0; !advance by one year
     max_step = max_step + 10; !unless it takes too long, in which case give up
     write(*,*)   max_time, max_step 
     call f_bisicles_advance(instance_id, start_time,  max_time, max_step)

     !ice <- snow conversion (a nasty UKESM feature)
     call f_bisicles_pop_thin_ice(instance_id, snow, snow2ice, dx, dims, boxlo, boxhi)
     write(*,*) 'ice to snow sum,max ', sum(snow), maxval(snow)
     
  end do

  call f_bisicles_get_2d_data(instance_id, melange, BISICLES_FIELD_MELANGE_THICKNESS, dx, dims, boxlo, boxhi)
  
  !free any memory allocated on the C++ side
  call f_bisicles_free_instance(instance_id)
#ifdef CH_MPI
  call MPI_Finalize ( ierr )
#endif

end program fwrapper
