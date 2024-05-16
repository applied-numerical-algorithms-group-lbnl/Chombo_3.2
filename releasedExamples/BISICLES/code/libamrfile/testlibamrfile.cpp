 #ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

/*===========================================================================
 testlibamrfile.c
 test the libamrfile functions
==========================================================================*/
#include <iostream>
#include <cstdlib>
#include "libamrfile.H"
#ifdef CH_USE_PETSC
#include "petsc.h"
#endif 
#ifdef CH_MPI
#include <mpi.h>
#endif


int main(int argc, char* argv[]) {

  int status;

#ifdef CH_USE_PETSC
  status = PetscInitialize(&argc, &argv,PETSC_NULL,PETSC_NULL); CHKERRQ(status);
#else
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif 
#endif // end petsc conditional

  int amr_id;
  char file[64] = "plot.amundsen.2d.hdf5";
  char outfile[64] = "plot.mod.amundsen.2d.hdf5";
  amr_read_file(&status, &amr_id, file);

  if (status != 0)
    {
      exit(status);
    }
  
   int n_level;
  amr_query_n_level(&status, &n_level, &amr_id);
  if (status != 0)
    {
      exit(status);
    }

  std::cout << " n_level = " <<  n_level << std::endl;

  for (int lev = 0; lev < n_level; lev++)
    {

      int n_fab;
      amr_query_n_fab(&status, &n_fab, &amr_id, &lev);
      if (status != 0)
	{
	  exit(status);
	}
  
      std::cout << " level = " << lev << " n_fab = " <<  n_fab << std::endl;

      for (int i = 0; i < n_fab; i++)
	{
	  int nx, ny, ncomp;
	  amr_query_fab_dimensions_2d(&status, &nx, &ny, &ncomp, &amr_id, &lev, &i);

	  if (status != 0)
	    {
	      exit(status);
	    } 

	  
	  int nghost = 1;
	  int nxg = nx + 2*nghost;
	  int nyg = ny + 2*nghost;
	  double *fab_data = new double[nxg*nyg];
	  double *x_data = new double[nxg];
	  double *y_data = new double[nyg];
	  int comp = 0;
	  
	  amr_read_fab_data_2d(&status, fab_data , x_data, y_data, &amr_id, &lev, &i, &comp, &nghost);


	  std::cout << " level = " << lev << " fab = " <<  i 
		    << " nx = " << nx << " ny = " << ny << " ncomp = " << ncomp
	            << " " << x_data[0] << " <= x <= " << x_data[nx-1] 
		    << " "  << y_data[0] << " <= y <= " << y_data[ny-1] 
		    << std::endl;
	  
	
	  if (status != 0)
	    {
	      exit(status);
	    } 
	  
	  for (int ix = 0; ix < nxg; ix++)
	    {
	      for (int iy = 0; iy < nyg; iy++)
		{
		  fab_data[iy*(nxg) + ix] = double(iy);
		}
	    }

	  
	  amr_write_fab_data_2d(&status, fab_data,  &nx, &ny , &amr_id, &lev, &i, &comp, &nghost);

	  delete[] fab_data; 
	  delete[] x_data;
	  delete[] y_data;
	}
    }

  {
    int lo[2]; 
    int hi[2]; 
    int lev = 0;
    amr_query_domain_corners(&status, lo, hi, &amr_id, &lev);
    std::cout << "level = " << lev 
	      << "lo = ("  << lo[0] << "," << lo[1] << ")"
	      << "hi = ("  << hi[0] << "," << hi[1] << ")"
	      << std::endl;
    lev++;
    amr_query_domain_corners(&status, lo, hi, &amr_id, &lev);
    std::cout << "level = " << lev 
	      << "lo = ("  << lo[0] << "," << lo[1] << ")"
	      << "hi = ("  << hi[0] << "," << hi[1] << ")"
	      << std::endl;
    
    lev++;
    amr_query_domain_corners(&status, lo, hi, &amr_id, &lev);
    std::cout << "level = " << lev 
	      << "lo = ("  << lo[0] << "," << lo[1] << ")"
	      << "hi = ("  << hi[0] << "," << hi[1] << ")"
	      << std::endl;
    
    int nx = hi[0] - lo[0] + 1;
    int ny = hi[1] - lo[1] + 1;
    double *fab_data = new double[nx*ny];
    double *x_data = new double[nx];
    double *y_data = new double[ny];
    int comp = 0;
    int order = 1;
    amr_read_box_data_2d(&status, fab_data , x_data, y_data, &amr_id, &lev, lo, hi, &comp, &order);
    
    delete[] fab_data; 
    delete[] x_data;
    delete[] y_data;
    
  }
  
  
  amr_write_file(&status, &amr_id, outfile);


  amr_free(&status, &amr_id);
  if (status != 0)
    {
      exit(status);
    } 

#ifdef CH_USE_PETSC
  status = PetscFinalize(); CHKERRQ(status);
#else
#ifdef CH_MPI
  MPI_Finalize();
#endif // mpi conditional
#endif // petsc conditional


  return 0;
}

