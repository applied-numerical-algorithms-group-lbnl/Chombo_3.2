#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef CWRAPPER_H
#define CWRAPPER_H
#include "NamespaceHeader.H"

/** 
    \file 
    \brief  Provides functions inside a extern "C" block that can be called
    by C (or FORTRAN) program to control BISICLES, pass data to it,
    and retrieve data from it. 

 */
extern "C"
{
#ifndef CH_MPI
#define MPI_Comm int
#endif
  ///create a new BISICLES instance, return an integer key
  /**
     bisicles_new_instance (or its fortran compatible wrapper f_bisicles_new_instance)
     will normally be the first called by program making use of BISICLES. It creates
     a new BISICLES instance in memory and provides a unique key by which
     the calling program can identify it. Calls to  bisicles_new_instance will
     typically be followed by calls to bisicles_set_2d_data, one call to 
     bisicles_init_instance, several calls to bisicles_advance and bisicles_get_2d_data,
     and a final call to f_bisicles_free_instance
     
     @param[out] instance_key integer key, used to look up the instance
     @param[in]  input_fname configuration file name
     @param[in]  comm for MPI compiled versions: allows the MIP communicator to be set by the calling program
   */
  void bisicles_new_instance(int *instance_key, const char *input_fname, MPI_Comm comm);

  ///free the BISICLES instance specified by instance_key
  void bisicles_free_instance(int *instance_key);
  
  /// Instruct BISICLES to read 2d data from an array at data_ptr.
  /**
      @param[in] instance_key pointer to integer key, used to look up the instance
      @param[in] data_ptr  array data managed by the caller.
      @param[in] field     purpose of the data
      @param[in] dx        mesh spacing 
      @param[in] dims      grid (i,j) domain dimensions, ie domain size is dx*(i,j)
      @param[in] boxlo     lower left (i,j) corner of the data supplied by __this__ rank
      @param[in] boxhi     upper right (i,j) corner of the data supplied by __this__ rank

      BISICLES will attempt to read FORTRAN ordered 2d data from data_ptr[boxlo[0]:boxlo[0],boxhi[0]:boxhi[0]]
      when computing certain fields, for example the surface mass balance. The calling program
      is free to decompose the domain over multiple MPI ranks: BISICLES will then transfer data between 
      ranks as it requires (which in inevitable in all but he simplest case where the calling program 
      and BISICLES choose identical domain decompositions). Taken in total, the data should
      cover the whole of a uniform grid with mesh spacing dx and extent dims.

      What BISICLES actually does with the data depends on the field parameter. Supported
      values are defined in cdriverconstants.h and include 
       - BISICLES_FIELD_SURFACE_FLUX,
       - BISICLES_FIELD_BASAL_FLUX, 
       - BISICLES_FIELD_FLOATING_ICE_BASAL_FLUX 
       - BISICLES_FIELD_GROUNDED_ICE_BASAL_FLUX 
       - BISICLES_FIELD_SURFACE_TEMPERATURE 
       - BISICLES_FIELD_SURFACE_HEAT_FLUX (for surface energy balance)
       - BISICLES_FIELD_BASAL_HEAT_FLUX (for geothermal heat) 
       - BISICLES_FIELD_TOPOGRAPHY_FLUX (for GIA, sea level change)
      

      
  */
  void bisicles_set_2d_data(int *instance_key,  double *data_ptr, const int *field, 
			  const double *dx, const int *dims, 
			  const int *boxlo, const int *boxhi);
  /**
     \brief Instruct BISICLES to read 2d geometry (ice thickness and bedrock eleveation) data from an arrays at thck_data_ptr, topg_data_ptr
     

     BISICLES will attempt to read FORTRAN ordered 2d data from 
     thck_data_ptr[boxlo[0]:boxlo[0],boxhi[0]:boxhi[0]] and  topg_data_ptr[boxlo[0]:boxlo[0],boxhi[0]:boxhi[0]]
     when computing the initial thickness and bedrock elevation, and (optionally) to recompute the
     bedrock following a regrid.  The calling program
     is free to decompose the domain over multiple MPI ranks: BISICLES will then transfer data between 
     ranks as it requires (which in inevitable in all but he simplest case where the calling program 
     and BISICLES choose identical domain decompositions). Taken in total, the data should
     cover the whole of a uniform grid with mesh spacing dx and extent dims.

     This special case exists because thck_data_ptr and *topg_data_ptr must be defined on the same BoxLayout,
     ruling out the use of bisicles_set_2d_data.

     __Do not__ attempt to use this function to modify the bedrock elevation due to e.g. GIA or sea level
     change. BISICLES only reads topg_data_ptr to set the initial condition and when the mesh is
     recomputed. Use bisicles_set_2d_data to set a rate of change of bedrock elevation.
     
     
     @param[in] instance_key pointer to integer key, used to look up the instance
     @param[in] thck_data_ptr  ice thickness array data, managed by the caller.
     @param[in] topg_data_ptr  bedrock elevation array data, managed by the caller.
     @param[in] dx        mesh spacing 
     @param[in] dims      grid (i,j) domain dimensions, ie domain size is dx*(i,j)
     @param[in] boxlo     lower left (i,j) corner of the data supplied by __this__ rank
     @param[in] boxhi     upper right (i,j) corner of the data supplied by __this__ rank


  */
  void bisicles_set_2d_geometry(int *instance_key,  double *thck_data_ptr, double *topg_data_ptr,
				const double *dx, const int *dims, 
				const int *boxlo, const int *boxhi);

  /**
     \brief Fill a rectangular array with 2D data from a BISICLES instance
     
     BISICLES will attempt to write FORTRAN ordered 2d data into data_ptr[boxlo[0]:boxlo[0],boxhi[0]:boxhi[0]].     
     Data will be averaged, interpolated or copied as needed: all three operations are conservative unless
     othrewise specified.  The calling program
     is free to decompose the domain over multiple MPI ranks: BISICLES will then transfer data between 
     ranks as it requires (which in inevitable in all but he simplest case where the calling program 
     and BISICLES choose identical domain decompositions). Taken in total, the data should
     cover the whole of a uniform grid with mesh spacing dx and extent dims.

     @param[in] instance_key pointer to integer key, used to look up the instance
     @param[in] data_ptr  array data, managed by the caller (BISICLES will not change the data_ptr, but will modify *data_ptr).
     @param[in] dx        uniform mesh spacing 
     @param[in] dims      grid (i,j) domain dimensions, ie domain size is dx*(i,j)
     @param[in] boxlo     lower left (i,j) corner of the data supplied by __this__ rank
     @param[in] boxhi     upper right (i,j) corner of the data supplied by __this__ rank

     The data writtem depends on the field parameter. Supported values include 
      - BISICLES_FIELD_SURFACE_ELEVATION
      - BISICLES_FIELD_BEDROCK_ELEVATION
      - BISICLES_FIELD_ICE_THICKNESS (ice thickness, arithmetic coarse average)
      - BISICLES_FIELD_ICE_THICKNESS_HARMONIC (ice thickness, harmonic coarse average, not conservative)
      - BISICLES_FIELD_CALVED_THICKNESS
      - BISICLES_FIELD_SURFACE_TEMPERATURE
      - BISICLES_FIELD_SURFACE_HEAT_FLUX
      
  */
  void bisicles_get_2d_data(int *intance_id, double *data_ptr, const int *field,
			    const double *dx, const int *dims, 
			    const int *boxlo, const int *boxhi);


  
  /** 
      \brief push thin ice onto the ice sheet where thickness < thin_ice_limit, record the quantity needed in data_ptr
      
      Exists to support thin ice -> snow conversion in UniCicles. Adds ice 
      according to the thickness on BISICLES grid:
      fills data_ptr with the change in thickness in that region, coarsened/refined as needed

  */
  
  void bisicles_push_thin_ice(int *intance_id, double *data_ptr, const double *thin_ice_limit,
				  const double *dx, const int *dims, 
				  const int *boxlo, const int *boxhi);
  /** 
      \brief pop ice from the top of the sheet where 
      ice thinner is than thin_ice_limit, record the quantity needed in data_ptr
      
      Exists to support thin ivce -> snow conversion in UniCicles. Removes the
      thin ice according to the thickness on BISICLES grid:
      fills data_ptr with the thickness in that region, coarsened/refined as needed

  */
  void bisicles_pop_thin_ice(int *intance_id, double *data_ptr, const double *thin_ice_limit,
				  const double *dx, const int *dims, 
				  const int *boxlo, const int *boxhi);
  
  /// initialize bisicles 
  void bisicles_init_instance(int *instance_key);
  ///step BISICLES forward in time by up to max_time or max_steps (whichever is sooner)
  void bisicles_advance(int *instance_key, double *start_time, double *max_time, int *max_step);
  
  ///write a checkpoint file 
  void bisicles_write_checkpoint(int *instance_key);
  //write a plot file
  void bisicles_write_plot(int *instance_key);
  ///read a checkpoint file 
  void bisicles_read_checkpoint(int *instance_key, const char *checkpoint_fname);
  
  void bisicles_set_header_int(int *instance_key, const char *attr_key,  const int *val);
  void bisicles_set_header_dble(int *instance_key, const char *attr_key, const double *val);
  void bisicles_set_header_char(int *instance_key, const char *attr_key, const char *val);
  void bisicles_get_header_int(int *instance_key, const char *attr_key, int *val);
  void bisicles_get_header_dble(int *instance_key, const char *attr_key, double *val);
  void bisicles_get_header_char(int *instance_key, const char *attr_key, char *val);

  // FORTRAN compatibility stuff : function foo gets a f_foo_ and f_foo
  void f_bisicles_new_instance_(int *instance_key,  char *input_fname, const int *len_fname, const int *mpi_comm);
  void f_bisicles_free_instance_(int *instance_key);
  void f_bisicles_set_2d_data_(int *instance_key,  double *data_ptr, const int *field, 
			  const double *dx, const int *dims, 
			     const int *boxlo, const int *boxhi);
  void f_bisicles_set_2d_geometry_(int *instance_key,  double *thck_data_ptr, double *topg_data_ptr,
				 const double *dx, const int *dims, 
				 const int *boxlo, const int *boxhi);
  void f_bisicles_get_2d_data_(int *intance_id, double *data_ptr, const int *field,
			    const double *dx, const int *dims, 
			     const int *boxlo, const int *boxhi);
  
  void f_bisicles_push_thin_ice_(int *intance_id, double *data_ptr, const double *thin_ice_limit,
				  const double *dx, const int *dims, 
				  const int *boxlo, const int *boxhi);
  
  void f_bisicles_pop_thin_ice_(int *intance_id, double *data_ptr, const double *thin_ice_limit,
				  const double *dx, const int *dims, 
				  const int *boxlo, const int *boxhi);
  
  void f_bisicles_init_instance_(int *instance_key);
  void f_bisicles_advance_(int *instance_key, double *start_time, double *max_time, int *max_step);
  void f_bisicles_write_checkpoint_(int *instance_key);
  void f_bisicles_write_plot_(int *instance_key);
  void f_bisicles_read_checkpoint_(int *instance_key, char *checkpoint_fname, const int *len_fname);
  void f_bisicles_set_header_int_(int *instance_key, char *attr_key, const int *attr_key_len, const int *val);
  void f_bisicles_set_header_dble_(int *instance_key, char *attr_key, const int *attr_key_len, const double *val);
  void f_bisicles_set_header_char_(int *instance_key, char *attr_key, const int *attr_key_len,  char *val, const int *len_val);
  void f_bisicles_get_header_int_(int *instance_key, char *attr_key, const int *attr_key_len,  int *val);
  void f_bisicles_get_header_dble_(int *instance_key, char *attr_key, const int *attr_key_len,  double *val);
  void f_bisicles_get_header_char_(int *instance_key, char *attr_key, const int *attr_key_len,  char *val, const int *len_val);


  void f_bisicles_new_instance(int *instance_key,  char *input_fname, const int *len_fname, const int *mpi_comm);
  void f_bisicles_free_instance(int *instance_key);
  void f_bisicles_set_2d_data(int *instance_key,  double *data_ptr, const int *field, 
			  const double *dx, const int *dims, 
			     const int *boxlo, const int *boxhi);
  void f_bisicles_set_2d_geometry(int *instance_key,  double *thck_data_ptr, double *topg_data_ptr,
				 const double *dx, const int *dims, 
				 const int *boxlo, const int *boxhi);
  void f_bisicles_get_2d_data(int *intance_id, double *data_ptr, const int *field,
			    const double *dx, const int *dims, 
			     const int *boxlo, const int *boxhi);

  void f_bisicles_push_thin_ice(int *intance_id, double *data_ptr, const double *thin_ice_limit,
				  const double *dx, const int *dims, 
				  const int *boxlo, const int *boxhi);
  
  void f_bisicles_pop_thin_ice(int *intance_id, double *data_ptr, const double *thin_ice_limit,
				  const double *dx, const int *dims, 
				  const int *boxlo, const int *boxhi);
  
  void f_bisicles_init_instance(int *instance_key);
  void f_bisicles_advance(int *instance_key, double *start_time, double *max_time, int *max_step);
  void f_bisicles_write_checkpoint(int *instance_key);
  void f_bisicles_write_plot(int *instance_key);
  void f_bisicles_read_checkpoint(int *instance_key, char *checkpoint_fname, const int *len_fname);
  void f_bisicles_set_header_int(int *instance_key, char *attr_key, const int *attr_key_len, const int *val);
  void f_bisicles_set_header_dble(int *instance_key, char *attr_key, const int *attr_key_len, const double *val);
  void f_bisicles_set_header_char(int *instance_key, char *attr_key, const int *attr_key_len,  char *val, const int *len_val);
  void f_bisicles_get_header_int(int *instance_key, char *attr_key, const int *attr_key_len,  int *val);
  void f_bisicles_get_header_dble(int *instance_key, char *attr_key, const int *attr_key_len,  double *val);
  void f_bisicles_get_header_char(int *instance_key, char *attr_key, const int *attr_key_len,  char *val, const int *len_val);
}



//set header data
template<typename T>
void bisicles_set_header(int *instance_key, const char *attr_key, const T *val);
//retrieve header data
template<typename T>
void bisicles_get_header(int *instance_key, const char *attr_key, T *val);
void bisicles_get_header(int *instance_key, const char *key, char *val);
#include "NamespaceFooter.H"
#endif
