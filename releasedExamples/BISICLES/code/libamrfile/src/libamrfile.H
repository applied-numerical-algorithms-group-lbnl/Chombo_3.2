
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
 libamrfile.cpp
 C-interface functions to read Chombo's amr files and access their data
==========================================================================*/


#ifdef __cplusplus
extern "C"
  {
#endif

#define LIBAMRFILE_ERR_UNSPECIFIED   1
#define LIBAMRFILE_ERR_READ_FAILED   2
#define LIBAMRFILE_ERR_NULL_POINTER   3
#define LIBAMRFILE_ERR_NO_SUCH_AMR_ID    4
#define LIBAMRFILE_ERR_BAD_AMR_HIERARCHY 5
#define LIBAMRFILE_ERR_NO_SUCH_LEVEL     6
#define LIBAMRFILE_ERR_NO_SUCH_FAB       7
#define LIBAMRFILE_ERR_WRITE_FAILED  8
#define LIBAMRFILE_ERR_BAD_NGHOST 9
#define LIBAMRFILE_ERR_NO_SUCH_COMP 10
#define LIBAMRFILE_ERR_CREATE_FAILED 11
#define LIBAMRFILE_ERR_TOO_MANY_AMR_HIERARCHIES  12
#define LIBAMRFILE_ERR_BAD_BOX 13
#define LIBAMRFILE_ERR_BAD_REFINEMENT_RATIO 14
#define LIBAMRFILE_ERR_BAD_INTERPOLATION_ORDER 15
#define LIBAMRFILE_MAX_AMR_HIERARCHIES 3

    /// read in amr file into memory, provide an integer id (amr_id). Although
    /// amr_id *could* be a simple pointer, an int has been chosen for portability
    void amr_read_file(int *status, int *amr_id, const char *file);
    void amr_read_file_R(int *status, int *amr_id, char **file);

    /// create a new amr structure with a coarse level only and n_components
    void amr_create_coarse_2d(int *status, int *amr_id, 
			   const int *nx, const int *ny, 
			   const double* dx, const int* n_comp, const int* n_ghost);

    /// write an amr file with data already in memory, lookup up via an integer id (amr_id).
    void amr_write_file(int *status, int *amr_id, const char *file);
    void amr_write_file_R(int *status, int *amr_id, char **file);

    /// free amr data in memory, lookup up via an integer id
    void amr_free(int *status, int *amr_id);

    /// free all amr data in memory
    void amr_free_all();

    /// how many levels are there in the amr data of amr_id ?
    void amr_query_n_level(int *status, int *n_level, const int *amr_id);

    /// how many components are there in the amr data of amr_id ? 
    void amr_query_n_comp(int *status, int* n_comp, const int* amr_id);

    /// return domain corners for a given level grid
    void amr_query_domain_corners(int *status, int *lo, int* hi, const int *amr_id, const int *level);

    /// return current time
    void amr_query_time(int *status, double *time, const int *amr_id);

    /// retrieve the first buflen characters of the name of component comp?
    void amr_query_comp_name(int *status, char *name, const int* amr_id, const int* comp, const int* buflen);
    void amr_query_comp_name_R(int *status, char **name, const int* amr_id, const int* comp, const int* buflen);

    /// look up component number (int) given a name (char*)
    void amr_query_comp_id(int *status, int *comp, const int* amr_id, const char* name , const int* namelen);
    void amr_query_comp_id_R(int *status, int *comp, const int* amr_id, const char** name, const int* namelen);

    /// set the name of component comp
    void amr_set_comp_name(int *status, const char *name, const int* amr_id, const int* comp);

    void amr_set_comp_name_R(int *status, const char **name, const int* amr_id, const int* comp);

    /// how many fabs are there in data amr_id level level_id ?
    void amr_query_n_fab(int *status, int *n_fab, const int *amr_id, const int *level_id);

    /// find dimensions and number of components given (amr_id, level_id, fab_id) 
    void amr_query_fab_dimensions_2d(int *status, int *nx, int *ny,  int *ncomp, 
				     const int *amr_id, 
				     const int *level_id, 
				     const int* fab_id);
    /// read x,y,z and components data given (amr_id, level_id, fab_id, comp_id)?
    void amr_read_fab_data_2d(int *status, 
			      double *fab_data, 
			      double *x_data, 
			      double *y_data,
			      const int *amr_id, 
			      const int *level_id, 
			      const int* fab_id,
			      const int* comp_id,
			      const int* nghost);
    /// read x,y and components data over a level lev box bounded by lo and hi, averaging and interpolating as needed
    void amr_read_box_data_2d(int *status, 
			      double *comp_data, 
			      double *x_data, 
			      double *y_data,
			      const int *amr_id, 
			      const int *level_id, 
			      const int *lo,
			      const int *hi,
			      const int* comp_id,
			      const int* interp_order);

    /// write component data given (amr_id, level_id, fab_id, comp_id)?
    void amr_write_fab_data_2d(int *status, 
			       double *fab_data, 
			       int *nx, int *ny,
			       const int *amr_id, 
			       const int *level_id, 
			       const int* fab_id,
			       const int* comp_id,
			       const int* nghost);

#ifdef __cplusplus
  }
#endif

