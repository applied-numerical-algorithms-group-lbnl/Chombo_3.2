#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CHOMBOSUNDIALS_H
#define _CHOMBOSUNDIALS_H

#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>

#include "ChomboSundialsAdaptor.H"

/// Macros and definitions for SUNDIALS nvector interface

#ifdef __cplusplus
extern "C" {
#endif

//! Forward declare pointer-to type
typedef struct _N_VectorContent_Ch *N_VectorContent_Ch;

//! content in SUNDIALS nvector, pointer to Chombo adaptor
struct _N_VectorContent_Ch {
  ChomboSundialsAdaptor* adaptor;
  bool own_data;
};

//! public API
N_Vector N_VNew_Ch(SUNContext sunctx, ChomboSundialsAdaptor* adaptor);
void N_VDataExchange(N_Vector v);
bool N_VEquate(N_Vector x, N_Vector z);

// accessor macros
#define NV_CONTENT_CH(v) ( static_cast<N_VectorContent_Ch>(v->content) )
#define NV_ADAP_CH(v)    ( (*NV_CONTENT_CH(v)->adaptor) )
#define NV_OWN_DATA_CH(v)   ( NV_CONTENT_CH(v)->own_data )

#ifdef __cplusplus
} // extern "C"
#endif

#endif
