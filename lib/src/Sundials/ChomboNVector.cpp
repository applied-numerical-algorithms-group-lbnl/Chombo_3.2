#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "sundials/sundials_types.h"

#include "ChomboNVector.H"

///// Operations /////
// forward declaration
N_Vector N_VCreate_Ch(SUNContext sunctx);

// USER FUNCTION: Associate an NVector wrapper with an adaptor
// Only called by the user
N_Vector N_VNew_Ch(SUNContext sunctx, ChomboSundialsAdaptor* adaptor) {
  N_Vector v = N_VCreate_Ch(sunctx);
  NV_CONTENT_CH(v)->adaptor = adaptor;
  NV_OWN_DATA_CH(v) = false;
  return v;
}

// REQUIRED: Returns the global length of the NVector
/* Called by SUNDIALS to get the global vector length to convert between the
   linear solver tolerance (L2 norm) and integrator tolerance (WRMS norm).
   Alternatively, the factor for converting between norms may be changed by
   calling <pkg>SetLSNormFactor. Additionally, when using the SUNDIALS
   ManyVector, this function is used to compute the overall vector length from
   the subvectors. */
sunindextype N_VGetLength_Ch(N_Vector v) {
  return NV_ADAP_CH(v).getLength();
}

// REQUIRED: Create a new NVector wrapper using the input vector as a template
// Called by SUNDIALS to create internal workspace vectors
N_Vector N_VClone_Ch(N_Vector v) {
  N_Vector retv = N_VCreate_Ch(v->sunctx);
  NV_CONTENT_CH(retv)->adaptor = NV_ADAP_CH(v).newAdaptor();
  NV_OWN_DATA_CH(retv) = true;
  return retv;
}

// REQUIRED: Destroy an NVector wrapper and the data (if owned)
// Called by SUNDIALS to free workspace vectors
void N_VDestroy_Ch(N_Vector v) {
  if (NV_OWN_DATA_CH(v))
  {
    delete NV_CONTENT_CH(v)->adaptor;
  }
  delete NV_CONTENT_CH(v); v->content = NULL;
  N_VFreeEmpty(v);
  v = NULL;
}

N_Vector_ID N_VGetVectorID_Ch(N_Vector w) {
  return SUNDIALS_NVEC_CUSTOM;
}

// OPTIONAL: Return the real and integer workspace sizes
// Only used by SUNDIALS when the user asks for the workspace size
void N_VSpace_Ch(N_Vector nvSpec, sunindextype *lrw, sunindextype *liw) {
 // dummy
}

// REQUIRED: Vector linear sum: z_i = a * x_i + b * y_i
void N_VLinearSum_Ch(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y, N_Vector z) {
  NV_ADAP_CH(z).linearSum(NV_ADAP_CH(x), NV_ADAP_CH(y), a, b);
}

// REQUIRED: Set all vector entries to a constant: z_i = c
void N_VConst_Ch(sunrealtype c, N_Vector z) {
  NV_ADAP_CH(z).setConst(c);
}

// OPTIONAL: Component-wise vector product: z_i = x_i .* y_i
// Required when using component-wise constrain handling
void N_VProd_Ch(N_Vector x, N_Vector y, N_Vector z) {
  NV_ADAP_CH(z).prod(NV_ADAP_CH(x), NV_ADAP_CH(y));
}

// REQUIRED: Component-wise vector quotient: z_i = x_i ./ y_i
void N_VDiv_Ch(N_Vector x, N_Vector y, N_Vector z) {
  NV_ADAP_CH(z).div(NV_ADAP_CH(x), NV_ADAP_CH(y));
}

// REQUIRED: Scale a vector: z_i = c * x_i
void N_VScale_Ch(sunrealtype c, N_Vector x, N_Vector z) {
  NV_ADAP_CH(z).scale(NV_ADAP_CH(x), c);
}

// REQUIRED: Component-wise absolute value: z_i = |x_i|
void N_VAbs_Ch(N_Vector x, N_Vector z) {
  NV_ADAP_CH(z).abs(NV_ADAP_CH(x));
}

// REQUIRED: Component-wise inverse: z_i = 1 / x_i
void N_VInv_Ch(N_Vector x, N_Vector z) {
  NV_ADAP_CH(z).inv(NV_ADAP_CH(x));
}

// REQUIRED: Add constant to all vector components: z_i = x_i + b
void N_VAddConst_Ch(N_Vector x, sunrealtype b, N_Vector z) {
  NV_ADAP_CH(z).addConst(NV_ADAP_CH(x), b);
}

// REQUIRED: Euclidean dot product: d = sum_{i=0}^{n-1} x_i * y_i
// Optional when supplying a direct linear solver
sunrealtype N_VDotProd_Ch(N_Vector x, N_Vector y) {
  return NV_ADAP_CH(x).dotProd(NV_ADAP_CH(y)); // MPI sum across ranks
}

// REQUIRED: Max norm: m = max_{0 <= i <= n-1} |x_i|
sunrealtype N_VMaxNorm_Ch(N_Vector x) {
  return NV_ADAP_CH(x).maxNorm(); // MPI comm
}

// REQUIRED: Weighted root-mean-square norm:
//   m = sqrt( (sum_{i=0}^{n-1} (x_i w_i)^2) / n)
// Used to compute norm of time step error estimate
// Used to compute norms in SUNDIALS nonlinear solvers
sunrealtype N_VWrmsNorm_Ch(N_Vector x, N_Vector w) {
  return NV_ADAP_CH(x).wRMSNorm(NV_ADAP_CH(w)); // MPI reduce
}

// OPTIONAL: Masked weighted root-mean-square norm:
//   m = sqrt( (sum_{i=0}^{n-1} (x_i w_i H(id_i))^2) / n)
//   where H(id_i) = 1 for id_i > 0 and = 0 for id_i <= 0
// Required when using IDA(S) and suppressing algebraic variables in the time
// step local error test
sunrealtype N_VWrmsNormMask_Ch(N_Vector x, N_Vector w, N_Vector id) {
  return NV_ADAP_CH(x).wRMSNormMask(NV_ADAP_CH(w), NV_ADAP_CH(id)); // MPI reduce
}

// REQUIRED: Min entry value: m = min_{0 <= i <= n-1} x_i
sunrealtype N_VMin_Ch(N_Vector x) {
  return NV_ADAP_CH(x).min(); // MPI reduce
}

// OPTIONAL: Weighted Euclidean l2 norm: m = sqrt((sum_{i=0}^{n-1} (x_i w_i)^2))
// Required with KINSOL
sunrealtype N_VWL2Norm_Ch(N_Vector x, N_Vector w) {
  return NV_ADAP_CH(x).wL2Norm(NV_ADAP_CH(w)); // MPI reduce
}

// OPTIONAL: l1 norm: m = sum_{i=0}^{n-1} |x_i|
// Required with KINSOL when using the internal finite difference Jacobian
sunrealtype N_VL1Norm_Ch(N_Vector x) {
  return NV_ADAP_CH(x).l1Norm(); // MPI reduce
}

// OPTIONAL: Compare components to a value: z_i = 1 if |x_i| >= c, = 0 otherwise
// Required with CVODE(S) or IDA(S) when using component-wise constrain handling
// or when using the diagonal linear solver in CVODE(S)
void N_VCompare_Ch(sunrealtype c, N_Vector x, N_Vector z) {
  NV_ADAP_CH(z).compare(NV_ADAP_CH(x), c);
}

// OPTIONAL: Component-wise inverse: z_i = 1 / x_i with check of x_i = 0
// Required with CVODE(S) when using the diagonal linear solver
sunbooleantype N_VInvTest_Ch(N_Vector x, N_Vector z) {
  return NV_ADAP_CH(z).invTest(NV_ADAP_CH(x));
}

// OPTIONAL: Checks if the components of x violate the constrains c and returns
// the mask vector m_i = 1 if check fails and m_i = 0 if passed (the function
// returns SUNTRUE if any constraint was violated). The checks are as follows:
//   if c_i =  2 then x_i must be >  0
//   if c_i =  1 then x_i must be >= 0
//   if c_i =  0 then no check is performed
//   if c_i = -1 then x_i must be <= 0
//   if c_i = -2 then x_i must be <  0
// Required when using component-wise constrain handling
sunbooleantype N_VConstrMask_Ch(N_Vector c, N_Vector x, N_Vector m) {
  return NV_ADAP_CH(m).constrMask(NV_ADAP_CH(c), NV_ADAP_CH(x));
}

// OPTIONAL: Minimum component-wise quotient:
//   m = min_{0 <= i <= n-1} num_i / denom_i
// Required when using component-wise constrain handling
sunrealtype N_VMinQuotient_Ch(N_Vector num, N_Vector denom) {
  return NV_ADAP_CH(num).minQuotient(NV_ADAP_CH(denom));
}

// OPTIONAL: Print to stdout
void N_VPrint_Serial_Ch(N_Vector x) {
  NV_ADAP_CH(x).print();
}

// OPTIONAL: Print to C-style file 
void N_VPrintFile_Serial_Ch(N_Vector x, FILE* outfile) {
  NV_ADAP_CH(x).printFile(outfile);
}

// USER FUNCTION: Create a new NVector wrapper
N_Vector N_VCreate_Ch(SUNContext sunctx) {
  // create an empty vector (member data structure and vtable pointers are
  // initialized to NULL)
  N_Vector v = N_VNewEmpty(sunctx);
  if (v == NULL) return NULL;

  // Allocate the vector member data structure
  v->content = (void *)new _N_VectorContent_Ch;

  // Assign the communicator in the member data structure
  // NV_CONTENT_CH(v)->comm = Chombo_MPI::comm;

  // Set the function pointers in the vtable (skipped some optional methods)
  v->ops->nvgetlength     = N_VGetLength_Ch;
  v->ops->nvclone         = N_VClone_Ch;
  // v->ops->nvcloneempty    = N_VCloneEmpty_Ch;
  v->ops->nvdestroy       = N_VDestroy_Ch;
  v->ops->nvgetvectorid   = N_VGetVectorID_Ch;
  v->ops->nvspace         = N_VSpace_Ch;
  v->ops->nvlinearsum     = N_VLinearSum_Ch;
  v->ops->nvconst         = N_VConst_Ch;
  v->ops->nvprod          = N_VProd_Ch;
  v->ops->nvdiv           = N_VDiv_Ch;
  v->ops->nvscale         = N_VScale_Ch;
  v->ops->nvabs           = N_VAbs_Ch;
  v->ops->nvinv           = N_VInv_Ch;
  v->ops->nvaddconst      = N_VAddConst_Ch;
  v->ops->nvdotprod       = N_VDotProd_Ch;
  v->ops->nvmaxnorm       = N_VMaxNorm_Ch;
  v->ops->nvwrmsnorm      = N_VWrmsNorm_Ch;
  v->ops->nvwrmsnormmask  = N_VWrmsNormMask_Ch;
  v->ops->nvmin           = N_VMin_Ch;
  v->ops->nvwl2norm       = N_VWL2Norm_Ch;
  v->ops->nvl1norm        = N_VL1Norm_Ch;
  v->ops->nvcompare       = N_VCompare_Ch;
  v->ops->nvinvtest       = N_VInvTest_Ch;
  v->ops->nvconstrmask    = N_VConstrMask_Ch;
  v->ops->nvminquotient   = N_VMinQuotient_Ch;
  /* debugging functions */
  // FIXME - should these be parallel / hdf5?
  v->ops->nvprint     = N_VPrint_Serial_Ch; 
  v->ops->nvprintfile = N_VPrintFile_Serial_Ch; 

  return v;
}

// USER FUNCTION: Exchange data
void N_VDataExchange(N_Vector v) {
  // using Chombo's exchange feature
  NV_ADAP_CH(v).exchange();
}

// USER FUNCTION: Copy data
bool N_VEquate(N_Vector x, N_Vector z) {
  return NV_ADAP_CH(x).copyTo(NV_ADAP_CH(x));
}

