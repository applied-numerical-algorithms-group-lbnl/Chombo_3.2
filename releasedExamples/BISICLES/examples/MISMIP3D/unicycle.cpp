#include <valarray>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>

using namespace std;

//#include "ArrayTools.hh"

#define RHO 900.0
#define RHOW 1000.0
#define GRAVITY 9.8
#define EPS 1.0e-9
#define NLAYER 11
bool gm_stop;



//#define EXPLICITMETH
//#define PERTMETH
//#define SIAMETH

#ifndef PERTMETH
#define NEGLECTMETH
#endif

enum model_type {SSA,L1L2}; // supported models
enum boundary_condition {NATURAL,PERIODIC,DIVIDE,MARINE}; // supported boundary conditions

#define MY_ASSERT(expression) assert(expression);
inline 
double harmonicMean(double a, double b){
  double d = a+b;
  return (abs(d) > EPS)?2.0*a*b/d:0.0;
} 


inline
double L2norm(const valarray<double>& x){

  return std::sqrt((x*x).sum());

}


inline
double maxnorm(const valarray<double>& x){

  return std::max(x.max(),-x.min());

}

void thicken(valarray<double>& cH,
	      valarray<double>& cR,
	      int nadv)
{

  int n = 0;
  double f = (1.0-RHO/RHOW);
  int nx = cH.size();
  for (int i = 0; i < nx; ++i)
    {
      double sg= cR[i]+cH[i];
      double sf = cH[i]*f;
      if (n < nadv && sf > sg + EPS)
	{
	  n++;
	  cH[i] = -cR[i]/(1-f) + 1.0;
	}
      
    }

}


//solve 3 DoF problem Ax = b
void solve3(valarray<double>& A,
	    valarray<double>& b,
	    valarray<double>& x)
{

  assert(A.size()>8);
  assert(b.size()>2);
  assert(x.size()>2);

  b[1] = A[0]*b[1] - A[3]*b[0];
  b[2] = A[0]*b[2] - A[6]*b[0];
  for (int i = 2; i >= 0; --i){
    A[3+i] = A[0] * A[3+i] - A[3]*A[0+i];
    A[6+i] = A[0] * A[6+i] - A[6]*A[0+i];
  }
 
  b[2] = A[4]*b[2] - A[7]*b[1];
  for (int i = 2; i > 0; --i){
    A[6+i] = A[4] * A[6+i] - A[7]*A[3+i];
  }
 
  x[2] = b[2]/A[8];
  x[1] = (b[1] - A[5]*x[2])/A[4];
  x[0] = (b[0] - A[1]*x[1] - A[2]*x[2])/A[0];

}


inline void simpleSmooth(valarray<double>& x)
{
  int n = x.size();
  valarray<double> y(n+2);
  y[slice(1,n,1)] = x[slice(0,n,1)];
  y[0] = x[0];
  y[n]=x[n-1];

  for (int i = 0; i < n; ++i)
    {
      x[i] = 0.25 * (2.0 * y[i] + y[i+1] + y[i-1]); 
    }

}


inline
void relaxGSRB(double r, bool red,
	       const valarray<double>& ad,
	       const valarray<double>& au,
	       const valarray<double>& al,
	       const valarray<double>& b, 
	       valarray<double>& x)
{

  int n = x.size();

  for (int i = (red?1:2) ; i < n-1; i+=2)
    {
      x[i] = r * x[i] 
	+ (1.0-r)* (b[i] - al[i]*x[i-1] - au[i]*x[i+1])/ad[i];
    }

}
void smoothGSRB
(int n,
 const valarray<double>& ad,
 const valarray<double>& au,
 const valarray<double>& al,
 const valarray<double>& b, 
 valarray<double>& x)
{
  //GSRB smoother
  double r = 1.0;
  for (int i = 0; i < n; ++i)
    {
      relaxGSRB(r,true,ad,au,al,b,x);
      relaxGSRB(r,false,ad,au,al,b,x);
    }
}

extern "C" {
  void dgbsv_(int *n, int *kl, int *ku, int *nrhs,
	      double *ab, int *ldab, int *ipiv,
	      double *b, int *ldb, int* info);

}



void prolongCell
(valarray<double>& xf, 
 const valarray<double>& xc)
{
  // cell-cente prolongation : linear interpolation
  int nf = xf.size();
  int nc = xc.size();
  assert( nc%2 == 0);// even arrays only
  assert(2*nc == nf);
  for (int j = 1; j < nc; ++j)
    {
      xf[2*j-1] = .25 * (3.0 * xc[j-1] + xc[j]);
      xf[2*j] = .25 * (xc[j-1] + 3.0 * xc[j]);
    }
  xf[0] = .25 * (5.0*xc[0]-xc[1]);
  xf[nf-1] = .25 * (5.0*xc[nc-1]-xc[nc-2]);

}

void restrictCell
(valarray<double>& xc, 
 const valarray<double>& xf)
{

  int nf = xf.size();
  int nc = xc.size();
  assert( nc%2 == 0);// even arrays only
  assert(2*nc == nf);
  for (int j = 1; j < nc + 1; ++j)
    {
      xc[j-1] = .5 * (xf[2*j-1] + xf[2*j]);
    }

}




// virtual base for classes which solve the
// linear (discrete) Poisson equation
// d/dx (mu du/dx) - C u = r, on a uniform
// mesh with spacing h. 
class PoissonSolver
{
public:
  virtual ~PoissonSolver(){};
  virtual void setCoeffs
  (double h, const valarray<double>& mu,
   const valarray<double>& C,const valarray<double>& r) = 0;
  
  virtual void computeResidual
  (const valarray<double>& u,valarray<double>& r) = 0;
  
  virtual double residualNorm
  (const valarray<double>& u) = 0;

  virtual void solve
  (valarray<double>& u) = 0;

  virtual void imiCorrect
  ( int i, double h, double d, double mu, double mux,
    double Cw, double Ce, double jf, double jC ) = 0;

};

//Poisson solver using the TDMA method
class TDMAPoissonSolver : public PoissonSolver
{

  int n;
  valarray<double> ad, au, al, b, res, p, q;
  bool periodic;
  TDMAPoissonSolver();

public:
  TDMAPoissonSolver
  (int n, bool periodic) 
  : n(n), ad(n), au(n), al(n), b(n), 
    res(n), p(n), q(n), periodic(periodic)
  {;}

  virtual void setCoeffs
  (double h,const valarray<double>& mu,
   const valarray<double>& C,const valarray<double>& r);
  
  virtual void computeResidual
  (const valarray<double>& u,valarray<double>& r);
  
  virtual double residualNorm
  (const valarray<double>& u)
  {
    computeResidual(u,res);
    //return L2norm(res);
    return maxnorm(res);
  }

  virtual void solve
  (valarray<double>& u);
  
  virtual void imiCorrect
  ( int i, double h, double d, double mu, double mux,
    double Cw, double Ce, double jf, double jC );


};


#define GL_RHS_NONE     0
#define GL_RHS_BISICLES 1
#define GL_MU_NONE     0
#define GL_MU_DU_SEAWARD 1
#define GL_MU_DU_LANDWARD 2
#define GL_MU_DU_SEALIMITED 3
#define GL_MU_DU_LANDLIMITED 4
#define GL_MU_DU_MEAN 5
#define GL_MU_DU_HARM 6
#define GL_MU_DU_LAST 7
#define GL_MU_MU_SEAWARD 11
#define GL_MU_MU_LANDWARD 12
#define GL_MU_MU_SEALIMITED 13
#define GL_MU_MU_LANDLIMITED 14
#define GL_MU_MU_MEAN 15
#define GL_MU_MU_HARM 16
#define GL_MU_MU_LAST 17

class GLCorrection
{
  int _M_rhs, _M_mu;

 public:
  GLCorrection(int rhs, int mu)
    : _M_rhs(rhs), _M_mu(mu)
    {;}

  void correctRhs(const valarray<double>& cS,
		  const valarray<double>& cH,
		  const vector<bool> floating,
		  double coeff,
		  valarray<double>& rhs) const
  {

    if (_M_rhs == GL_RHS_BISICLES){
      int nx = cS.size();
      
      for (int i = 1; i < nx-1; ++i){
	
	if (floating[i]){
	  double a = 0.5;
	   if (!floating[i-1] && floating[i+1] ){
	    
	     rhs[i] = -(cS[i+1] - cS[i]) *  
	       ((1.0-a) * cH[i+1] + (a)* cH[i]) * coeff;
	     
	     int j = i - 1;
	     
	     double testrhs = -(cS[j] - cS[j-1]) *
	       ((1.0-a) * cH[j-1] + (a)* cH[j]) * coeff;

	     //if (std::abs(rhs[j]) > std::abs( testrhs))
	       rhs[j] =  testrhs;
	     
	     
	   }
	   if (!floating[i+1] && floating[i-1]){
	      
	     rhs[i] = -(cS[i] - cS[i-1]) *  
	       ((1.0-a) * cH[i-1] + (a)* cH[i]) * coeff;
	     
	     int j = i + 1;
	     
	     
	     double testrhs = -(cS[j+1] - cS[j]) *
	       ((1.0-a) * cH[j+1] + (a)* cH[j]) * coeff;
	     //if (std::abs(rhs[j]) > std::abs( testrhs))
	       rhs[j] = testrhs;
	     
 
	   }
	}

      }
    }

  }

  double glDu(double landGrad, double glGrad, double seaGrad) const
  {
    switch(_M_mu){
      
    case GL_MU_DU_SEAWARD:
      return seaGrad;
      break;
      
    case GL_MU_DU_LANDWARD:
      return landGrad;
      break;  
      
    case GL_MU_DU_SEALIMITED:
      return (abs(seaGrad)<abs(glGrad))?seaGrad:glGrad;
      break;  
    
    case GL_MU_DU_LANDLIMITED:
      return (abs(landGrad)<abs(glGrad))?landGrad:glGrad;
      break;  
      
    case GL_MU_DU_MEAN:
      return 0.5 * (landGrad + seaGrad);
      break;  
      
    case GL_MU_DU_HARM:
      return (landGrad*seaGrad < EPS)?EPS
	:(2.0 * landGrad*seaGrad / (landGrad + seaGrad));
      break;  
     
    default:
      return glGrad;
      break;
    }
    
  }


  void correctGradU
  (valarray<double>& fGradU, 
   const vector<bool>& floating) const
  {
    if (_M_mu > GL_MU_NONE && _M_mu < GL_MU_DU_LAST)
      {
	int nx = fGradU.size() - 1;
	for (int i = 1; i < nx-1; ++i){
	  if (floating[i]){
	    if (!floating[i-1] && floating[i+1]) 
	      {
		//gl between i and i - 1,
		fGradU[i] = glDu(fGradU[i-1] , fGradU[i], fGradU[i + 1]);
	      }
	    if (!floating[i+1] && floating[i-1])
	      { 
		//gl between i and i + 1, 
		fGradU[i+1] = glDu(fGradU[i+1] , fGradU[i], fGradU[i-1]);
	      }
	  }
	}
      }
  }


  double glMu(double landMu, double glMu, double seaMu) const
  {
    switch(_M_mu){
      
    case GL_MU_MU_SEAWARD:
      return seaMu;
      break;
      
    case GL_MU_MU_LANDWARD:
      return landMu;
      //return std::max(landMu,seaMu);
      break;  
      
    case GL_MU_MU_SEALIMITED:
      return (abs(seaMu)<abs(glMu))?seaMu:glMu;
      break;  
    
    case GL_MU_MU_LANDLIMITED:
      return (abs(landMu)<abs(glMu))?landMu:glMu;
      break;  
      
    case GL_MU_MU_MEAN:
      return 0.5 * (landMu + seaMu);
      break;  
      
    case GL_MU_MU_HARM:
      return (landMu*seaMu < EPS)?EPS
	:(2.0 * landMu*seaMu / (landMu + seaMu));
      break;  
     
    default:
      return glMu;
      break;
    }
    
  }



  void correctMuH(valarray<double>& fMuH, const vector<bool>& floating) const
  {
     if (_M_mu >= GL_MU_MU_SEAWARD && _M_mu < GL_MU_MU_LAST)
       {
	 int nx = fMuH.size() - 1;
	 for (int i = 1; i < nx-1; ++i){
	   if (floating[i]){
	    if (!floating[i-1]  && floating[i+1]) 
	      {
		//gl between i and i - 1,
		fMuH[i] = glMu(fMuH[i-1] , fMuH[i], fMuH[i + 1]);
	      }
	    if (!floating[i+1] && floating[i-1] )
	      { 
	    	//gl between i and i + 1, 
	    	fMuH[i+1] = glMu(fMuH[i+1] , fMuH[i], fMuH[i-1]);
	      }
	  }
	}
      }

  }
      
  
};





//virtual base which implements specialist parts
//of the velocity model: calculation of the 
//effective viscosity , surface velocity,
//and vertically averaged velocity. 
class VelocityModel 
{
public:
  virtual ~VelocityModel(){};
  
  virtual void computeFaceMu
  (const valarray<double>& fGradUb,
   valarray<double>& fMu) const = 0;

  //compute face velocity fUa and face diffision fD
  virtual void computeFaceVelocity
  (const valarray<double>& cU,
   valarray<double>& fUa,
   valarray<double>& fD
   ) const = 0;



  
};

inline
void faceAverage(const valarray<double>& cU,
		 valarray<double>& fU)
{
  
  int nx = cU.size();
  assert(fU.size() == nx + 1);
  
  for (int i = 1 ; i < nx ; ++i){
      fU[i] =  0.5 * (cU[i-1]  + cU[i]);
    }

  fU[0] = cU[0];
  fU[nx] = cU[nx-1];

}

inline
void cellAverage(const valarray<double>& fU,
		 valarray<double>& cU)
{
  
  int nx = cU.size();
  assert(fU.size() == nx + 1);
  
  for (int i = 0 ; i < nx ; ++i){
    cU[i] = 0.5 * (fU[i]+fU[i+1]);
    }
}


inline 
void faceGradient(const valarray<double>& cU,
		  double dx,
		  valarray<double>& gU,
		  bool periodic= false)
{

  int nx = cU.size();
  assert(gU.size() == nx + 1);
  
  for (int i = 1 ; i < nx ; ++i)
      gU[i] =  (cU[i]  - cU[i-1])/dx;

  if (periodic)
    {
      gU[0] = cU[0] - cU[nx-2];
      gU[nx] = cU[nx] - cU[1];
    }
  else
    {
      //hmmm... should bring the BC in here.q
      //OTOH, [0] and [nx] are outside the domain
      gU[0] = gU[1];
      gU[nx] = gU[nx-1];
    }

  

}


class SSAVelocityModel : public VelocityModel
{
  double dx, n;
  const valarray<double> fH;
  const valarray<double> fA;
  SSAVelocityModel();

public:
  virtual ~SSAVelocityModel(){};
  
  SSAVelocityModel
  (const valarray<double> fH,
   const valarray<double> fA,
   double dx, double n)
    : fH(fH), fA(fA), n(n), dx(dx)
  {;}

  virtual void computeFaceVelocity
  (const valarray<double>& cU,
   valarray<double>& fUa,
   valarray<double>& fD) const
  {
    int nx = cU.size();
    assert(fUa.size() == nx + 1);
    assert(fD.size() == nx + 1);

    faceAverage(cU, fUa);
  
    fD = 0.0;
  }
 

  virtual void computeFaceMu
  (const valarray<double>& fGradU,
   valarray<double>& fMu) const;
  
  
};

void SSAVelocityModel::computeFaceMu
(const valarray<double>& fGradU,
 valarray<double>& fMu) const
{

  assert(fMu.size() == fA.size());
  assert(fMu.size() == fGradU.size());
  
  if ( std::abs(n - 1.0) < 1.0e-6)
    {
      //linear rheology
      fMu = 0.5 / fA;
    }
  else
    {
      double q = 1.0/n;
      double p = (1.0-n) / (2.0*n);
      int nx = fGradU.size() - 1;


      for (unsigned int i = 1; i < nx ; i++)
	{
	  double gu = abs(fGradU[i]);
	  double es = gu*gu;
	  double es0 = 1.0e-16;
	  fMu[i] = 0.5 * pow(fA[i],-q) * pow (es + es0 , p);
	}
    }
  fMu[0] = fMu[fMu.size()-1] = 0.0;
  fMu*=fH;
  fMu*=4.0;
}


class L1L2VelocityModel : public VelocityModel
{
  double dx, n;
  const valarray<double> fH;
  const valarray<double> cH;
  const valarray<double> fA;
  const valarray<double> cA;
  mutable valarray<double> fGradS;
  valarray<double> cGradS;
  const vector<bool> floating;
  mutable  valarray<double> phiTildeSqr;
  

  L1L2VelocityModel();

  //residual of the L1L2 constutive equation
  double res
  (double mu, double A, double es,
   double ts, double p) const
  {
    return 1.0 - 2.0 * A * mu * pow(4.0 * pow(mu,2.0) * es + ts, p);
  } 

  //solve the L1L2 constutive equation via the secant method
  void secantSolve(double& mu, const double mua, const double A,
		   const double es, const double ts, const double p,
		   const int niter, const double tol, const bool dbg) const
  {

    double fa, fb, t, tmua;
    tmua = mua;
    fa = res(tmua, A, es, ts, p);
    for (int i = 0; i < niter; ++i)
      {
	fb = res(mu, A, es, ts, p);
	if (abs(fb) < tol) break;
	t = mu;
	mu = tmua - (tmua - mu)/(fa - fb)*fa;
	tmua = t;
	fa = fb;
      }
  }

  void facePerturb(const double amplitude) const
  {
    
    valarray<double> cHp(cH);
    int nx = cHp.size();
    for (int i = 0; i < nx; i+=2)
      {
	cHp[i] += amplitude;
      }
    for (int i = 0; i < nx; i+=2)
      {
	cHp[i] -= amplitude;
      }
    
    valarray<double> cSp(0.0,nx);
    for (int i = 0; i < nx; i++)
      {
	 if (floating[i])
	   {
	     cSp[i] = (1-RHO/RHOW)*(cHp[i] - cH[i]);
	   }
	 else
	   {
	     cSp[i] = cHp[i] - cH[i];
	   }
      }
    for (int i = 1; i < nx; i++)
      {
	fGradS += (cSp[i] - cSp[i-1])/dx;
      }

    phiTildeSqr = RHO*GRAVITY*fH*fGradS;
    phiTildeSqr *= phiTildeSqr;

  }

 public:
  virtual ~L1L2VelocityModel(){};
  
  L1L2VelocityModel
  (const valarray<double>& fH,
   const valarray<double>& cH,
   const valarray<double>& fA,
   const valarray<double>& cA,
   const valarray<double>& afGradS,
   const valarray<double>& acGradS,
   const vector<bool>& floating,
   double dx, double n)
    : fH(fH), cH(cH),  fA(fA), cA(cA), 
      fGradS(afGradS), cGradS(acGradS),
      floating(floating),
      dx(dx), n(n),
      phiTildeSqr(fGradS.size())
  {


#ifdef SMOOTHA
    int nx = cH.size();

    valarray<double> mu(1.0e+6,nx+1);
    mu[0] = 0.0;
    mu[nx] = 0.0;
    valarray<double> C(1.0,nx+1);
    TDMAPoissonSolver ps(nx+1,false);
    ps.setCoeffs(dx,mu,C,fGradS);
    ps.solve(fGradS);
    TDMAPoissonSolver psc(nx,false);
    mu[nx-1] = 0.0;
    psc.setCoeffs(dx,mu,C,fGradS);
    psc.solve(cGradS);
#endif
#ifdef LIMIT
    double limit = 1.0e-5;
    int nx = cH.size();
    for (int i = 0; i < nx; ++i)
      {
	if (fGradS[i] > 0)
	  fGradS[i] = max(fGradS[i],limit);
	else
	  fGradS[i] = min(fGradS[i],-limit);

 	if (cGradS[i] > 0)
	  cGradS[i] = max(cGradS[i],limit);
	else
	  cGradS[i] = min(cGradS[i],-limit);
      }


#endif
    // {
    //   int nx = cH.size();
      
    //   valarray<double> mu(1.0e+3,nx+2);
    //   mu[0] = 0.0;
    //   mu[nx] = 0.0;
    //   TDMAPoissonSolver ps(nx+1,false);
    //   valarray<double> C(1.0,nx+1);
    //   ps.setCoeffs(dx,mu,C,fGradS);
    //   ps.solve(fGradS);
    // }
    phiTildeSqr = RHO*GRAVITY*fH*fGradS;
    phiTildeSqr *= phiTildeSqr;
    
  }
  
  virtual void computeFaceMu
  (double sigma,
   const valarray<double>& epsSpr,
   valarray<double>& fMu, 
   const valarray<double>& fMua) const;
  
  virtual void computeFaceMu
  (const valarray<double>& sigma,
   const valarray<double>& fGradU,
   vector<valarray<double> >& fMu) const;
  
  virtual void computeFaceMu
  (const valarray<double>& fGradU,
   valarray<double>& fMu) const;

  //velocity in each layer
  virtual void  computeFaceVelocity
  (const valarray<double>& cUb,
   const valarray<double>& fGradUb,
   const valarray<double>& sigma, 
   vector<valarray<double> >& fU) const;

    //velocity in each layer
  virtual void  computeFaceVelocity
  (const valarray<double>& cUb,
   const valarray<double>& sigma, 
   vector<valarray<double> >& fU) 
  {
    valarray<double> fGradUb(0.0, cUb.size() + 1);
    faceGradient(cUb,dx,fGradUb);
    computeFaceVelocity(cUb, fGradUb, sigma, fU);
  }


  //vertically averaged/ surface velocities and diffusion
  virtual void computeFaceVelocity
  (const valarray<double>& cU,
   valarray<double>& fUa,
   valarray<double>& fD) const;


};

inline
void raise(valarray<double> v, double power)
{

  for (int i = 0; i < v.size(); ++i)
    v[i] = pow(v[i],power);

}

// void L1L2VelocityModel::computeFaceDiffusion
//   (valarray<double>& fD) const
// {

//   fD = 2.0/(n+2) * std::pow(RHO*GRAVITY,n) * fA
//     * std::pow(fGradS,n-1) *  std::pow(fH,n+2);
//   //fD *= (1.0 - RHO/RHOW);
//   fD*=0.0;
// }


void L1L2VelocityModel::computeFaceVelocity
(const valarray<double>& cUb,
 const valarray<double>& sigma,
 const valarray<double>& fGradUb,
 vector<valarray<double> >& fU) const
{

  vector<valarray<double> > fMu(sigma.size(), fGradUb);
  computeFaceMu(sigma, fGradUb, fMu);
    
  //sigma-derivative of non-gravitational phizx at the i cell face
  vector<valarray<double> > fdphizx
    (sigma.size(), valarray<double> (fGradUb.size()));

  valarray<double> cGradUb(cUb.size());
  cellAverage(fGradUb,cGradUb);

 
  
  for (int l = 0; l < sigma.size(); ++l)
    {
      valarray<double> cMu(cUb.size());
      cellAverage(fMu[l],cMu);
      for (int i = 2; i < cUb.size() -1; ++i)
	{
	  fdphizx[l][i] = 4.0 / dx * 
	    (cH[i]*cMu[i]*cGradUb[i] 
	     - cH[i-1]*cMu[i-1]*cGradUb[i-1]);

	  //assert(abs(fdphizx[l][i]) < 1.0e+2);

	}
    }

  
  vector<valarray<double> > fphizx
     (sigma.size(), valarray<double> (0.0,fGradUb.size()));

  //non-gravitational phizx at the i cell face
  for (int l = 1; l < sigma.size(); ++l)
    {
      double s = 0.5 * (sigma[l] - sigma[l-1]);
      fphizx[l] =  fphizx[l-1] + s * (fdphizx[l] + fdphizx[l-1]);
    }


  // add gravational stress to phizx
  for (int l = 0; l < sigma.size(); ++l)
    {
      fphizx[l] -= RHO * GRAVITY * sigma[l] * fGradS * fH;
    }

  //compute u(sigma) - ub = int( 
  vector<valarray<double> > fdU
    (sigma.size(), valarray<double> (fGradUb.size()));
  
  double p = 0.5 * (n - 1.0);

  for (int l = sigma.size() - 1; l >= 0; --l)
    {
      fdU[l] = (4.0 * fMu[l] * fMu[l] * fGradUb * fGradUb 
      		+ fphizx[l]*fphizx[l]);
      fdU[l] = fphizx[l]*fphizx[l];
      raise(fdU[l],p);
      fdU[l] *= 2.0 * fphizx[l] * fH * fA;
      
    }

  fU[sigma.size() - 1] = 0.0;
  for (int l = sigma.size() - 2; l >=0 ; --l)
    {
      double s = 0.5 * (sigma[l+1]-sigma[l]); 
      fU[l] = s * (fdU[l] + fdU[l+1]) + fU[l+1];
#ifdef LIMITFVERT
      limit(fU[l],LIMITFVERTVAL)
#endif
    }

  
}

void L1L2VelocityModel::computeFaceVelocity
(const valarray<double>& cUb,
 valarray<double>& fUa,
 valarray<double>& fD) const
{
  int nLayer = NLAYER; // arbitrary for now
  valarray<double> sigma(nLayer);
  double ds = 1.0 / double(nLayer);
  for (int l = 1; l < sigma.size(); ++l)
    sigma[l] = ds + sigma[l-1];

  
 
  valarray<double> fGradUb(0.0,fUa.size());
  for (int i = 1; i < cUb.size() ; ++i)
    {
      fGradUb[i] = (cUb[i] - cUb[i-1])/dx;
    }


#ifdef EXPLICITMETH

  vector<valarray<double> > fU (nLayer, valarray<double>(0.0,fUa.size()));
  computeFaceVelocity(cUb, sigma,  fGradUb, fU);
  //midpoint rule inegration over sigma to find vertically averaged velocity
  faceAverage(cUb, fUa);
  for (int l = 0; l < nLayer; l++){
    fUa += ds*(fU[l]);
  }
  fD = 0.0;
#endif

#ifdef PERTMETH
  faceAverage(cUb, fUa);
  //velocity with perturbed thickness field
  double P = 1.0e-1; //amplitude of perturbation
  //perturbation with phase 0
  facePerturb(P);
  vector<valarray<double> > fUp (nLayer, valarray<double>(0.0,fUa.size()));
  computeFaceVelocity(cUb, sigma,  fGradUb, fUp);
  //perturbation with phase pi
  facePerturb(-2.0 * P);
  vector<valarray<double> > fUm (nLayer, valarray<double>(0.0,fUa.size()));
  computeFaceVelocity(cUb, sigma,  fGradUb, fUm);
  //back to unperturbed state
  facePerturb(P);
  int nx = cUb.size();
  //midpoint rule inegration over sigma to find vertically averaged velocity
  valarray<double> up(0.0,nx+1);
  valarray<double> um(0.0,nx+1);
  for (int l = 0; l < nLayer; l++){
     up += ds*(fUp[l]);
     um += ds*(fUm[l]);
  }
  for (int i = 1; i < nx; i++)
    {

      double dP = (i%2 == 0)?(2.0*P):(-2.0*P);
      fD[i] = fH[i]*std::max(0.0 , 0.5*(um[i]-up[i]) * dx / dP) ;

    }
  
  
#define FUDGE
#ifdef FUDGE
  double fudgeu = 200.0;
  for (int i = 1; i < nx; i++)
    {
      double fudge = std::exp( - std::abs(fUa[i]/fudgeu));
      up[i]*=fudge;  um[i]*=fudge;  
      fD[i] *= fudge;
    }
#endif

  //fUa += .5 * (up + um)
    //int k = 1;
    for (int i = 1; i < nx; i++)
    {
    fUa[i] += (.5 * (up[i] + um[i]) + fD[i]/fH[i]*(cH[i] - cH[i-1])/dx);
    }
    //k = 2;
#endif

#ifdef SIAMETH
  int nx = cH.size();
  vector<valarray<double> > fU (nLayer, valarray<double>(0.0,fUa.size()));
  computeFaceVelocity(cUb, sigma,  fGradUb, fU);
  valarray<double> u(0.0,nx+1);

  double rn = 2.0 * std::pow(RHO * GRAVITY,n) / (n+2.0);
  fD = rn * std::pow(fH,n+2.0) * fA * std::pow ( std::abs(fGradS), n-1.0);
  for (int i = 0; i < nx; i++)
    {
      if (floating[i])
	fD[i] *= (1-RHO/RHOW);
    }
  for (int l = 0; l < nLayer; l++){
    u += ds*(fU[l]);
  }

  // u*=0.25; fD *= 0.25;
  faceAverage(cUb, fUa);
  
 
  for (int i = 1; i < nx; i++)
    {
      fUa[i] += (u[i] + fD[i]/fH[i]*(cH[i] - cH[i-1])/dx);
    }
  
#endif

#ifdef NEGLECTMETH
 faceAverage(cUb, fUa);
 fD = 0.0;
#endif

  // for (int i = 1; i < nx; i++)
  //   {
  //     if (!floating[i] && !floating[i-1])
  // 	fD[i] += fH[i];
  //   }
}


void L1L2VelocityModel::computeFaceMu
(double sigma, 
 const valarray<double>& epsSqr,
 valarray<double>& fMu, 
 const valarray<double>& fMua) const
{

  double p1 = 0.5 * (n - 1.0);
  double p2 = - p1 / n;
  double p3 = - 1.0 / n;
  const double es0 = 1.0e-16;

  for (int i = 0; i < fMu.size(); ++i)
    {
      
      const double& es = epsSqr[i];//pow(fGradU[i],2.0) + es0;
      const double ts = phiTildeSqr[i] * sigma * sigma;
      //initialise with SSA Glen's law 
      //double mua = 0.5 * pow(fA[i],p3)* pow(es,p2);
      
      if (ts > es0)
	{
	  const int niter = 15;
	  const double tol = 1.0e-9;
	  double mua =  fMua[i];
	  if (std::abs(mua-fMu[i]) < tol)
	    mua = .95 * fMu[i];
	  secantSolve(fMu[i], mua, fA[i], epsSqr[i], ts, p1, niter, tol, i < 2);

	} 
      else 
	{
	  //mua will do in this case
	  fMu[i] =  0.5 * pow(fA[i],p3)* pow(es,p2);
	}
    }

}


void L1L2VelocityModel::computeFaceMu
(const valarray<double>& sigma,
 const valarray<double>& fGradU,
 vector<valarray<double> >& fMu) const
{
  const double es0 = 1.0e-16;
  valarray<double> epsSqr(fGradU);
  epsSqr*=epsSqr; epsSqr+=es0;
  
  const double p1 = 0.5 * (n - 1.0);
  const double p2 = - p1 / n;
  const double p3 = - 1.0 / n;

 
  //top layer : use glen's law mu_G and 0.95* mu_G to initialize
  fMu[1] = 0.5 * pow(fA,p3)* pow(epsSqr,p2);
  fMu[0] = 0.95 * fMu[1];
  computeFaceMu(sigma[0], epsSqr, fMu[0], fMu[1]);

  //second layer : use mu[0] and 0.95*mu[0]
  fMu[1] = 0.95 * fMu[0];
  computeFaceMu(sigma[1], epsSqr, fMu[1], fMu[0]);

  //remaining layers: use mu[l-1] and mu[l-2] to initialize
  for (int l = 2; l < sigma.size(); ++l)
    {
      fMu[l] = fMu[l-1];
      computeFaceMu(sigma[l], epsSqr, fMu[l], fMu[l-2]);
    }
}

void L1L2VelocityModel::computeFaceMu
(const valarray<double>& fGradU,
 valarray<double>& fMu) const
{

  int nLayer = NLAYER; // arbitrary for now
  valarray<double> sigma(nLayer);
  double ds = 1.0 / double(nLayer);
  sigma[0] = ds/2.0;
  for (int l = 1; l < sigma.size(); ++l)
    sigma[l] = ds + sigma[l-1]; 

  vector<valarray<double> > layerFMu(nLayer, fMu);
  computeFaceMu(sigma,fGradU,layerFMu);
  
  //midpoint rule inegration over sigma to find integral(mu,dz)/H
  fMu = 0.0;
  for (int l = 0; l < nLayer; l++){
    fMu += ds*layerFMu[l];
  }

  fMu*=fH;
  fMu*=4.0;

}

// virtual base for classes which solve the
// non linear (discrete) Poisson equation
// d/dx (mu du/dx) - (C |u|^m-1 + C0) u = r, on a uniform
// mesh with spacing h.
// mu is provided by a subclass of velocityModel
class VelocitySolver
{
  public:
  virtual ~VelocitySolver(){};

  virtual void solve
  (const VelocityModel& vModel,
   const valarray<double>& fA,
   const valarray<double>& cH,
   const valarray<double>& cRhs,
   const valarray<double>& cC0,
   const valarray<double>& cC,
   const valarray<double>& fH,
   const valarray<double>& cS,
   const valarray<double>& cR,
   const vector<bool>& floating, 
   const double gn, const double gm,
   const double h, const GLCorrection& gc,
   valarray<double>& u, 
   valarray<double>& res)=0;
  
};






void  TDMAPoissonSolver::imiCorrect
(int i, double h, double d, double mu, double mux, 
 double Cw, double Ce, double jf, double jC )
{
  
  if (i > 0 && i < n - 1){

    double dp1 = d + h;
    double dm1 = d - h;
    double dp2 = dp1 + h; 
    //correction on the left
    { 
      valarray<double> J(9);
      J[0] = 1.0; J[1] = 1.0; J[2] = 1.0 + (dp1*dp1)/(2.0*mu) * jC;
      J[3] =  dm1 ; J[4] = d ; J[5] = dp1;
      J[6] = J[3]*J[3]/2.0; J[7] = J[4]*J[4]/2.0 ; J[8] = J[5]*J[5]/2.0;
      valarray<double> r(3);
      r[0] = 0.0 ;
      r[1] = mux;
      r[2] = mu;
      valarray<double> x(3);
      solve3(J,r,x);
      x*=-h;
      
      al[i] = x[0];
      ad[i] = x[1] + Cw*h;
      au[i] = x[2];
      
      b[i] += au[i] * (-(dp1*dp1)/(2.0*mu)) * jf;
    }
    //correction on the right
    {
      valarray<double> J(9);
      J[2] = 1.0; J[1] = 1.0; J[0] = 1.0 - (d*d)/(2.0*mu) * jC;
      J[3] = d ; J[4] = dp1 ; J[5] = dp2;
      J[6] = J[3]*J[3]/2.0; J[7] = J[4]*J[4]/2.0 ; J[8] = J[5]*J[5]/2.0;
      valarray<double> r(3);
      r[0] = 0.0 ;
      r[1] = mux ;
      r[2] = mu ;
      valarray<double> x(3);
      solve3(J,r,x);
      x*=-h;

      al[i+1] = x[0];
      ad[i+1] = x[1] + Ce*h;
      au[i+1] = x[2];

      b[i+1] += al[i+1] * (-((d)*(d))/(2.0*mu)) * jf;
    }
    
  }

}


void  TDMAPoissonSolver::solve
  (valarray<double>& u)
{

  assert(u.size() == n);

  p[0] = - au[0]/ad[0];
  q[0] = b[0]/ad[0];
  
  for (unsigned int i = 1; i < n ; i++)
    {
      double f = 1.0 / (ad[i] + al[i] * p[i-1]);
      p[i] = -au[i] * f;
      q[i] = (b[i] - al[i]*q[i-1]) * f;
    }
  
  u[n-1] = q[n-1];
  
  for (int i = n-2; i >=0; i--)
    {
      u[i] = p[i]*u[i+1] + q[i];
    }

  
}

void  TDMAPoissonSolver::computeResidual
(const valarray<double>& u,valarray<double>& r)
{
  assert(r.size() == n);
  assert(u.size() == n);
  
  r[0] = ad[0]*u[0] +  au[0]*u[1] - b[0];
  for (unsigned int i = 1; i < n-1 ; i++)
    {
      r[i] = ad[i]*u[i] + al[i]*u[i-1] + au[i]*u[i+1] - b[i];
    }
  r[n-1] = ad[n-1]*u[n-1] +  al[n-1]*u[n-2] - b[n-1];
}

void TDMAPoissonSolver::setCoeffs
(double h,const valarray<double>& mu,
 const valarray<double>& C, const valarray<double>& r)
{
  //assert(mu.size() == n + 1);
  //assert(C.size() == n);
  //assert(r.size() == n);
  //bulk PDE
  for (int i = 1; i < n -1; ++i)
    {
      au[i] = - mu[i+1]/h;
      al[i] = - mu[i]/h ;
      ad[i] = C[i]*h - al[i] - au[i];
      b[i] = r[i]*h;
    }

  if (periodic)
    {
      //not done yet...
    }
  else {
    //other boundary conditions,
    //imposed through sources
    {
      // x = 0
      int i = 0;
      al[i] = 0.0;
      au[i] = - mu[i+1]/h;
      ad[i] = C[i]*h - al[i] - au[i];
      b[i] = r[i]*h;
    }
    {
      // x = L
      int i = n -1 ;
      al[i] = - mu[i]/h;
      au[i] = 0.0;
      ad[i] = C[i]*h - al[i] - au[i];
      b[i] = r[i]*h;
      
    }
  }


  // for (int i = 0; i < n; ++i)
  //   {
  //     au[i] = au[i] / ad[i];
  //     al[i] = al[i] / ad[i];
  //     b[i] = b[i] / ad[i];
  //     ad[i] = 1.0;
  //   }

  double scale = 1.0 / (EPS + abs(b.max()));
  if (scale < 1.0){
    au *= scale;
    al *= scale;
    ad *= scale;
    b *= scale;
  }
}



// compute Cu = C |u|^(m - 1)
void computeCu (const valarray<double>& C,
		const valarray<double>& u,
		const double m,
		valarray<double>& Cu)
{

  for (int i = 0; i < C.size(); ++i)
    Cu[i] = C[i] * pow ( std::abs(u[i]) + 1.0e-6 , m - 1.0); 
} 



class PicardSolver : public VelocitySolver
{
  valarray<double> ad, au, al, b, mu, res, cCu;
  bool periodic;
  int maxIter;
  double tol, thresh;
  PicardSolver(); 


public:

  PicardSolver(int n, bool periodic = false, int maxIter = 100, 
	       double tol = 1.0e-8, double thresh = 1.0e-12) 
    
    : ad(n), au(n), al(n), b(n), mu(n+1), res(n), cCu(n), 
      periodic(periodic), maxIter(maxIter), tol(tol), thresh(thresh)
  {;}
  
  void solve(const VelocityModel& vModel,
	     const valarray<double>& fA,
	     const valarray<double>& cH,
	     const valarray<double>& cRhs,
	     const valarray<double>& cC0,
	     const valarray<double>& cC,
	     const valarray<double>& fH,
	     const valarray<double>& cS,
	     const valarray<double>& cR,
	     const vector<bool>& floating, 
	     const double n, const double m,
	     const double dx, const GLCorrection& glCorrection,
	     valarray<double>& u, valarray<double>& res);
};


double glGradient(double landGrad, double glGrad, double seaGrad)
{
  return (landGrad + glGrad + seaGrad)/3.0; // gradUmeanOfThree
  //return (abs(seaGrad) < abs(glGrad))?seaGrad:glGrad;
}


//solve the nonlinear velocity problem with Picard iterations
void PicardSolver::solve(const VelocityModel& vModel,
			 const valarray<double>& fA,
			 const valarray<double>& cH,
			 const valarray<double>& cRhs,
			 const valarray<double>& cC0,
			 const valarray<double>& cC,
			 const valarray<double>& fH,
			 const valarray<double>& cS,
			 const valarray<double>& cR,
			 const vector<bool>& floating, 
			 const double n, const double m,
			 const double dx, const GLCorrection& glCorrection,
			 valarray<double>& cU, valarray<double>& res)
  
{
  int iter = 0;
  bool done = false;
  double res0Norm = 0.0;
  int nx = cU.size();
  valarray<double> fMu(0.0,nx+1);

  TDMAPoissonSolver poissonSolver(nx, periodic);
  valarray<double> fGradU(0.0,nx+1);
  //std::cout << "residual norms:"; 

  do {
    
    


    faceGradient(cU,dx,fGradU);
    //glCorrection.correctGradU(fGradU,floating); 
    vModel.computeFaceMu(fGradU,  fMu);
    //glCorrection.correctMuH(fMu,floating); 
 
    // for (int i = 1; i < nx-1; ++i){
    // 	  if (floating[i]){
	   
    // 	    if (!floating[i-1] || !floating[i+1])
    // 	      {	
		
    // 		assert(abs(fMu[i] - std::max(fMu[i+1],fMu[i-1])) < 1.0) ;
		
    // 	      }
    // 	  }
    // }
	
    computeCu(cC, cU , m , cCu);
    cCu += cC0;

    poissonSolver.setCoeffs(dx, fMu, cCu , cRhs);


    double resNorm = poissonSolver.residualNorm(cU);
    //double resNorm = poissonSolver.maxNorm(cU);
    //std::cout << resNorm << ":" ; 

    MY_ASSERT(isfinite(resNorm))

    if (iter == 0)
      res0Norm = resNorm;

    res[iter] = resNorm;
    done = iter > maxIter || resNorm < tol * res0Norm || resNorm < thresh;

    if (!done){
      poissonSolver.solve(cU);
    }
    
    iter++;

  } while(!done);
  //std::cout << std::endl;
}


double qSchoof(double H, double A, double C, double n, double m)
{
  return std::pow(H, (n+m+3)/(m+1.0)) 
    * std::pow( A * std::pow(RHO*GRAVITY,n+1.0) * std::pow(1.0-RHO/RHOW,n) 
		/ (std::pow(4.0,n)*C), 1.0/(m+1.0));

}

void computeVelocity
(double n, double m, double L,
 const valarray<double> &cR,
 const valarray<double> &fA,
 valarray<double> &cC,
 const valarray<double> &cH,
 const valarray<double> &fH,
 valarray<double>& cUb,
 valarray<double>& fUs,
 valarray<double>& fUa,
 valarray<double>& fD,
 valarray<double>& cS,
 boundary_condition bc0,
 boundary_condition bcL,
 model_type modelType,
 int maxIter, double relTol, double absTol,
 const GLCorrection& glCorrection, double regPar,
 valarray<double>& nlRes)
{
 
  int nx = cH.size();
  double dx = L/double(nx);
  
  //cell centred surface elevation, and basal traction coeff.
  //also bool to indicate cell-centred flotation
  valarray<double>& C = cC;
  vector<bool> floating(nx,false);


  bool first = true;
  //surface gradient
  for (int i = 0; i < nx; ++i)
    {
      double sg= cR[i]+cH[i];
      double sf = cH[i]*(1.0-RHO/RHOW);
      cS[i] = max(sg,sf);
      floating[i] = sf > sg + EPS;
      if (floating[i])
	{
	  if (first)
	    {
	      //std::cout << " first floating cell is " << i << std::endl;
	      first = false;
	    }

	  C[i] = 0.0;
	}
    }
  // bulk rhs of elliptic equation rho * g * H * s'
  valarray<double> rhs(0.0,nx);
  double twoDx = 2.0 * dx;
  for (int i = 1; i < nx-1; ++i){
    rhs[i] = -(cS[i+1] - cS[i-1])/twoDx * cH[i] * RHO * GRAVITY;
  }
  rhs[0] = -(cS[1] - cS[0])/dx * cH[0] * RHO * GRAVITY;
  rhs[nx-1] = -(cS[nx-1] - cS[nx-2])/dx * cH[nx-1] * RHO * GRAVITY;


  //grounding line correction of rhs
  //glCorrection.correctRhs(cS, cH , floating,  RHO*GRAVITY/dx, rhs);
  for (int i = 1; i < nx-1; ++i){
    if (floating[i] && !floating[i-1]){
      rhs[i] = -(cS[i+1] - cS[i])/dx * 0.5 * (cH[i+1] + cH[i]) * RHO * GRAVITY;
      int j = i-1;
      rhs[j] = -(cS[j] - cS[j-1])/dx * 0.5*  (cH[j] + cH[j-1]) * RHO * GRAVITY;
    }
  }
  
  //regularization ; set C0 and rhs to lambda * Hf and lambda * qSchoof
  //                 at the grounding line
  valarray<double> C0(cC);
   C0 = 0.0;
  // for (int i = 1; i < nx-1; ++i){
  //   //fD[i] = 0.0;
  //   if (floating[i])
  //     {
	
  // 	if (!floating[i-1] && floating[i+1])
  // 	  {
  // 	    double lambda = std::abs(regPar);
  // 	    double Hf = - cR[i-1] * RHOW / RHO;
  // 	    double q = qSchoof( Hf,fA[i],cC[i-1], n, m);
  // 	    C0[i] = lambda * Hf;
  // 	    rhs[i] += lambda * q;
  // 	  }

  // 	if (!floating[i+1] && floating[i-1])
  // 	  {
  // 	    double lambda = 1.0e-1;
  // 	    double Hf = - cR[i+1] * RHOW / RHO;
  // 	    double q = qSchoof( Hf,fA[i+1],cC[i+1], n, m);
  // 	    C0[i] = lambda * Hf;
  // 	    rhs[i] += lambda * q;
  // 	  }

  //     }
  //   else
  //     {
  // 	;
  //     }
  // }
  
  // boundary condition rhs
  switch (bc0)
    {
    case DIVIDE:
      C[0] += 1.0/ EPS;
      //rhs[0] = 0.0;
      break;
    case MARINE:
      //rhs[0] += pow(fH[0],2.0) * 0.5 * GRAVITY*RHO/RHOW * (RHOW-RHO) / dx;
      rhs[0] = 0.5 * cH[0] * (cS[1] - 0.0) * GRAVITY*RHO / dx;
      break;
    case PERIODIC:
      rhs[0] = -(cS[1] - cS[nx-2])/twoDx * cH[0] * RHO * GRAVITY;
      break;
    default:
      //doing nothing gives a natural boundary condition
      break;
    }

  switch (bcL)
    {
    case DIVIDE:
      C[nx-1] += 1.0/ EPS;
      //rhs[nx-1] = 0.0;
      break;
    case MARINE:
      rhs[nx-1] = -0.5 * cH[nx-1] * (0.0 - cS[nx-2]) * GRAVITY*RHO / dx;
      break;
    case PERIODIC:
      rhs[nx-1] = rhs[0];
      break;
    default:
      //doing nothing gives a natural boundary condition
      break;
    }


  VelocityModel* vModel_ptr;
  

  switch (modelType){
  case SSA:
    vModel_ptr = new SSAVelocityModel(fH,fA,dx,n);
  
    break;
  case L1L2:
    valarray<double> fGradS(0.0,nx+1);
    faceGradient(cS, dx, fGradS, bc0 == PERIODIC );

    valarray<double> cGradS(0.0,nx);
    for (int i = 1; i < nx - 1; ++i)
      cGradS[i] = 0.5 * (cS[i+1] - cS[i-1])/dx;
    
    if (bc0 == PERIODIC && bcL == PERIODIC)
      {
	cGradS[0] = cGradS[nx - 1] = 0.5 * (cS[1] - cS[nx-2])/dx;
      }
    else 
      {
	;//other boundary conditions have grad(s) = 0
      }
    
    valarray<double> cA(0.0,nx);
    for (int i = 0; i < nx; ++i)
      cA[i] = 0.5 * (fA[i] + fA[i+1]);


    vModel_ptr = new L1L2VelocityModel(fH,cH,fA,cA,fGradS,cGradS,floating,dx,n);
   
    break;
  }

  VelocityModel& vModel = *vModel_ptr;

  PicardSolver p(nx, false, maxIter, relTol, absTol);
  p.solve(vModel, fA, cH, rhs , C0, C, fH, cS, cR,  
	  floating, n, m, dx , glCorrection, cUb, nlRes);


  vModel.computeFaceVelocity(cUb, fUa, fD);
 

  delete vModel_ptr;

}

//advance thickness, Collala PPM, 
//given a face-centred velocity profile
void advanceThicknessPPM(double L, double dt,
		      const valarray<double>& cHold,
		      const valarray<double>& ca,
		      const valarray<double>& fU,
		      const valarray<double>& fD,
		      valarray<double>& cHnew)
{

  int nx = cHold.size();
  assert(cHnew.size() == nx);
  assert(ca.size() == nx);
  assert(fU.size() == nx + 1);
  double dx = L / double(nx);  
  
  valarray<double> adflux(0.0,nx+1); //advective flux
  valarray<double> Hhalf(0.0,nx+1); //face H at half times
  
  const valarray<double>& W = cHold;

  // undivided slopes of H, left face, right face, centre
  valarray<double> dW(0.0,nx); 
  valarray<double> dWp(0.0,nx); 
  valarray<double> dWm(0.0,nx); 
  for (int i = 1; i < nx ; ++i){
    dWm[i] = dWp[i-1] = W[i] - W[i-1];
  }
  for (int i = 1; i < nx-1 ; ++i){
    dW[i] = 0.5*(dWm[i]+dWp[i]);
  }
  dW[0] = dWp[0];
  dW[nx-1] = dWm[nx-1];

  //van leer limiter
  for (int i = 0; i < nx ; ++i)
    {
      double lim = 2.0*std::min(std::abs(dWm[i]),std::abs(dWp[i]));
      lim = std::min(lim,std::abs(dW[i]));
      if (dWm[i]*dWp[i] < 0.0)
	{
	  dW[i] = 0.0;
	} 
      else
	{
	  double s = (dWm[i]>0.0)?1.0:(-1.0);
	  dW[i] = lim*s;
	}
      
    }

  // ppm face values
  valarray<double> Wface(0.0,nx);
  for (int i = 1; i < nx ; ++i)
    {
      Wface[i] = 0.5 * ( W[i-1] + dW[i-1]/3.0 + W[i] - dW[i]/3.0);
    }

  valarray<double> Wm(W);
  valarray<double> Wp(W);
  Wm*=-1;Wp*=-1;
  for (int i = 0; i < nx ; ++i){
    Wm[i] +=  Wface[i];
    Wp[i] +=  Wface[i+1];
  }

  // ppm limiter
  for (int i = 0; i < nx ; ++i){
    double q = Wm[i]*Wp[i];
    if (q < 0.0)
      {
	double p = (Wp[i]+Wm[i])*0.5;
	double m = (Wp[i]-Wm[i]);
	double s = (p > 0)?1.0:(-1.0);
	if (p*m > 0.0)
	  {
	    Wp[i] = s*std::min(-2.0*s*Wm[i],s*Wp[i]);
	  }
	else
	  {
	    Wm[i] = s*std::min(-2.0*s*Wp[i],s*Wm[i]);
	  }
      }
    else
      {
	Wm[i]=Wp[i]=0.0;
      }

  }

  
  valarray<double> u(0.0,nx);
  cellAverage(fU,u);
  // ppm predictor
  for (int i = 0; i < nx ; ++i){
    double sigmin = -std::min(u[i]*dt/dx,0.0);
    double sigmax = std::max(u[i]*dt/dx,0.0);
    double sigm=0.0;
    double sigp=0.0;
    if (u[i] > 0.0)
      {
	sigm = sigmin;
	sigp = u[i]*dt/dx;
      } 
    else
      {
	sigm = -u[i]*dt/dx;
	sigp = sigmax;
      }

    double dm = Wm[i];
    double dp = Wp[i];

    Wm[i] = dm + sigm * ( (dp-dm) - (dp+dm)*(3.0-2.0*sigm)) * 0.5;
    Wp[i] = dp + sigp * ( (dm-dp) - (dp+dm)*(3.0-2.0*sigp)) * 0.5;

  }

  // post predictor normal source
  for (int i = 0; i < nx ; ++i){
    double dudx = (fU[i+1]-fU[i]) * dt /dx * 0.5;
    Wp[i] -= dudx * W[i];
    Wm[i] -= dudx * W[i];
  }


  Wm+=W; Wp+=W;

  Wm+=0.5*ca*dt;
  Wp+=0.5*ca*dt;

  //diffusion source
  valarray<double> faceGradH(0.0,nx+1);
  for (int i = 1; i < nx; ++i)
    {
      faceGradH[i] = (cHold[i]-cHold[i-1])/dx;
    }
  for (int i = 0; i < nx; ++i)
    {
      double s = 0.5*dt*(fD[i+1]*faceGradH[i+1]-fD[i]*faceGradH[i])/dx;
      Wm[i]+=s;
      Wp[i]+=s;
    }

  //Reimman problem
  for (int i = 1; i < nx ; ++i){
    Hhalf[i] = (fU[i] > 0.0)?Wp[i-1]:Wm[i];
  }

  adflux = fU * Hhalf;

  //  // subtract diffusion
  // valarray<double> faceGradH(0.0,nx+1);
  // for (int i = 1; i < nx; ++i)
  //   {
  //     faceGradH[i] = (cHold[i]-cHold[i-1])/dx;
  //   } 
  
  // adflux += faceGradH * fD;

  //upstream boundary
  adflux[0] = 0.0;
  //downstream boundary
  adflux[nx] = fU[nx] * cHold[nx-1];

  TDMAPoissonSolver hsolver(nx,false);
  valarray<double> dtD(dt,nx+1);dtD*=fD;
  valarray<double> one(1.0,nx);
  valarray<double> rhs(0.0,nx);
  for (int i = 0; i < nx ; ++i){
    rhs[i] = cHold[i] + dt * (ca[i]-(adflux[i+1]-adflux[i])/dx);
    //cHnew[i] = cHold[i] + dt * (-(adflux[i+1]-adflux[i])/dx + ca[i]);
  }
  hsolver.setCoeffs(dx,dtD,one,rhs);
  hsolver.solve(cHnew);

}		     
		     
//advance thickness, first order upwind, 
//given a face-centred velocity profile
void advanceThicknessFO(double L, double dt,
		      const valarray<double>& cHold,
		      const valarray<double>& ca,
		      const valarray<double>& fU,
		      const valarray<double>& fD,
		      valarray<double>& cHnew)
{

  int nx = cHold.size();
  assert(cHnew.size() == nx);
  assert(ca.size() == nx);
  assert(fU.size() == nx + 1);
    
  double dx = L / double(nx);

  valarray<double> adflux(0.0,nx+1); //advective flux
  for (int i = 1; i < nx ; ++i){
    //simple upwind
    adflux[i] = max(fU[i],0.0) * cHold[i-1] 
      + min(fU[i],0.0) * cHold[i];
  }
  //upstream boundary
  adflux[0] = 0.0;
  //downstream boundary
  adflux[nx] = fU[nx] * cHold[nx-1]; 

  
  TDMAPoissonSolver hsolver(nx,false);
  valarray<double> dtD(dt,nx+1);dtD*=fD;
  valarray<double> one(1.0,nx);
  valarray<double> rhs(0.0,nx);
  for (int i = 0; i < nx ; ++i){
    rhs[i] = cHold[i] + dt * (ca[i]-(adflux[i+1]-adflux[i])/dx);
    //cHnew[i] = cHold[i] + dt * (-(adflux[i+1]-adflux[i])/dx + ca[i]);
  }
  hsolver.setCoeffs(dx,dtD,one,rhs);
  hsolver.solve(cHnew);
}



//advance thickness for time-time0 years, solving the momentum equation
//to determine the velocity field

void advanceThickness(double n, double m, double L, 
		      double time0, double time,
		      double cfl,
		      const valarray<double>& cR,
		      const valarray<double>& fA,
		      const valarray<double>& ca,
		      valarray<double>& cC,
		      const valarray<double>& cC0,
		      const valarray<double>& cHstart,
		      valarray<double>& cHend,
		      valarray<double>& cUb,
		      valarray<double>& fUs,
		      valarray<double>& fUa,
		      valarray<double>& cS,
		      boundary_condition bc0,
		      boundary_condition bcL,
		      model_type modelType,
		      int maxIter, double relTol, 
		      double absTol, 
		      const GLCorrection& glCorrection,
		      double regPar)
{
  
  int nx = cHstart.size();
  double dx = L / double(nx); 
  valarray<double> cH(cHstart);
  valarray<double> fH(0.0, nx + 1);
  
  valarray<double> nlRes(2*maxIter); 

  double t = time0;
  while (t <= time){

    for (int i = 1; i < nx; ++i){
      fH[i] = 0.5 * (cH[i] + cH[i-1]);
    } 

    fH[0] = cH[0];
    fH[nx] = cH[nx-1];
    valarray<double> fD(0.0,nx+1);
    cC = cC0;
    computeVelocity(n,m,L,cR,fA,cC,cH,fH,cUb,fUs,fUa,fD,cS,bc0,bcL,
		    modelType, maxIter, relTol, absTol, glCorrection,
		    regPar,nlRes);

    double eps = 1.0e-9;
    double umax  = maxnorm(fUa);
    
    double dt = min( cfl * dx / (abs(fUa).max() + 1.0), eps + time -t);


    dt = std::pow(2.0, std::floor(std::log(dt)/std::log(2.0)));
    t += dt;
    //std::cout << " t = " << t << " dt = " << dt  << " umax = "  << umax <<  std::endl;
    MY_ASSERT(umax < 10000);
    {
      int nx = cUb.size();
      fUa[nx] = 2*cUb[nx-1] - fUa[nx-1];
    }
    advanceThicknessPPM(L,dt,cH,ca,fUa,fD,cHend);
    cH = cHend;

    

  }
}


#if 0

//C interface, designed to be accessed from R, Fortran, etc
extern "C" {



  // advance the cell-centred thickness field, explicitly, 
  // though an advection equation
  void advanceThickness
  (int *nx_ptr, 
   double *n_ptr, 
   double *m_ptr, 
   double *L_ptr, 
   double *dt_ptr,
   double *cR_ptr,
   double *fA_ptr, 
   double *cC_ptr,
   double *ca_ptr,
   double *cHstart_ptr,
   double *cHend_ptr, 
   double *cUb_ptr,
   double *fUs_ptr, 
   double *fUa_ptr, 
   double *cS_ptr, 
   int *bc_ptr, 
   int *model_ptr, 
   int *maxIter_ptr, 
   double *relTol_ptr, 
   double *absTol_ptr, 
   int *gl_rhs_c_ptr,
   int *gl_mu_c_ptr,
   double *regPar_ptr
   )
  {

    if (nx_ptr && n_ptr && m_ptr && L_ptr && dt_ptr
	&& cR_ptr && cHstart_ptr && cHend_ptr
	&& fA_ptr && cC_ptr && cUb_ptr
	&& fUs_ptr && fUa_ptr && cS_ptr 
	&& bc_ptr && model_ptr && maxIter_ptr 
	&& relTol_ptr && absTol_ptr && 
	gl_rhs_c_ptr && gl_mu_c_ptr && regPar_ptr)
      {

	int nx = *nx_ptr;
	
	//copy in
	valarray<double> cHstart(cHstart_ptr, nx);
	valarray<double> cHend(cHstart_ptr, nx);
	valarray<double> fA(fA_ptr, nx + 1);
	valarray<double> cC(cC_ptr, nx);
	valarray<double> ca(ca_ptr, nx);
	valarray<double> cR(cR_ptr, nx);
	valarray<double> cS(cS_ptr, nx);
	valarray<double> cUb(cUb_ptr, nx);
	valarray<double> fUs(fUs_ptr, nx + 1);
	valarray<double> fUa(fUa_ptr, nx + 1);

	GLCorrection glCorrection(*gl_rhs_c_ptr , *gl_mu_c_ptr);

	
	model_type modelType = SSA;
	if (*model_ptr != 0)
	  modelType = L1L2;
	
	boundary_condition bc[2] = {NATURAL,NATURAL};
	for (int i = 0; i < 2; ++i)
	  {
	    switch (bc_ptr[i])
	      {
	      case 0: 
		bc[i] = NATURAL; break;
	      case 1:
		bc[i] = PERIODIC; break;
	      case 2: 
		bc[i] = DIVIDE; break;
	      case 3:
		bc[i] = MARINE; break;
	      default:
		break;
	      }
	  }
	

	advanceThickness(*n_ptr,*m_ptr,*L_ptr, 0.0, *dt_ptr, 0.5,
			 cR, fA, ca, cC, cHstart, cHend,
			 cUb, fUs, fUa, cS, bc[0], bc[1],
			 modelType, *maxIter_ptr, *relTol_ptr, 
			 *absTol_ptr, glCorrection, *regPar_ptr);

	

	//copy out
	for (int i =0; i < nx; ++i)
	  cHend_ptr[i] =cHend[i];
	for (int i =0; i < nx; ++i)
	  cUb_ptr[i] = cUb[i];
	for (int i =0; i < nx; ++i)
	  cS_ptr[i] = cS[i];      
	for (int i =0; i < nx + 1; ++i)
	  fUs_ptr[i] = fUs[i];
	for (int i =0; i < nx + 1; ++i)
	  fUa_ptr[i] = fUa[i];
	


      }
   
  }



  //compute cell centred base velocity ub and face centred surface and
  //vertically averaged velocities us and ua on a regular grid
  //of nx poinuts with x[1] = 0, x[i] = i - h, x[n] = L
  //
  //on a regular mesh of length(cell_H).
  //an elliptic equation
  //
  // 1. ( 4 * mu * H * ub')' - C * |ub|^m-1 * ub = rho * g * H * s'
  //
  //is solved for ub
  //
  // boundary conditions are either
  //  0 "divide" - u[1] (or u[n]) = 0
  //  1 "marine" - a Neumman condition with a flux depending on A and H
  //  2 "periodic" u[0] = u[n], u'[0] = u'[n]
  //
  // mu is either computed from Glens's law (in the SSA case)
  //
  // 2. mu = 1/2 A^{-1/n} |du/dx|^{1/n-1} 
  //
  // or through L1L2 model, with 2. as the underlying relation
  //
  // In the SSA case us and ua are identical and are simply face
  // centred averages of ua. In the L1L2 case, they include additional
  // contributions
  void computeVelocity
  (int *nx_ptr, 
   double *n_ptr, 
   double *m_ptr, 
   double *L_ptr,
   double *cR_ptr,
   double *fA_ptr,
   double *cC_ptr,
   double *cH_ptr, 
   double *fH_ptr, 
   double *cUb_ptr,
   double *fUs_ptr, 
   double *fUa_ptr, 
   double *cS_ptr, 
   int *bc_ptr, 
   int *model_ptr, 
   int *maxIter_ptr, 
   double *relTol_ptr, 
   double *absTol_ptr, 
   int *gl_rhs_c_ptr,
   int *gl_mu_c_ptr,
   double *reg_par_ptr,
   double *nlres_ptr)
  {
 

   if (nx_ptr && n_ptr && m_ptr && L_ptr
	&& cR_ptr && cH_ptr && fH_ptr
	&& fA_ptr && cC_ptr && cUb_ptr
	&& fUs_ptr && fUa_ptr && cS_ptr 
	&& bc_ptr && model_ptr && maxIter_ptr 
	&& relTol_ptr && absTol_ptr 
       && gl_rhs_c_ptr && gl_mu_c_ptr
        && nlres_ptr && reg_par_ptr)
      {
      int nx = *nx_ptr;

    //copy in
    valarray<double> cH(cH_ptr, nx);
    valarray<double> fH(fH_ptr, nx + 1);
    valarray<double> fA(fA_ptr, nx + 1);
    valarray<double> cC(cC_ptr, nx);
    valarray<double> cR(cR_ptr, nx);
    valarray<double> cS(cS_ptr, nx);
    valarray<double> cUb(cUb_ptr, nx);
    valarray<double> fUs(fUs_ptr, nx + 1);
    valarray<double> fUa(fUa_ptr, nx + 1);
    valarray<double> nlRes(0.0,*maxIter_ptr + 100);

   
    GLCorrection glCorrection(*gl_rhs_c_ptr,*gl_mu_c_ptr);
	
    model_type modelType = SSA;
    if (*model_ptr != 0)
      modelType = L1L2;

    boundary_condition bc[2] = {NATURAL,NATURAL};
    for (int i = 0; i < 2; ++i)
      {
        switch (bc_ptr[i])
	  {
	  case 0: 
	    bc[i] = NATURAL; break;
	  case 1:
	    bc[i] = PERIODIC; break;
	  case 2: 
	    bc[i] = DIVIDE; break;
	  case 3:
	    bc[i] = MARINE; break;
	  default:
	    break;
	  }
      }

    valarray<double> fD(0.0,nx+1);
    computeVelocity(*n_ptr, *m_ptr, *L_ptr, cR, fA,cC, cH,  fH,
		    cUb, fUs, fUa, fD, cS, bc[0], bc[1],
		    modelType, *maxIter_ptr, *relTol_ptr, *absTol_ptr, 
		    glCorrection, *reg_par_ptr, nlRes);

    //copy out
    for (int i =0; i < nx; ++i)
      cUb_ptr[i] = cUb[i];
    for (int i =0; i < nx; ++i)
      cS_ptr[i] = cS[i];      
    for (int i =0; i < nx + 1; ++i)
      fUs_ptr[i] = fUs[i];
    for (int i =0; i < nx + 1; ++i)
      fUa_ptr[i] = fUa[i];
    for (int i =0; i < *maxIter_ptr; ++i)
      nlres_ptr[i] = nlRes[i];
      }
  }

}

#endif

#ifdef MAINA


#ifdef HDF
#include "hdf5.h"
//write each valarray x[i] to set s[i] in hdf5 file filename
void fout(const string& filename,
	    const vector<string>& s, 
	    const vector<valarray<double>* >& x)
{

  MY_ASSERT(s.size() == x.size());
  herr_t      status;

  hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, 
			    H5P_DEFAULT, H5P_DEFAULT);
  
  hsize_t n = (*x[0]).size();
  hid_t  space_id = H5Screate_simple(1,&n,NULL);
  
  for (int i = 0; i < s.size(); ++i)
    {
      hid_t dset_id = H5Dcreate(file_id, s[i].c_str(), 
				H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);

      status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			&(*x[i])[0]);


      status = H5Dclose(dset_id);
    }
  status = H5Sclose(space_id);
  status = H5Fclose(file_id); 


}

void fin (const string& filename,
	  const vector<string>& s, 
	  const vector<valarray<double>* >& x)
{

  
  MY_ASSERT(s.size() == x.size());
  herr_t      status;
  hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY,  H5P_DEFAULT);

  hsize_t n = (*x[0]).size();
  hid_t  space_id = H5Screate_simple(1,&n,NULL);

  for (int i = 0; i < s.size(); ++i)
    {
      hid_t dset_id = H5Dopen(file_id,s[i].c_str());
      
      status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			&(*x[i])[0]);
      
      status = H5Dclose(dset_id);
    }

}

#else

void fout(const string& fnam,
	  const vector<string>& s, 
	  const vector<valarray<double>* >& x)
{

  MY_ASSERT(s.size() == x.size());
  int n = (*x[0]).size();

  ofstream os; os.open(fnam.c_str());

  for (int j = 0; j < s.size(); j++)
    {
      os << s[j] << " ";
    }
  os << endl;
  for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < s.size(); j++)
	{
	  os << (*x[j])[i] << " ";
	}
      os << endl;
    }

  os.close();
}

#endif

int main(int argc, char* argv[]){


  if (argc != 5 && argc != 6)
    {
      cerr << "usage: " << argv[0] << " nx model accumulation (mm) output [input]" << endl;
      exit(EXIT_FAILURE);
    }
  int nx = atoi(argv[1]);
  string modelname(argv[2]);
  model_type model;
  if (modelname == "l1l2")
    {
      model = L1L2;
    }
  else if (modelname == "ssa")
    {
      model = SSA;
    }
  else 
    {
      cerr << "unknown model: " << modelname << endl;
      exit(EXIT_FAILURE);
    }
 
  int acabmm = atoi(argv[3]);
#ifdef MISMIP3D
  double L =  800e+3;
  double A = 3.1536e-18;
  double C = 31651.76;
  double origin = -100.0;
  double slope = -0.001;
  double H0 = 150;
#endif
#ifdef MISMIP
 double L = 1800e+3;
 double origin = 710.0;
 double slope = -0.001038;
 double C = 24125.96;
 double A = 1.464788e-16;
 double H0 = 80;
#endif
  double n = 3.0;
  double m = 0.333333333333333;
  double h = L/double(nx);
 
  double cfl = 0.5;

  valarray<double>fX(0.0,nx+1);
  for (int i = 1; i < nx + 1; ++i){
    fX[i] = fX[i-1] + h;
 }
 
 valarray<double>x(0.0,nx);
 cellAverage(fX,x);

 valarray<double>fR(origin,nx+1);
 fR +=  slope * fX;

 valarray<double>cR(0.0,nx);
 cellAverage(fR,cR);

 valarray<double>cHstart(H0,nx);
 valarray<double>cHend(H0,nx);
 valarray<double>cS(0.0,nx);
 valarray<double>fA(A,nx+1);
 valarray<double>fUs(A,nx+1);
 valarray<double>fUa(A,nx+1);
 valarray<double>cC0(C,nx);
 valarray<double>cC(C,nx);
 valarray<double>cUb(0.0,nx);
 valarray<double>ca(double(acabmm)*1.0e-3,nx);
 GLCorrection glCorrection(1 , 0);
 gm_stop = false;
 vector<string> ccs(5);
 vector<valarray<double>* > cc(ccs.size());
 int j = 0;
 cc[j] = &x;ccs[j++]="x";
 cc[j] = &cHend;ccs[j++]="thck";
 cc[j] = &cR;ccs[j++]="topg";
 cc[j] = &cS;ccs[j++]="usrf";
 cc[j] = &cUb;ccs[j++]="uvelb";
#ifdef HDF
#define ext "hdf5"

 string fnamout(argv[4]);
 //input
 if (argc > 5)
     {
       string fnamin(argv[5]);
       
       fin(fnamin,ccs,cc);

       cHstart = cHend;
     }

   //extra output fields
   ccs.push_back("C");
   cc.push_back(&cC);

#else
#define ext "txt"
#endif
   double dt = 1000.0;
   double t = 0;
   for (int i =0; i < 50; ++i)
   {
     advanceThickness(n,m,L,t,t+dt, cfl, cR , fA, ca, cC, cC0, 
		      cHstart , cHend, cUb, fUs, fUa, 
		      cS, DIVIDE, MARINE, model, 5, 1e-6, 1e-6,
		      glCorrection,0.0);
     t += dt;
     std::cout << " nx = " << nx << " t = " << t << std::endl; 
     // stringstream ss; ss << modelname << "." << nx << "." << acabmm << "mm." << ext;;
     fout(fnamout,ccs,cc);

     cHstart = cHend;
   }

 {
 

  
 }

#if 0
 { 
   ofstream os; 
   os.open("thickness.dump");
   os.precision(16);
   for (int i =0; i < cHend.size(); ++i)
     os << cHend[i] << " ";
   os.close();
   
   os.open("velocity.dump");
   os.precision(16);
   for (int i =0; i < cUb.size(); ++i)
     os << cUb[i] << " ";
   os.close();
 }
#endif
 //  gm_stop = true;
 // // //this one is to compare with the BISICLES when it has read thickness.dump
 //  advanceThickness(n,m,L,t,t+dt, cfl, cR , fA, ca, cC, 
 //  		    cHstart , cHend, cUb, fUs, fUa, 
 //  		  cS, DIVIDE, MARINE, model, 3, 1e-2, 1e-6,
 //  		    glCorrection,0.0);

}


#endif
#ifdef  MAINB
#include "netcdf.h"
int wnc(const std::string filename,
	int ny, double y0, double y1,
	const valarray<double>& x,
	const valarray<double>& topg,
	const valarray<double>& usrf,
	const valarray<double>& lsrf,
	const valarray<double>& thck,
	const valarray<double>& beta)
{

  int nx = x.size();

  valarray<double> y(y0,ny);

 
  valarray<double> topg2(0.0,ny*nx);
  valarray<double> thck2(0.0,ny*nx);
  valarray<double> usrf2(0.0,ny*nx);
  valarray<double> lsrf2(0.0,ny*nx);
  valarray<double> beta2(0.0,(ny-1)*(nx-1));
  double dy = (y1-y0) / double(ny-1);
  for (int j = 1; j < ny; j++)
    {
      y[j] = y[j-1] + dy;
    }

 valarray<double> cX(0.0,nx-1);
  cellAverage(x,cX);
  valarray<double> cY(0.0,ny-1);
  cellAverage(y,cY);


  for (int j = 0; j < ny; j++)
    {
    
      for (int i = 0; i < nx; ++i)
	{
	  topg2[j*nx + i] = topg[i];
	  thck2[j*nx + i] = thck[i];
	  usrf2[j*nx + i] = usrf[i];
	  lsrf2[j*nx + i] = lsrf[i];
	}
    }

  for (int j = 0; j < ny-1; j++)
    {
      double f = std::abs(cY[j])/10.0e+3;
      f = 1.0 - 0.5 * exp(-f*f);
      f = 1.0;
      for (int i = 0; i < nx-1; ++i)
	{
	 
	  beta2[j*(nx-1) + i] = f* beta[i];

	}
    }

  int nc_id, var_id, rc;
  int node_dim_ids[3], cell_dim_ids[3];
 
  if ((rc = nc_create(filename.c_str(),NC_CLOBBER,&nc_id)))
    return 1;

  if ((rc = nc_def_dim(nc_id, "time", 1 , &node_dim_ids[0])))
    return 1;
  
  if ((rc = nc_def_dim(nc_id, "x1", nx, &node_dim_ids[2])))
    return 1;

  if ((rc = nc_def_dim(nc_id, "y1", ny, &node_dim_ids[1])))
    return 1;

  if ((rc = nc_def_dim(nc_id, "x0", nx-1, &cell_dim_ids[2])))
    return 1;

  if ((rc = nc_def_dim(nc_id, "y0", ny-1, &cell_dim_ids[1])))
    return 1;
  
  cell_dim_ids[0] = node_dim_ids[0];

  int x1_id;
  if ((rc = nc_def_var(nc_id, "x1", NC_DOUBLE,1,  &node_dim_ids[2], &x1_id)))
    return 1;
  
  int y1_id;
  if ((rc = nc_def_var(nc_id, "y1", NC_DOUBLE,1,  &node_dim_ids[1], &y1_id)))
    return 1;

  int x0_id;
  if ((rc = nc_def_var(nc_id, "x0", NC_DOUBLE,1,  &cell_dim_ids[2], &x0_id)))
    return 1;
  
  int y0_id;
  if ((rc = nc_def_var(nc_id, "y0", NC_DOUBLE,1,  &cell_dim_ids[1], &y0_id)))
    return 1;
  
  int time_id;
  if ((rc = nc_def_var(nc_id, "time", NC_DOUBLE,1,  &node_dim_ids[0], &time_id)))
    return 1;

  
  int thck_id;
  if ((rc = nc_def_var(nc_id, "thk", NC_DOUBLE,3,  node_dim_ids, &thck_id)))
    return 1;
  
  int topg_id;
  if ((rc = nc_def_var(nc_id, "topg", NC_DOUBLE,3,  node_dim_ids, &topg_id)))
    return 1;
 
  int usrf_id;
  if ((rc = nc_def_var(nc_id, "usrf", NC_DOUBLE,3,  node_dim_ids, &usrf_id)))
    return 1;

  int lsrf_id;
  if ((rc = nc_def_var(nc_id, "lsrf", NC_DOUBLE, 3, node_dim_ids, &lsrf_id)))
    return 1; 
 
  int beta_id;
  if ((rc = nc_def_var(nc_id, "beta", NC_DOUBLE,3,  cell_dim_ids, &beta_id)))
    return 1;


  if ((rc = nc_enddef(nc_id)))
    return 1;

  if ((rc = nc_put_var_double(nc_id, x1_id, &x[0])))
    return 1;
  
  if ((rc = nc_put_var_double(nc_id, y1_id, &y[0])))
    return 1;
  
  double time = 0.0;
  if ((rc = nc_put_var_double(nc_id, time_id, &time)))
    return 1;


  if ((rc = nc_put_var_double(nc_id, x0_id, &cX[0])))
    return 1;
  
  if ((rc = nc_put_var_double(nc_id, y0_id, &cY[0])))
    return 1;


  if ((rc = nc_put_var_double(nc_id, thck_id, &thck2[0])))
    return 1;
  if ((rc = nc_put_var_double(nc_id, topg_id, &topg2[0])))
    return 1;
  if ((rc = nc_put_var_double(nc_id, usrf_id, &usrf2[0])))
    return 1;
  if ((rc = nc_put_var_double(nc_id, lsrf_id, &lsrf2[0])))
    return 1;
  
  if ((rc = nc_put_var_double(nc_id, beta_id, &beta2[0])))
    return 1;
  


  if ((rc = nc_close(nc_id)))
    return 1;

}


int main(int argc, char* argv[]){

 double L =  800e+3;
 double n = 3.0;
 double m = 0.33333333333333;
 int nx = 256;
 nx = 512;;
 int ny = nx / 8 + 1;
 double h = L/double(nx);
 
 valarray<double>fX(0.0,nx+1);
 for (int i = 1; i < nx + 1; ++i){
   fX[i] = fX[i-1] + h;
 }
 
 valarray<double>x(0.0,nx);
 cellAverage(fX,x);

 valarray<double>fR(-100.0,nx+1);
 fR -= 0.001 * fX;

 valarray<double>cR(0.0,nx);
 cellAverage(fR,cR);

 valarray<double>cHstart(150.0,nx);
 valarray<double>cHend(150.0,nx);
 valarray<double>cS(0.0,nx);
 valarray<double>fA(3.1536e-18,nx+1);
 valarray<double>fUs(3.1536e-18,nx+1);
 valarray<double>fUa(3.1536e-18,nx+1);
 valarray<double>cC(31651.76,nx);
 valarray<double>cUb(0.0,nx);
 valarray<double>ca(0.49,nx);
 GLCorrection glCorrection(1 , 0);

 for (int i =0; i < 300; ++i)
   {

     //face centred data needed by Glimmer
     valarray<double>fS(0.0,nx+1);
     faceAverage(cS,fS);
     valarray<double>fH(0.0,nx+1);
     faceAverage(cHend,fH);
     valarray<double>fL(0.0,nx+1);
     fL = fS - fH;
     int r = wnc("ssa.a049.nc",ny,-50.0e+3,50.0e+3, fX, fR, fS, fL, fH,  cC);
     if (r != 0)
       //std::cout << "error writing netdcf file" << std::endl;
     
     cHstart = cHend;
     advanceThickness(n,m,L,100.0, cR , fA, ca, cC, 
 		  cHstart , cHend, cUb, fUs, fUa, 
 		  cS, DIVIDE, MARINE, SSA, 3, 1e-2, 1e-6,
 		  glCorrection,0.0);
   }
 cHend = 150.0;
 for (int i =0; i < 300; ++i)
   {
     //face centred data needed by Glimmer
     valarray<double>fS(0.0,nx+1);
     faceAverage(cS,fS);
     valarray<double>fH(0.0,nx+1);
     faceAverage(cHend,fH);
     valarray<double>fL(0.0,nx+1);
     fL = fS - fH;
     int r = wnc("l1l2.a049.nc",ny,-50.0e+3,50.0e+3, fX, fR, fS, fL, fH,  cC);
     if (r != 0)
       //std::cout << "error writing netdcf file" << std::endl;
     
     cHstart = cHend;
     advanceThickness(n,m,L,100, cR , fA, ca, cC, 
		      cHstart , cHend, cUb, fUs, fUa, 
		      cS, DIVIDE, MARINE, L1L2, 3, 1e-2, 1e-6,
		      glCorrection,0.0);
 
   }
 

}

#endif
