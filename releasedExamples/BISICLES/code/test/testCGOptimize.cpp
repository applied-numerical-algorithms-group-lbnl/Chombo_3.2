#include "CGOptimize.H"
#include "MayDay.H"
#include "AMRIO.H"
//simple test of CG optimize :: minimise \sum _n (x_n-n)^2

class RealVector
{
  Real* m_x;
  int m_n;
public:
  RealVector()
  {
    m_n = 0;
    m_x = NULL;
  }

  void define(int n)
  {
    if (m_x!=NULL)
      {
	delete[] m_x;
	m_x = NULL;
      }
    m_n = n;
    m_x = new Real[m_n];
  }

  RealVector(int n)
  {
    m_n = n;
    m_x = new Real[m_n];
  }
  ~RealVector()
  {
    delete[] m_x;
    m_x = NULL;
  }
  const int n() const
  {
    return m_n;
  }

  const Real& operator()(int i) const
   {
   
     if ( (i >= m_n) || (i < 0) || !m_x)
       MayDay::Error("i out of bounds");
       return m_x[i];
   }

  Real& operator()(int i)
  {
    if ( (i >= m_n) || (i < 0) || !m_x)
      MayDay::Error("i out of bounds");
    return m_x[i];
  }

};


class QuarticObjective : public ObjectiveWithGradient<RealVector>
{
public:

 
  void create(RealVector& a, const  RealVector& b)
  {
    a.define(b.n());
  }
  void free(RealVector& a)
  {
    
  }
  Real operator()(const  RealVector& x)
  {
    Real f = 0.0;
    for (int i =0; i < x.n(); i++)
      {
	Real z = x(i) - Real(i);
	f += z*z*z*z  / Real(i+1);
      }
    return f;
  }
  
  void restart(){}

  void computeObjectiveAndGradient(Real& f, Real& r, RealVector& g, const  RealVector& x, bool a_inner)
  {
    f = 0.0;
    r = 0.0;
    for (int i =0; i < x.n(); i++)
      {
	Real z = x(i) - Real(i);
	g(i) = 4.0 * z * z * z / Real(i+1);
	f+= z * z * z * z / Real(i+1);
      }
  }

  void scale(RealVector& x, const  Real s)
  {
    for (int i =0; i < x.n(); i++)
      x(i) *= s;
  }

  void assign(RealVector& y, const RealVector& x)
  {
    for (int i =0; i < x.n(); i++)
      y(i) = x(i);
  }

  void preCond(RealVector& s, const RealVector& r)
  {
    assign(s,r);
  }

  //y = y + x * s
  void incr(RealVector& y, const RealVector& x, Real s)
  {
    for (int i =0; i < x.n(); i++)
      y(i) = y(i) + s*x(i);
  }


  Real dotProduct(RealVector& y, const RealVector& x)
  {
    Real r = 0.0;
    for (int i =0; i < x.n(); i++)
      r+=x(i)*y(i);
    return r;
  }

  int nDoF(const RealVector& x)
  {
    return x.n();
  }

 
};


int main(int argc, char* argv[])
{

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  { // Begin nested scope

#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif
    int rank, number_procs;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
    MPI_Comm_size(Chombo_MPI::comm, &number_procs);
#else
    rank=0;
    number_procs=1;
#endif

  }

  int n = 300;
  RealVector x(n);
  QuarticObjective F;

  pout() << "initial F(x) = " << F(x) << std::endl;

  CGOptimize(F,x,100,1.0e-10,0.999,1.0e-3,-1.0,100,1.0e-10);

  Real final = F(x);

  int status = 0;
  if (Abs(final) > 5.0e-5) status = -1;

  pout() << "final F(x) = " << final << std::endl;

  if( status == 0 )
    pout() << argv[0] << " passed." << endl ;
  else
    pout() << argv[0] << " failed with return code " << status << endl ;


#ifdef CH_MPI
  MPI_Finalize();
#endif
  return status;
}
