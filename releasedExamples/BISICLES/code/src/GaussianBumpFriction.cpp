#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
#include "GaussianBumpFriction.H"
#include "NamespaceHeader.H"

GaussianBumpFriction::GaussianBumpFriction(const Vector<Real>& a_t,
					   const Vector<Real>& a_C0, 
					   const Vector<Real>& a_a,
					   const Vector<RealVect>& a_b,
					   const Vector<RealVect>& a_c)
  :m_t(a_t), m_C0(a_C0), m_a(a_a), m_b(a_b), m_c(a_c)
  {
    CH_assert(m_C0.size() >= m_t.size() + 1);
    CH_assert(m_a.size() >= m_t.size() + 1);
    CH_assert(m_b.size() >= m_t.size() + 1);
    CH_assert(m_c.size() >= m_t.size() + 1);
  } 


BasalFriction* GaussianBumpFriction::new_basalFriction() const
{
  GaussianBumpFriction* ptr = new GaussianBumpFriction(m_t,m_C0,m_a,m_b,m_c);
  return static_cast<BasalFriction*>(ptr);
}

void GaussianBumpFriction::setBasalFriction(LevelData<FArrayBox>& a_C,
					    LevelSigmaCS& a_coordSys,
					    Real a_time,
					    Real a_dt)
{

  //inefficient search for the correct time interval. I'm 
  //assuming there will only be a small number.
  int i = 0;
  while (i < m_t.size() && a_time >= m_t[i])
    i++;

  RealVect dx = a_coordSys.dx(); 
  RealVect d(RealVect::Unit);
  d /= m_c[i];
  d /= sqrt(2.0);
  for (DataIterator dit = a_C.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& C = a_C[dit];
      for (BoxIterator bit(C.box()); bit.ok(); ++bit)
	{
	  IntVect iv = bit();
	      RealVect x(iv);
	      x += 0.5*RealVect::Unit;
	      x *= dx;
	      x -= m_b[i];
	      x *= d;

	      C(iv,0) = m_C0[i] * (1 - m_a[i] * exp(- x.dotProduct(x)));
	}
    }
  }

#include "NamespaceFooter.H"
