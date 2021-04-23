#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PARTICLEBC_H_
#define _PARTICLEBC_H_

#include "ParticleData.H"
#include "DataIterator.H"
#include "REAL.H"
#include "NamespaceHeader.H"

class ParticleBC
{
public:

  /// enforces periodic boundary conditions on the particles in 
  /// a_particleData. Particles that are outside the range given 
  /// by a_leftEdge and a_rightEdge will have their positions 
  /// periodiocally shifted until that is no longer true
  template <class P>
  static void enforcePeriodic(ParticleData<P>& a_particleData,
                              const RealVect& a_leftEdge,
                              const RealVect& a_rightEdge);
};

#include "NamespaceFooter.H"

// implementation
#include "ParticleBCI.H"

#endif // include guard
