#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <fstream>
#include <iostream>
using std::cout;
using std::endl;
using std::ofstream;

#include "defineLevelSigmaCS.H"
#include "BoxIterator.H"


#include "NamespaceHeader.H"

// utility function for defining coordinate systems for testing
void defineLevelSigmaCS(LevelSigmaCS& a_sigmaCoords,
                        const RealVect& a_domainSize,
                        int a_thicknessType, 
                        int a_basalType,
			const RealVect& a_basalSlope,
                        Real a_thicknessScale)
{
  const DisjointBoxLayout& grids = a_sigmaCoords.grids();
  DataIterator dit = grids.dataIterator();

  const LevelData<FArrayBox>& zBrefLevel = a_sigmaCoords.getTopography();
  LevelData<FArrayBox> zBLevel(grids,1, zBrefLevel.ghostVect());

  const RealVect& dx = a_sigmaCoords.dx();

  //Real Cbase = 0.1;
  // set to 1.5 so that ice thickness doesn't actually go to zero
  Real Cbase = 1.5;
  
  for (dit.begin(); dit.ok(); ++dit)
    { 
       FArrayBox& zB = zBLevel[dit];

      // first test -- zb - C*sin(4*pi*x)*sin(4*pi*y)
      
       {
	 BoxIterator bit(zB.box());
	 for (bit.begin(); bit.ok(); ++bit)
	   {
	     IntVect iv = bit();
	     RealVect x(iv);
	     x += 0.5*RealVect::Unit;
              x *= dx;
              
              //Real baseHeight = x.dotProduct(a_basalSlope);
	      Real baseHeight = 0.0;
	      if (a_basalType == constantZb)
		{
		  baseHeight += Cbase;
		}
              else if (a_basalType == sinusoidalZb)
                {
                  if (SpaceDim == 2)
                    {
                      baseHeight += Cbase+sin(4*Pi*x[0]/a_domainSize[0])
                        *sin(4.0*Pi*x[1]/a_domainSize[1]);
                    }
                  else if (SpaceDim == 3)
                    {
                      baseHeight += D_TERM(Cbase,+sin(4*Pi*x[1]/a_domainSize[1]),
                                          *sin(4.0*Pi*x[2]/a_domainSize[2]));
                    }
                }
              else if (a_basalType == sinusoidalYZb)
                {
                  if (SpaceDim == 2)
                    {
                      baseHeight += 250.0+250.0*sin(2.0*Pi*x[1]/a_domainSize[1]);
                        
                    }
                  else if (SpaceDim == 3)
                    {
                      baseHeight += D_TERM(250.0,+250.0,
                                          *sin(2.0*Pi*x[2]/a_domainSize[2]));
                    }
                }
	      else if (a_basalType == xInclineZb)
		{
		  baseHeight += Cbase + x[SpaceDim-2] 
		    / a_domainSize[SpaceDim-2]; 
		}
	       else if (a_basalType == yInclineZb)
		{
		  baseHeight += Cbase + x[SpaceDim-1] 
		    / a_domainSize[SpaceDim-1]; 
		}
	       else if (a_basalType == pattynAZb)
		{
		  
		  
		  if (SpaceDim == 2)
                    {
		      Real Cslope = 1000.0;
		      baseHeight += Cslope  + a_thicknessScale * (sin(2.0*Pi*x[0]/a_domainSize[0])
						    *(sin(2.0*Pi*x[1]/a_domainSize[1])));
		    } 
		  else if (SpaceDim == 3)
		    {
		      Real Cslope = 1000.0 ;// - 0.008726868 * x[1];
		      baseHeight += D_TERM(Cslope,+  a_thicknessScale * sin(2*Pi*x[1]/a_domainSize[1]),
                                          * sin(2.0*Pi*x[2]/a_domainSize[2]));
		    }
		}
	       else if (a_basalType == pattynBZb)
		{
		  
		  
		  if (SpaceDim == 2)
                    {
		      Real Cslope = 1000.0;
		      baseHeight += Cslope  + a_thicknessScale * sin(2.0*Pi*x[0]/a_domainSize[0]);
		    } 
		  else if (SpaceDim == 3)
		    {
		      Real Cslope = 1000.0 ;// - 0.008726868 * x[1];
		      baseHeight += D_TERM(Cslope,+  a_thicknessScale * sin(2*Pi*x[1]/a_domainSize[1]), * 1.0);
		    }
		}
              else
                {
                  MayDay::Error("basalType undefined in defineSigmaCS");
                }
              zB(iv,0) = baseHeight;
            }
        } // end if not constant zB
    } // end loop over grids
  
  a_sigmaCoords.setTopography(zBLevel);
  a_sigmaCoords.setBackgroundSlope(a_basalSlope);
  a_sigmaCoords.exchangeTopography();
  // set thickness
  LevelData<FArrayBox>& levelH = a_sigmaCoords.getH();
  
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& H = levelH[dit];
      FArrayBox& zB = zBLevel[dit];
      if (a_thicknessType == constantThickness)
        {
          // constant thickness
          H.setVal(a_thicknessScale);
        } 
      else if (a_thicknessType == constantThickness1km)
        {
          // constant thickness
          H.setVal(1000.0);
        } 
      else if (a_thicknessType == pattynAH)
	{

	  BoxIterator bit(H.box());
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              RealVect x(iv);
              x += 0.5*RealVect::Unit;
              x *= dx;
	      if (SpaceDim == 2)
		H(iv,0) = 1000 - a_thicknessScale  * (sin(2.0*Pi*x[0]/a_domainSize[0])
					 *(sin(2.0*Pi*x[1]/a_domainSize[1])));
	      else if (SpaceDim ==3)
		H(iv,0) = D_TERM(1000,- a_thicknessScale * sin(2*Pi*x[1]/a_domainSize[1]),
                                          * sin(2.0*Pi*x[2]/a_domainSize[2]));
	    } 
	}
      else if (a_thicknessType == pattynBH)
	{

	  BoxIterator bit(H.box());
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              RealVect x(iv);
              x += 0.5*RealVect::Unit;
              x *= dx;
	      if (SpaceDim == 2)
		H(iv,0) = 1000 - a_thicknessScale  * sin(2.0*Pi*x[0]/a_domainSize[0]);
	      else if (SpaceDim ==3)
		H(iv,0) = D_TERM(1000,- a_thicknessScale * sin(2*Pi*x[1]/a_domainSize[1]),
                                          * 1.0);
	    } 
	}
      else if (a_thicknessType == constantZs1km)
        {
          //const FArrayBox& zB = a_sigmaCoords[dit].getTopography();
          BoxIterator bit(H.box());
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              Real thickness = 1000.0 - zB(iv,0);
              H(iv,0) = thickness*a_thicknessScale;
              
            }
        }
      else if (a_thicknessType == sinusoidalH)
        {
          BoxIterator bit(H.box());
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              RealVect x(iv);
              x += 0.5*RealVect::Unit;
              x *= dx;
              Real thickness;
              if (SpaceDim == 1)
                {
                  thickness = Cbase+sin(2.0*Pi*x[0]/a_domainSize[0]);
                }
              else if (SpaceDim == 2)
                {
                  thickness = Cbase+sin(2.0*Pi*x[0]/a_domainSize[0])
                    *sin(2.0*Pi*x[1]/a_domainSize[1]);
                }
              else if (SpaceDim == 3)
                {
                  thickness = D_TERM(Cbase,+sin(4*Pi*x[1]/a_domainSize[1]),
                                      *sin(4.0*Pi*x[2]/a_domainSize[2]));
                }
              H(iv,0) = a_thicknessScale*thickness;
              
            }
        }
      else if (a_thicknessType == sinusoidalHx)
        {
          BoxIterator bit(H.box());
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              RealVect x(iv);
              x += 0.5*RealVect::Unit;
              x *= dx;
              Real thickness;
              if (SpaceDim == 1)
                {
                  thickness = Cbase+sin(2.0*Pi*x[0]/a_domainSize[0]);
                }
              if (SpaceDim == 2)
                {
                  thickness = Cbase+sin(2.0*Pi*x[0]/a_domainSize[0]);
                }
              else if (SpaceDim == 3)
                {
                  thickness = Cbase+sin(4*Pi*x[1]/a_domainSize[1]);
                }
              H(iv,0) = a_thicknessScale*thickness;
              
            }
        }
      else if (a_thicknessType == singleSinBump)
        {
          BoxIterator bit(H.box());
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              RealVect x(iv);
              x += 0.5*RealVect::Unit;
              x *= dx;
              Real thickness;
              if (SpaceDim == 2)
                {
                  thickness = Cbase+sin(Pi*x[0]/a_domainSize[0])
                    *sin(Pi*x[1]/a_domainSize[1]);
                }
              else if (SpaceDim == 3)
                {
                  thickness = D_TERM(Cbase,+sin(Pi*x[1]/a_domainSize[1]),
                                      *sin(Pi*x[2]/a_domainSize[2]));
                }
              H(iv,0) = a_thicknessScale*thickness;
              
            }
        }
      else if (a_thicknessType == doubleZb)
        {
          // thickness == zB, so the surface z is 2*Zb
          H.copy(zB);
        }
      else if (a_thicknessType == circle)
        {
          // center at middle of domain
          RealVect center(a_domainSize);
          center *= 0.5;
          Real circleRadSqr = (0.3*a_domainSize[1]*0.3*a_domainSize[1]);
          
          BoxIterator bit(H.box());
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              RealVect x(iv);
              x += 0.5*RealVect::Unit;
              x *= dx;
              x -= center;
              Real radSqr;
              if (SpaceDim == 2)
                {
                  radSqr = x[0]*x[0] + x[1]*x[1];
                }
              else if (SpaceDim == 3)
                {
                  radSqr = x[1]*x[1] + x[2]*x[2];
                }
              else 
                {
                  MayDay::Error("circle not defined if not in 1 or 2D");
                }
              Real localThickness = 0;
              // polynomial thickness in r
              Real A = -1.0;
              Real C = 1.0;
              radSqr /= circleRadSqr;
              if (radSqr < 1.0)
                {
                  localThickness = A*radSqr +C;
                }
              H(iv,0) = a_thicknessScale*localThickness;
            } // end loop over cells
        }
      else
        {
          MayDay::Error("bad thickness type") ;
        }
       
    } // end loop over boxes
  //a_sigmaCoords.recomputeGeometry();
}



#include "NamespaceFooter.H"

