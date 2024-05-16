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

#include "LayoutData.H"
#include "LevelSigmaCS.H"
#include "RealVect.H"
#include "CONSTANTS.H"

#ifndef _DEFINELEVELSIGMACS_H_
#define _DEFINELEVELSIGMACS_H_

#include "NamespaceHeader.H"

enum basalTypes{constantZb=0,
                sinusoidalZb,
                sinusoidalYZb,
                xInclineZb,
                yInclineZb,
		pattynAZb,
		pattynBZb,
                num_basaltype};

// thickness distribution enum
enum thicknessEnum{constantThickness = 0,
                   constantThickness1km,
                   constantZs1km,
                   doubleZb, 
                   sinusoidalH,
                   sinusoidalHx,
                   singleSinBump,
                   circle,
		   pattynAH,
		   pattynBH,
                   num_thicknesstype};


// utility function for defining coordinate systems for testing
void defineLevelSigmaCS(LevelSigmaCS& a_sigmaCoords,
                        const RealVect& a_domainSize,
                        int a_thicknesstype,
                        int a_basalType,
			const RealVect& a_basalSlope,
                        Real a_thicknessScale = 1.0);


#include "NamespaceFooter.H"

#endif
