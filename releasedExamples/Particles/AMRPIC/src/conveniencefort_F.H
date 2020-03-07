#ifndef _CONVENIENCEFORT_F_H_
#define _CONVENIENCEFORT_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDINTERPOLATEINTIME 
#define GUARDINTERPOLATEINTIME 
// Prototype for Fortran procedure INTERPOLATEINTIME ...
//
void FORTRAN_NAME( INTERPOLATEINTIME ,interpolateintime )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(phiNew)
      ,CHFp_CONST_FRA(phiOld)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_INT(srcComp)
      ,CHFp_CONST_INT(dstComp)
      ,CHFp_CONST_INT(sizeComp)
      ,CHFp_BOX(box) );

#define FORT_INTERPOLATEINTIME FORTRAN_NAME( inlineINTERPOLATEINTIME, inlineINTERPOLATEINTIME)
#define FORTNT_INTERPOLATEINTIME FORTRAN_NAME( INTERPOLATEINTIME, interpolateintime)

inline void FORTRAN_NAME(inlineINTERPOLATEINTIME, inlineINTERPOLATEINTIME)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(phiNew)
      ,CHFp_CONST_FRA(phiOld)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_INT(srcComp)
      ,CHFp_CONST_INT(dstComp)
      ,CHFp_CONST_INT(sizeComp)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_INTERPOLATEINTIME");
 FORTRAN_NAME( INTERPOLATEINTIME ,interpolateintime )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(phiNew)
      ,CHFt_CONST_FRA(phiOld)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_INT(srcComp)
      ,CHFt_CONST_INT(dstComp)
      ,CHFt_CONST_INT(sizeComp)
      ,CHFt_BOX(box) );
}
#endif  // GUARDINTERPOLATEINTIME 

}

#endif
