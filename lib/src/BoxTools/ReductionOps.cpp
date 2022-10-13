#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ReductionOps.H"
#include "BoxIterator.H"
#include "CellToEdge.H"
#include "NamespaceHeader.H"


void computeFaceReductionWeights(LevelData<FluxBox>& a_weights)
{
  // Before calling op in FaceSumOp, one must
  //a)Construct a LevelData<FluxBox> weights using the same layout as the data
  //b)Call  computeFaceReductionWeights
  //c)Pointwise multiply the data by the weights

//pout()<<"Beginning computeFaceReductionWeights"<<endl;
  const DisjointBoxLayout grids = a_weights.getBoxes();
  IntVect ghostVect = IntVect::Unit;
  // cell-centered temp data
  LevelData<FArrayBox> ccTemp(grids, 1, ghostVect);
  
  // first, set ghosts to 1 and interiors to zero
  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& validBox = grids[dit];
      FArrayBox& thisTemp = ccTemp[dit];
      
      // first, set *everything* to 1
      thisTemp.setVal(1.0);
      // now set valid cells to 0
      thisTemp.setVal(0.0, validBox, 0);
    }
  
  // now call exchange -- this will set ghosts at shared faces to be zero
  ccTemp.exchange();

  // now set interiors back to 1
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& validBox = grids[dit];
      ccTemp[dit].setVal(1.0, validBox, 0);
    }
   
  // now all interior and ghost cells will be 1, except for ghost cells at shared faces.
  // averaging cell->faces should result in 1 on all faces except for shared faces, which will be 0.5
   
  CellToEdge(ccTemp, a_weights);
   
       
  for (dit.begin(); dit.ok(); ++dit)
    {
      FluxBox& weightFluxBox = a_weights[dit()];
      for (int idir = 0; idir <SpaceDim; ++idir)
        {
          FArrayBox& weightFab = weightFluxBox[idir]; 
       // pout()<<"box = "<<weightFab.box()<<endl;
       // pout()<<"idir = "<<idir<<endl;
          for (BoxIterator bit(weightFab.box()); bit.ok(); ++bit)
            {
              IntVect iv = bit();
//pout()<<"weightFab("<<iv<<","<<0<<") = "<<weightFab(iv,0)<<endl;
            }
        }
    }
//pout()<<"Finishing computeFaceReductionWeights"<<endl;
}

SumOp::SumOp():scale(1.0)
{
}

SumOp::SumOp(int a_summingDir ):scale(1.0)
{
  m_summingDir.resize(1);
  m_summingDir[0] = a_summingDir;
}

SumOp::SumOp(const Vector<int>& a_summingDir ):scale(1.0)
{
  m_summingDir = a_summingDir;
 
}

void
SumOp::linearOut(const FArrayBox& arg, void* buf, const Box& R,
                 const Interval& comps) const
{
  Real* buffer = (Real*)buf;

  Box flattenedBox(R);
  for (int n=0; n<m_summingDir.size(); n++)
    {
      flattenedBox.setBig(m_summingDir[n],R.smallEnd(m_summingDir[n]));
    }

  // single summing direction
  if (m_summingDir.size() == 1)
    {
      int transverseLo = R.smallEnd(m_summingDir[0]);
      int transverseHi = R.bigEnd(m_summingDir[0]);
      // don't apply scale here, since we'll do it in the linearIn side of
      // things.
      ForAllXCBNN(Real, arg, flattenedBox, comps.begin(), comps.size())
        {
          // this gets rid of the "unused variable" warning
          // while we also initialize the buffer
          //*buffer = 0;
          *buffer = 0*argR;

          IntVect iv(D_DECL6(iR, jR, kR,
                             _iv[3], _iv[4], _iv[5]));
          for (int itrans=transverseLo; itrans<=transverseHi; itrans++)
            {
              iv[m_summingDir[0]] = itrans;
              *buffer += arg(iv,_n);
            }
          buffer++;
        } EndFor
    }
  else
    {
      // multiple summing directions -- need to do this a bit differently
      Box transverseBox(IntVect::Zero, IntVect::Zero);
      for (int n=0; n<m_summingDir.size(); n++)
        {
          transverseBox.setBig(m_summingDir[n],
                               (R.bigEnd(m_summingDir[n])-R.smallEnd(m_summingDir[n])));
          //transverseBox.setSmall(m_summingDir[n],R.smallEnd(m_summingDir[n]));
        }
      // don't apply scale here, since we'll do it in the linearIn side of
      // things.
      ForAllXCBNN(Real, arg, flattenedBox, comps.begin(), comps.size())
        {
          // this gets rid of the "unused variable" warning
          // while we also initialize the buffer
          //*buffer = 0;
          *buffer = 0*argR;

          IntVect iv(D_DECL6(iR, jR, kR,
                             _iv[3], _iv[4], _iv[5]));

          BoxIterator transverseIt(transverseBox);
          for (transverseIt.begin(); transverseIt.ok(); ++transverseIt)
                {
                  IntVect shift = transverseIt();
                  IntVect srcLoc = iv + shift;
                  *buffer += arg(srcLoc,_n);
                } // end loop over transverse directions
              buffer++;
        } EndFor
    }
  
}

void
SumOp::linearIn(FArrayBox& arg,  void* buf, const Box& R,
                const Interval& comps) const
{

  Real* buffer = (Real*)buf;
  if (scale != 1.0)
    {
      ForAllXBNNnoindx(Real, arg, R, comps.begin(), comps.size())
        {
          argR+=*buffer * scale;
          buffer++;
        } EndFor
            }
  else
    {
      ForAllXBNNnoindx(Real, arg, R, comps.begin(), comps.size())
        {
          argR+=*buffer;
          buffer++;
        } EndFor
            }
}

void
SumOp::op(FArrayBox& dest,
          const Box& RegionFrom,
          const Interval& Cdest,
          const Box& RegionTo,
          const FArrayBox& src,
          const Interval& Csrc) const
{
//int numComp = Cdest.size();
  CH_assert(Cdest.size() == Csrc.size());
  CH_assert(dest.nComp() > Cdest.end());
  CH_assert(src.nComp() > Csrc.end());

  int offset = Cdest.begin() - Csrc.begin();

  // start by doing this with a boxIterator -- can make this more efficient
  // later if necessary
  if (m_summingDir.size() == 1)
    {
      int toBoxLo = RegionTo.smallEnd(m_summingDir[0]);
      int toBoxHi = RegionTo.bigEnd(m_summingDir[0]);

      BoxIterator fromBit(RegionFrom);
      for (fromBit.begin(); fromBit.ok(); ++fromBit)
        {
          IntVect fromIV = fromBit();
          //int summingIndex = fromIV[m_summingDir];

          if (toBoxLo == toBoxHi)
            {
              IntVect toIV = fromIV;
              toIV[m_summingDir[0]] = toBoxLo;
              for (int comp=Csrc.begin(); comp<=Csrc.end(); comp++)
                {
                  dest(toIV, comp+offset) += scale*src(fromIV, comp);
                } // end loop over components
            }
          else
            {
              // this is in case the toBox is more than one cell wide in
              // the summing direction.  This is most likely unnecessary,
              // to be honest.
              Box toBox(fromIV, fromIV);
              toBox.setSmall(m_summingDir[0], toBoxLo);
              toBox.setBig(m_summingDir[0], toBoxHi);

              BoxIterator toBit(toBox);
              for (toBit.begin(); toBit.ok(); ++toBit)
                {
                  IntVect toIV = toBit();

                  for (int comp=Csrc.begin(); comp<=Csrc.end(); comp++)
                    {
                      dest(toIV, comp+offset) += scale*src(fromIV, comp);
                    } // end loop over components
                } // end loop over toBox
            }
        } // end loop over from cells
    } // end if single summing direction
  else
    {
      int width = 1;
      for (int n=0; n<m_summingDir.size(); n++)
        {
          int toBoxLo = RegionTo.smallEnd(m_summingDir[n]);
          int toBoxHi = RegionTo.bigEnd(m_summingDir[n]);
          width *= (toBoxHi - toBoxLo+1);
        }

      BoxIterator fromBit(RegionFrom);
      for (fromBit.begin(); fromBit.ok(); ++fromBit)
        {
          IntVect fromIV = fromBit();
          //int summingIndex = fromIV[m_summingDir];

          if (width == 1)
            {
              IntVect toIV = fromIV;
              for (int n=0; n<m_summingDir.size(); n++)
                {
                  toIV[m_summingDir[n]] = RegionTo.smallEnd(m_summingDir[n]);
                }
              for (int comp=Csrc.begin(); comp<=Csrc.end(); comp++)
                {
                  dest(toIV, comp+offset) += scale*src(fromIV, comp);
                } // end loop over components
            }
          else
            {
              // this is in case the toBox is more than one cell wide in
              // the summing direction.  This is most likely unnecessary,
              // to be honest.
              Box toBox(fromIV, fromIV);
              for (int n=0; n<m_summingDir.size(); n++)
                {
                  toBox.setSmall(m_summingDir[n], RegionTo.smallEnd(m_summingDir[n]));
                  toBox.setBig(m_summingDir[n], RegionTo.bigEnd(m_summingDir[n]));
                }

              BoxIterator toBit(toBox);
              for (toBit.begin(); toBit.ok(); ++toBit)
                {
                  IntVect toIV = toBit();

                  for (int comp=Csrc.begin(); comp<=Csrc.end(); comp++)
                    {
                      dest(toIV, comp+offset) += scale*src(fromIV, comp);
                    } // end loop over components
                } // end loop over toBox
            } // end if toRegion is more than one cell "wide"
        } // end loop over from cells
    } // end if more than one summing direction
}

// FaceSumOp

FaceSumOp::FaceSumOp():m_scale(1.0), m_projectionDirections(-1*IntVect::Unit)
{
}

FaceSumOp::FaceSumOp(const int& a_summingDir ):m_scale(1.0)
{
  m_projectionDirections = BASISV(a_summingDir);
}

FaceSumOp::FaceSumOp(const IntVect& a_summingDir ):m_scale(1.0)
{
  m_projectionDirections = a_summingDir;
}

void FaceSumOp::setReductionDir(const IntVect& a_sumDir)
{
  m_projectionDirections = a_sumDir;
}

void FaceSumOp::setScale(const Real& a_scale)
{
  m_scale = a_scale;
}

void FaceSumOp::linearIn(FluxBox       & a_arg  ,  
                         void*           a_buf  , 
                         const Box     & a_R    ,
                         const Interval& a_comps) const
{
  //pout()<<"Beginning LinearIn"<<endl;
  //pout()<<endl;
    
  //construct a temporary FluxBox, which will be used to accumulate the values into a_arg
  FluxBox temp(a_R, a_arg.nComp());
  temp.setVal(99999999);
  
  temp.linearIn(a_buf  , 
                a_R    , 
                a_comps);
    
  for (int idir = 0; idir<SpaceDim; ++idir)
    {
      Box intersectBox = temp[idir].box();
      //intersectBox.surroundingNodes(idir);
      intersectBox    &= a_arg[idir].box();
      
      if (!intersectBox.isEmpty())
        {
          // comp loop
          for (int iComp = a_comps.begin(); iComp <= a_comps.end(); ++iComp)
            {
              //boxIterator
              for(BoxIterator bit(intersectBox); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  a_arg[idir](iv, iComp) += temp[idir](iv,iComp);
                }
            }
        }
    }

  //  pout()<<"Finishing LinearIn"<<endl;
}

void FaceSumOp::linearOut(const FluxBox       & a_arg  ,  
                          void*                 a_buf  , 
                          const Box           & a_R    ,
                          const Interval      & a_comps) const
{
  //  pout()<<"Beginning LinearOut"<<endl;
  //pout()<<endl;

  IntVect big   = a_R.bigEnd  ();
  IntVect small = a_R.smallEnd();
  for (int idir= 0; idir < SpaceDim; ++idir)
    {
      if (m_projectionDirections[idir] == 1)
        {
          big  [idir] = 0;
          small[idir] = 0;
        }
    }

  Box reducedBox(small,big);
  //  pout()<<"reducedBox = "<< reducedBox<<endl;
  //pout()<<endl;

  // use temp to call op before writing to buffer
  FluxBox temp(reducedBox, a_arg.nComp());
  
  temp.setVal(0.0);
  
  // reduce a face
  op(temp      ,
     a_R       ,
     a_comps   ,
     reducedBox,
     a_arg     ,
     a_comps   );
  
  // call flux box linear out
  temp.linearOut(a_buf      , 
                 reducedBox , 
                 a_comps   );

  //  pout()<<"finished LinearOut"<<endl;
}

void FaceSumOp::op(FluxBox       & a_dataTo    ,
                   const Box     & a_regionFrom,
                   const Interval& a_compTo    ,
                   const Box     & a_regionTo  ,
                   const FluxBox & a_dataFrom  ,
                   const Interval& a_compFrom  ) const
{
  //  pout()<<"Entering op"<<endl;

  // Before calling this function, one must
  //a)Construct a LevelData<FluxBox> weights using the same layout as a_dataFrom
  //b)Call  computeFaceReductionWeights
  //c)Pointwise multiply the a_dataFrom by the weights
  
  CH_assert(a_compTo  .size () == a_compFrom .size());
  CH_assert(a_dataTo  .nComp()  > a_compTo   .end ());
  CH_assert(a_dataFrom.nComp()  > a_compFrom .end ());
  
  // cout<<"entering FaceSumOp::op"<<endl;
  
  int offset = a_compTo.begin() - a_compFrom.begin();

  // iDir = 0 corresponds to a face that is constant in the 0 direction
  for (int iDir = 0; iDir < SpaceDim; iDir++)
    {
      Box regionToDir  (a_regionTo  );
      Box regionFromDir(a_regionFrom);
                   
      regionToDir  .surroundingNodes(iDir);
      regionFromDir.surroundingNodes(iDir);

      // src and destination for corresponding elements of two flux boxes
      FArrayBox      & dataToDir   = a_dataTo  [iDir];
      const FArrayBox& dataFromDir = a_dataFrom[iDir];

      const Box toFabBox   (dataToDir  .box());
      const Box fromFabBox (dataFromDir.box());           

      // use these for projecting and injecting
      setDestBox   (toFabBox  );
      setSourceBox (fromFabBox);

      // iterate over the source region      
      for (BoxIterator fromBit(regionFromDir); fromBit.ok(); ++fromBit)
        {
          // convention of flux box data: each cell centered intvect corresponds to the low face
          IntVect fromIv = fromBit();

          //unmodified modify destination intvects
          IntVect toIv   = project(fromIv);
          IntVect toIvHi = toIv + BASISV(iDir);

          if (m_projectionDirections[iDir] == 0)
            {
              if (toFabBox.contains(toIv) && fromFabBox.contains(fromIv))   
                {   
                  // for all components
                  for (int comp = a_compFrom.begin(); comp <= a_compFrom.end(); comp++)
                    {
                      dataToDir(toIv, comp + offset) += m_scale*dataFromDir(fromIv, comp);
                    } 
                }
            }
          else
            {
              // m_projectionDir[iDir] == 1
              if ( toFabBox.contains(toIv)  && fromFabBox.contains(fromIv) )   
                {
                  // write the (scaled) data to both faces of the one dimensional box (iDir != non-summing directions)
                  for (int comp = a_compFrom.begin(); comp <= a_compFrom.end(); comp++)
                    {
                      CH_assert(dataToDir(toIv, comp+offset) == dataToDir(toIvHi, comp+offset));
                                         
                      dataToDir(toIv  , comp+offset) += m_scale*dataFromDir(fromIv, comp);
                      dataToDir(toIvHi, comp+offset) += m_scale*dataFromDir(fromIv, comp);
                    }
                }
            } 
        } 
    }
  //  pout()<<"Leaving op"<<endl;
}

IntVect  FaceSumOp::project(const IntVect & a_iv )const
{
  IntVect retval = m_destBox.smallEnd(); 
  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      if (m_projectionDirections[idir] != 1)
        {
          retval[idir] = a_iv[idir];
        }
    }
  
  return retval;
}

Box FaceSumOp::project (const Box & a_box)const
{
  IntVect smallEnd = m_destBox.smallEnd();
  IntVect bigEnd   = m_destBox.bigEnd();

  IntVect smallEndIn = a_box.smallEnd();
  IntVect bigEndIn   = a_box.bigEnd();

  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      if (m_projectionDirections[idir] != 1)
        {
          smallEnd[idir] = smallEndIn[idir];
          bigEnd  [idir] = bigEndIn  [idir];
        }
    }
  
  Box retval(smallEnd,bigEnd);
  return retval;
}

// expand a destination box
Box FaceSumOp::inject (const Box     & a_box)const 
{
  IntVect smallEnd = m_sourceBox.smallEnd();
  IntVect bigEnd   = m_sourceBox.bigEnd();

  IntVect smallEndIn = a_box.smallEnd();
  IntVect bigEndIn   = a_box.bigEnd();

  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      if (m_projectionDirections[idir] != 1)
        {
          smallEnd[idir] = smallEndIn[idir];
          bigEnd  [idir] = bigEndIn  [idir];
        }
    }
  
  Box retval(smallEnd,bigEnd);
  return retval;
  return retval;
}

void  FaceSumOp::setDestBox(const Box & a_destBox   )const
{
   m_destBox = a_destBox;
}
void  FaceSumOp::setSourceBox    (const Box & a_sourceBox )const
{
  m_sourceBox = a_sourceBox;
}

Box  FaceSumOp::getDestinationBox() const
{
  return m_destBox;
}
Box  FaceSumOp::getSourceBox     ()const
{
  return m_sourceBox;
}

void  FaceSumOp::setProjectionDirections(const IntVect    & a_directions)const
{
  m_projectionDirections = a_directions;
}

IntVect FaceSumOp::getProjectionDirections(const IntVect    & a_directions)const
{
  IntVect retval;
  return m_projectionDirections;
}

int FaceSumOp::size(const FluxBox & a_fluxBox,
                    const Box     & a_bx     , 
                    const Interval& a_comps  ) const
{
  //pout()<<"Using FaceSumOp::size()"<<endl;
  int totalSize = 0;
  Box reducedBox = a_bx;
  for (int idir = 0; idir<SpaceDim; ++idir)
    {
      if (m_projectionDirections[idir] == 1)
        {
          int lo = reducedBox.smallEnd(idir);  
          reducedBox.setSmall(idir,lo);
          reducedBox.setBig  (idir,lo);
          // pout()<<"reducedBox = "<< reducedBox<<endl;
        }
    }
  return totalSize = a_fluxBox.size(reducedBox,
                                    a_comps   );
}

// --------------------------------------------
// spreading operator
// --------------------------------------------

SpreadingOp::SpreadingOp():scale(1.0)
{
}

SpreadingOp::SpreadingOp(int a_spreadingDir ):scale(1.0)
{
  m_spreadingDir.resize(1);
  m_spreadingDir[0] = a_spreadingDir;
}

SpreadingOp::SpreadingOp(const Vector<int>& a_spreadingDir ):scale(1.0)
{
  m_spreadingDir = a_spreadingDir;
}

void
SpreadingOp::linearIn(FArrayBox& arg,  void* buf, const Box& R,
                      const Interval& comps) const
{
  Real* buffer = (Real*)buf;
  Box flattenedBox(R);
  for (int n=0; n<m_spreadingDir.size(); n++)
    {
      flattenedBox.setBig(m_spreadingDir[n],R.smallEnd(m_spreadingDir[n]));
    }

  if (scale != 1.0)
    {
      // this copies from buf to the low-end of arg in the spreadingDir
      ForAllXBNNnoindx(Real, arg, flattenedBox, comps.begin(), comps.size())
        {
          argR=*buffer * scale;
          buffer++;
        } EndFor
            }
  else
    {
      ForAllXBNNnoindx(Real, arg, flattenedBox, comps.begin(), comps.size())
        {
          argR=*buffer;
          buffer++;
        } EndFor
            }
  // now need to copy from low-end of R to fill the rest of the box
  // do this by calling SpreadingOp::applyOp with scale = 1 (since
  // we already applied the scale above)
  this->applyOp(arg, flattenedBox, comps, R, arg, comps, 1.0);

}

void
SpreadingOp::op(FArrayBox& dest,
                const Box& RegionFrom,
                const Interval& Cdest,
                const Box& RegionTo,
                const FArrayBox& src,
                const Interval& Csrc) const
{
  applyOp(dest, RegionFrom, Cdest, RegionTo, src, Csrc, scale);
}

void
SpreadingOp::applyOp(FArrayBox& dest,
                     const Box& RegionFrom,
                     const Interval& Cdest,
                     const Box& RegionTo,
                     const FArrayBox& src,
                     const Interval& Csrc,
                     Real a_scale) const
{
//int numComp = Cdest.size();
  CH_assert(Cdest.size() == Csrc.size());
  CH_assert(dest.nComp() > Cdest.end());
  CH_assert(src.nComp() > Csrc.end());

  int offset = Cdest.begin() - Csrc.begin();

  // start by doing this with a boxIterator -- can make this more efficient
  // later if necessary
  if (m_spreadingDir.size() == 1)
    {
      int fromBoxLo = RegionFrom.smallEnd(m_spreadingDir[0]);
      int fromBoxHi = RegionFrom.bigEnd(m_spreadingDir[0]);

      // this doesn't make any sense if regionFrom is more than one wide
      CH_assert(fromBoxLo == fromBoxHi);

      BoxIterator toBit(RegionTo);
      for (toBit.begin(); toBit.ok(); ++toBit)
        {
          IntVect toIV = toBit();

          IntVect fromIV = toIV;
          fromIV[m_spreadingDir[0]] = fromBoxLo;
          for (int comp=Csrc.begin(); comp<=Csrc.end(); comp++)
            {
              dest(toIV, comp+offset) = a_scale*src(fromIV, comp);
            } // end loop over components

        } // end loop over from cells
    } // end if only one spreading direction
  else
    {
      // multiple spreading directions
      // this doesn't make any sense if fromBox is more than
      // one cell wide in any of the transverse directions
      for (int n=0; n<m_spreadingDir.size(); n++)
        {
          CH_assert(RegionFrom.size(m_spreadingDir[n]) == 1);
        }

      BoxIterator toBit(RegionTo);
      for (toBit.begin(); toBit.ok(); ++toBit)
        {
          IntVect toIV = toBit();
          IntVect fromIV = toIV;
          for (int n=0; n<m_spreadingDir.size(); n++)
            {
              fromIV[m_spreadingDir[n]] =RegionFrom.smallEnd(m_spreadingDir[n]);
            }
          for (int comp=Csrc.begin(); comp<=Csrc.end(); comp++)
            {
              dest(toIV, comp+offset) = a_scale*src(fromIV, comp);
            } // end loop over components
        } // end loop over from cells
    } // end if more than one spreading direction
}

// --------------------------------------------
// face-centered spreading operator
// --------------------------------------------

FaceSpreadingOp::FaceSpreadingOp():scale(1.0)
{
}

FaceSpreadingOp::FaceSpreadingOp(int a_spreadingDir ):scale(1.0)
{
  m_spreadingDir.resize(1);
  m_spreadingDir[0] = a_spreadingDir;
}

FaceSpreadingOp::FaceSpreadingOp(const Vector<int>& a_spreadingDir ):scale(1.0)
{
  m_spreadingDir = a_spreadingDir;
}

void
FaceSpreadingOp::linearIn(FluxBox& arg,  void* buf, const Box& R,
                          const Interval& comps) const
{
  Real* buffer = (Real*)buf;

  for (int dir=0; dir<SpaceDim; dir++)
    {
      FArrayBox& argDir = arg[dir];
      Box faceR(R);
      Box flattenedBox(faceR);
      for (int n=0; n<m_spreadingDir.size(); n++)
        {
          flattenedBox.setBig(m_spreadingDir[n],faceR.smallEnd(m_spreadingDir[n]));
        }

      flattenedBox.surroundingNodes(dir);
      if (scale != 1.0)
        {
          ForAllXBNNnoindx(Real, argDir, flattenedBox, comps.begin(), comps.size())
            {
              argDirR=*buffer * scale;
              buffer++;
            } EndFor
                }
      else
        {
          ForAllXBNNnoindx(Real, argDir, flattenedBox, comps.begin(), comps.size())
            {
              argDirR=*buffer;
              buffer++;
            } EndFor
                }
    }

  // now need to copy from low-end of R to fill the rest of the box
  // do this by calling SpreadingOp::applyOp with scale = 1 (since
  // we already applied the scale above)
  Box flattenedBox(R);
  for (int n=0; n<m_spreadingDir.size(); n++)
    {
      flattenedBox.setBig(m_spreadingDir[n],R.smallEnd(m_spreadingDir[n]));
    }
  this->applyOp(arg, flattenedBox, comps, R, arg, comps, 1.0);

}

void
FaceSpreadingOp::op(FluxBox& dest,
                    const Box& RegionFrom,
                    const Interval& Cdest,
                    const Box& RegionTo,
                    const FluxBox& src,
                    const Interval& Csrc) const
{
  applyOp(dest, RegionFrom, Cdest, RegionTo, src, Csrc, scale);
}

void
FaceSpreadingOp::applyOp(FluxBox& dest,
                         const Box& RegionFrom,
                         const Interval& Cdest,
                         const Box& RegionTo,
                         const FluxBox& src,
                         const Interval& Csrc,
                         Real a_scale) const
{
//int numComp = Cdest.size();
  CH_assert(Cdest.size() == Csrc.size());
  CH_assert(dest.nComp() > Cdest.end());
  CH_assert(src.nComp() > Csrc.end());

  int offset = Cdest.begin() - Csrc.begin();

  if (m_spreadingDir.size() == 1)
    {
      // start by doing this with a boxIterator -- can make this more efficient
      // later if necessary
      for (int dir=0; dir<SpaceDim; dir++)
        {

          FArrayBox& destDir = dest[dir];
          const FArrayBox& srcDir = src[dir];
          Box RegionFromDir(RegionFrom);
          RegionFromDir.surroundingNodes(dir);

          // since a 1-wide cell-centered box will wind up with 2 faces
          // in the normal direction, make a design decision to only spread
          // from the lower face in the spreadingDir direction
          if (dir == m_spreadingDir[0]) RegionFromDir.growHi(dir,-1);

          Box RegionToDir(RegionTo);
          RegionToDir.surroundingNodes(dir);

          int fromBoxLo = RegionFromDir.smallEnd(m_spreadingDir[0]);
          int fromBoxHi = RegionFromDir.bigEnd(m_spreadingDir[0]);

          // this doesn't make any sense if regionFrom is more than one wide
          CH_assert(fromBoxLo == fromBoxHi);

          BoxIterator toBit(RegionToDir);
          for (toBit.begin(); toBit.ok(); ++toBit)
            {
              IntVect toIV = toBit();

              IntVect fromIV = toIV;
              fromIV[m_spreadingDir[0]] = fromBoxLo;
              for (int comp=Csrc.begin(); comp<=Csrc.end(); comp++)
                {
                  destDir(toIV, comp+offset) = a_scale*srcDir(fromIV, comp);
                } // end loop over components

            } // end loop over from cells
        } // end loop over face dirs
    } // end if only one spreading direction
  else
    {
      // multiple directions
      for (int dir=0; dir<SpaceDim; dir++)
        {
          FArrayBox& destDir = dest[dir];
          const FArrayBox& srcDir = src[dir];
          Box RegionToDir(RegionTo);
          RegionToDir.surroundingNodes(dir);

          // this doesn't make any sense if the FromRegion is more than
          // one cell wide in the spreading dir
          for (int n=0; n<m_spreadingDir.size(); n++)
            {
              CH_assert(RegionFrom.size(m_spreadingDir[n]) == 1);
            }

          BoxIterator toBit(RegionToDir);
          for (toBit.begin(); toBit.ok(); ++toBit)
            {
              IntVect toIV = toBit();
              IntVect fromIV = toIV;

              for (int n=0; n<m_spreadingDir.size(); n++)
                {
                  fromIV[m_spreadingDir[n]] = RegionFrom.smallEnd(m_spreadingDir[n]);
                }
              for (int comp=Csrc.begin(); comp<=Csrc.end(); comp++)
                {
                  destDir(toIV, comp+offset) = a_scale*srcDir(fromIV, comp);
                } // end loop over components
            } // end loop over "to" cells
        } // end loop over face directions
    } // end if more than one spreading direction
}

#include "NamespaceFooter.H"
