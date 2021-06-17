#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "LevelData.H"
#include "FluxBox.H"
#include "BoxIterator.H"

#include "PetscSolver.H"
#include "EBCellFAB.H"
#include "EBPoissonPetscSolver.H"

#include "NamespaceHeader.H"



#ifdef CH_USE_PETSC
//
void 
EBPoissonPetscSolver::
solve( LevelData<EBCellFAB>& a_phi, const LevelData<EBCellFAB>& a_rhs )
{
  CH_assert(m_isDefined);
  if(!m_homogeneous)
    {

      const DisjointBoxLayout&  dbl = a_phi.disjointBoxLayout();
      const ProblemDomain &  dom = dbl.physDomain();
      //  pout() << "ebpetscsolver solving for domain = " << dom.domainBox() << " and dbl = "  << dbl << endl;
      pout() << "ebpetscsolver solving for domain = " << dom.domainBox() << endl;
      LevelData<EBCellFAB> rhshomog;
      m_op->create(   rhshomog, a_rhs);
      m_op->setToZero(rhshomog);
      m_op->incr(     rhshomog, a_rhs  ,  1.0);

      // make sure covered cells are set to 0.0 because PetscSolver assumes all cells are defined
      EBLevelDataOps::setCoveredVal(rhshomog, 0.0);
      //call the base solver
      PetscSolver<LevelData<EBCellFAB> >::solve(a_phi, rhshomog);
    }
  else
    {
      LevelData<EBCellFAB>& rhscast = const_cast<LevelData<EBCellFAB> &>(a_rhs);
      EBLevelDataOps::setCoveredVal(rhscast, 0.0);
      PetscSolver<LevelData<EBCellFAB> >::solve(a_phi, rhscast);
    }
}
///
int 
EBPoissonPetscSolver::
formMatrix( Mat a_mat, 
            const LevelData<EBCellFAB> *a_dummy,
            PetscInt dummy_my0, PetscInt dummy_nloc,
            PetscInt *dummy_d, PetscInt *dummy_o )
{
  CH_assert(m_isDefined);

  PetscInt  ierr;
  Real nil = 1;
  DisjointBoxLayout dbl = m_gids.disjointBoxLayout();
  const EBLevelGrid& eblg = m_op->getEBLG();
  const EBISLayout & ebisl = eblg.getEBISL();
  for ( DataIterator dit = m_gids.dataIterator() ; dit.ok() ; ++dit )
    {
      const Box &box = dbl[dit];
      for(BoxIterator bit(box); bit.ok(); ++bit)
        {
          int irow  = m_gids[dit()](bit(), 0); //this will break with multicells
          if(ebisl[dit()].isCovered(bit()))
            {
              //have to add something everywhere because petscsolver is not too smart
              PetscInt irowpet = irow;
              ierr = MatSetValues(a_mat,1,&irowpet,1,&irowpet,&nil,INSERT_VALUES); CHKERRQ(ierr); 
            }
          else
            {
              Vector<VolIndex> vofs = ebisl[dit()].getVoFs(bit());
              for(int qvof = 0; qvof < vofs.size(); qvof++)
                {
                  const VolIndex& vof = vofs[qvof];
                  VoFStencil stenc;
                  //get the stencil for this point

                  m_op->getVoFStencil(stenc, vof, dit());
                  for(int ivof = 0; ivof < stenc.size(); ivof++)
                    {
                      const VolIndex& stenvof = stenc.vof(ivof);
                      const Real    & stenwei = stenc.weight(ivof);
                      int icol = m_gids[dit()](stenvof.gridIndex(), 0);  //this will also break with multicells
                      PetscInt irowpet = irow;
                      PetscInt icolpet = icol;
                      ierr = MatSetValues(a_mat,1,&irowpet,1,&icolpet,&stenwei,INSERT_VALUES);
                      CHKERRQ(ierr);
                    }//end loop over stencil
                } //end loop over vofs in cell
            }// end if(notcovered)
        } //box iterator over box because petscsolver is not too smart
    }//data iterator

  ierr = MatAssemblyBegin(a_mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(a_mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  return 0;

//  CH_assert(m_isDefined);
//
//  m_op->getPetscMatrix(a_mat);

  return 0;
}
#endif

#include "NamespaceFooter.H"
