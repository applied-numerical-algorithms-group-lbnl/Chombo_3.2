#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif
#ifdef CH_USE_PETSC

#include "PetscAMRSolver.H"
#include "JFNKSolver.H"
#include "NamespaceHeader.H"
#include "petscviewer.h"

//
PetscAMRSolver::PetscAMRSolver(int a_verb /* =0 */) : m_op_mfree(0),m_mfree_homogeneous(true),m_verbose(a_verb)
{
}

#undef __FUNCT__
#define __FUNCT__ "apply_mfree"
PetscErrorCode PetscAMRSolver::apply_mfree(Mat A, Vec x, Vec f)
{
  CH_TIME("PetscAMRSolver::apply_mfree");
  //Whenever I write any PETSc code, I look forward to pulling classes back from the void.
  PetscFunctionBeginUser;
  void *ctx;
  MatShellGetContext(A, &ctx);
  PetscAMRSolver *tthis = (PetscAMRSolver*)ctx;
  tthis->m_petscCompMat.putPetscInChombo(x, tthis->m_phi_mfree);
  tthis->m_op_mfree->applyOp(tthis->m_Lphi_mfree,tthis->m_phi_mfree,tthis->m_mfree_homogeneous);
  tthis->m_petscCompMat.putChomboInPetsc(tthis->m_Lphi_mfree,f);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "solve_mfree"
PetscErrorCode
PetscAMRSolver::solve_mfree( Vector<LevelData<FArrayBox>*>& a_phi, 
			     const Vector<LevelData<FArrayBox>*>& a_rhs, 
			     LinearizedVTOp *a_op )
{
  CH_TIME("PetscSolver::solve_mfree");
#ifdef CH_MPI
  MPI_Comm wcomm = Chombo_MPI::comm;
#else
  MPI_Comm wcomm = PETSC_COMM_SELF;
#endif
  PetscErrorCode ierr;
  Vec  x, b;      /* approx solution, RHS */
  Mat  A;         /* linear system matrix */
  KSP  ksp;       /* linear solver context */
  PetscFunctionBeginUser;

  m_petscCompMat.setVerbose(m_verbose);
  ierr = m_petscCompMat.createMatrix();CHKERRQ(ierr); 
  A = m_petscCompMat.getMatrix();
  ierr = MatViewFromOptions(A,NULL,"-mat_view");CHKERRQ(ierr);

  //create an operator matrix shell with same dimensions as m_mat
  PetscInt m, n, M, N;
  ierr = MatGetSize(A, &M, &N);CHKERRQ(ierr); CH_assert(M == N);
  ierr = MatGetLocalSize(A, &m, &n);CHKERRQ(ierr);
  Mat L; 
  ierr = MatCreateShell(wcomm,m,n,N,N,(void *)this,&L);CHKERRQ(ierr);
  ierr = MatShellSetOperation(L,MATOP_MULT,(void(*)(void))apply_mfree);
  m_op_mfree = a_op;
  //allocate space for a vector and a matrix-vector product in Chombo-land
  a_op->create( m_phi_mfree , a_phi);
  a_op->create( m_Lphi_mfree , a_rhs);

  ierr = MatCreateVecs(A,&x,&b); CHKERRQ(ierr);
  ierr = m_petscCompMat.putChomboInPetsc(a_rhs,b);CHKERRQ(ierr);
  ierr = KSPCreate(wcomm, &ksp); CHKERRQ(ierr);
#if PETSC_VERSION_LT(3,5,0)
  ierr = KSPSetOperators(ksp, L, A, SAME_NONZERO_PATTERN);CHKERRQ(ierr);
#else
  ierr = KSPSetOperators(ksp, L, A);CHKERRQ(ierr);
#endif
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

  PetscBool ism = PETSC_FALSE;

#if PETSC_VERSION_GE(3,7,0)
  PetscOptionsGetBool(PETSC_NULL,PETSC_NULL,"-ksp_monitor",&ism,PETSC_NULL);
#else
 PetscOptionsGetBool(PETSC_NULL,"-ksp_monitor",&ism,PETSC_NULL);
#endif

  if(ism)
    {
      ierr = KSPMonitorSet(ksp,ksp_monitor_pout,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
    }

  setCoords(ksp, a_rhs);
  ierr = KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);CHKERRQ(ierr);
  ierr = VecViewFromOptions(b,NULL,"-vec_view");CHKERRQ(ierr);
  ierr = KSPSolve(ksp, b, x); CHKERRQ(ierr);

  ierr = m_petscCompMat.putPetscInChombo(x,a_phi);CHKERRQ(ierr);
  
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = MatDestroy(&L);CHKERRQ(ierr);

  if(false){  
    int nlevel = a_phi.size();
    string fname = "PetscAMRSolver_solve_mfree.";      
    char suffix[30];
    sprintf(suffix, "%dd.%d.hdf5",SpaceDim,nlevel);
    fname += suffix;
    plot(fname,a_phi);
  }

  //clean up 
  a_op->clear(m_phi_mfree);
  a_op->clear(m_Lphi_mfree);
  PetscFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "solve"
PetscErrorCode
PetscAMRSolver::solve( Vector<LevelData<FArrayBox>*>& a_phi, 
		       const Vector<LevelData<FArrayBox>*>& a_rhs )
{
  CH_TIME("PetscAMRSolver::solve");
#ifdef CH_MPI
  MPI_Comm wcomm = Chombo_MPI::comm;
#else
  MPI_Comm wcomm = PETSC_COMM_SELF;
#endif
  PetscErrorCode ierr;
  Vec  x, b;      /* approx solution, RHS */
  Mat  A;         /* linear system matrix */
  KSP  ksp;       /* linear solver context */
  PetscFunctionBeginUser;

  ierr = m_petscCompMat.createMatrix(); CHKERRQ(ierr); 
  A = m_petscCompMat.getMatrix();
  
  //create an operator matrix shell with same dimensions as m_mat
  ierr = MatCreateVecs(A,&x,&b); CHKERRQ(ierr);
  ierr = m_petscCompMat.putChomboInPetsc(a_rhs,b); CHKERRQ(ierr);
  ierr = KSPCreate(wcomm, &ksp); CHKERRQ(ierr);
#if PETSC_VERSION_LT(3,5,0)
  ierr = KSPSetOperators(ksp, A, A, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
#else
  ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
#endif
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

  PetscBool ism = PETSC_FALSE;
#if PETSC_VERSION_GE(3,7,0)
  PetscOptionsGetBool(PETSC_NULL,PETSC_NULL,"-ksp_monitor",&ism,PETSC_NULL);
#else
  PetscOptionsGetBool(PETSC_NULL,"-ksp_monitor",&ism,PETSC_NULL);
#endif
  if(ism)
    {
      ierr = KSPMonitorSet(ksp,ksp_monitor_pout,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
    }

  setCoords(ksp, a_rhs);
  ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); CHKERRQ(ierr);
  ierr = m_petscCompMat.putChomboInPetsc(a_phi,x); CHKERRQ(ierr);
  ierr = KSPSolve(ksp, b, x); CHKERRQ(ierr);

  ierr = m_petscCompMat.putPetscInChombo(x,a_phi); CHKERRQ(ierr);

  if(false){  
    int nlevel = a_phi.size();
    string fname = "PetscAMRSolver_solve.";      
    char suffix[30];
    sprintf(suffix, "%dd.%d.hdf5",SpaceDim,nlevel);
    fname += suffix;
    plot(fname,a_phi);
  }

  ierr = VecDestroy(&x); CHKERRQ(ierr);
  ierr = VecDestroy(&b); CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

void PetscAMRSolver::plot(const string a_fname, const Vector<LevelData<FArrayBox>*>& a_phi)
{
  int nlevel = a_phi.size();
  Vector<string> varNames(2);
  varNames[0] = "vel-x";
  varNames[1] = "vel-y";
  Real bogusVal = 1.0;
  Vector<DisjointBoxLayout> grids(nlevel);
  Vector<int> refrat(nlevel-1);
  ProblemDomain dom = a_phi[0]->getBoxes().physDomain();
  for (int i=0;i<nlevel;i++) 
    {
      grids[i] = a_phi[i]->getBoxes();
      if(i!=0) refrat[i-1] = 2;	
    }
#ifdef CH_USE_HDF5
  WriteAMRHierarchyHDF5( a_fname,
			 grids,
			 a_phi,
			 varNames,
			 dom.domainBox(),
			 1.0,
			 bogusVal,
			 bogusVal,
			 refrat,
			 nlevel);
#endif
}

#undef __FUNCT__
#define __FUNCT__ "setCoords"
PetscErrorCode
PetscAMRSolver::setCoords(KSP a_ksp, const Vector<LevelData<FArrayBox>*>& a_rhs)
{
  PC pc; PetscInt gid,sz,bs,n,m,nGrids=a_rhs.size();
  const int my0eq=CH_SPACEDIM*m_petscCompMat.m_gid0;
#if PETSC_VERSION_LT(3,4,0) & PETSC_VERSION_RELEASE
  const PCType type;
#else
  PCType type;
#endif	
  PetscErrorCode ierr;
  Mat A;
  PetscFunctionBeginUser;

  A = m_petscCompMat.getMatrix();  
  ierr = KSPGetPC( a_ksp, &pc );     CHKERRQ(ierr);
  ierr = PCGetType( pc, &type );    CHKERRQ(ierr);
  ierr = MatGetBlockSize( A, &bs );               CHKERRQ( ierr );
  if ( strcmp(type,PCGAMG) == 0 && bs > 1 )
    {
      PetscReal    *coords;      
      ierr = MatGetLocalSize( A, &m, &n );  CHKERRQ(ierr);
      sz = CH_SPACEDIM*(m/bs);
      ierr = PetscMalloc( (sz+1)*sizeof(PetscReal), &coords ); CHKERRQ(ierr);
      for (int ilev=nGrids-1;ilev>=0;ilev--)
	{
	  const DisjointBoxLayout &dbl = a_rhs[ilev]->getBoxes();
	  for ( DataIterator dit(dbl) ; dit.ok() ; ++dit )
	    {
	      const Box &box = dbl[dit];
	      const BaseFab<PetscInt>& gidfab = (*m_petscCompMat.m_GIDs[ilev])[dit];
	      BoxIterator bit(box);
	      for (bit.begin(); bit.ok(); bit.next())
		{
		  IntVect iv = bit(); // coordinate in any scaled, shifted, rotated frame.		  
		  if( (gid=gidfab(iv,0)) >= 0)
		    {
		      for (int i=CH_SPACEDIM*gid-my0eq,n=0;n<CH_SPACEDIM;n++,i++) coords[i] = (PetscReal)iv[n];
		    }
		}
	    }
	}
      ierr = PCSetCoordinates( pc, CH_SPACEDIM, sz/CH_SPACEDIM, coords ); CHKERRQ(ierr);
      ierr = PetscFree( coords );  CHKERRQ(ierr);
    } 
  PetscFunctionReturn(0);
}

#include "NamespaceFooter.H"
#endif
