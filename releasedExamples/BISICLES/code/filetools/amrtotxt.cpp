 #ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

//===========================================================================
// amrtotxt.cpp
// read in data from an AMR hierarchy, and write to pout(), with
// space separated lines like
// amr level [j [k]] x [y [z]] time var1 var2 var3...
// *only* the valid regions are written. amr is a literal.
//===========================================================================

#include <iostream>
#include "ParmParse.H"
#include "AMRIO.H"
#include "LoadBalance.H"

int main(int argc, char* argv[]) {

#ifdef CH_MPI
#warning "Parallel amrtotxt might well not work";
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
   
    if(argc < 3) 
      { 
	std::cerr << " usage: " << argv[0] << " <input_file>  <var 1 [,var 2, ...] " << std::endl; 
	exit(1); 
      }

    char* in_file = argv[1];
    //char* out_file = argv[2];

    std::string exec(argv[0]);
    
    Vector<std::string> var;
    for (int i = 2; i < argc; i++)
      {
	var.push_back(std::string(argv[i]));
      } 

    if (exec.find("amrtotxt") < exec.size())
      {
	//pout() << "converting AMR file " << in_file << " to text file " << out_file << std::endl;
	Vector<std::string> names;
	Vector<LevelData<FArrayBox>* > data;
	Vector<DisjointBoxLayout> grids;
	Vector<int> ratio;
	Real crseDx = 0.0, dt = 0.0, time = 0.0;
	Box domBox;
	int numLevels;
	int status = ReadAMRHierarchyHDF5
	  (in_file,grids,data,names,domBox,crseDx,dt,time,
	   ratio,numLevels);
	if (status != 0)
	{
	  MayDay::Error("failed to read AMR hierarchy");
	}
	Vector<Real> dx(numLevels,crseDx);
	for (int lev = 1; lev < numLevels; lev++)
	  {
	    dx[lev] = dx[lev-1] / Real(ratio[lev-1]);
	  }
	
	//build map from column number (var index) to AMR component
	Vector<int> comp(var.size());
	for (int i = 0; i < var.size(); i++)
	  {
	    for (int j = 0; j < names.size(); j++)
	      {
		if (names[j] == var[i])
		  {
		    comp[i] = j;
		  }
	      }
	  }

	//build valid region masks
	Vector<LevelData<BaseFab<int> >* > mask(numLevels,NULL);
	for (int lev = numLevels - 1; lev >= 0; lev--)
	  {
	    mask[lev] = new LevelData<BaseFab<int> >(grids[lev],1,IntVect::Zero);
	    for (DataIterator dit(grids[lev]); dit.ok(); ++dit)
	      {
		(*mask[lev])[dit].setVal(1);

		if (lev < numLevels - 1)
		  {
		    for (DataIterator fit(grids[lev+1]); fit.ok(); ++fit)
		      {
			Box covered = grids[lev+1][fit];
			covered.coarsen(ratio[lev]);
			covered &= grids[lev][dit];
			if (!covered.isEmpty())
			  {
			    (*mask[lev])[dit].setVal(0,covered,0,1);
			  }
		      }
		  }
	      }
	  
	  }

	for (int lev = 0; lev < numLevels; lev++)
	  {
	    for (DataIterator dit(grids[lev]); dit.ok(); ++dit)
	      {
		const FArrayBox& fab = (*data[lev])[dit];

		for (BoxIterator bit(grids[lev][dit]);bit.ok();++bit)
		  {
		    const IntVect& iv = bit();
		    if ( (*mask[lev])[dit](iv) != 0)
		      {
			
			pout() << "amr "  << lev << " "; 
			for (int dir = 0; dir < SpaceDim; dir++)
			  pout() << iv[dir] << " ";
			for (int dir = 0; dir < SpaceDim; dir++)
			  pout() << (Real(iv[dir]) + 0.5)*dx[lev] << " ";
			pout() << time << " ";
			for (int ic = 0; ic < comp.size(); ic++)
			  {
			    pout() << fab(iv,comp[ic]) << " ";
			  }


			pout() << std::endl;
		      }
		  }
	      }
	  }	


	for (int lev = 0; lev < numLevels; lev++)
	  {
	    if (mask[lev] != NULL)
	      {
		delete mask[lev];
		mask[lev] = NULL;
	      }

            if (data[lev] != NULL)
              {
                delete data[lev];
                data[lev] = NULL;
              }
	  }


      }

  }  // end nested scope
  CH_TIMER_REPORT();
  
#ifdef CH_MPI
  MPI_Finalize();
#endif
  
  return 0;
}
