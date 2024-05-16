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
// amrtoplot.cpp
// read in data from an AMR hierarchy, and write to files.
// if the input file is name.hdf5, variable phi0 is written to 
// phi0.<level>.name.txt
// which has space separated lines like
//  x [y [z]]  var
// *only* the valid regions are written. 
//===========================================================================

#include <iostream>
#include <fstream>
#include "ParmParse.H"
#include "AMRIO.H"
#include "LoadBalance.H"

// if FINITE_VOLUME is defined, then plot finite-volume view of data
// (constant over cells)
// otherwise, only plot cell-center values
//#define FINITE_VOLUME


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
   
    if(argc < 2) 
      { 
	std::cerr << " usage: " << argv[0] << " <input_file> " << std::endl; 
	exit(1); 
      }

    char* in_file = argv[1];
    //char* out_file = argv[2];
    
    string file_root(in_file);
    
    // strip off "hdf5" suffix and replace with "txt"
    int position = file_root.find("hdf5");
    file_root.replace(position, 4, "txt");

    Vector<int> comp;
    

    std::string exec(argv[0]);
    
    if (exec.find("amrtoplot") < exec.size())
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

        pout() << "File = " << in_file << ": Time = " << time << endl;
        pout() << "writing: ";
        for (int var=0; var<names.size(); var++)
          {
            std::string compName = names[var];

            if (var > 0 ) pout() << "         ";
            pout() << compName << endl;

            for (int lev = 0; lev < numLevels; lev++)
              {
                char iter_str[100];
                sprintf(iter_str,"%s.%d.%s",compName.c_str(), lev,
                        file_root.c_str());
                
                ofstream os(iter_str, ios::out);
                
                for (DataIterator dit(grids[lev]); dit.ok(); ++dit)
                  {
                    const FArrayBox& fab = (*data[lev])[dit];
                    
                    for (BoxIterator bit(grids[lev][dit]);bit.ok();++bit)
                      {
                        const IntVect& iv = bit();
                        if ( (*mask[lev])[dit](iv) != 0)
                          {   
#ifdef FINITE_VOLUME                            
                            // low corner
                            for (int dir = 0; dir < SpaceDim; dir++)
                              os << (Real(iv[dir]))*dx[lev] << " ";
                            os << fab(iv,var) << std::endl;
#endif

                            // cell center
                            for (int dir = 0; dir < SpaceDim; dir++)
                              os << (Real(iv[dir]) + 0.5)*dx[lev] << " ";
                            os << fab(iv,var) << std::endl;

#ifdef FINITE_VOLUME                            
                            // high corner
                            for (int dir = 0; dir < SpaceDim; dir++)
                              os << (Real(iv[dir]) + 1.0)*dx[lev] << " ";
                            os << fab(iv,var) << std::endl;
#endif
			  }
                      }
                    os << std::endl;
		  }
                os.close();
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
