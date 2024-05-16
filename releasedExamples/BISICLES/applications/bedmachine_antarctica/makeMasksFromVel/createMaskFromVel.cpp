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
// 
// read an AMR hierarchy, loop over subsets of its domain and set mask to
// a single mask value wherever there is ice, based on a trigger variable.
// By default, "on" will be 10000, "off" wil be zero (both defaults are
// also settable below)
//
// arguments:
//  input_file  => file containing trigger
//  output_file => file to write mask to
//  trigger_vars => variables to trigger "on" or "off" (nonzero is "on")
//
//
//===========================================================================
#include <iostream>
#include <fstream>
#include "AMRIO.H"
#include "LoadBalance.H"
#include "ParmParse.H"

#include "NamespaceHeader.H"

Real mask_on = 10000.0;
Real mask_off = 0.0;

bool verbose = true;

int main(int argc, char* argv[]) {

#ifdef CH_MPI
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

    if(argc < 4) 
      { 
	std::cerr << " usage: " << argv[0] << " <input_file> <output_file> <var1> [<var2>, etc]" << std::endl; 
	exit(0); 
      }

    char* in_file = argv[1];
    char* out_file = argv[2];

    int i=3;
    int nvars = argc-3;
    Vector<string> vars(nvars);
    for (int n=0; n<nvars; n++)
      {
        vars[n] = argv[i];
        i++;
      }

    
    // read velocity file...
    if (verbose) 
      {
        pout() << "reading " << in_file << "..." << endl;
      }


    Vector<std::string> infile_names;
    Vector<LevelData<FArrayBox>* > infile_data;
    Vector<DisjointBoxLayout> infile_grids;
    Vector<int> infile_ratio;
    Real infile_crseDx = 0.0, infile_dt = 0.0, infile_time = 0.0;
    Box infile_crseBox;
    int infile_numLevels;
    int status = ReadAMRHierarchyHDF5
      (in_file,infile_grids,infile_data,infile_names,
       infile_crseBox,infile_crseDx,infile_dt,infile_time,
       infile_ratio,infile_numLevels);

    // will eventually need a Vector of dx for this
    Vector<RealVect> infile_dx(infile_numLevels);
    {
      RealVect levelDx = infile_crseDx*RealVect::Unit;
      for (int lev=0; lev < infile_numLevels; lev++)
        {
          infile_dx[lev] = levelDx;
          if (lev < infile_ratio.size())
            {
              levelDx /= infile_ratio[lev];
            }
        }
    }

    if (status != 0)
      {
	MayDay::Error("failed to read AMR hierarchy");
      }

    // find trigger variables in input file
    Vector<int> trigger_comps(vars.size(), -1);

    for (int n=0; n<vars.size(); n++)
      {
        for (int j=0; j<infile_names.size(); j++)
          {
            if (infile_names[j] == vars[n])
              {
                trigger_comps[n] = j;
              }
          } // end loop over datafile names
        if (trigger_comps[n] == -1)
          {
            pout() << "trigger variable " << vars[n] << " not found!" << endl;
            MayDay::Error("exiting");
          }
      }
    
    // 
    Vector<LevelData<FArrayBox>* > outfile_data(infile_data.size(), NULL);

    for (int lev=0; lev<outfile_data.size(); lev++)
      {
        const DisjointBoxLayout& grids = infile_grids[lev];
        IntVect ghostVect = IntVect::Zero;
        LevelData<FArrayBox>* mask = new LevelData<FArrayBox>(grids, 1,
                                                              ghostVect);
        outfile_data[lev] = mask;

        LevelData<FArrayBox>& data = *infile_data[lev];

        DataIterator dit = grids.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
          {
            const FArrayBox& dataFab = data[dit];
            FArrayBox& maskFab = (*mask)[dit];
            maskFab.setVal(mask_off);
            BoxIterator bit(grids[dit]);
            
            // now loop over cells and tag where discriminant isn't zero
            for (int i=0; i<trigger_comps.size(); i++)
              {
                int thisVar = trigger_comps[i];
                for (bit.begin(); bit.ok(); ++bit)
                  {
                    IntVect iv = bit();
                    if (dataFab(iv,thisVar) != 0.0)
                      {
                        maskFab(iv,0) = mask_on;
                      } // end if there is ice here
                  } // end loop over cells
              } // end loop over trigger variables
          } // end loop over grids
      } // end loop over levels

    
    
    {
      pout() << "... done." << endl;
      pout() << "writing " << out_file << "..." << endl;
    }

    Vector<string> outfile_names(1, "maskVal");
    
    WriteAMRHierarchyHDF5(out_file, infile_grids, outfile_data, 
                          outfile_names, infile_crseBox, infile_crseDx, 
                          infile_dt, infile_time, infile_ratio,
                          infile_numLevels);
			  
    

    if (verbose)
      {
        pout() << "... done." << endl;
      }
  
    for (int lev=0; lev<outfile_data.size(); lev++)
      {
        if (outfile_data[lev] != NULL)
          {
            delete outfile_data[lev];
            outfile_data[lev] = NULL;
          }
      }

    for (int lev=0; lev<infile_data.size(); lev++)
      {
        if (infile_data[lev] != NULL)
          {
            delete infile_data[lev];
            infile_data[lev] = NULL;
          }
      }

  }  // end nested scope
  CH_TIMER_REPORT();

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return 0;
}

#include "NamespaceFooter.H"
