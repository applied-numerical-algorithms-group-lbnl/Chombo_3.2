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
// subregion.cpp
// read an AMR hierarchy, select a subset of its domain (defined in the
// index space of the coarsest level), and write the result to a Chombo
//  hdf5 file
//===========================================================================
#include <iostream>
#include "AMRIO.H"
#include "LoadBalance.H"
#include "NamespaceHeader.H"

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
	std::cerr << " usage: " << argv[0] << " <input_file> <output_file> region_lo_x_index region_lo_y_index region_hi_x_index region_hi_y_index" << std::endl; 
	exit(0); 
      }

    char* in_file = argv[1];
    char* out_file = argv[2];

    Vector<int> subregion_indices;
    for (int i = 3; i < 3+SpaceDim*2; i++)
      {
	subregion_indices.push_back(atoi(argv[i]));
      }  

    IntVect loVect, hiVect;
    for (int i=0; i<SpaceDim; i++)
      {
        loVect[i] = subregion_indices[i];
        hiVect[i] = subregion_indices[SpaceDim+i];
      }

    pout() << "loVect = " << loVect << endl;
    pout() << "hiVect = " << hiVect << endl;

    Box subDomain(loVect,hiVect);

    if (verbose) 
      {
        pout() << "subDomain = " << subDomain << endl;
      }


    if (verbose) 
      {
        pout() << "reading " << in_file << "..." << endl;
      }

    Vector<std::string> names;
    Vector<LevelData<FArrayBox>* > data;
    Vector<DisjointBoxLayout> grids;
    Vector<int> ratio;
    Real crseDx = 0.0, dt = 0.0, time = 0.0;
    Box crseBox;
    int numLevels;
    int status = ReadAMRHierarchyHDF5
      (in_file,grids,data,names,crseBox,crseDx,dt,time,
       ratio,numLevels);
    if (status != 0)
      {
	MayDay::Error("failed to read AMR hierarchy");
      }
    
    if (verbose)
      {
        pout() << "... done." << endl;
        pout() << endl;
        pout() << "extracting subregion..." << endl;
      }
    
    Vector<DisjointBoxLayout> outGridsTemp(numLevels);
    int outNumLevels = 0;
    Box levelSubdomain = subDomain;
    for (int lev=0; lev<numLevels; lev++)
      {
        LayoutIterator lit = grids[lev].layoutIterator();
        Vector<Box> levelBoxes;
        for (lit.begin(); lit.ok(); ++lit)
          {
            Box intersectBox = grids[lev][lit];
            intersectBox &= levelSubdomain;
            if (verbose) 
              {
                pout() << "subBox = " << intersectBox << endl;
              }
            if (!intersectBox.isEmpty())
              {
                levelBoxes.push_back(intersectBox);
              }
          }

        if (lev < ratio.size()) 
          {
            levelSubdomain.refine(ratio[lev]);
          }
        if (levelBoxes.size() > 0)
          {
            Vector<int> procAssign(levelBoxes.size());
            LoadBalance(procAssign, levelBoxes);
            DisjointBoxLayout levelGrids(levelBoxes, procAssign);
            outGridsTemp[lev] = levelGrids;
            outNumLevels = lev+1;                
          }
      }
    Vector<DisjointBoxLayout> outGrids(outNumLevels);
    for (int i=0; i<outNumLevels; i++)
      {
        outGrids[i] = outGridsTemp[i];
      }
    
    if (verbose)
      {
        pout() << "... done." << endl;
        pout() << "constructing output data" << endl;
      }

    Vector<LevelData<FArrayBox>* > outData(outNumLevels,NULL);
    for (int lev = 0; lev < outNumLevels; lev++)
      {
	
	outData[lev] = new LevelData<FArrayBox>(outGrids[lev],data[lev]->nComp(),data[lev]->ghostVect());
        data[lev]->copyTo(data[lev]->interval(),*outData[lev], outData[lev]->interval());
      }

    if (verbose)
      {
        pout() << "writing file..." << endl;
      }
    WriteAMRHierarchyHDF5(out_file, outGrids, outData, names, 
			  subDomain, crseDx, dt, time, ratio,outNumLevels);
    

    if (verbose)
      {
        pout() << "... done." << endl;
      }
  
    for (int lev=0; lev<data.size(); lev++)
      {
        if (data[lev] != NULL)
          {
            delete data[lev];
            data[lev] = NULL;
          }
      }
    for (int lev=0; lev<outData.size(); lev++)
      {
        if (outData[lev] != NULL)
          {
            delete outData[lev];
            outData[lev] = NULL;
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
