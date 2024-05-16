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
// sum.cpp -- reads in two hdf5 files, sums the components 
// which are present in both. assumes same grid hierarchies.
//===========================================================================

#include <iostream>
#include "AMRIO.H"
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
	std::cerr << " usage: " << argv[0] << " <input_file1> <input_file2> <output_file> " << std::endl; 
	exit(0); 
      }

    char* in_file1 = argv[1];
    char* in_file2 = argv[2];
    char* out_file = argv[3];

 
    if (verbose)
      {
        pout ();
        pout() << "reading first AMR file ("<< in_file1 << ")..." << endl;
      }
    Vector<LevelData<FArrayBox>* > data1;
    Vector<DisjointBoxLayout> grids1;
    Vector<string> names1;
    Vector<int> ratio1;
    Real crseDx1 = 0.0, dt = 0.0, time = 0.0;
    Box crseBox1;
    int numLevels1;
    int status = ReadAMRHierarchyHDF5
      (in_file1,grids1,data1,names1,crseBox1,crseDx1,dt,time,
       ratio1,numLevels1);
    if (status != 0)
      {
        MayDay::Error("failed to read AMR hierarchy");
      }
    
    if (verbose)
      {
        pout() << "... done." << endl;
      }

    if (verbose)
      {
        pout ();
        pout() << "reading second AMR file ("<< in_file2 << ")..." << endl;
      }
    Vector<LevelData<FArrayBox>* > data2;
    Vector<DisjointBoxLayout> grids2;
    Vector<string> names2;
    Vector<int> ratio2;
    Real crseDx2 = 0.0;
    Box crseBox2;
    int numLevels2;
    status = ReadAMRHierarchyHDF5
      (in_file2,grids2,data2,names2,crseBox2,crseDx2,dt,time,
       ratio2,numLevels2);
    if (status != 0)
      {
        MayDay::Error("failed to read AMR hierarchy");
      }
    
    if (verbose)
      {
        pout() << "... done." << endl;
      }

    CH_assert(numLevels1 == numLevels2);

    for (int lev=0; lev<data1.size(); lev++)
      {
        LevelData<FArrayBox>& dataLev1 = *data1[lev];
        LevelData<FArrayBox>& dataLev2 = *data2[lev];
        LevelData<FArrayBox> tempData(dataLev1.getBoxes(),
                                     dataLev2.nComp(),
                                     dataLev1.ghostVect());
        dataLev2.copyTo(dataLev2.interval(), tempData,
                        tempData.interval());

        DataIterator dit = dataLev1.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
          {
            FArrayBox& thisData1 = dataLev1[dit];
            FArrayBox& thisTemp = tempData[dit];
            int comp2 = -1;
            for (int comp1=0; comp1<thisData1.nComp(); comp1++)
              {
                for (int tempComp=0; tempComp<tempData.nComp();tempComp++)
                  {
                    if (names1[comp1] == names2[tempComp])
                      {
                        comp2 = tempComp;
                      }
                  } // end loop over data2 names
                if (comp2 >= 0)
                  {
                    thisData1.plus(thisTemp,comp2, comp1, 1);
                  }
              } // end loop over comps in dest
          } // end loop over boxes
      } // end loop over levels

    WriteAMRHierarchyHDF5(out_file, grids1, data1, names1, 
                          crseBox1, crseDx1, dt, time, ratio1, numLevels1);
    
    // clean up memory
    for (int lev=0; lev<data1.size(); lev++)
      {
        if (data1[lev] != NULL)
          {
            delete data1[lev];
            data1[lev] = NULL;
          }

        if (data2[lev] != NULL)
          {
            delete data2[lev];
            data2[lev] = NULL;
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
