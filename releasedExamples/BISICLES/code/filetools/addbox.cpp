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
// addbox.cpp
// read a multi-level AMR hierarchy, add a new box to one of the higher
// levels, creating said higher levels if needed, write the result to a Chombo hdf5 file
//===========================================================================


#include <iostream>
#include "AMRIO.H"
#include "FineInterp.H"
#include "NamespaceHeader.H"

bool verbose = true;

//given a std::string a,b,c create an IntVect with compeonets a,b,c
IntVect split(const std::string a_s, const char a_c)
{
  IntVect r;
  int i = 0;
  int k = 0;
  
  while (a_s[i] == a_c)
    i++;
  do {
    int j = std::min( a_s.find(a_c,i), a_s.length());
    r[k] = atoi( (a_s.substr(i,j-i)).c_str() ); k++;
    i = j + 1;
  } while (k < SpaceDim && i < a_s.length()); 
  if (k != SpaceDim)
    {
      cerr << "failed to parse " << a_s << " as an IntVect" << endl;
      MayDay::Error("failed to parse IntVect");
    }
  return r;
}

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
    if(argc < 6) 
      { 
	std::cerr << " usage: " << argv[0] << " <input_file> <output_file> nRef ilo[,jlo[,klo]] ihi[,jhi[,khi]] " << std::endl; exit(1); 
      }
	char* in_file = argv[1];
	char* out_file = argv[2];
	int nRef =  atoi(argv[3]);
	IntVect lo = split(std::string(argv[4]),',');
	IntVect hi = split(std::string(argv[5]),',');
	Box b(lo,hi);

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

	if (ratio[0] == 1 && numLevels == 1)
	  {
	    ratio[0] = nRef;
	  }

	int lev = 1;
	int levelRef = ratio[0];
	while (nRef > levelRef && lev < numLevels )
	  {
	    lev++;
	    levelRef *= ratio[lev-1];
	  }
	
	if (nRef < levelRef)
	  {
	    cerr << "nRef = " << nRef << "falls between levels";
	    MayDay::Error("invalid nRef");
	  }

	
	int outNumLevels = numLevels;

	Vector<int> outRatio(ratio);
	if (numLevels == 1)
	  {
	    outNumLevels += 1;
	  }
	else if (nRef > levelRef)
	  { 
	    outNumLevels += 1;
	    outRatio.resize(outNumLevels);
	    int lastRef = nRef/levelRef;
	    CH_assert (lastRef == 2 || lastRef == 4 || lastRef == 8 || lastRef == 16);
	    CH_assert (lastRef*levelRef == nRef);
	    outRatio[outNumLevels-1] = nRef/levelRef;
	  }

	
	pout() << "adding box " << b << std::endl;
	
	Vector<DisjointBoxLayout> outGrids(grids);
	Vector<LevelData<FArrayBox>* > outData(data);
	if (outNumLevels > numLevels)
	  {
	    //in this case we need to create a new level with a single box
	    Vector<Box> boxes(1,b);
	    Vector<int> procIDs(1,0);
	    //LoadBalance(procIDs, boxes);
	    ProblemDomain domain = grids[numLevels-1].physDomain();
	    domain.refine(outRatio[outNumLevels-2]);
	    const DisjointBoxLayout dbl(boxes,procIDs,domain);
	    outGrids.resize(outNumLevels,dbl);
	    outData.resize(outNumLevels,new LevelData<FArrayBox>(dbl,data[0]->nComp(),data[0]->ghostVect()));
	    FineInterp fi(dbl,data[0]->nComp(),outRatio[outNumLevels-2],dbl.physDomain());
	    fi.interpToFine(*outData[outNumLevels-1],*data[numLevels-1]);
	    
	  }
	else 
	  {
	    //in this case we need to add a single box to the top leve
	    const DisjointBoxLayout& levelDBL = grids[lev];
	    
	    Vector<Box> boxes(levelDBL.size() + 1);
	    LayoutIterator lit = levelDBL.layoutIterator();
	    int i = 0;
	    for (lit.begin(); lit.ok(); ++lit, ++i) 
	      {
		boxes[i] = levelDBL[lit()];
	      }
	    boxes[i] = b; 
	    Vector<int> procIDs(levelDBL.size() + 1,0);
	    ProblemDomain domain = grids[lev].physDomain();
	    const DisjointBoxLayout dbl(boxes,procIDs,domain);
	    outData[lev] = new LevelData<FArrayBox>(dbl,data[lev]->nComp(),data[lev]->ghostVect());
	    FineInterp fi(dbl,data[lev]->nComp(),outRatio[outNumLevels-2],dbl.physDomain());
	    fi.interpToFine(*outData[lev],*data[lev-1]);
	    int nComp = data[lev]->nComp();
	    data[lev]->copyTo(Interval(0,nComp-1),*outData[lev],Interval(0,nComp-1));

	  }
	

	WriteAMRHierarchyHDF5(out_file, outGrids, outData, names, 
			      crseBox, crseDx, dt, time, outRatio,outNumLevels);
	
	
	
	for (int lev = 0; lev < numLevels; lev++)
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
