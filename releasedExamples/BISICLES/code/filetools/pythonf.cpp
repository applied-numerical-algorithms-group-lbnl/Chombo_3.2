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
// pythonf.cpp
//
// pythonf A.hdf5 B.hdf5 s.py f a,b,c,... p,q,r,s,...
//
// read a multi-level AMR hierarchy A with variables (a,b,c...), 
// and writre multi-level AMR hierarchy B with variables (p,q,r,s...)
// given a python script s.py and function name f such that
// (p,q,r,s,...) = f(a,b,c,...)
//===========================================================================

#include <iostream>
#include "AMRIO.H"
#include "BoxIterator.H"
#include <Python.h>
#include "NamespaceHeader.H"

//given a std::string a<a_c>b<a_c>... , create a Vector<std::string>
//whose elements are a,b,...
Vector<std::string> split(const std::string a_s, const char a_c)
{
  Vector<std::string> r;
  int i =0;
  while (a_s[i] == a_c)
    i++;
  do {
    int j = std::min( a_s.find(a_c,i), a_s.length());
    r.push_back( a_s.substr(i,j-i));
    i = j + 1;
  } while (i < a_s.length()); 
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

    if(argc != 7) 
      { 
	std::cerr << " usage: " << argv[0] << " <input_file> <output_file> <python script> <python function> <input tuple> <output tuple>" << std::endl; 
	exit(0); 
      }

    std::string inFile(argv[1]);
    std::string outFile(argv[2]);
    char* pyScript = (argv[3]);
    char* pyFunc = (argv[4]);
    Vector<std::string> inTuple = split(std::string(argv[5]),',');
    Vector<std::string> outTuple = split(std::string(argv[6]),',');

    //read the input data;
    Vector<std::string> inNames;
    Vector<LevelData<FArrayBox>* > inData;
    Vector<DisjointBoxLayout> grids;
    Vector<int> ratio;
    Real crseDx = 0.0, dt = 0.0, time = 0.0;
    Box crseBox;
    int numLevels;
    int status = ReadAMRHierarchyHDF5
      (inFile,grids,inData,inNames,crseBox,crseDx,dt,time,
       ratio,numLevels);
    if (status != 0)
      {
	MayDay::Error("failed to read AMR hierarchy");
      }

    Vector<Real> dx(1,crseDx);
    for (int lev = 1; lev < numLevels; lev++)
      {
	dx.push_back(dx[lev-1]/Real(ratio[lev-1]));
      }


    //lookup table inTuple -> inNames
    Vector<int> inComps(inTuple.size());
    for (int j =0; j < inTuple.size(); j++)
      {
	bool found = false;
	for (int i =0; i < inNames.size(); i++)
	  {
	    if (inTuple[j] == inNames[i])
	      {
		inComps[j] = i;
		found = true;
		break;
	      }
	  }
	if (!found)
	  {
	    pout() << inTuple[j] << " not found in input ";
	    MayDay::Error("missing variable in input file");
	  }
      }


    //allocate memory for the output data
    Vector<LevelData<FArrayBox>* > outData(numLevels,NULL);
    for (int lev = 0; lev < numLevels; lev++)
      {
	outData[lev] = new LevelData<FArrayBox>(grids[lev],int(outTuple.size()),inData[lev]->ghostVect());
      }

    //pythonate
    {
      PyObject *pName, *pModule,  *pFunc;
      PyObject *pArgs, *pValue;

      Py_Initialize();

      pName = PyUnicode_FromString(pyScript);
      pModule = PyImport_Import(pName); // imports the script pName into the python interpreter
      Py_DECREF(pName);//presumably, python can throw away this object now

      if (pModule != NULL) 
	{
	  pFunc = PyObject_GetAttrString(pModule, pyFunc);
	  if (pFunc && PyCallable_Check(pFunc)) 
	    {
	      for (int lev = 0; lev < numLevels; lev++)
		{
		  const LevelData<FArrayBox>& inLevelData = *inData[lev];
		  LevelData<FArrayBox>& outLevelData = *outData[lev];
		  for (DataIterator dit(grids[lev]);dit.ok();++dit)
		    {
		      const FArrayBox& inFab = inLevelData[dit];
		      FArrayBox& outFab = outLevelData[dit];
		      for (BoxIterator bit(grids[lev][dit]);bit.ok();++bit)
			{
			  const IntVect& iv = bit();

			  //construct a python tuple of args, built from the input tuple
                          // and ( x[,y[,z]] i[,j[,k]] level dx_level)
			  int nExtraArgs = SpaceDim + SpaceDim + 2;
			  pArgs = PyTuple_New(inTuple.size() + nExtraArgs);
			  int i;
			  for (i = 0; i < inComps.size(); i++)
			    {
			      pValue = PyFloat_FromDouble(inFab(iv,inComps[i]));
			      if (!pValue)
				{
				  Py_DECREF(pArgs);
				  Py_DECREF(pModule);
				  MayDay::Error("fab(iv) -> python failed");
				}
			      PyTuple_SetItem(pArgs, i, pValue);
			    }
			  
			  RealVect x = dx[lev]*(RealVect(iv) + 0.5*RealVect::Unit);
			  for (int dir = 0; dir < SpaceDim; dir++)
			    {
			      
			      PyTuple_SetItem(pArgs, i, PyFloat_FromDouble(x[dir]));i++;
			    }
			  for (int dir = 0; dir < SpaceDim; dir++)
			    {
			      
			      PyTuple_SetItem(pArgs, i, PyFloat_FromDouble(Real(iv[dir])));i++;
			    }
			  
			  
			  PyTuple_SetItem(pArgs, i, PyFloat_FromDouble(dx[lev]));i++;

			  
			  PyTuple_SetItem(pArgs, i, PyFloat_FromDouble(Real(lev)));i++;

			  //call the python function
			  pValue = PyObject_CallObject(pFunc, pArgs);
			  Py_DECREF(pArgs);
			  if (pValue != NULL) 
			    {
			      //retrive data from python
			      if (PyTuple_CheckExact(pValue))
				{
				  for (int i = 0; i < outFab.nComp() ; ++i)
				    {
				      PyObject *pComp = PyTuple_GetItem(pValue,i);
				      if (pComp != NULL)
					{
					  outFab(iv,i) = PyFloat_AS_DOUBLE(pComp);
					}
				      else
					{
					  Py_DECREF(pFunc);
					  Py_DECREF(pModule);
					  Py_DECREF(pValue);
					  MayDay::Error("could not extract component");
					}
				    }
				}
			      else if (outFab.nComp() == 1 && PyFloat_CheckExact(pValue))
				{
				  outFab(iv) = PyFloat_AS_DOUBLE(pValue);
				}
			      else
				{
				  Py_DECREF(pFunc);
				  Py_DECREF(pModule);
				  Py_DECREF(pValue);
				  MayDay::Error("return value is an unsupported object");
				}

			      Py_DECREF(pValue);
			    }
			  else
			    {
			      Py_DECREF(pFunc);
			      Py_DECREF(pModule);
			      PyErr_Print(); 
			      MayDay::Error("python call failed");
			    }
			}
		    }
		}
	    }
	  else
	    {
	      if (PyErr_Occurred())
                PyErr_Print();
	      MayDay::Error("python cannot find function");
	    }
	} 
      else
	{
	  PyErr_Print();
	  MayDay::Error("python import failed");
	}


      Py_Finalize();
    }

    //write output
    WriteAMRHierarchyHDF5(outFile,grids,outData,outTuple,crseBox,crseDx,dt,time,
			  ratio,numLevels);


    for (int lev = 0; lev < numLevels; lev++)
      {
	if (outData[lev] != NULL)
	  {
	    delete outData[lev];
	    outData[lev] = NULL;
	  }
	if (inData[lev] != NULL)
	  {
	    delete inData[lev];
	    inData[lev] = NULL;
	  }
      }
    

    int debug = 0; debug++;

 }  // end nested scope
  CH_TIMER_REPORT();

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return 0;
}

#include "NamespaceFooter.H"
