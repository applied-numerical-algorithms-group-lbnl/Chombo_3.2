#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SPMD_H_
#define _SPMD_H_
///
/**
  This is a grab bag of useful stuff used in Chombo MPI code.
  This includes utilities for Chombo serialization routines used in more advanced 
  parallel data structures.
*/
#include "REAL.H"
#include <string>
#include "Vector.H"
#include "MayDay.H"
#include "Misc.H"

#ifdef CH_MPI
#include "mpi.h"
#endif

#include "memtrack.H"
#include "BaseNamespaceHeader.H"

struct Chombo_MPI
{
  
#ifdef CH_MPI
  static MPI_Comm comm;
#endif
  
};


extern unsigned long long CH_MAX_MPI_MESSAGE_SIZE;
extern unsigned long long CH_MaxMPISendSize;
extern unsigned long long CH_MaxMPIRecvSize;

/// local process ID
/**
    Returns the ID of the locally running process in the range 0 <=
    procID() < numProc().  This has no relation to the operating system
    pid.  There is always a procID() == 0.  */
inline int
procID(
#ifdef CH_MPI  
  MPI_Comm a_comm = Chombo_MPI::comm
#endif  
  )
{
#ifdef CH_MPI  
  int iproc =  4586;
  MPI_Comm_rank(a_comm,  &iproc);
  return iproc;
#else
  return 0;
#endif  
}

/// number of parallel processes
/**
    Returns the number of parallel processes running.
    Always returns at least 1. */
inline unsigned int
numProc(
#ifdef CH_MPI  
  MPI_Comm a_comm = Chombo_MPI::comm
#endif  
  )
{
#ifdef CH_MPI
  static int ret = -1;
  if (ret == -1)
  {
    MPI_Comm_size(a_comm, &ret);
  }
  return ret;
#else
  return 1;
#endif  
}

///all ranks wait here to sync-up
/**
   All MPI ranks wait here to sync-up.  Calls MPI_Barrier(comm)  */
inline void barrier(
#ifdef CH_MPI            
       MPI_Comm a_comm = Chombo_MPI::comm
#endif
  )
{
#ifdef CH_MPI            
  MPI_Barrier(a_comm);
#endif
}

template <class T>
int linearSize(const T& inputT);

template <class T>
void linearIn(T& a_outputT, const void* const inBuf);

template <class T>
void linearOut(void* const a_outBuf, const T& inputT);

///gather from multiple processes
/**
  Gathers a_input of every process onto Vector<T> a_outVec on process
  number a_dest.  It is an error if a_dest < 0 or a_dest >= numProc().
  a_outVec is numProc() long with the value of a_input on every
  process in its elements.\\

T must have have the following functions:
  \\

  int linearSize<T>(const T& inputT); \\
  void linearIn<T>(T& a_outputT, const void* const inBuf); \\
  void linearOut<T>(void* a_outBuf, const T& inputT); \\

*/
template <class T>
inline void
gather(Vector<T>& a_outVec,
       const T  & a_input,
       int        a_dest
#ifdef CH_MPI
            ,MPI_Comm a_comm = Chombo_MPI::comm
#endif            
  )
{
#ifndef CH_MPI
  a_outVec.resize(1);
  a_outVec[0] = a_input;
#else
  CH_assert (a_dest >= 0);
  CH_assert(a_dest <  numProc());
  //now THIS size lives on THIS processor
  int isize = linearSize(a_input);

  //make stuff for linearout
  void* loclBuf = mallocMT(isize);
  if (loclBuf == NULL)
    MayDay::Error("out of memory in gather 1");

  //put linearized T into its proper buffer
  linearOut(loclBuf, a_input);

  int nProcess = numProc();
  int sendCount = 1;
  int recdCount = 1;

  //need to gather isizes onto processor a_dest
  int* vectSize = NULL;
  int* vectDisp = NULL;
  void* sendBuf = static_cast<void*>(&isize);
  //allocate received buffer
  if (procID(a_comm) == a_dest)
    {
      vectSize = new int[nProcess];
      vectDisp = new int[nProcess];
    }

  int result1 = MPI_Gather(sendBuf, sendCount, MPI_INT,
                           vectSize,recdCount, MPI_INT,
                           a_dest,  a_comm);

  if (result1 != MPI_SUCCESS)
    MayDay::Error("Gather<T> failed in MPI_Gather 1");

  //make memory for gather, linearin
  void* recdBuf = NULL;
  if (procID(a_comm) == a_dest)
    {
      size_t itotsize=0;
      for (int iproc = 0; iproc < nProcess; iproc++)
        {
          vectDisp[iproc] = itotsize;
          itotsize += vectSize[iproc];
        }
      recdBuf = mallocMT(itotsize);
      if (recdBuf == NULL)
        {
          MayDay::Error("out of memory in gather 2");
        }
    }

  //gather data
  int result2 = MPI_Gatherv(loclBuf, isize, MPI_BYTE,
                            recdBuf, vectSize, vectDisp, MPI_BYTE,
                            a_dest, a_comm);
  if (result2 != MPI_SUCCESS)
    MayDay::Error("Gather<T> failed in MPI_Gather 2");

  if (procID(a_comm) == a_dest)
    {
      //calculate offset into array for current processor
      int ioffset = 0;
      a_outVec.resize(nProcess);
      //need to cast to char* to do pointer arithmetic
      char* arithPtr = (char*)recdBuf;
      for (int iproc = 0; iproc < nProcess; iproc++)
        {
          ioffset = vectDisp[iproc];
          char* thisProcBuf = arithPtr + ioffset;
          linearIn(a_outVec[iproc], thisProcBuf);
        }

      //delete memory for dest-specific arrays
      delete[] vectSize;
      delete[] vectDisp;
      freeMT(recdBuf);
    }

  //delete memory for local buffer
  freeMT(loclBuf);
#endif  
}

///broadcast to every process
/**
  Broadcasts a_inAndOut to every process from process number a_src.
  It is an error if a_src < 0 or a_src >= numProc().  \\

  T must have have the following functions: \\

  int linearSize<T>(const T& inputT); \\
  void linearIn<T>(T& a_outputT, const void* const inBuf); \\
  void linearOut<T>(void* a_outBuf, const T& inputT); \\ */
template <class T>
inline void
broadcast(T& a_inAndOut,
          int a_src
#ifdef CH_MPI               
          ,MPI_Comm a_comm = Chombo_MPI::comm
#endif
  )
{
  //only does anything if in parallel 
#ifdef CH_MPI
  CH_assert (a_src >= 0);
  CH_assert(a_src <  numProc());
  int isize;
  if (procID(a_comm) == a_src)
  {
    isize = linearSize(a_inAndOut);
  }

  MPI_Bcast(&isize, 1, MPI_INT, a_src, a_comm);

  void* broadBuf = mallocMT(isize);

  if (broadBuf == NULL)
  {
    MayDay::Error("out of memory in broadcast");
  }

  //take inAndOut from src and put it into broadBuf
  if (procID(a_comm) == a_src)
  {
    linearOut(broadBuf, a_inAndOut);
  }

  //broadcast broadBuf to all procs
  MPI_Bcast(broadBuf, isize, MPI_BYTE, a_src, a_comm);

  if (procID(a_comm)==a_src)
  {
    CH_MaxMPISendSize = Max<long long>(CH_MaxMPISendSize, isize);
  }
  else
  {
    CH_MaxMPIRecvSize = Max<long long>(CH_MaxMPIRecvSize, isize);
  }
  //take broadBuf and put back into inAndOut if not src
  if (procID(a_comm) != a_src)
  {
    linearIn(a_inAndOut, broadBuf);
  }

  //delete memory for buffer
  freeMT(broadBuf);
#endif  
}
//u l l specialization of linearSize
template < >
int linearSize(const long long& a_input);

//u l l specialization of linearIn
template < >
void linearIn(long long& a_outputT, const void* const inBuf);

//unsigned long long specialization of linearOut
template < >
void linearOut(void* const a_outBuf, const long long& a_inputT);

//integer specialization of linearSize
template < >
int linearSize(const int& a_input);

//u l l specialization of linearSize
template < >
int linearSize(const unsigned long long& a_input);

//integer specialization of linearIn
template < >
void linearIn(int& a_outputT, const void* const inBuf);

//u l l specialization of linearIn
template < >
void linearIn(unsigned long long& a_outputT, const void* const inBuf);

//integer specialization of linearOut
template < >
void linearOut(void* const a_outBuf, const int& a_inputT);

//unsigned long long specialization of linearOut
template < >
void linearOut(void* const a_outBuf, const unsigned long long& a_inputT);

//long integer specialization of linearSize
template < >
int linearSize(const long& a_input);

//long integer specialization of linearIn
template < >
void linearIn(long& a_outputT, const void* const inBuf);

//long integer specialization of linearOut
template < >
void linearOut(void* const a_outBuf, const long& a_inputT);

//unsigned long integer specialization of linearSize
template < >
int linearSize(const unsigned long& a_input);

//unsigned long integer specialization of linearIn
template < >
void linearIn(unsigned long& a_outputT, const void* const inBuf);

//unsigned long integer specialization of linearOut
template < >
void linearOut(void* const a_outBuf, const unsigned long& a_inputT);

//Real specialization of linearSize
template < >
int linearSize(const float& a_input);

template < >
int linearSize(const double& a_input);

//Real specialization of linearIn
template < >
void linearIn(float& a_outputT, const void* const a_inBuf);

template < >
void linearIn(double& a_outputT, const void* const a_inBuf);

//Real specialization of linearOut
template < >
void linearOut(void* const a_outBuf, const float& a_inputT);

template < >
void linearOut(void* const a_outBuf, const double& a_inputT);

// std::string specialization.
template <>
int linearSize(const std::string& a_input);
template <>
void linearIn(std::string& a_outputT, const void* const a_inBuf);
template <>
void linearOut(void* const a_outBuf, const std::string& a_inputT);

//Vector<int>  specialization
template < >
int linearSize(const Vector<int>& a_input);
template < >
void linearIn(Vector<int>& a_outputT, const void* const inBuf);
template < >
void linearOut(void* const a_outBuf, const Vector<int>& a_inputT);

//Vector<unsigned long long>  specialization
template < >
int linearSize(const Vector<unsigned long long>& a_input);
template < >
void linearIn(Vector<unsigned long long>& a_outputT, const void* const inBuf);
template < >
void linearOut(void* const a_outBuf, const Vector<unsigned long long>& a_inputT);

//Vector<long>  specialization
template < >
int linearSize(const Vector<long>& a_input);
template < >
void linearIn(Vector<long>& a_outputT, const void* const inBuf);
template < >
void linearOut(void* const a_outBuf, const Vector<long>& a_inputT);

//Vector<Real>  specialization
template < >
int linearSize(const Vector<float>& a_input);
template < >
void linearIn(Vector<float>& a_outputT, const void* const inBuf);
template < >
void linearOut(void* const a_outBuf, const Vector<float>& a_inputT);

template < >
int linearSize(const Vector<double>& a_input);
template < >
void linearIn(Vector<double>& a_outputT, const void* const inBuf);
template < >
void linearOut(void* const a_outBuf, const Vector<double>& a_inputT);

//Vector<std::string>  specialization
template < >
int linearSize(const Vector<std::string>& a_input);
template < >
void linearIn(Vector<std::string>& a_outputT, const void* const inBuf);
template < >
void linearOut(void* const a_outBuf, const Vector<std::string>& a_inputT);

//Vector<Vector<int> > specialization
template < >
int linearSize(const Vector<Vector<int> >& a_input);
template < >
void linearIn(Vector<Vector<int> >& a_outputT, const void* const inBuf);
template < >
void linearOut(void* const a_outBuf, const Vector<Vector<int> >& a_inputT);

//Vector<T> specialization of linearSize
template <class T>
int linearListSize(const Vector<T>& a_input);

//Vector<T> specialization of linearIn
template <class T>
void linearListIn(Vector<T>& a_outputT, const void* const a_inBuf);

//Vector<T> specialization of linearOut
template <class T>
void linearListOut(void* const a_outBuf, const Vector<T>& a_inputT);

class SerialTask
{
public:
    enum task
    {
      compute=0
    };
};

int GetPID(int rank);
int GetRank(int pid);

// return id of unique processor for special serial tasks
int
uniqueProc(const SerialTask::task& a_task);

#include "BaseNamespaceFooter.H"

#include "SPMDI.H"

#endif
