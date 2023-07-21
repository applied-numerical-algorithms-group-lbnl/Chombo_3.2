#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SPMDI_H_
#define _SPMDI_H_

#include "memtrack.H"
#include "parstream.H"
#include "BaseNamespaceHeader.H"

///
/**
 Default implementations for linearization routines.
 These should work independent of MPI
**/

template <class T>
int linearSize(const T& inputT)
{
  return inputT.linearSize();
}

template <class T>
void linearIn(T& a_outputT, const void* const inBuf)
{
  a_outputT.linearIn(inBuf);
}

template <class T>
void linearOut(void* const a_outBuf, const T& inputT)
{
  inputT.linearOut(a_outBuf);
}



//Vector<T> specialization of linearIn
template <class T>
void
linearListIn(Vector<T>& a_outputT, const void* const a_inBuf)
{
  //first entry is the size of the vector
  const int* const intBuf = (int*)a_inBuf;
  int vecsize = intBuf[0];
  Vector<int> vecOffset(vecsize);
  //next vecsize entries are offsets of data into buffer
  for (int ivec = 0; ivec < vecsize; ivec++)
    {
      vecOffset[ivec] = intBuf[ivec+1];
    }
  //next vecsize entries are the actual data
  //yes I could do this in one loop but that would
  // either
  // a) make it less symmetric with linearOut
  // and/or
  // b) make both of them far less readable
  a_outputT.resize(vecsize);
  const char* const charbuf = (char*)a_inBuf;
  for (int ivec = 0; ivec < vecsize; ivec++)
    {
      const char* const dataLoc = charbuf + vecOffset[ivec];
      linearIn(a_outputT[ivec], dataLoc);
    }
}

//Vector<T> specialization of linearOut
template <class T>
void
linearListOut(void* const a_outBuf, const Vector<T>& a_input)
{
  //first entry is the size of the vector
  int* const intBuf = (int*)a_outBuf;
  intBuf[0] = a_input.size();
  int vecsize = intBuf[0];
  Vector<int> vecOffset(vecsize);
  //next vecsize entries are offsets of data into buffer
  //next vecsize entries are the actual data
  int ioffset = (vecsize+1)*sizeof(int);
  for (int ivec = 0; ivec < vecsize; ivec++)
    {
      intBuf[ivec+1] = ioffset;
      vecOffset[ivec] = ioffset;
      ioffset += linearSize(a_input[ivec]);
    }
  //yes I could do this in one loop but that would
  // either
  // a) make it less symmetric with linearIn
  // and/or
  // b) make both of them far less readable
  char* const charBuf = (char*)a_outBuf;
  for (int ivec = 0; ivec < vecsize; ivec++)
    {
      char* const dataLoc = charBuf + vecOffset[ivec];
      linearOut(dataLoc, a_input[ivec]);
    }
}

//Vector<T> specialization of linearSize
template <class T>
int
linearListSize(const Vector<T>& a_input)
{
  //first entry is the size of the vector (int)
  //next vecsize entries are offsets of data into buffer (int)
  //next vecsize entries are the actual data
  int itotsize = (a_input.size() + 1)*sizeof(int);
  for (unsigned int ivec = 0; ivec < a_input.size(); ivec++)
    {
      itotsize += linearSize(a_input[ivec]);
    }
  return itotsize;
}

#include "BaseNamespaceFooter.H"

#endif
