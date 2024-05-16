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
// ReadLevelData.cpp
//===========================================================================


#include "ReadLevelData.H"
#include "LoadBalance.H"
#include "BRMeshRefine.H"
#include "CoarseAverage.H"
#include "NamespaceHeader.H"

//fill a LevelData<FArrayBox> from data stored in an AMR file
//for each name in a_names look for either 
//(a) a single component LevelData named "name"
//(b) a multi-component LevelData, stored consecutively and beginning with name
void readLevelData(Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_data,
		   Real& a_dx,
		   const std::string a_file,
		   const Vector<std::string>& a_names, 
		   int a_nComp)
{

  Vector<LevelData<FArrayBox>* > vectData;
  Vector<DisjointBoxLayout> vectGrids;
  Vector<int> vectRatio;
  Vector<std::string> names;
  Real dt = 0.0,time = 0.0;
  Box domBox;
  int numLevels;

  pout() << " attempting to open file  " << a_file << std::endl;
#ifdef CH_USE_HDF5
  int status = ReadAMRHierarchyHDF5
    (a_file,vectGrids,vectData,names,domBox,a_dx,dt,time,
     vectRatio,numLevels);
  CH_assert(status == 0);
 if (status != 0)
   MayDay::Error("failed to read file");
#endif
 CH_assert(vectData.size() == 1);
 if (vectData.size() != 1)
   MayDay::Error("bad data");
 
 
 Vector<Box> boxes;
 int max_box_size = 64;
 int block_factor = 8;
 domainSplit(domBox, boxes, max_box_size, block_factor);
 Vector<int> procAssign(boxes.size());
 LoadBalance(procAssign,boxes);

 DisjointBoxLayout grids(boxes, procAssign, domBox);

 for (int i =0; i < a_data.size(); ++i)
   a_data[i]->define(grids,a_nComp,IntVect::Zero);

 int read = 0;
 
 for (int j = 0; j < a_names.size(); j++)
   {
     pout() << " looking for variable " << a_names[j] << std::endl;
     for (int i = 0; i < names.size(); i++)
       {
	 if (names[i] == a_names[j])
	   {
	     pout() << " found variable " << names[i] << std::endl;
	     vectData[0]->copyTo(Interval(i,i+a_nComp-1),*a_data[j],Interval(0,a_nComp-1));
	     read++;
	   }
       }
   }

 CH_assert(read == a_names.size());
 if (!read)
   MayDay::Error("no data in AMR Hierarchy");

 for (int i = 0; i < vectData.size(); i++)
   {
     if (vectData[i] != NULL)
       {
	 delete vectData[i]; vectData[i] = NULL;
       }
   }

}


///create one or more AMR Hieracrhies from data stored in an AMR file
/**
   
   \param Vector<Vector<RefCountedPtr<LevelData<FArrayBox> > > > a_data : each 
   Vector<RefCountedPtr<LevelData<FArrayBox> > > is an AMR hierarchy.
   \param a_crseDx : mesh spacing on level 0
   \param a_refRatio : refinement ratio between meshes.
   \param a_file : file name to read from
   \param a_name : for each name in a_names look for either 
   (a) a single component LevelData named "name"
   (b) a multi-component LevelData, stored consecutively and beginning with name
   \param a_nComp : the number of components in each heirachy
*/
void readMultiLevelData
(Vector< Vector <RefCountedPtr<LevelData<FArrayBox> > > >& a_data,
 Real& a_crseDx,
 Vector<int>& a_refRatio,
 const std::string a_file,
 const Vector<std::string>& a_names, 
 int a_nComp)
{

  Vector<LevelData<FArrayBox>* > vectData;
  Vector<DisjointBoxLayout> vectGrids;
  Vector<std::string> names;
  Real dt = 0.0,time = 0.0;
  Box domBox;
  int numLevels;

  pout() << " attempting to open file  " << a_file << std::endl;
#ifdef CH_USE_HDF5
  int status = ReadAMRHierarchyHDF5
    (a_file,vectGrids,vectData,names,domBox,a_crseDx,dt,time,
     a_refRatio,numLevels);
 CH_assert(status == 0);

 if (status != 0)
   MayDay::Error("failed to read file");
#endif 
 //We probably don't want to use the load balancing provided
 //given in the file, e.g that might lead to all the data being
 //stored on a single processor.

 //decompose level 0
 Vector<Box> boxes;
 int maxBoxSize = 64;
 // start with a blockfactor of 8, but reduce if needed
 int blockFactor = 8;
 Box tempBox(domBox);
 tempBox.coarsen(blockFactor);
 tempBox.refine(blockFactor);
 if (domBox != tempBox)
   {
     // need to reduce blockFactor
     bool ok = false;
     while (!ok)
       {
         blockFactor /= 2;
         tempBox = domBox;
         tempBox.coarsen(blockFactor);
         tempBox.refine(blockFactor);
         if (tempBox == domBox) ok = true;
       }
   } 
 domainSplit(domBox, boxes, maxBoxSize, blockFactor);
 Vector<int> procAssign(boxes.size());
 LoadBalance(procAssign,boxes);
 Vector<DisjointBoxLayout> grids(1,DisjointBoxLayout(boxes, procAssign, domBox));

 //decompose the higher levels, if the exist
 for (int lev = 1; lev < vectData.size(); lev++)
   {
     const DisjointBoxLayout& levelDBL =  vectGrids[lev];
     Vector<Box> levelBoxes(levelDBL.size());
     int index = 0;
     for (LayoutIterator lit = levelDBL.layoutIterator();lit.ok(); ++lit, ++index)
       {
	 levelBoxes[index] =  levelDBL[lit()];
       }
     for (int dir = 0; dir < SpaceDim;  dir++)
       {
	 breakBoxes(levelBoxes,  maxBoxSize, dir);
       }
     procAssign.resize(levelBoxes.size());
     LoadBalance(procAssign,levelBoxes);
     grids.push_back(DisjointBoxLayout(levelBoxes, procAssign, levelDBL.physDomain()));
   }


 //define all of the LevelData

 a_data.resize(a_names.size());
 for (int i =0; i < a_data.size(); ++i)
   { 
     for (int lev = 0; lev < grids.size(); lev++)
       {
	 a_data[i].push_back( RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>));
	 a_data[i][lev]->define(grids[lev],a_nComp,IntVect::Zero);       
       }
   }
 
 int read = 0;
 
 for (int j = 0; j < a_names.size(); j++)
   {
     pout() << " looking for variable " << a_names[j] << std::endl;
     for (int i = 0; i < names.size(); i++)
       {
	 if (names[i] == a_names[j])
	   {
	     pout() << " found variable " << names[i] << std::endl;
	     for (int lev = 0; lev <  grids.size(); lev++)
	       {
		 vectData[lev]->copyTo(Interval(i,i+a_nComp-1),*a_data[j][lev],Interval(0,a_nComp-1));
	       }
	     read++;
	   }
       }

     //SLC : possibly this coarse average should be optional, but it seems to me that
     //in the most common use cases (reading DEMs etc) we will always want to 
     //ensure that coarse levels approximate fine levels.
     for (int lev = a_data[j].size() - 1; lev > 0; lev--)
       {
	 const LevelData<FArrayBox>& fine = *a_data[j][lev];
	 LevelData<FArrayBox>& crse = *a_data[j][lev-1];
	 CoarseAverage av(fine.disjointBoxLayout(), fine.nComp(), a_refRatio[lev-1] );
	 av.averageToCoarse(crse,fine);
       }

   }

 CH_assert(read == a_names.size());
 if (!read)
   MayDay::Error("no data in AMR Hierarchy");

 for (int i = 0; i < vectData.size(); i++)
   {
     if (vectData[i] != NULL)
       {
	 delete vectData[i]; vectData[i] = NULL;
       }
   }

}


#include "NamespaceFooter.H"
