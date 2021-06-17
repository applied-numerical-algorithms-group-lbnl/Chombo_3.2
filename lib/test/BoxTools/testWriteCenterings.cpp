#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "REAL.H"
#include "parstream.H"

#include "LevelData.H"
#include "LoadBalance.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "EdgeDataBox.H"
#include "IntVectSet.H"
#include "AMRIO.H"
#include "BoxIterator.H"
#include "BRMeshRefine.H" // contains domainSplit function
#include "CH_Attach.H"
#include "CH_HDF5.H"
#include "UsingNamespace.H"

using std::endl;

/// Global variables for handling output:
static const char *pgmname = "testWriteCenterings" ;
static const char *indent = "   ", *indent2 = "      " ;
static bool verbose = true ;

static bool write_cell = true;
static bool write_face = true;
static bool write_edge = true;


/***************/
/***************/
void dumpmemoryatexit();
/***************/
/***************/


void
parseTestOptions( int argc ,char* argv[] ) ;

int writeLevelTest();

int
main(int argc,char **argv)
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
  //  registerDebugger();
#endif
  parseTestOptions( argc ,argv ) ;

  int retval = 0;

  retval = writeLevelTest();

  if (retval == 0)
    {
      pout() << "testWriteCenterings passed all tests!" << endl;
    }
  else
    {
      pout() << "testWriteCenterings failed at least one test, return value = "
             << retval << endl;
    }

#ifdef CH_MPI
  MPI_Finalize();
#endif

  return retval;
}


int writeLevelTest()
{
  int status = 0;
  
  int domLength =  128;
  int maxBoxSize =   64;

  if (SpaceDim > 3)
    {
      domLength =  4;
      maxBoxSize =  2;  
    }

  Box domainBox(IntVect::Zero, (domLength-1)*IntVect::Unit);
  
  // start with non-periodic case
  ProblemDomain domain(domainBox);

  pout() << "calling domainSplit" << std::endl;
  Vector<Box> gridBoxes;
  domainSplit(domainBox, gridBoxes, maxBoxSize);
  int numBoxes = gridBoxes.size();
  Vector<int> procAssign(numBoxes,0);

  pout() << "calling loadbalance" << std::endl;
  // for mpi case, distribute boxes among processors
  LoadBalance(procAssign,gridBoxes);
  pout() << "calling dbl define" << std::endl;
  DisjointBoxLayout grids(gridBoxes, procAssign, domain);

  
  // cell-centered data
  //const int nComp = SpaceDim;
  const int nComp = 1;
  const Interval comps(0,nComp-1);
  //IntVect ghostVect = IntVect::Unit;
  IntVect ghostVect = IntVect::Zero;
  LevelData<FArrayBox> cellData(grids, nComp, ghostVect);
  
  DataIterator dit = grids.dataIterator();
  Real n = 0.0;
  
  for (dit.begin(); dit.ok(); ++dit)
    {
      cellData[dit].setVal(n);
      n += 1;
    }
  
  // might as well..
  //writeLevel(&cellData);

  
  // face-centered data
  LevelData<FluxBox> faceData(grids, nComp, ghostVect);
    
  n = 0.0;
  
  for (dit.begin(); dit.ok(); ++dit)
    {
      for (int dir=0; dir<SpaceDim; dir++)
        {
          Real faceN = n + 0.1*dir;
          faceData[dit][dir].setVal(faceN);          
        }
      n += 1;                
    }
    
  // edge-centered data
  LevelData<EdgeDataBox> edgeData(grids, nComp, ghostVect);
  
  n = 0.0;
    
  for (dit.begin(); dit.ok(); ++dit)
    {
      for (int dir=0; dir<SpaceDim; dir++)
        {
          Real edgeN = n + 0.1*dir;
          edgeData[dit][dir].setVal(edgeN);          
        }
      n += 1;                
    }
  
#ifdef CH_USE_HDF5

  const std::string plotFileName = "CenteringTest.hdf5";
  
  //  test writing
  {
    HDF5Handle handle( plotFileName.c_str(), HDF5Handle::CREATE );
    
    handle.setGroup("/");
    
    //
    // write the default header
    //
    
    // set the component names
    HDF5HeaderData header;
    
    Vector<string> vectNames;
    char field_name[50];
    char comp_name[50];
    char coords[SpaceDim];
    coords[0] = 'x';
    if(SpaceDim==2) coords[1] = 'z';
    if(SpaceDim==3) {
      coords[1] = 'y';
      coords[2] = 'z';
    }
    
    int numMeshComps = SpaceDim;
    for (int dir=0; dir<SpaceDim; dir++)
      {
        sprintf(field_name, "grid_%c", coords[dir]);
        vectNames.push_back(field_name);
      }
    
    for (int i = 0; i < numMeshComps; ++i)
      {
        sprintf(comp_name, "component_%d", i);
        header.m_string[comp_name] = vectNames[i];
      }
    header.m_int["num_components"] = numMeshComps;
    
    header.writeToFile(handle);
    
    //
    // write the cell centered grid data
    //
    if (write_cell)
      {     
        const std::string group1Name= std::string("cell_centered_grid");
        handle.setGroup(group1Name);
        
        header.clear();
        header.m_int["num_components"] = cellData.nComp();
        header.m_box["prob_domain"] = domain.domainBox(); 
        header.writeToFile(handle);
        
        
        write(handle, cellData.boxLayout());
        write(handle, cellData, "data", cellData.ghostVect());
      }
    
    //
    // write the face centered grid data
    //
    if (write_face)
      {
        //const LevelData<FluxBox>& gridOnFaces(m_mesh.getXfc());
        
        const std::string group2Name = std::string("face_centered_grid");
        handle.setGroup(group2Name);
        
        header.clear();
        header.m_int["is_fluxbox"] = 1;
        header.m_int["num_components"] = faceData.nComp();
        header.m_box["prob_domain"] = domain.domainBox(); 
        header.writeToFile(handle);
        
        write(handle, faceData.boxLayout());
        write(handle, faceData, "data", faceData.ghostVect());
      }
    
    // edge data
    if (write_edge)
      {
        //const LevelData<EdgeDataBox>& gridOnEdges(m_mesh.getXec());
        
        const std::string group3Name = std::string("edge_centered_grid");
        handle.setGroup(group3Name);
        
        header.clear();
        header.m_int["is_edgebox"] = 1;
        header.m_int["num_components"] = edgeData.nComp();
        header.m_box["prob_domain"] = domain.domainBox(); 
        header.writeToFile(handle);
        
        write(handle, edgeData.boxLayout());
        write(handle, edgeData, "edgeData", edgeData.ghostVect());
      }
    
    // close the handle
    handle.close();
  }

  // now test reading
  {
    HDF5Handle handle;
    status = handle.open( plotFileName.c_str(), HDF5Handle::OPEN_RDONLY );
    if (status != 0)
      {
        if ( verbose )
          pout() << "File open failed "<<status<<endl;
        return status;
      }        

    handle.setGroup("/");

    //
    // read the cell centered grid data
    //
    if (write_cell)
      {     
        const std::string group1Name= std::string("cell_centered_grid");
        handle.setGroup(group1Name);
        
        Vector<Box> readBoxes;
        
        status = read(handle, readBoxes);
        
        if (status != 0)
          {
            if ( verbose )
              pout() << indent2 << "box read failed "<<status<<endl;
            return status;
          }
        readBoxes.sort();
        Vector<int> assign;
        status = LoadBalance(assign, readBoxes);
        if (status != 0)
          {
            if ( verbose )
              pout() << indent2 << "BoxLayout LoadBalance failed "<<status<<endl;
            return status;
          }
        DisjointBoxLayout readDBL(readBoxes, assign, domain);

        if (verbose)
          {
            pout() << "read in cell-centered boxes" << endl;
            DataIterator readDit = readDBL.dataIterator();        
            for (readDit.begin(); readDit.ok(); ++readDit)
              {
                pout() << readDBL[readDit] << endl;
              }
          }

        LevelData<FArrayBox> readCellData;
        // start with assumption that grids are correct
        read(handle, readCellData, "data", grids, comps, true);

        // now check that solutions match
        bool matchedData = true;
        Real tol = 1e-10;
        DataIterator dit = grids.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
          {
            readCellData[dit].minus(cellData[dit]);
            BoxIterator bit(readCellData[dit].box());
            for (bit.begin(); bit.ok(); ++bit)
              {
                IntVect iv = bit();
                if (fabs(readCellData[dit](iv,0)) > tol)
                  {
                    if ( verbose )
                      pout() << "cell Data read-write mismatch at cell "
                             << iv << endl;
                    matchedData = false;
                  }
              } // end loop over cells
          } // end loop over grids

        if (!matchedData)
          {
            status = 10;
            if ( verbose )
              pout() << indent2 << "Cell-centered data read test failed "<<endl;
            return status;
          }
      } // end cell-centered data test
    
    //
    // read the face centered grid data
    //
    if (write_face)
      {
        //const LevelData<FluxBox>& gridOnFaces(m_mesh.getXfc());
        
        const std::string group2Name = std::string("face_centered_grid");
        handle.setGroup(group2Name);

        Vector<Box> readBoxes;
        
        status = read(handle, readBoxes);
        
        if (status != 0)
          {
            if ( verbose )
              pout() << indent2 << "box read failed "<<status<<endl;
            return status;
          }
        readBoxes.sort();
        Vector<int> assign;
        status = LoadBalance(assign, readBoxes);
        if (status != 0)
          {
            if ( verbose )
              pout() << indent2 << "BoxLayout LoadBalance failed "<<status<<endl;
            return status;
          }
        DisjointBoxLayout readDBL(readBoxes, assign, domain);

        if (verbose)
          {
            pout() << "read in face-centered boxes" << endl;
            DataIterator readDit = readDBL.dataIterator();        
            for (readDit.begin(); readDit.ok(); ++readDit)
              {
                pout() << readDBL[readDit] << endl;
              }
          }

        LevelData<FluxBox> readFaceData;
        // start with assumption that grids are correct
        read(handle, readFaceData, "data", grids, comps, true);

        // now check that solutions match
        bool matchedData = true;
        Real tol = 1e-10;
        DataIterator dit = grids.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
          {
            for (int dir=0; dir<SpaceDim; dir++)
              {
                (readFaceData[dit])[dir].minus((faceData[dit])[dir]);
                BoxIterator bit((readFaceData[dit])[dir].box());
                for (bit.begin(); bit.ok(); ++bit)
                  {
                    IntVect iv = bit();
                    if (fabs((readFaceData[dit])[dir](iv,0)) > tol)
                      {
                        if ( verbose )
                          pout() << dir << "-face Data read-write mismatch at face "
                                 << iv << endl;
                        matchedData = false;
                      } // end if no match
                  } // end loop over faces
              } // end loop over face directions
          } // end loop over grids
        
        if (!matchedData)
          {
            status = 10;
            if ( verbose )
              pout() << indent2 << "face-centered data read test failed "<<endl;
            return status;
          }
      }

    // edge data
    if (write_edge)
      {
        //const LevelData<EdgeDataBox>& gridOnEdges(m_mesh.getXec());
        
        const std::string group3Name = std::string("edge_centered_grid");
        handle.setGroup(group3Name);

        Vector<Box> readBoxes;
        
        status = read(handle, readBoxes);

        if (status != 0)
          {
            if ( verbose )
              pout() << indent2 << "box read failed "<<status<<endl;
            return status;
          }
        readBoxes.sort();
        Vector<int> assign;
        status = LoadBalance(assign, readBoxes);
        if (status != 0)
          {
            if ( verbose )
              pout() << indent2 << "BoxLayout LoadBalance failed "<<status<<endl;
            return status;
          }
        DisjointBoxLayout readDBL(readBoxes, assign, domain);

        if (verbose)
          {
            pout() << "read in edge-centered boxes" << endl;
            DataIterator readDit = readDBL.dataIterator();        
            for (readDit.begin(); readDit.ok(); ++readDit)
              {
                pout() << readDBL[readDit] << endl;
              }
          }

        LevelData<EdgeDataBox> readEdgeData;
        // start with assumption that grids are correct
        read(handle, readEdgeData, "edgeData", grids, comps, true);

        // now check that solutions match
        bool matchedData = true;
        Real tol = 1e-10;
        DataIterator dit = grids.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
          {
            for (int dir=0; dir<SpaceDim; dir++)
              {
                FArrayBox& thisOrigData = (edgeData[dit])[dir];                
                FArrayBox& thisReadData = (readEdgeData[dit])[dir];
                thisReadData.minus(thisOrigData);
                BoxIterator bit((readEdgeData[dit])[dir].box());
                for (bit.begin(); bit.ok(); ++bit)
                  {
                    IntVect iv = bit();
                    if (fabs(thisReadData(iv,0)) > tol)
                      {
                        if ( verbose )
                          pout() << dir << "-edge Data read-write mismatch at edge "
                                 << iv << endl;
                        matchedData = false;
                      } // end if no match
                  } // end loop over edges
              } // end loop over edge directions
          } // end loop over grids
        
        if (!matchedData)
          {
            status = 10;
            if ( verbose )
              pout() << indent2 << "edge-centered data read test failed "<<endl;
            return status;
          }

      }

    // close the handle

    handle.close();
    
  }
  
#else
  // this test doesn't make any sense w/o HDF5, so don't do anything 
  //MayDay::Error( "plotting only defined with hdf5" );
#endif
  return status;
}


///
// Parse the standard test options (-v -q) out of the command line.
// Stop parsing when a non-option argument is found.
///
void
parseTestOptions( int argc ,char* argv[] )
{
  for ( int i = 1 ; i < argc ; ++i )
  {
    if ( argv[i][0] == '-' ) //if it is an option
    {
      // compare 3 chars to differentiate -x from -xx
      if ( strncmp( argv[i] ,"-v" ,3 ) == 0 )
      {
        verbose = true ;
        // argv[i] = "" ;
      }
      else if ( strncmp( argv[i] ,"-q" ,3 ) == 0 )
      {
        verbose = false ;
        // argv[i] = "" ;
      }
      else
      {
        break ;
      }
    }
  }
  return ;
}
