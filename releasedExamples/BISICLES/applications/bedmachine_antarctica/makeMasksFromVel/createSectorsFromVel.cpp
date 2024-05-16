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
// read an AMR hierarchy, loop over subsets of its domain and set "sector"
// values based on the velocity based on relevance for a given catchment,
// and write the result
// to a Chombo hdf5 file.
//
// input file needs the following:
//
// main.subregion_file -- contains information about which regions to test 
//                        and the velocity-field parameters to test against
// main.velocity_file -- contains directionally-resolved velocity field
// (to decide if flow is into a given basin or not)
// main.xVel_name -- name of x-velocity in velocity_file
// main.yVel_name -- name of y-velocity in velocity_file
// main.mask_file -- contains mask to be modified
// main.mask_name -- name of mask value in data files
// main.output_file -- file into which modified data is written
//
//===========================================================================
#include <iostream>
#include <fstream>
#include "AMRIO.H"
#include "LoadBalance.H"
#include "ParmParse.H"
#include "FillFromReference.H"

#include "NamespaceHeader.H"

bool verbose = true;
// number to use as a test to see of sector has already been set
int maskTestVal = 100;


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

    if(argc < 2) 
      { 
	std::cerr << " usage: " << argv[0] << " <input_file>" << std::endl; 
	exit(0); 
      }

    char* in_file = argv[1];
    ParmParse pp(argc-2,argv+2,NULL,in_file);
    ParmParse ppMain("main");

    // parse inputs. 

    string subregions_file;
    string velocity_file;
    string xVel_name;
    string yVel_name;
    string data_file;
    string mask_name;
    int sector_val;
    string output_file;


    ppMain.get("subregion_file", subregions_file);
    ppMain.get("velocity_file",  velocity_file);
    ppMain.get("xVel_name",  xVel_name);
    ppMain.get("yVel_name",  yVel_name);
    ppMain.get("data_file",  data_file);
    ppMain.get("mask_name",  mask_name);
    ppMain.get("sector_val",  sector_val);
    ppMain.get("output_file",  output_file);

    // read velocity file...
    if (verbose) 
      {
        pout() << "reading " << velocity_file << "..." << endl;
      }


    Vector<std::string> velfile_names;
    Vector<LevelData<FArrayBox>* > velfile_data;
    Vector<DisjointBoxLayout> velfile_grids;
    Vector<int> velfile_ratio;
    Real velfile_crseDx = 0.0, velfile_dt = 0.0, velfile_time = 0.0;
    Box velfile_crseBox;
    int velfile_numLevels;
    int status = ReadAMRHierarchyHDF5
      (velocity_file,velfile_grids,velfile_data,velfile_names,
       velfile_crseBox,velfile_crseDx,velfile_dt,velfile_time,
       velfile_ratio,velfile_numLevels);

    // will eventually need a Vector of dx for this
    Vector<RealVect> velfile_dx(velfile_numLevels);
    {
      RealVect levelDx = velfile_crseDx*RealVect::Unit;
      for (int lev=0; lev < velfile_numLevels; lev++)
        {
          velfile_dx[lev] = levelDx;
          if (lev < velfile_ratio.size())
            {
              levelDx /= velfile_ratio[lev];
            }
        }
    }

    if (status != 0)
      {
	MayDay::Error("failed to read AMR hierarchy");
      }
    
    int xVelComp = -1;
    int yVelComp = -1;




    for (int comp=0; comp<velfile_names.size(); comp++)
      {
        if (velfile_names[comp] == xVel_name) 
          {
            xVelComp = comp;
          }
        else if (velfile_names[comp] == yVel_name)
          {
            yVelComp = comp;
          }
      }
    CH_assert (xVelComp >= 0);
    CH_assert (yVelComp >= 0);

    if (verbose)
      {
        pout() << "... done." << endl;
        pout() << endl;
      }

    // read data file...

    if (verbose) 
      {
        pout() << "reading " << data_file << "..." << endl;
      }


    Vector<std::string> datafile_names;
    Vector<LevelData<FArrayBox>* > datafile_data;
    Vector<DisjointBoxLayout> datafile_grids;
    Vector<int> datafile_ratio;
    Real datafile_crseDx = 0.0, datafile_dt = 0.0, datafile_time = 0.0;
    Box datafile_crseBox;
    int datafile_numLevels;
    status = ReadAMRHierarchyHDF5
      (data_file,datafile_grids,datafile_data,datafile_names,
       datafile_crseBox,datafile_crseDx,datafile_dt,datafile_time,
       datafile_ratio,datafile_numLevels);

    if (status != 0)
      {
	MayDay::Error("failed to read AMR hierarchy");
      }
    
    int data_maskComp = -1;

    for (int comp=0; comp<datafile_names.size(); comp++)
      {
        if (datafile_names[comp] == mask_name) 
          {
            data_maskComp = comp;
          }
      }
    
    CH_assert (data_maskComp >= 0);

    if (verbose)
      {
        pout() << "... done." << endl;
        pout() << endl;
        pout() << "modifying sectors..." << endl;
      }

    
    // set up test regions from file
    // format of file: 
    // number of subregions
    // then for each subregion (on its own line):
    // low x index, low y index, high x index, high y index, velocity test direction, sign
    // where the subregion (on the coarsest level) is 
    // (low_x, low_y)->(high_x,high_y)
    // velocity test direction is 0 if we're testing x-vel, 
    //                            1 if we're testing y-vel
    //                            2 if we just want to grab the entire block
    //                            3 if we want the entire block and want to force it (so clobber existing tags)
    //                            4 if like 2, but upper- or lower-triangular
    // sign is +1 if we're masking positive velocities,
    //         -1 if we're masking negative velocities
    //

    int numTestRegions = -1;
    Vector<Box> vectSubregion;
    // this will indicate the velocity test direction
    // if 0, we're testing based on xvel, if 1, then yvel
    Vector<int> vectTestDir;
    Vector<int> vectSign;

    if (procID() == uniqueProc(SerialTask::compute))
      {
        ifstream is(subregions_file.c_str(), ios::in);
        if (is.fail())
          {
            MayDay::Error("Cannot open subregion file");
          }
        is >> numTestRegions;
        
        // advance pointer in file
        while (is.get() != '\n');
        
        vectSubregion.resize(numTestRegions);
        // this will indicate the velocity test direction
        // if 0, we're testing based on xvel, if 1, then yvel
        vectTestDir.resize(numTestRegions);
        vectSign.resize(numTestRegions);
        
        for (int i=0; i< numTestRegions; i++)
          {
            int lo_x, lo_y, hi_x, hi_y;
            is >> lo_x;
            is >> lo_y;
            is >> hi_x;
            is >> hi_y;
            is >> vectTestDir[i];
            is >> vectSign[i];

            // advance pointer in file
            while (is.get() != '\n');
            
            // set up subregion box 
            // (note that this is 2D-specific)
            IntVect loVect(lo_x,lo_y);
            IntVect hiVect(hi_x,hi_y);
            vectSubregion[i] = Box(loVect,hiVect);
          } // end loop over subregion read

      } // end if serial proc

    // broadcast results
    broadcast(vectSubregion, uniqueProc(SerialTask::compute));
    broadcast(vectTestDir, uniqueProc(SerialTask::compute));
    broadcast(vectSign, uniqueProc(SerialTask::compute));


    for (int i=0; i<numTestRegions; i++)
      {
        Box subregion = vectSubregion[i];
        int testVelComp = -1;
        if (vectTestDir[i] == 0)
          {
            testVelComp = 0;
          }
        else if (vectTestDir[i] == 1)
          {
            testVelComp =1;            
          }
        else if (vectTestDir[i] == 2)
          {
            testVelComp =2;            
          }
        else 
          {
            testVelComp =vectTestDir[i];            
          }        
          
        RealVect dxLev = datafile_crseDx*RealVect::Unit;
        for (int lev=0; lev<datafile_numLevels; lev++)
          {
            LevelData<FArrayBox>& dataLev = *datafile_data[lev];
            
            // this is going to work a whole lot better if everything is 
            // on the same grids and can use the same iterators...
            // in an ideal world, we could assume that the velocity and 
            // data are on the same AMR hierarchies. 
            // Sadly, we don't live in such a world, and so we need
            // to define the velocity on the same grids as the data using 
            // FillFromReference.

            LevelData<FArrayBox> velTemp(datafile_grids[lev],
                                         velfile_data[0]->nComp());

            flattenCellData(velTemp, dxLev, velfile_data,
                            velfile_dx, false);
            

            
            LevelData<FArrayBox> velData(datafile_grids[lev],
                                         SpaceDim);
            { // gratuitous use of scoping to keep things neat
              
              // first, make a version of the entire vel

              // don't assume that yVelComp = xVelComp+1
              Interval xVelSrcInt(xVelComp,xVelComp);
              Interval xVelDestInt(0,0);
              velTemp.copyTo(xVelSrcInt, velData, xVelDestInt);
              
              Interval yVelSrcInt(yVelComp,yVelComp);
              Interval yVelDestInt(1,1);
              velTemp.copyTo(yVelSrcInt, velData, yVelDestInt);
            }

            
            DataIterator dit = dataLev.dataIterator();
            for (dit.begin(); dit.ok(); ++dit)
              {
                FArrayBox& thisData = dataLev[dit];
                FArrayBox& thisVel = velData[dit];
                Box intersectBox = thisData.box();
                intersectBox &= subregion;
                
                if (!intersectBox.isEmpty())
                  {
                    BoxIterator bit(intersectBox);
                    for (bit.begin(); bit.ok(); ++bit)
                      {
                        IntVect iv = bit();
                        if (testVelComp == 3)
                          {
                            // mark as part of this sector regardless of
                            // what's already there
                            thisData(iv,data_maskComp) = sector_val;
                          }                            
                        // don't poach from other sectors
                        else if (thisData(iv, data_maskComp) > maskTestVal)
                          {
                            if (testVelComp == 2)
                              {
                                // mark as part of this sector
                                thisData(iv,data_maskComp) = sector_val;
                              }
                            else if (testVelComp == 4)
                              {
                                int ilo = subregion.smallEnd()[0];
                                int ihi = subregion.bigEnd()[0];
                                int jlo = subregion.smallEnd()[1];
                                int jhi = subregion.bigEnd()[1];
                                if (vectSign[i]*((jlo-iv[1])*(ihi-ilo)+(jhi-jlo)*(iv[0]-ilo)) > 0)
                                  {
                                    thisData(iv,data_maskComp) = sector_val;
                                  }
                              }
                            else if (testVelComp == 5)
                              {
                                int ilo = subregion.smallEnd()[0];
                                int ihi = subregion.bigEnd()[0];
                                int jlo = subregion.smallEnd()[1];
                                int jhi = subregion.bigEnd()[1];
                                if (vectSign[i]*((jhi-iv[1])*(ihi-ilo)+(jlo-jhi)*(iv[0]-ilo)) > 0)
                                  {
                                    thisData(iv,data_maskComp) = sector_val;
                                  }
                              }                            
                            else if (vectSign[i]*thisVel(iv,testVelComp) > 0.0)
                              {
                                // mark as part of this sector
                                thisData(iv,data_maskComp) = sector_val;
                              }
                          } // end if this cell hasn't been claimed already
                      } // end loop over cells in subregion
                  } // end if intersect box not empty
              } // end loop over boxes on this level
            
            if (lev < datafile_ratio.size()) 
              {
                subregion.refine(datafile_ratio[lev]);
                dxLev /= datafile_ratio[lev];
              }
            
          } // end loop over levels
      } // end loop over xvel tests

    
    {
      pout() << "... done." << endl;
      pout() << "writing " << output_file << "..." << endl;
    }

    
    WriteAMRHierarchyHDF5(output_file, datafile_grids, datafile_data, 
                          datafile_names, datafile_crseBox, datafile_crseDx, 
                          datafile_dt, datafile_time, datafile_ratio,
                          datafile_numLevels);
			  
    

    if (verbose)
      {
        pout() << "... done." << endl;
      }
  
    for (int lev=0; lev<datafile_data.size(); lev++)
      {
        if (datafile_data[lev] != NULL)
          {
            delete datafile_data[lev];
            datafile_data[lev] = NULL;
          }
      }

    for (int lev=0; lev<velfile_data.size(); lev++)
      {
        if (velfile_data[lev] != NULL)
          {
            delete velfile_data[lev];
            velfile_data[lev] = NULL;
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
