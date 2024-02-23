#ifndef ____Debugging_Tools____
#define ____Debugging_Tools____
#define DBG_LOX -1
#define DBG_HIX 3
#define DBG_LOY -1
#define DBG_HIY 3
#define DBG_LOZ 0
#define DBG_HIZ 0
#define VALID_ONLY false

/**
//   for dumpFracPr and dumpFracCh
//   Makes box = defined fraction of the domain.
// set the numerators and denominators of fraction you want displayed
// This:
//#define DBG_NUMER_LO 1
//#define DBG_DENOM_LO 8
//#define DBG_NUMER_HI 1
//#define DBG_DENOM_HI 4
//  Will display the fraction of the domain that is between 1/8 and 1/4 of
//   of the domain area.  (so near the lower left corner)
//  This may go down in history as the worst API ever but it is helpful when you
//  need to keep an eye on a particular region of the domain when AMR is involved.
 **/
#define DBG_NUMER_LO_X 1
#define DBG_DENOM_LO_X 8
#define DBG_NUMER_HI_X 1
#define DBG_DENOM_HI_X 4

#define DBG_NUMER_LO_Y 0
#define DBG_DENOM_LO_Y 1
#define DBG_NUMER_HI_Y 1
#define DBG_DENOM_HI_Y 8

#define DBG_NUMER_LO_Z 0 0
#define DBG_DENOM_LO_Z 01
#define DBG_NUMER_HI_Z 00
#define DBG_DENOM_HI_Z 0 1

#include "Proto.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "DisjointBoxLayout.H"



/// Old school Debugging utilities 
/**
  // If you want to print areas of data to screen during a debug session,
  // his is the way you do it.  All of these functions can be called from gdb.
  // TLDR:
  // This is a set  of crude but effective tools to  print out stuff from gdb.   
  // Send pointers to the data to these functions and the data gets printed out neatly.    
  // The user is meant to hack the the bounds according to her needs.   
  // The macros above are the primary buttons to modify.    
  // 
  // Set the bounds (DBG_LOX and so on) to customize which area gets printed.
  // Set VALID_ONLY to TRUE if you only want valid data printed when a 
  //     level gets printed out (sometimes it is help ful to see ghost data, other times not so much).
  // 
  // This command:
  // (gdb) p PrChDebug::dumpLDAreaPr(&a_cor)  
  // yeilds:
  // ibox = 0:
  // data_box = [(-2,-2),( 3, 3)], area_box = [( 0, 0),( 3, 3)]
  // data for variable 0:
  // ( 0, 3): 4.86993301e-05  ( 1, 3): 3.13317884e-04  ( 2, 3): 1.72362925e-04  ( 3, 3): 6.44824962e-05  
  // ( 0, 2): 3.13287135e-04  ( 1, 2): 5.52541534e-04  ( 2, 2): 8.82170603e-04  ( 3, 2): 1.72362925e-04  
  // ( 0, 1): 1.72411307e-04  ( 1, 1): 8.82090868e-04  ( 2, 1): 5.52541534e-04  ( 3, 1): 3.13317884e-04  
  // ( 0, 0): 6.46298245e-05  ( 1, 0): 1.72411307e-04  ( 2, 0): 3.13287135e-04  ( 3, 0): 4.86993301e-05  
  // and so on.
**/
namespace PrChDebug
{
  typedef   Chombo::IntVect                      ch_iv;
  typedef   Chombo::Box                          ch_box;
  typedef   Chombo::FArrayBox                    ch_fab;
  typedef   Chombo::DataIterator                 ch_dit;
  typedef   Chombo::DisjointBoxLayout            ch_dbl;
  typedef   Chombo::LevelData<ch_fab>            ch_ldf;
                                                  
  typedef    Proto::Point                        pr_pt;
  typedef    Proto::Box                          pr_box;
  typedef    Proto::BoxData<double, 1>           pr_fab_1;
  typedef    Proto::BoxData<double, DIM>         pr_fab_vec;
  typedef    Proto::BoxData<double, 3>           pr_fab_three;
  typedef    Proto::BoxData<double, 3*DIM>       pr_fab_threedim;
  typedef    Proto::BoxData<double, DIM*DIM>     pr_fab_ten;
  typedef    Proto::DisjointBoxLayout            pr_dbl;
  typedef    Proto::LevelBoxData<double,1>       pr_lbd_sca;
  typedef    Proto::LevelBoxData<double,DIM>     pr_lbd_vec;
  typedef    Proto::LevelBoxData<double,3  >     pr_lbd_three;
  typedef    Proto::LevelBoxData<double,3*DIM>   pr_lbd_threedim;
  typedef    Proto::LevelBoxData<double,DIM*DIM> pr_lbd_ten;
  typedef    Proto::Stencil<double>              pr_sten;
  
  void
  dumpSten(pr_sten* a_sten_ptr)
  {
    if(a_sten_ptr != NULL)
    {
      a_sten_ptr->print();
    }
  }

  void 
  dumpCh(ch_fab* a_dataPtr, ch_box a_areaBox)
  {

    if(a_dataPtr != NULL)
    {
      using std::cout;
      using std::end;
      cout    << setprecision(8)
              << setiosflags(ios::showpoint)
              << setiosflags(ios::scientific);

      ch_box area = a_areaBox;

      ch_fab & data    =*a_dataPtr;
      ch_box databox   = a_dataPtr->box();
      ch_iv  hidata    =   databox.bigEnd();
      ch_iv  lodata    =   databox.smallEnd();

      ///to get around weird box types
      ch_box printBox(lodata, hidata);
      printBox &= a_areaBox;
      ch_iv lodeb = printBox.smallEnd();
      ch_iv hideb = printBox.bigEnd();
      
      for(int ivar = 0; ivar < data.nComp(); ivar++)
      {
        cout << "data_box = " << databox << ", area_box = " << a_areaBox <<  endl;
        cout << "data for variable "<< ivar << ":" << endl;
#if DIM==3        
        for(int k = hideb[2]; k >= lodeb[2]; k--)
        {
#endif          
        for(int j = hideb[1]; j >= lodeb[1]; j--)
        {
          for(int i = lodeb[0]; i <= hideb[0]; i++)
          {
            ch_iv pt(D_DECL(i,j,k));
            cout << pt << ":";
            double data_value = data(pt, ivar);
            //keeps positive and negative numbers lined up
            if(    data_value >= 0)
            {
              cout << " "; 
            }
            // 
            cout  << data_value << "  ";
          } //end loop in x dir
          cout << endl;
        } //end loop in  y dir
#if DIM==3
        } //end loop in z direction
#endif
      }// end loop over components
    } //end if not null
  } //end function dumpAreaCh

  ///
  void 
  dumpAreaCh(ch_fab* dataPtr)
  {
    static const ch_iv       s_areaLo(D_DECL( DBG_LOX, DBG_LOY, DBG_LOZ));
    static const ch_iv       s_areaHi(D_DECL( DBG_HIX, DBG_HIY, DBG_HIZ));
    static const ch_box     s_areaBox(s_areaLo, s_areaHi);
    dumpCh(dataPtr, s_areaBox);
  }
  
  ///
  void 
  dumpPr(pr_fab_1* dataPtr, pr_box a_areaBox)
  {
    if(dataPtr != NULL)
    {
      using std::cout;
      using std::end;

      cout    << setprecision(8)
              << setiosflags(ios::showpoint)
              << setiosflags(ios::scientific);

      pr_box area = a_areaBox;

      auto & data   = *dataPtr;
      pr_box databox  = dataPtr->box();
      pr_box printBox = databox & a_areaBox;
      pr_pt  hideb    = printBox.high();
      pr_pt  lodeb    = printBox.low();
      
      int ivar = 0;
      {
        cout << "data_box = " << databox << ", area_box = " << a_areaBox <<  endl;
        cout << "data for variable "<< ivar << ":" << endl;
#if DIM==3        
        for(int k = hideb[2]; k >= lodeb[2]; k--)
        {
#endif          
        for(int j = hideb[1]; j >= lodeb[1]; j--)
        {
          for(int i = lodeb[0]; i <= hideb[0]; i++)
          {
            pr_pt pt(D_DECL(i,j,k));
            cout << pt << ":";
            double data_value = data(pt, ivar);
            //keeps positive and negative numbers lined up
            if(    data_value >= 0)
            {
              cout << " "; 
            }
            // 
            cout  << data_value << "  ";
          } //end loop in x dir
          cout << endl;
        } //end loop in  y dir
#if DIM==3
        } //end loop in z direction
#endif
      }// end loop over components
    } //end if not null
  } //end function dumpPr

  void 
  dumpAreaPr(pr_fab_1* dataPtr)
  {
    static const pr_pt       s_areaLo(D_DECL( DBG_LOX, DBG_LOY, DBG_LOZ));
    static const pr_pt       s_areaHi(D_DECL( DBG_HIX, DBG_HIY, DBG_HIZ));
    static const pr_box     s_areaBox(s_areaLo, s_areaHi);
    dumpPr(dataPtr, s_areaBox);
  }

  void 
  dumpAreaVecPr(pr_fab_vec* dataPtr)
  {
    static const pr_pt       s_areaLo(D_DECL( DBG_LOX, DBG_LOY, DBG_LOZ));
    static const pr_pt       s_areaHi(D_DECL( DBG_HIX, DBG_HIY, DBG_HIZ));
    static const pr_box     s_areaBox(s_areaLo, s_areaHi);

    if(dataPtr != NULL)
    {
      using std::cout;
      using std::end;

      cout    << setprecision(8)
              << setiosflags(ios::showpoint)
              << setiosflags(ios::scientific);

      pr_box area = s_areaBox;

      auto & data   = *dataPtr;
      pr_box databox  = dataPtr->box();
      pr_box printBox = databox & s_areaBox;
      pr_pt  hideb    = printBox.high();
      pr_pt  lodeb    = printBox.low();
      
      for(int ivar = 0; ivar < DIM; ivar++)
      {
        cout << "data_box = " << databox << ", area_box = " << s_areaBox <<  endl;
        cout << "data for variable "<< ivar << ":" << endl;
#if DIM==3        
        for(int k = hideb[2]; k >= lodeb[2]; k--)
        {
#endif          
          for(int j = hideb[1]; j >= lodeb[1]; j--)
          {
            for(int i = lodeb[0]; i <= hideb[0]; i++)
            {
              pr_pt pt(D_DECL(i,j,k));
              cout << pt << ":";
              double data_value = data(pt, ivar);
              //keeps positive and negative numbers lined up
              if(    data_value >= 0)
              {
                cout << " "; 
              }
              // 
              cout  << data_value << "  ";
            } //end loop in x dir
            cout << endl;
          } //end loop in  y dir
#if DIM==3
        } //end loop in z direction
#endif
      }// end loop over components
      
    } //end if not null
  } //end function dumpAreaVecPr

  void 
  dumpAreaThreePr(pr_fab_three* dataPtr)
  {
    static const pr_pt       s_areaLo(D_DECL( DBG_LOX, DBG_LOY, DBG_LOZ));
    static const pr_pt       s_areaHi(D_DECL( DBG_HIX, DBG_HIY, DBG_HIZ));
    static const pr_box     s_areaBox(s_areaLo, s_areaHi);

    if(dataPtr != NULL)
    {
      using std::cout;
      using std::end;

      cout    << setprecision(8)
              << setiosflags(ios::showpoint)
              << setiosflags(ios::scientific);

      pr_box area = s_areaBox;

      auto & data   = *dataPtr;
      pr_box databox  = dataPtr->box();
      pr_box printBox = databox & s_areaBox;
      pr_pt  hideb    = printBox.high();
      pr_pt  lodeb    = printBox.low();
      
      for(int ivar = 0; ivar < 3; ivar++)
      {
        cout << "data_box = " << databox << ", area_box = " << s_areaBox <<  endl;
        cout << "data for variable "<< ivar << ":" << endl;
#if DIM==3        
        for(int k = hideb[2]; k >= lodeb[2]; k--)
        {
#endif          
          for(int j = hideb[1]; j >= lodeb[1]; j--)
          {
            for(int i = lodeb[0]; i <= hideb[0]; i++)
            {
              pr_pt pt(D_DECL(i,j,k));
              cout << pt << ":";
              double data_value = data(pt, ivar);
              //keeps positive and negative numbers lined up
              if(    data_value >= 0)
              {
                cout << " "; 
              }
              // 
              cout  << data_value << "  ";
            } //end loop in x dir
            cout << endl;
          } //end loop in  y dir
#if DIM==3
        } //end loop in z direction
#endif
      }// end loop over components
      
    } //end if not null
  } //end function dumpAreaThreePr

  void 
  dumpAreaThreeDimPr(pr_fab_threedim* dataPtr)
  {
    static const pr_pt       s_areaLo(D_DECL( DBG_LOX, DBG_LOY, DBG_LOZ));
    static const pr_pt       s_areaHi(D_DECL( DBG_HIX, DBG_HIY, DBG_HIZ));
    static const pr_box     s_areaBox(s_areaLo, s_areaHi);

    if(dataPtr != NULL)
    {
      using std::cout;
      using std::end;

      cout    << setprecision(8)
              << setiosflags(ios::showpoint)
              << setiosflags(ios::scientific);

      pr_box area = s_areaBox;

      auto & data   = *dataPtr;
      pr_box databox  = dataPtr->box();
      pr_box printBox = databox & s_areaBox;
      pr_pt  hideb    = printBox.high();
      pr_pt  lodeb    = printBox.low();
      
      for(int ivar = 0; ivar < 3*DIM; ivar++)
      {
        cout << "data_box = " << databox << ", area_box = " << s_areaBox <<  endl;
        cout << "data for variable "<< ivar << ":" << endl;
#if DIM==3        
        for(int k = hideb[2]; k >= lodeb[2]; k--)
        {
#endif          
          for(int j = hideb[1]; j >= lodeb[1]; j--)
          {
            for(int i = lodeb[0]; i <= hideb[0]; i++)
            {
              pr_pt pt(D_DECL(i,j,k));
              cout << pt << ":";
              double data_value = data(pt, ivar);
              //keeps positive and negative numbers lined up
              if(    data_value >= 0)
              {
                cout << " "; 
              }
              // 
              cout  << data_value << "  ";
            } //end loop in x dir
            cout << endl;
          } //end loop in  y dir
#if DIM==3
        } //end loop in z direction
#endif
      }// end loop over components
      
    } //end if not null
  } //end function dumpAreaThreeDimPr
  
  void 
  dumpAreaTenPr(pr_fab_ten* dataPtr)
  {
    static const pr_pt       s_areaLo(D_DECL( DBG_LOX, DBG_LOY, DBG_LOZ));
    static const pr_pt       s_areaHi(D_DECL( DBG_HIX, DBG_HIY, DBG_HIZ));
    static const pr_box     s_areaBox(s_areaLo, s_areaHi);

    if(dataPtr != NULL)
    {
      using std::cout;
      using std::end;

      cout    << setprecision(8)
              << setiosflags(ios::showpoint)
              << setiosflags(ios::scientific);

      pr_box area = s_areaBox;

      auto & data   = *dataPtr;
      pr_box databox  = dataPtr->box();
      pr_box printBox = databox & s_areaBox;
      pr_pt  hideb    = printBox.high();
      pr_pt  lodeb    = printBox.low();
      
      for(int ivar = 0; ivar < DIM*DIM; ivar++)
      {
        cout << "data_box = " << databox << ", area_box = " << s_areaBox <<  endl;
        cout << "data for variable "<< ivar << ":" << endl;
#if DIM==3        
        for(int k = hideb[2]; k >= lodeb[2]; k--)
        {
#endif          
          for(int j = hideb[1]; j >= lodeb[1]; j--)
          {
            for(int i = lodeb[0]; i <= hideb[0]; i++)
            {
              pr_pt pt(D_DECL(i,j,k));
              cout << pt << ":";
              double data_value = data(pt, ivar);
              //keeps positive and negative numbers lined up
              if(    data_value >= 0)
              {
                cout << " "; 
              }
              // 
              cout  << data_value << "  ";
            } //end loop in x dir
            cout << endl;
          } //end loop in  y dir
#if DIM==3
        } //end loop in z direction
#endif
      }// end loop over components
      
    } //end if not null
  } //end function dumpAreaTenPr
  
  /// print to screen (very neatly) data from the area whereever it is found in chombo data
  void 
  dumpLDAreaCh(ch_ldf* dataPtr)
  {
    if(dataPtr != NULL)
    {
      static const ch_iv       s_areaLo(D_DECL( DBG_LOX, DBG_LOY, DBG_LOZ));
      static const ch_iv       s_areaHi(D_DECL( DBG_HIX, DBG_HIY, DBG_HIZ));
      static const ch_box     s_areaBox(s_areaLo, s_areaHi);

      ch_ldf & data = *dataPtr;
      ch_dbl   grid = data.disjointBoxLayout();
      ch_dit  dit = data.dataIterator();
      for(int ibox = 0; ibox < dit.size(); ibox++)
      {
        ch_box valid = grid[dit[ibox]];
        ch_box valid_area = valid & s_areaBox;
        cout << "ibox = " << ibox << ", valid = " << valid << ", valid_area = :" << valid_area << endl;
        bool valid_only = VALID_ONLY;

        if(valid_only)
        {
          dumpCh(&data[dit[ibox]], valid_area);
        }
        else if(!valid_area.isEmpty())
        {
          dumpCh(&data[dit[ibox]], s_areaBox);
        }
      }
    }
  }

  void 
  dumpFracCh(ch_ldf* dataPtr)
  {
    if(dataPtr != NULL)
    {

      ch_ldf & data = *dataPtr;
      ch_dbl   grid = data.disjointBoxLayout();
      ch_box dom_box = grid.physDomain().domainBox();
      ch_iv areaLo, areaHi;
      ch_iv size_iv = dom_box.bigEnd() + ch_iv::Unit;

      int lo_dir_x =     (DBG_NUMER_LO_X*size_iv[0])/(DBG_DENOM_LO_X) - 1;
      int hi_dir_x = 1 + (DBG_NUMER_HI_X*size_iv[0])/(DBG_DENOM_HI_X);
      areaLo[0] = lo_dir_x;
      areaHi[0] = hi_dir_x;
      
      int lo_dir_y =     (DBG_NUMER_LO_Y*size_iv[1])/(DBG_DENOM_LO_Y) - 1;
      int hi_dir_y = 1 + (DBG_NUMER_HI_Y*size_iv[1])/(DBG_DENOM_HI_Y);
      areaLo[1] = lo_dir_y;
      areaHi[1] = hi_dir_y;
#if DIM==3
      int lo_dir_z =     (DBG_NUMER_LO_Z*size_iv[2])/(DBG_DENOM_LO_Z) - 1;
      int hi_dir_z = 1 + (DBG_NUMER_HI_Z*size_iv[2])/(DBG_DENOM_HI_Z);
      areaLo[2] = lo_dir_z;
      areaHi[2] = hi_dir_z;
#endif      
      ch_box     areaBox(areaLo, areaHi);
      ch_dit  dit = data.dataIterator();
      for(int ibox = 0; ibox < dit.size(); ibox++)
      {
        ch_box valid = grid[dit[ibox]];
        ch_box valid_area = valid & areaBox;
        cout << "ibox = " << ibox << ", valid = " << valid << ", valid_area = :" << valid_area << endl;
        bool valid_only = VALID_ONLY;

        if(valid_only)
        {
          dumpCh(&data[dit[ibox]], valid_area);
        }
        else if(!valid_area.isEmpty())
        {
          dumpCh(&data[dit[ibox]], areaBox);
        }
      }
    }
  }
  
  /// print to screen (very neatly) data from the area whereever it is found in proto data
  void 
  dumpFracPr(pr_lbd_sca* dataPtr)
  {
    if(dataPtr != NULL)
    {
      pr_lbd_sca & data = *dataPtr;
      auto  dit = data.dataIterator();
      pr_pt      areaLo, areaHi;
      pr_pt dom_hi = data.layout().domain().box().high();
      pr_pt size_pt = dom_hi + pr_pt::Ones(1);
      
      int lo_dir_x =     (DBG_NUMER_LO_X*size_pt[0])/(DBG_DENOM_LO_X) - 1;
      int hi_dir_x = 1 + (DBG_NUMER_HI_X*size_pt[0])/(DBG_DENOM_HI_X);
      areaLo[0] = lo_dir_x;
      areaHi[0] = hi_dir_x;
      
      int lo_dir_y =     (DBG_NUMER_LO_Y*size_pt[1])/(DBG_DENOM_LO_Y) - 1;
      int hi_dir_y = 1 + (DBG_NUMER_HI_Y*size_pt[1])/(DBG_DENOM_HI_Y);
      areaLo[1] = lo_dir_y;
      areaHi[1] = hi_dir_y;
#if DIM==3
      int lo_dir_z =     (DBG_NUMER_LO_Z*size_pt[2])/(DBG_DENOM_LO_Z) - 1;
      int hi_dir_z = 1 + (DBG_NUMER_HI_Z*size_pt[2])/(DBG_DENOM_HI_Z);
      areaLo[2] = lo_dir_z;
      areaHi[2] = hi_dir_z;
#endif      
      pr_box     areaBox(areaLo, areaHi);
      for(int ibox = 0; ibox < dit.localSize(); ibox++)
      {
        pr_box valid = data.layout()[dit[ibox]];
        pr_box valid_area  = valid & areaBox;
        cout << "ibox = " << ibox << ", valid = " << valid << ", valid_area = " << valid_area << ":" << endl;
        bool valid_only = VALID_ONLY;
        if(valid_only)
        {
          dumpPr(&data[dit[ibox]], valid_area);
        }
        else if(!valid_area.empty())
        {
          dumpPr(&data[dit[ibox]], areaBox);
        }

      }
    }
  }
  /// print to screen (very neatly) data from the area whereever it is found in proto data
  void 
  dumpLDAreaPr(pr_lbd_sca* dataPtr)
  {
    if(dataPtr != NULL)
    {
      pr_lbd_sca & data = *dataPtr;
      auto  dit = data.dataIterator();
      for(int ibox = 0; ibox < dit.localSize(); ibox++)
      {
        static const pr_pt       s_areaLo(D_DECL( DBG_LOX, DBG_LOY, DBG_LOZ));
        static const pr_pt       s_areaHi(D_DECL( DBG_HIX, DBG_HIY, DBG_HIZ));
        static const pr_box     s_areaBox(s_areaLo, s_areaHi);
        pr_box valid = data.layout()[dit[ibox]];
        pr_box valid_area  = valid & s_areaBox;
        cout << "ibox = " << ibox << ", valid = " << valid << ", valid_area = " << valid_area << ":" << endl;
        bool valid_only = VALID_ONLY;
        if(valid_only)
        {
          dumpPr(&data[dit[ibox]], valid_area);
        }
        else if(!valid_area.empty())
        {
          dumpPr(&data[dit[ibox]], s_areaBox);
        }

      }
    }
  }

  /// print to screen (very neatly) data from the area whereever it is found in proto data
  void 
  dumpLDAreaTenPr(pr_lbd_ten* dataPtr)
  {
    if(dataPtr != NULL)
    {
      pr_lbd_ten & data = *dataPtr;
      auto  dit = data.dataIterator();
      for(int ibox = 0; ibox < dit.localSize(); ibox++)
      {
        cout << "ibox = " << ibox << ":" << endl;
        dumpAreaTenPr(&data[dit[ibox]]);
      }
    }
  }
  /// print to screen (very neatly) data from the area whereever it is found in proto data
  void 
  dumpLDAreaVecPr(pr_lbd_vec* dataPtr)
  {
    if(dataPtr != NULL)
    {
      pr_lbd_vec & data = *dataPtr;
      auto  dit = data.dataIterator();
      for(int ibox = 0; ibox < dit.localSize(); ibox++)
      {
        cout << "ibox = " << ibox << ":" << endl;
        dumpAreaVecPr(&data[dit[ibox]]);
      }
    }
  }
  void 
  dumpLDAreaThreePr(pr_lbd_three* dataPtr)
  {
    if(dataPtr != NULL)
    {
      pr_lbd_three & data = *dataPtr;
      auto  dit = data.dataIterator();
      for(int ibox = 0; ibox < dit.localSize(); ibox++)
      {
        cout << "ibox = " << ibox << ":" << endl;
        dumpAreaThreePr(&data[dit[ibox]]);
      }
    }
  }
  void 
  dumpLDAreaThreeDimPr(pr_lbd_threedim* dataPtr)
  {
    if(dataPtr != NULL)
    {
      pr_lbd_threedim & data = *dataPtr;
      auto  dit = data.dataIterator();
      for(int ibox = 0; ibox < dit.localSize(); ibox++)
      {
        cout << "ibox = " << ibox << ":" << endl;
        dumpAreaThreeDimPr(&data[dit[ibox]]);
      }
    }
  }

  void dumpLayoutCh(ch_dbl* a_data_ptr)
  {
    if(a_data_ptr != NULL)
    {
      a_data_ptr->print();
    }
  }
  void dumpLayoutPr(pr_dbl* a_data_ptr)
  {
    if(a_data_ptr != NULL)
    {
      cout << "pr_dbl size = " << a_data_ptr->size() << endl; 
      if(a_data_ptr != NULL)
      {
        auto dit = a_data_ptr->begin();
        auto dbl = *a_data_ptr;
        for(int ibox =  0; ibox < dit.localSize(); ibox++)
        {
          cout << dbl[dit[ibox]] << endl;
        }
      }
    }
  }
  void dumpLDLayoutPr(pr_lbd_sca* a_data_ptr)
  {
    if(a_data_ptr != NULL)
    {
      auto dit = a_data_ptr->begin();
      auto dbl = a_data_ptr->layout();
      for(int ibox =  0; ibox < dit.localSize(); ibox++)
      {
        cout << dbl[dit[ibox]] << endl;
      }
    }
  }
  void dumpLDLayoutCh(ch_ldf* a_data_ptr)
  {
    if(a_data_ptr != NULL)
    {
      a_data_ptr->disjointBoxLayout().print();
    }
  }
}//end namespace PrChDebug

///
inline void dummy_function_to_get_functions_into_symbol_table()
{
  PrChDebug::dumpLayoutCh(NULL);
  PrChDebug::dumpLayoutPr(NULL);
  PrChDebug::dumpLDLayoutCh(NULL);
  PrChDebug::dumpLDLayoutPr(NULL);
  PrChDebug::dumpAreaCh(NULL);
  PrChDebug::dumpLDAreaCh(NULL);
  PrChDebug::dumpLDAreaPr(NULL);
  PrChDebug::dumpAreaVecPr(NULL);
  PrChDebug::dumpAreaThreePr(NULL);
  PrChDebug::dumpAreaThreeDimPr(NULL);
  PrChDebug::dumpAreaTenPr(NULL);
  PrChDebug::dumpLDAreaVecPr(NULL);
  PrChDebug::dumpLDAreaThreePr(NULL);
  PrChDebug::dumpLDAreaThreeDimPr(NULL);
  PrChDebug::dumpLDAreaTenPr(NULL);
  PrChDebug::dumpSten(NULL);
  PrChDebug::dumpFracCh(NULL);
  PrChDebug::dumpFracPr(NULL);
}

#undef DBG_LOX 
#undef DBG_LOY 
#undef DBG_LOZ 
#undef DBG_HIX 
#undef DBG_HIY 
#undef DBG_HIZ 
#undef VALID_ONLY
#undef DBG_NUMER_LO 
#undef DBG_DENOM_LO 
#undef DBG_NUMER_HI 
#undef DBG_DENOM_HI 


#endif
