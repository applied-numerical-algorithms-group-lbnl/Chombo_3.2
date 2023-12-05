#ifndef ____Debugging_Tools____
#define ____Debugging_Tools____
///set this stuff so you can customize your area
//box bounds of area
#define DBG_LOX 0
#define DBG_LOY 0
#define DBG_LOZ 0
#define DBG_HIX 4
#define DBG_HIY 7
#define DBG_HIZ 0

#include "Proto.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "DisjointBoxLayout.H"



/// Old school Debugging utilities 
/**
  If you want to print areas of data to screen during a debug session,
  this is the way you do it.  All of these functions can be called from gdb.
**/
namespace PrChDebug
{
  typedef   Chombo::IntVect                  ch_iv;
  typedef   Chombo::Box                      ch_box;
  typedef   Chombo::FArrayBox                ch_fab;
  typedef   Chombo::DataIterator             ch_dit;
  typedef   Chombo::DisjointBoxLayout        ch_dbl;
  typedef   Chombo::LevelData<ch_fab>        ch_ldf;
                                              
  typedef    Proto::Point                    pr_pt;
  typedef    Proto::Box                      pr_box;
  typedef    Proto::BoxData<double, 1>       pr_fab_1;
  typedef    Proto::BoxData<double, DIM>     pr_fab_Dim;
  typedef    Proto::DisjointBoxLayout        pr_dbl;
  typedef    Proto::LevelBoxData<double,1>   pr_lbd_1;
  typedef    Proto::LevelBoxData<double,DIM> pr_lbd_DIM;
  typedef    Proto::Stencil<double>          pr_sten;
  
  void
  dumpSten(pr_sten* a_sten_ptr)
  {
    if(a_sten_ptr != NULL)
    {
      a_sten_ptr->print();
    }
  }
  ///
  void 
  dumpAreaCh(ch_fab* dataPtr)
  {
    static const ch_iv       s_areaLo(D_DECL( DBG_LOX, DBG_LOY, DBG_LOZ));
    static const ch_iv       s_areaHi(D_DECL( DBG_HIX, DBG_HIY, DBG_HIZ));
    static const ch_box     s_areaBox(s_areaLo, s_areaHi);

    if(dataPtr != NULL)
    {
      using std::cout;
      using std::end;

      cout    << setprecision(8)
              << setiosflags(ios::showpoint)
              << setiosflags(ios::scientific);

      ch_box area = s_areaBox;

      ch_fab & data   = *dataPtr;
      ch_box databox  = dataPtr->box();
      ch_box printBox = databox & s_areaBox;
      ch_iv  hideb    = printBox.bigEnd();
      ch_iv  lodeb    = printBox.smallEnd();
      
      for(int ivar = 0; ivar < data.nComp(); ivar++)
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
  dumpAreaPr(pr_fab_1* dataPtr)
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
      
      int ivar = 0;
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
  } //end function dumpAreaPr
  
  /// print to screen (very neatly) data from the area whereever it is found in chombo data
  void 
  dumpLDAreaCh(ch_ldf* dataPtr)
  {
    if(dataPtr != NULL)
    {
      ch_ldf & data = *dataPtr;
      ch_dit  dit = data.dataIterator();
      for(int ibox = 0; ibox < dit.size(); ibox++)
      {
        cout << "ibox = " << ibox << ":" << endl;
        dumpAreaCh(&data[dit[ibox]]);
      }
    }
  }

  /// print to screen (very neatly) data from the area whereever it is found in proto data
  void 
  dumpLDAreaPr(pr_lbd_1* dataPtr)
  {
    if(dataPtr != NULL)
    {
      pr_lbd_1 & data = *dataPtr;
      auto  dit = data.dataIterator();
      for(int ibox = 0; ibox < dit.localSize(); ibox++)
      {
        cout << "ibox = " << ibox << ":" << endl;
        dumpAreaPr(&data[dit[ibox]]);
      }
    }
  }
}//end namespace PrChDebug

inline void func_to_get_functions_into_namespace()
{
  PrChDebug::dumpAreaCh(NULL);
  PrChDebug::dumpAreaPr(NULL);
  PrChDebug::dumpLDAreaCh(NULL);
  PrChDebug::dumpLDAreaPr(NULL);
  PrChDebug::dumpSten(NULL);
}

#undef DBG_LOX 
#undef DBG_LOY 
#undef DBG_LOZ 
#undef DBG_HIX 
#undef DBG_HIY 
#undef DBG_HIZ 


#endif
