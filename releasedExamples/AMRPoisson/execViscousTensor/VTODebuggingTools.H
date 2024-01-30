#ifndef ___VTO_Debugging_Tools____
#define ___VTO_Debugging_Tools____
///set this stuff so you can customize your area
//box bounds of area
#define DBG_LOX -1
#define DBG_HIX 4
#define DBG_LOY -1
#define DBG_HIY 8
#define DBG_LOZ 0
#define DBG_HIZ 0

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
  typedef   IntVect                      ch_iv;
  typedef   Box                          ch_box;
  typedef   FArrayBox                    ch_fab;
  typedef   DataIterator                 ch_dit;
  typedef   DisjointBoxLayout            ch_dbl;
  typedef   LevelData<ch_fab>            ch_ldf;
                                                  
  
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
      ch_iv  hidata    = databox.bigEnd();
      ch_iv  lodata    = databox.smallEnd();

      ///to get around weird box types
      ch_box printBox(lodata, hidata);
      printBox &= s_areaBox;
      ch_iv lodeb = printBox.smallEnd();
      ch_iv hideb = printBox.bigEnd();
      
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

  void dumpLayoutCh(ch_dbl* a_data_ptr)
  {
    if(a_data_ptr != NULL)
    {
      a_data_ptr->print();
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
  PrChDebug::dumpLDLayoutCh(NULL);
  PrChDebug::dumpAreaCh(NULL);
  PrChDebug::dumpLDAreaCh(NULL);
}

#undef DBG_LOX 
#undef DBG_LOY 
#undef DBG_LOZ 
#undef DBG_HIX 
#undef DBG_HIY 
#undef DBG_HIZ 


#endif
