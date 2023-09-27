#ifndef __DebugFunctionssss__
#define __DebugFunctionssss__
#define PI 3.141592653589793
#include "NamespaceHeader.H"
#if CH_SPACEDIM==2
//for area functions
//cell centered box
IntVect areaLo( 0, 0);
IntVect areaHi( 3, 3);
//which side of cell gets printed in face dumps
Side::LoHiSide faceSideX = Side::Lo;
Side::LoHiSide faceSideY = Side::Hi;
Box areaBox(areaLo, areaHi);

void 
dumpFaceArea(EBFaceFAB* dataPtr)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);


    EBFaceFAB & data = *dataPtr;
    EBGraph graph = data.getEBISBox().getEBGraph();
    Box databox = data.getCellRegion();
    Box area = areaBox & graph.getDomain();
    cout << "data region contains:" << endl;
    IntVect lodeb = area.smallEnd();
    IntVect hideb = area.bigEnd();
    Side::LoHiSide faceSide = faceSideX;
    if(data.direction() == 1)
    {
      faceSide = faceSideY;
    }
    for(int icomp = 0; icomp < data.nComp(); icomp++)
    {
      cout << "data for component = " << icomp << endl;
      for(int j = hideb[1]; j >= lodeb[1]; j--)
      {
        for(int i = lodeb[0]; i <= hideb[0]; i++)
        {
          IntVect pt(i,j);
          if(databox.contains(pt) && areaBox.contains(pt))
          {
            Vector<VolIndex> vofs = graph.getVoFs(pt);
            cout << pt << ":";
            bool printLines = false;
            if(vofs.size() >  0)
            {
              Vector< FaceIndex>faces = graph.getFaces(vofs[0], data.direction(), faceSide);
              if(faces.size() == 0)
              {
                printLines = true;
              }
              else
              {
                cout  << data(faces[0], icomp) << "  ";
              }
            }
            else
            {
              printLines = true;
            }
            if(printLines)
            {
              cout  << "------------" << "  ";
            }
          }
        }//i loop
        cout << endl;
      } //j loop
    }
  } //if dataptr != null
}
void 
dumpFluxArea(EBFluxFAB* dataPtr)
{
  if(dataPtr != NULL)
  {
    auto& fab = *dataPtr;
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      cout <<  "data for face direction = " << idir << endl;
      dumpFaceArea(&fab[idir]);
    }
  }
}
void 
fullDumpFace(EBFaceFAB* dataPtr)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);


    EBFaceFAB & data = *dataPtr;
    EBGraph graph = data.getEBISBox().getEBGraph();
    Box databox = data.getCellRegion();
    Box area = databox & graph.getDomain();
    cout << "data region contains:" << endl;
    IntVect lodeb = area.smallEnd();
    IntVect hideb = area.bigEnd();
    Side::LoHiSide faceSide = faceSideX;
    if(data.direction() == 1)
    {
      faceSide = faceSideY;
    }
    for(int j = hideb[1]; j >= lodeb[1]; j--)
    {
      for(int i = lodeb[0]; i <= hideb[0]; i++)
      {
        IntVect pt(i,j);
        if(databox.contains(pt))
        {
          Vector<VolIndex> vofs = graph.getVoFs(pt);
          cout << pt << ":";
          bool printLines = false;
          if(vofs.size() >  0)
          {
            Vector< FaceIndex>faces = graph.getFaces(vofs[0], data.direction(), faceSide);
            if(faces.size() == 0)
            {
              printLines = true;
            }
            else
            {
              cout  << data(faces[0], 0) << "  ";
            }
          }
          else
          {
            printLines = true;
          }
          if(printLines)
          {
            cout  << "------------" << "  ";
          }
        }
      } //i loop
      cout << endl;
    } //j loop
  } //dataptr != null
}

void 
fullDumpCell(EBCellFAB* dataPtr)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    EBGraph graph = dataPtr->getEBISBox().getEBGraph();
    Box area = dataPtr->box();
    area &= graph.getDomain();

    EBCellFAB & data = *dataPtr;
    Box databox = dataPtr->box();
    IntVect hideb = databox.bigEnd();
    IntVect lodeb = databox.smallEnd();
    for(int ivar = 0; ivar < data.nComp(); ivar++)
    {
      cout << "data region contains, for variable "<< ivar << endl;
      for(int j = hideb[1]; j >= lodeb[1]; j--)
      {
        for(int i = lodeb[0]; i <= hideb[0]; i++)
        {
          IntVect pt(i,j);
          if(databox.contains(pt))
          {
            Vector<VolIndex > vofs = graph.getVoFs(pt);
            cout << pt << ":";
            if(vofs.size() >  0)
            {
              cout  << data(vofs[ivar], 0) << "  ";
              if(data(vofs[0], ivar)  >= 0)
              {
                cout << " ";
              }
            }
            else
            {
              cout  << "---------------" << "  ";
            }
          }
        }
        cout << endl;
      }
    }
  }
}
void 
dumpCellArea(EBCellFAB* dataPtr)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(8)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    EBGraph graph = dataPtr->getEBISBox().getEBGraph();
    Box area = areaBox;
    area &= graph.getDomain();

    EBCellFAB & data = *dataPtr;
    Box databox = dataPtr->box();
    Box printBox = databox & areaBox;
    IntVect hideb = printBox.bigEnd();
    IntVect lodeb = printBox.smallEnd();
    for(int ivar = 0; ivar < data.nComp(); ivar++)
    {
      cout << "data region contains, for variable "<< ivar << endl;
      for(int j = hideb[1]; j >= lodeb[1]; j--)
      {
        for(int i = lodeb[0]; i <= hideb[0]; i++)
        {
          IntVect pt(i,j);
          Vector<VolIndex > vofs = graph.getVoFs(pt);
          cout << pt << ":";
          if(vofs.size() >  0)
          {
            cout  << data(vofs[0], ivar) << "  ";
            if(data(vofs[0], ivar)  >= 0)
            {
              cout << " ";
            }// 
          } //end  if there is data
          else
          {
            cout  << "---------------" << "  ";
          }
        }
        cout << endl;
      }
    }
  }
}
void 
dumpLDCellArea(LevelData<EBCellFAB>* dataPtr)
{
  if(dataPtr != NULL)
  {
    LevelData<EBCellFAB>& data = *dataPtr;
    DataIterator dit = data.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      cout << "ibox = " << ibox << ":" << endl;
      dumpCellArea(&data[dit[ibox]]);
    }
  }
}
void 
dumpLDFaceArea(LevelData<EBFaceFAB>* dataPtr)
{
  if(dataPtr != NULL)
  {
    LevelData<EBFaceFAB>& data = *dataPtr;
    DataIterator dit = data.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      cout << "ibox = " << ibox << ":" << endl;
      dumpFaceArea(&data[dit[ibox]]);
    }
  }
}
void 
dumpLDFluxArea(LevelData<EBFluxFAB>* dataPtr)
{
  if(dataPtr != NULL)
  {
    LevelData<EBFluxFAB>& data = *dataPtr;
    DataIterator dit = data.dataIterator();
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      cout << "ibox = " << ibox << ":" << endl;
      dumpFluxArea(&data[dit[ibox]]);
    }
  }
}

inline void func_to_get_functions_into_namespace()
{
  fullDumpCell(NULL);
  fullDumpFace(NULL);
  dumpCellArea(NULL);
  dumpFaceArea(NULL);
  dumpFluxArea(NULL);
  dumpLDCellArea(NULL);
  dumpLDFaceArea(NULL);
  dumpLDFluxArea(NULL);
}
#endif

#include "NamespaceFooter.H"

#endif
