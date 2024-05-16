#include <stdlib.h>
#include <iostream>
#include <istream> 
#include <fstream> 

using std::cin; 
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::cerr;
using std::ios;

int main(int argc, char* argv[]) 
{
  std::string inFile, outFile;

  inFile = argv[1];
  outFile = argv[2];

  float scale = 1.0;
  if (argc > 3) 
    {
      //cout << argv[3] << endl;
      scale = atoi(argv[3]);
    }


  //cout << "filename" << endl;
  //cin >> inFile;

  //cout << "outfile" << endl;
  //cin >> outFile;

  
  cout << " infile = " << inFile << ", outfile = " << outFile << endl;

  ifstream is(inFile.c_str(), ios::in);
  ofstream os(outFile.c_str(), ios::out);

  if (is.fail())
    {
      cerr << "Cannot open input file";
      return 1;
    }

  int globalIter = 0;
  int iter = -100; 
  double val;
  int lastIter = -1;

  bool done = false;
  while (is.good())
    {
      while (is.get() != '\n');
      is >> iter;
      is >> val;
       
      globalIter++;

      if (iter == 0)
        {
          if (lastIter != 0)
            {
              globalIter--;
            }
          os << endl;
        }
      
      double factor = 1;
      os << globalIter << "    "  << val << endl;

      //cout << "iter = " << iter << ", val = " << val << endl;
      
      lastIter = iter;
    }

}
