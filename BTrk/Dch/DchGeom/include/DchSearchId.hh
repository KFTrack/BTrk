#ifndef DCHSEARCHID_HH
#define DCHSEARCHID_HH

#include <string>

struct DchSearchId{
  int*    idLim;
  std::string*   tag; // string expression to match set/element name
  DchSearchId(int* idlim, std::string* stag) :
    idLim(idlim), tag(stag)
  {;}
};

#endif
