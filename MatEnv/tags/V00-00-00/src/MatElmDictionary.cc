//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: MatElmDictionary.cc 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class MatElmDictionary (transient version)
//      Source file
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Mossadek Talby  (SLAC - CPPM/IN2P3 University of Aix-Marseille II)
//
// Modification History:
//   October 15, 1998 - Talby : created
//-----------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

//----------------------
// C++ Headers --
//----------------------
#include <fstream>
#include <assert.h>

#include <string>
#include <vector>
#include <algorithm>
#include "BaBar/BbrCollectionUtils.hh"
using babar::Collection::DeleteObject;

//----------------------
// Base Class Headers --
//----------------------
//#include "AbsEnv/AbsEnv.hh"
//#include "GenEnv/GenEnv.hh"
//#include "EidData/EidCondKeyTriplet.hh"
#include "MatEnv/MatElmDictionary.hh"
//#include "BdbTime/BdbTime.hh"
//#include "ProxyDict/AbsArgVal.hh"
//#include "ProxyDict/Ifd.hh"
//#include "ProxyDict/IfdStrKey.hh"
#include "ErrLogger/ErrLog.hh"
using std::fstream;

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

// Create Constructor 

MatElmDictionary::MatElmDictionary()
{
// Get the transient MatElementList

  /* MatElementList* elmList = 0;
  BdbTime* toUse;
  if ( gblEnv->getGen() !=0 &&
       gblEnv->getGen()->primaryCondKey() != 0 ) {
    toUse = new BdbTime( gblEnv->getGen()->primaryCondKey()->key() );
    ErrMsg(routine) << "Elements being fetched with BdbTime from Env: "
		    << *toUse << endmsg;
  } else {
    toUse = new BdbTime; // current program time
    ErrMsg(error) 
      << "BdbTime not in Env. Elements are being fetched using BdbTime 'now' ("
      << *toUse << ")" << endmsg;
  }

  AbsArgVal<BdbTime> aarg(*toUse);
  elmList = Ifd< MatElementList >::get( gblPEnv, aarg);

  if ( elmList == 0 ) {
    ErrMsg(fatal)
      << "MatEnv/MatElmDictionary: No access to the list of elements"
      << endmsg; 
  }
  
  FillElmDict(elmList);
  delete toUse;*/
  MatElementList* elmList = new MatElementList("BaBar/MatEnv/ElementsList.data");
  FillElmDict(elmList);
}

void MatElmDictionary::FillElmDict(MatElementList* elmList)
{
  std::vector<MatElementObj*>* elmVec = elmList->getElementVector();
  int nelement = elmVec->size();        
  for (size_t ie=0; ie<nelement; ie++) {
    //
    // copy the object into the dictionary. The disctionary now has
    // ownership of the copied objects. Use clearAndDestroy in dtor to remove
    // copies from memory.
    MatElementObj* Obj = new MatElementObj(*(*elmVec)[ie]);
    std::string* key = new std::string(Obj->getName());
    (*this)[key] = Obj;
  }   
}
MatElmDictionary::~MatElmDictionary()
{
  std::map<std::string*, MatElementObj*, babar::Collection::PtrLess>::iterator iter;
  for (iter = begin(); iter != end(); ++iter) {
    delete iter->first;
    delete iter->second;
  }
  clear();
}






