//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: MatIsoDictionary.cc 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class MatIsoDictionary (transient version)
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

#include "BTrk/BaBar/BaBar.hh"
//----------------------
// C++ Headers --
//----------------------
#include <fstream>
#include <assert.h>

#include <vector>
#include <algorithm>

//----------------------
// Base Class Headers --
//----------------------
//#include "AbsEnv/AbsEnv.hh"
//#include "GenEnv/GenEnv.hh"
//#include "EidData/EidCondKeyTriplet.hh"
#include "BTrk/MatEnv/MatIsoDictionary.hh"
//#include "BdbTime/BdbTime.hh"
//#include "ProxyDict/AbsArgVal.hh"
//#include "ProxyDict/Ifd.hh"
//#include "ProxyDict/IfdStrKey.hh"
#include "BTrk/BaBar/ExternalInfo.hh"
#include "BTrk/BaBar/FileFinderInterface.hh"
#include "BTrk/BaBar/ErrLog.hh"
using std::fstream;

// Create Constructor

MatIsoDictionary::MatIsoDictionary()
{
// Get the transient MatIsotopeList

  /*MatIsotopeList* isoList = 0;
  BdbTime* toUse;
  if ( gblEnv->getGen() !=0 &&
       gblEnv->getGen()->primaryCondKey() != 0 ) {
    toUse = new BdbTime( gblEnv->getGen()->primaryCondKey()->key() );
    ErrMsg(routine) << "Isotopes being fetched with BdbTime from Env: "
		    << *toUse << endmsg;
  } else {
    toUse = new BdbTime; // current program time
    ErrMsg(error)
      << "BdbTime not in Env. Isotopes are being fetched using BdbTime 'now' ("
      << *toUse << ")" << endmsg;
  }

  AbsArgVal<BdbTime> aarg(*toUse);
  isoList = Ifd< MatIsotopeList >::get( gblPEnv, aarg);

  if ( isoList == 0 ) {
    ErrMsg(fatal)
      << "BTrk/MatEnv/MatIsoDictionary: No access to the list of isotopes"
      << endmsg;
  }

  FillIsoDict(isoList);
  delete toUse;*/
  std::string fullPath = ExternalInfo::fileFinderInstance()->matIsoDictionaryFileName();
  MatIsotopeList* mtrList = new MatIsotopeList(fullPath);
  FillIsoDict(mtrList);

}

void MatIsoDictionary::FillIsoDict(MatIsotopeList* isoList)
{
  std::vector<MatIsotopeObj*>* isoVec = isoList->getIsotopeVector();
  size_t nisotope = isoVec->size();
  for (size_t is=0; is<nisotope; is++)
    {
      // copy the object into the dictionary. The disctionary now has
      // ownership of the copied objects. Use clearAndDestroy in dtor to remove
      // copies from memory.
      MatIsotopeObj* Obj = new MatIsotopeObj(*(*isoVec)[is]);
      std::string* key = new std::string(Obj->getName());
      (*this)[key] = Obj;
    }
}

MatIsoDictionary::~MatIsoDictionary()
{
  std::map<std::string*, MatIsotopeObj*, babar::Collection::PtrLess>::iterator iter;
  for (iter = begin(); iter != end(); ++iter) {
    delete iter->first;
    delete iter->second;
  }
  clear();
}






