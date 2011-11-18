//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: MatMtrDictionary.cc 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class MatMtrDictionary (transient version)
//      Source file
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Mossadek Talby  (SLAC - CPPM/IN2P3 University of Aix-Marseille II)
//      MatMaterialObj* theMaterial;

// Modification History:
//   October 15, 1998 - Talby : created
//-----------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
//----------------------
// C++ Headers --
//----------------------
#include <fstream>
#include <assert.h>

#include <vector>
#include <algorithm>
#include "BaBar/BbrCollectionUtils.hh"
#include <map>

//----------------------
// Base Class Headers --
//----------------------
//#include "AbsEnv/AbsEnv.hh"
//#include "GenEnv/GenEnv.hh"
//#include "EidData/EidCondKeyTriplet.hh"
#include "MatEnv/MatMtrDictionary.hh"
//#include "BdbTime/BdbTime.hh"
//#include "ProxyDict/AbsArgVal.hh"
//#include "ProxyDict/Ifd.hh"
//#include "ProxyDict/IfdStrKey.hh"
#include "ErrLogger/ErrLog.hh"

using std::fstream;

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "Mu2eUtilities/inc/ConfigFileLookupPolicy.hh"

// Create Constructor 

MatMtrDictionary::MatMtrDictionary()
{
// Get the transient MatMaterialList

  /* MatMaterialList* mtrList = 0;
  BdbTime* toUse;
  if ( gblEnv->getGen() !=0 &&
       gblEnv->getGen()->primaryCondKey() != 0 ) {
    toUse = new BdbTime( gblEnv->getGen()->primaryCondKey()->key() );
    ErrMsg(routine) << "MatMtrDictionary: Materials being fetched with BdbTime from Env: "
		    << *toUse << endmsg;
  } else {
    toUse = new BdbTime; // current program time
    ErrMsg(error) 
      << "MatMtrDictionary: BdbTime not in Env. Materials are being fetched using BdbTime 'now' ("
      << *toUse << ")" << endmsg;
  }

  AbsArgVal<BdbTime> aarg(*toUse);
  mtrList = Ifd< MatMaterialList >::get( gblPEnv, aarg);
  
  if ( mtrList == 0 ) {
    ErrMsg(fatal)
      << "MatMtrDictionary: No access to the list of materials"
      << endmsg; 
  }
  */
  mu2e::ConfigFileLookupPolicy findFile;
  std::string fullPath = findFile("BaBar/MatEnv/MaterialsList.data");
  MatMaterialList* mtrList = new MatMaterialList(fullPath);
  FillMtrDict(mtrList);
}

void MatMtrDictionary::FillMtrDict(MatMaterialList* mtrList)
{
  std::vector<MatMaterialObj*>* mtrVec = mtrList->getMaterialVector();
  size_t nmaterial = mtrVec->size();        
  for (size_t im=0; im<nmaterial; im++){
    //
    // copy the object into the dictionary. The disctionary now has
    // ownership of the copied objects. Use clearAndDestroy in dtor to remove
    // copies from memory.
    MatMaterialObj* Obj = new MatMaterialObj(*(*mtrVec)[im]);
    std::string* key = new std::string(Obj->getName());
    (*this)[key] = Obj;
    ErrMsg(routine) << "MatMtrDictionary: Inserted Material " << *key << endmsg;
    }   
}

MatMtrDictionary::~MatMtrDictionary()
{
  std::map<std::string*, MatMaterialObj*, babar::Collection::PtrLess>::iterator 
    iter = begin();
  for (; iter != end(); ++iter) {
    delete iter->first;
    delete iter->second;
  }
  clear();
}





