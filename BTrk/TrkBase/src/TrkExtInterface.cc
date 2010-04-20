//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkExtInterface.cc,v 1.8 2004/05/03 20:16:38 brownd Exp $
//
// Description:
//     
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner
//
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
#include "TrkBase/TrkExtInterface.hh"
#include "TrkBase/TrkRep.hh"
#include "ErrLogger/ErrLog.hh"

//------------------------------------------------------------------------
TrkExtInterface::TrkExtInterface() : _myRep(0), _nonconst(false){
//------------------------------------------------------------------------
}
 
//------------------------------------------------------------------------
TrkExtInterface::~TrkExtInterface() {
//------------------------------------------------------------------------
}


//------------------------------------------------------------------------
bool
TrkExtInterface::attach(TrkRep* newRep) {
//------------------------------------------------------------------------
  bool retval(false);
  if (myKey() == newRep->myKey()) {
    setRep(newRep);
    retval = true;
  }
  return retval;
}

//------------------------------------------------------------------------
bool
TrkExtInterface::attach(const TrkRep* newRep) {
//------------------------------------------------------------------------
  bool retval(false);
  if (myKey() == newRep->myKey()) {
    setRep(newRep);
    retval = true;
  }
  return retval;
}

void
TrkExtInterface::setRep(const TrkRep* newRep) {
  _myRep = (TrkRep*) newRep;   
  _nonconst = false;
}

void
TrkExtInterface::setRep(TrkRep* newRep) {
  _myRep = newRep;   
  _nonconst = true;
}

//------------------------------------------------------------------------
TrkRep* 
TrkExtInterface::myRep() {
//------------------------------------------------------------------------
  if(_myRep != 0 && _nonconst)
    return _myRep;
  else if (0 == _myRep) {
    ErrMsg(error) << "Cannot use interface without attached rep." 
		  << endmsg;
    return 0;
  } else
    ErrMsg(error) << "Cannot return non-const rep after const attachment" 
		  << endmsg;
    return 0;
}

//------------------------------------------------------------------------
const TrkRep* 
TrkExtInterface::myConstRep() const {
//------------------------------------------------------------------------
  if (0 == _myRep) {
    ErrMsg(error) << "Cannot use interface without attached rep." 
		  << endmsg;
    return 0;
  }
  return _myRep;
}

