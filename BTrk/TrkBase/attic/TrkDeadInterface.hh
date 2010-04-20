//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkDeadInterface.hh,v 1.1 2004/04/16 22:47:44 brownd Exp $
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

#ifndef TRKDEADINTERFACE_HH
#define TRKDEADINTERFACE_HH
#include "TrkBase/TrkExtInterface.hh"
#include "ProxyDict/IfdIntKey.hh"
// Class interface //

class TrkDeadInterface : public TrkExtInterface {
public:
  TrkDeadInterface() : _mykey(0) {;}
  virtual ~TrkDeadInterface() {;}
  virtual const IfdKey& myKey() const { return _mykey; }
private:
  IfdIntKey _mykey;
};

#endif







