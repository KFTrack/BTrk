//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalMiniInterface.hh,v 1.2 2000/12/05 19:56:18 brownd Exp $
//
// Description:
//     class KalMiniInterface.  Implementation of TrkExtInterface for
//     KalMiniRep
//     
// Copyright Information:
//	Copyright (C) 2000	Lawrence Berkeley Laboratory
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Dave Brown, 11/14/00
//
//------------------------------------------------------------------------

#ifndef KALMINIINTERFACE_HH
#define KALMINIINTERFACE_HH
#include "TrkBase/TrkExtInterface.hh"

class TrkErrCode;
class TrkVolume;
class KalMiniRep;

// Class interface //
class KalMiniInterface : public TrkExtInterface {

public:
  KalMiniInterface();
  virtual ~KalMiniInterface();
  virtual const IfdKey& myKey() const;
// allow acess to the rep
  const KalMiniRep* kalmanMiniRep() const;
  KalMiniRep* kalMiniRep();
private:
  // Preempt 
  KalMiniInterface&   operator= (const KalMiniInterface&);
  KalMiniInterface(const KalMiniInterface &);
};

#endif
