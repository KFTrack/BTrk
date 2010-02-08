//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalInterface.hh,v 1.11 2002/04/23 17:28:08 brownd Exp $
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

#ifndef KALINTERFACE_HH
#define KALINTERFACE_HH
#include "TrkBase/TrkExtInterface.hh"
#include "TrkBase/TrkDirection.hh"
#include "TrkBase/TrkSimpTraj.hh"
#include "KalmanTrack/KalMiniInterface.hh"
#include "KalmanTrack/KalPairInterface.hh"

class TrkErrCode;
class TrkVolume;
class KalRep;

// Class interface //
class KalInterface : public TrkExtInterface {

public:
  KalInterface();
  virtual ~KalInterface();
  virtual const IfdKey& myKey() const;
// I have to override attach, to allow attaching to either a KalRep or a KalMiniRep
  virtual bool attach(TrkRep*);
  virtual bool attach(const TrkRep*);
  TrkErrCode extendThrough(const TrkVolume&, trkDirection trkDir = trkOut);
  const TrkSimpTraj* seed() const;
  double refMomentum()const;
  double refMomFltLen()const;
  const KalRep* kalmanRep() const;
  KalRep* kalRep();
private:
  KalMiniInterface _miniiface;
  KalPairInterface _pairiface;
  bool _mini;
  // Preempt 
  KalInterface&   operator= (const KalInterface&);
  KalInterface(const KalInterface &);
  friend class KalMaker;
};

#endif







