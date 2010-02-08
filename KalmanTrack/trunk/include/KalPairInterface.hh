//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalPairInterface.hh,v 1.3 2002/04/24 00:41:18 brownd Exp $
//
// Description:
//     
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Doug Roberts
//
//------------------------------------------------------------------------

#ifndef KALPAIRINTERFACE_HH
#define KALPAIRINTERFACE_HH
#include "TrkBase/TrkExtInterface.hh"
#include "TrkBase/TrkDirection.hh"
#include "TrkBase/TrkSimpTraj.hh"

class TrkErrCode;
class TrkVolume;
class KalPairRep;

// Class interface //
class KalPairInterface : public TrkExtInterface {

public:
  KalPairInterface() {}
  virtual ~KalPairInterface();
  virtual const IfdKey& myKey() const;

  //  TrkErrCode extendThrough(const TrkVolume&, trkDirection trkDir = trkOut);
  const TrkSimpTraj* seed() const;
  double refMomentum()const;
  double refMomFltLen()const;
  const KalPairRep* kalmanPairRep() const;
  KalPairRep* kalPairRep();

private:	
  // Preempt 
  KalPairInterface&   operator= (const KalPairInterface&);
  KalPairInterface(const KalPairInterface &);
};

#endif







