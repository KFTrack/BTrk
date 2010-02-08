//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalPairMaker.hh,v 1.5 2000/07/17 18:27:03 roberts Exp $
//
// Description:
//   Creates tracks with KalPairReps.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 1998	Lawrence Berkeley Laboratory
//
// Author(s): David Brown, 3/21/98
//
//------------------------------------------------------------------------

#ifndef KALPAIRMAKER_HH
#define KALPAIRMAKER_HH

#include "KalmanTrack/KalPairConduit.hh"
#include "TrkBase/TrkFitMaker.hh"

// Class interface //
class KalPairMaker : public TrkFitMaker {
public:
  KalPairMaker();
  virtual ~KalPairMaker();
//
//  The input tracks should be 'const', but attach doesn't allow this
//
  TrkRecoTrk* makePairTrack(TrkRecoTrk& oldtrk, KalPairConduit* conduit);

private:	
  // Preempt
  KalPairMaker& operator= (const KalPairMaker&);
  KalPairMaker(const KalPairMaker &);
};
#endif
