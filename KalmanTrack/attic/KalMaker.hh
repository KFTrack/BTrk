//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalMaker.hh,v 1.26 2003/03/28 17:24:16 brownd Exp $
//
// Description:
//   Creates tracks with KalReps.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 1997	Lawrence Berkeley Laboratory
//
// Author(s): David Brown, 4/16/97
//
//------------------------------------------------------------------------

#ifndef KALMAKER_HH
#define KALMAKER_HH
#include "TrkBase/TrkFitMaker.hh"
#include "PDT/PdtPid.hh" 
#include "KalmanTrack/KalRep.hh"

class TrkRep;
class TrkExchangePar;
class TrkRecoTrk;
class TrkHitOnTrk;
class TrkSimpTraj;
class TrkErrCode;
class TrkContext;
class KalStub;
class KalContext;

// Class interface //
class KalMaker : public TrkFitMaker {
public:
  KalMaker(const KalContext& context);
  virtual ~KalMaker();
// changefit functions will create a KalRep for the mass hypo defined
// by the KalContext, deleting all previous reps.  This function
// assumes an existing fit (else it will return an error)
  TrkErrCode changeFit(TrkRecoTrk& theTrack) const;
// Similair to the above, but this function creates a KalRep
// _constrained_ to a subset of its seed parameters.  The constraint
// is applied with correct statistical weights for those parameters marked
// as 'true' in the bool vector, at the specified flight length.
  TrkErrCode changeFit(TrkRecoTrk& theTrack,
		       const TrkSimpTraj&,
		       double constraintfltlen,
		       bool* constrainparams) const;
// create a new track with a KalRep.  By default, this will create the
// hypo specified by the KalContext object, which becomes the tracks default
// type.  The KalContext hypo can be overridden by setting a non-null argument
// in the call (which will become the tracks default hypo).
  TrkRecoTrk* makeTrack(const TrkExchangePar&,
			const TrkContext&, double trackT0,
			PdtPid::PidType hypo=PdtPid::null) const;
// same, seeding the fit with a SimpTraj.  Charge and momentum will be obtained
// from TrkMomCalculator.  Optionally constrain the fit according to a subset
// of the input parameters
  TrkRecoTrk* makeTrack(const TrkSimpTraj&, 
			TrkHotList* hots,
			const TrkContext&, 
			double trackT0,
			PdtPid::PidType hypo=PdtPid::null,
			double constraintfltlen=0,
			bool* constrainparams=0,
			bool stealhots=false) const;

// create a new track from a KalStub
  TrkRecoTrk* makeTrack(KalStub& stub) const;
// This is used to reconstitute the default rep from persistent tracks
  KalRep* makeRep(TrkRecoTrk& theTrack,
		  const TrkExchangePar& helix,
		  const TrkHotList* hotlist) const;
// Add a new mass hypothesis to an existing track.  The default hypothesis
// MUST ALREADY BE A KALREP for this function to succeed (the default hypo
// need not have been fit). Optionally make this new hypo the default
// for this track.
  TrkErrCode addHypo(TrkRecoTrk& theTrack,
		     PdtPid::PidType hypo=PdtPid::null,
		     bool makeDefault=false) const;
// fit a hypo (by default, the tracks default hypo.  this does all the ugly
// machinations.  It will only work if the rep for this hypo exists and is
// a KalRep
  TrkErrCode fitHypo(TrkRecoTrk& theTrack,
		     PdtPid::PidType hypo=PdtPid::null) const;
private:	
// Preempt
  KalMaker&   operator= (const KalMaker&);
  KalMaker(const KalMaker &);
  const KalContext& _kalcon;
};
#endif
