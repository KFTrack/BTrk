//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalMiniPromoter.hh,v 1.4 2007/07/09 21:54:59 brownd Exp $
//
// Description:
//   Promote KalMiniRep-based track fits to refit mode
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 2004	Lawrence Berkeley Laboratory
//
// Author(s): David Brown, Brian Aagaard
//
//------------------------------------------------------------------------

#ifndef KALMINIPROMOTER_HH
#define KALMINIPROMOTER_HH
#include "PDT/PdtPid.hh"
#include "TrkBase/TrkErrCode.hh"
class TrkHotSelector;
class TrkRecoTrk;
class KalMiniRep;
class KalRep;

namespace KalMiniPromoter {
  
// promote the requested fit to 'hots' (refit) mode.  Optionally filter the hots before the fit
// the return code is either the status of the last successful operation, or the first
// failing operation.  If the promotion fails, the fit is returned to its original state.
// default fit hypo is used if none is provided.
  TrkErrCode promoteToHots(TrkRecoTrk* track,
			   PdtPid::PidType hypo=PdtPid::null,
			   const TrkHotSelector* selector=0,
			   bool restorefail=true);
// optionally refine hotlist before refitting
  void selectHots(KalMiniRep*,const TrkHotSelector* selector);
// find the KalRep underneath the track rep for a given hypo.  Works on mini or reconstruction events
// this function DOES NOT fit the rep!!
//
  KalRep* findKalRep(TrkRecoTrk* trk,
                     PdtPid::PidType hypo=PdtPid::null);
};
#endif
