//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id
//
// Description:
//	namespace for utility functions on Kalman tracks
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//     Dave Brown, LBNL April 1, 2005
//
// Copyright Information:
//	Copyright (C) 2005		Lawrence Berkeley Laboratory
//
//------------------------------------------------------------------------

#ifndef KalUtils_HH
#define KalUtils_HH

//----------------------
// Base Class Headers --
//----------------------
#include "PDT/PdtPid.hh"
#include <map>
class TrkRecoTrk;
class TrkHotList;

//		---------------------
// 		-- Class Interface --
//		---------------------
namespace KalUtils {
// compute the layer overlap fraction for 2 tracks
  double layerOverlap(const TrkRecoTrk* trk1,
		      const TrkRecoTrk* trk2);
// compare 2 track fits at a particular flightlength, optionally
// masking which parameters to compare
  double paramChisq(const TrkRecoTrk* trk1,double flt1,
		    const TrkRecoTrk* trk2,double flt2,
		    PdtPid::PidType hypo=PdtPid::null,
		    bool* tparams=0);
// get full hot list from a track, even mini in cache mode (assuming esd was read)
  const TrkHotList* fullHotList(const TrkRecoTrk*);
// a simple map of which layers are hit on a track; phi and z are counted separate, and Dch
// is offset
  void fillLayerMap(std::map<unsigned, unsigned >& layermap,
		    const TrkHotList* hots,
		    int itrk);

  void countHits(const TrkRecoTrk* trk,double fltlen,
		 int& nbefore,int& nafter);

};
#endif

