//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id
//
// Description:
//	namespace for functions that compare 2 tracks.
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

#include "BaBar/BaBar.hh"
#include "KalmanTrack/KalUtils.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "TrkBase/TrkHotList.hh"
#include "TrkBase/TrkHitOnTrk.hh"
#include "KalmanTrack/KalInterface.hh"
#include "KalmanTrack/KalMiniInterface.hh"
#include "KalmanTrack/KalMiniRep.hh"
#include "KalmanTrack/KalRep.hh"

double 
KalUtils::layerOverlap(const TrkRecoTrk* trk1,
			 const TrkRecoTrk* trk2) {
  double retval(0.0);
// get the default hot lists for each track
  const TrkHotList* hots1 = fullHotList(trk1);
  const TrkHotList* hots2 = fullHotList(trk2);
  if(hots1 != 0 && hots2 != 0 ){
// build a map of hit layers on each track
    std::map<unsigned, unsigned > layermap;
// initialize to 0
    for(unsigned ilay=0;ilay<65;++ilay)
      layermap[ilay]=0;
// fill with hits
    fillLayerMap(layermap,hots1,0);
    fillLayerMap(layermap,hots2,1);
    double nlayer(0.0);
    double nover(0.0);
    std::map<unsigned,unsigned>::iterator ilay = layermap.begin();
    while(ilay != layermap.end()){
      unsigned layer = (ilay->second)&0x3;
// bit-packed
      if( layer == 0x3)
	nover += 1.0;
      if( layer != 0)
	nlayer += 1.0;
      ++ilay;
    }
    retval = nover/nlayer;
  }
  return retval;
}

void
KalUtils::fillLayerMap(std::map<unsigned,unsigned >& layermap,
			 const TrkHotList* hots,
			 int itrk) {
  TrkHotList::hot_iterator ihot = hots->begin();
  while(ihot != hots->end()){
    if(ihot->isActive()){
      int layer = ihot->layerNumber();
      if(ihot->dchHitOnTrack() != 0)
	layer += 22;
      else if(ihot->svtHitOnTrack() != 0){
	layer  *= 2;
	if(ihot->whatView() == TrkEnums::zView)
	  layer += 1;
      }
      layermap[layer] |= (1<<itrk);
    }
    ++ihot;
  }
}

const TrkHotList*
KalUtils::fullHotList(const TrkRecoTrk* trk) {
  static KalMiniInterface kalmi;
  const TrkHotList* retval(0);
  if(trk->hots() != 0 &&
     trk->hots()->hitCapable())
    retval = trk->hots();
  else if(trk->attach(kalmi,trk->defaultType()) &&
	  kalmi.kalmanMiniRep() != 0)
// mini track: attach the mini interface and get the full hot list
    retval = kalmi.kalmanMiniRep()->fullHotList();
  return retval;
}

double
KalUtils::paramChisq(const TrkRecoTrk* trk1,double flt1,
		     const TrkRecoTrk* trk2,double flt2,
		     PdtPid::PidType hypo,
		     bool* tparams) {
  double retval(-1.0);
  const TrkFit* fit1 = trk1->fitResult();
  const TrkFit* fit2 = trk2->fitResult();
  if(fit1 !=0 && fit2 != 0){
    double lflt1,lflt2;
    const TrkSimpTraj* ltraj1 = fit1->traj().localTrajectory(flt1,lflt1);
    const TrkSimpTraj* ltraj2 = fit2->traj().localTrajectory(flt2,lflt2);    
// find fit result at the respective flightlengths
// convert these parameters to KalParams object
    KalParams kp1(*ltraj1->parameters());
    KalParams kp2(*ltraj2->parameters());
// take their chisq-space difference, including masking if requested
    retval = kp2.chisq(kp1,tparams);
  }
  return retval;
}

void
KalUtils::countHits(const TrkRecoTrk* trk,
		    double fltlen,
		    int& nbefore,
		    int& nafter) {
  nbefore = 0;
  nafter = 0;
// find hots and loop
  const TrkHotList* thots = KalUtils::fullHotList(trk);
  if(thots != 0){
    for(TrkHotList::hot_iterator ihot(thots->begin());ihot!=thots->end();ihot++){
      const TrkHitOnTrk* hot = ihot.get();
      if(hot->isActive()){
        if(hot->fltLen() < fltlen)
          nbefore++;
        else
          nafter++;
      }
    }
  } else
    ErrMsg(error) << "Full hots list not found!  cannot count hots" << endmsg;
}
