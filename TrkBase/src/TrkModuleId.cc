//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkModuleId.cc,v 1.3 2007/02/05 22:16:38 brownd Exp $
//
//
// Copyright Information:
//	Copyright (C) 2004	Lawrence Berkeley Laboratory
//
// Author(s): David Brown, 11/22/04
//
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
#include "TrkBase/TrkModuleId.hh"
#include "TrkBase/TrkHistory.hh"
#include "TrkBase/TrkFitStatus.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include <assert.h>

TrkModuleId* TrkModuleId::_instance(0);

const TrkModuleId&
TrkModuleId::instance() {
  if(_instance == 0){
    _instance = new TrkModuleId();
    assert(_instance != 0);
  }
  return *_instance;
}

TrkModuleId::TrkModuleId() {
// fill the list of modules
  _finders["Unknown"] = nofinder;
  _finders["DchL3TrkConverter"] = dchl3trkconverter;
  _finders["DchTrackFinder"] = dchtrackfinder;
  _finders["DchRadTrkFinder"] = dchradtrkfinder;
  _finders["DcxTrackFinder"] =  dcxtrackfinder;
  _finders["DcxSparseFinder"] =  dcxsparsefinder;
  _finders["SvtTrackFinder"] =  svttrackfinder;
  _finders["SvtCircle2Helix"] =  svtcircle2helix;
//
  _modifiers["Unknown"] = nomodifier;
  _modifiers["DcxHitAdder"] =  dcxhitadder;
  _modifiers["DchTrkFitUpdater"] =  dchtrkfitupdater;
  _modifiers["TrkSetTrackT0"] =  trksettrackt0;
  _modifiers["TrkSetTrackT1"] =  trksettrackt1;
  _modifiers["DchKal1DFit"] =  dchkal1dfit;
  _modifiers["DchKalRX"] =  dchkalrx;
  _modifiers["SvtKal1DFit"] =  svtkal1dfit;
  _modifiers["SvtKalRX"] =  svtkalrx;
  _modifiers["DchKalFinalFit"] =  dchkalfinalfit;
  _modifiers["SvtKalFinalFit"] =  svtkalfinalfit;
  _modifiers["TrkSvtHitAdder"] =  trksvthitadder;
  _modifiers["TrackMerge"] =  trackmerge;
  _modifiers["TrkDchHitAdder"] = trkdchhitadder;
  _modifiers["TrkDchRadHitAdder"] =  trkdchradhitadder;
  _modifiers["DefaultKalRX"] =  defaultkalrx;
  _modifiers["TrkHitFix"] =  trkhitfix;
  _modifiers["TrkLoopFix"] =  trkloopfix;
  _modifiers["TrkSvtHitAdderFixup"] =  trksvthafix;
  _modifiers["TrkFailedRecovery"] =  trkfailedfix;
  _modifiers["TrkdEdxMomConstrain"] = trkmomfix;
// translate these to vectors
  _findernames = std::vector<std::string>(_finders.size());
  std::map<std::string,int>::iterator ifnd = _finders.begin();
  while(ifnd != _finders.end()){
    _findernames[ifnd->second] = ifnd->first;
    ++ifnd;
  }
  _modifiernames = std::vector<std::string>(_modifiers.size());
  std::map<std::string,int>::iterator imod = _modifiers.begin();
  while(imod != _modifiers.end()){
    _modifiernames[imod->second] = imod->first;
    ++imod;
  }
// test consistency
  for(unsigned i=0;i<_findernames.size();i++)
    assert(i == _finders[_findernames[i]]);
  for(unsigned j=0;j<_modifiernames.size();j++)
    assert(j == _modifiers[_modifiernames[j]]);
}

TrkModuleId::~TrkModuleId()
{}



int
TrkModuleId::finder(const std::string& fndmod) const {
  std::map<std::string,int>::const_iterator ifnd =
    _finders.find(fndmod);
  if(ifnd != _finders.end())
    return (TrkModuleId::trkFinders)ifnd->second;
  else
    return nofinder;
}


int
TrkModuleId::modifier(const std::string& modmod) const {
  std::map<std::string,int>::const_iterator imod =
    _modifiers.find(modmod);
  if(imod != _modifiers.end())
    return (TrkModuleId::trkModifiers)imod->second;
  else
    return nomodifier;
}

const std::string&
TrkModuleId::finder(int ifnd) const {
  if(ifnd < _findernames.size())
    return _findernames[ifnd];
  else {
    return _findernames[nofinder];
  }
}

const std::string&
TrkModuleId::modifier(int imod) const {
  if(imod < _modifiernames.size())
    return _modifiernames[imod];
  else {
    return _modifiernames[nomodifier];
  }
}

int
TrkModuleId::finder(const TrkHistory& hist) const {
  return finder(hist.module());
}

int
TrkModuleId::modifier(const TrkHistory& hist) const {
  return modifier(hist.module());
}


void
TrkModuleId::addFinder(const TrkHistory& hist) {
  const std::string& fndmod = hist.module();
// make sure this finder isn't already there
  std::map<std::string,int>::const_iterator ifnd =
    _finders.find(fndmod);
  if(ifnd == _finders.end()){
    _findernames.push_back(fndmod);
    _finders[fndmod] = _findernames.size();
  }
}

void
TrkModuleId::addModifier(const TrkHistory& hist) {
  const std::string& modmod = hist.module();
// make sure this finder isn't already there
  std::map<std::string,int>::const_iterator imod =
    _modifiers.find(modmod);
  if(imod == _modifiers.end()){
    _modifiernames.push_back(modmod);
    _modifiers[modmod] = _modifiernames.size();
  }
}

unsigned
TrkModuleId::finderMap(const TrkRecoTrk* trk) const {
  unsigned finders(0);
  const TrkFitStatus* status = trk->status(); // default rep status
  if(status != 0){
// loop over history
    std::vector<TrkHistory>::const_iterator ihist = status->beginHistory();
    while( ihist != status->endHistory()){
      int ifnd = finder(*ihist);
      if(ifnd != TrkModuleId::nofinder)
        finders |= (1 << ifnd);
      ++ihist;
    }
  }
  return finders;
}

unsigned
TrkModuleId::modifierMap(const TrkRecoTrk* trk) const {
  unsigned modifiers(0);
  const TrkFitStatus* status = trk->status(); // default rep status
  if(status != 0){
// loop over history
    std::vector<TrkHistory>::const_iterator ihist = status->beginHistory();
    while( ihist != status->endHistory()){
      int imod = modifier(*ihist);
      if(imod != TrkModuleId::nomodifier)
        modifiers |= (1 << imod);
      ++ihist;
    }
  }
  return modifiers;
}

void
TrkModuleId::setFinders(TrkRecoTrk* trk,unsigned findermap) const {
  TrkFitStatus* status = trk->status();
  if(status != 0){
// loop over bits in finder and modifier words
    unsigned ibit;
    for(ibit=0;ibit<32;++ibit){
      if(findermap & (1 << ibit)){
// extend the vector if this bit hasn't been seen before
	status->addHistory(TrkErrCode(),finder(ibit).c_str());
      }
    }
  }
}

void
TrkModuleId::setModifiers(TrkRecoTrk* trk,unsigned modifiermap) const {
  TrkFitStatus* status = trk->status();
  if(status != 0){
// loop over bits in finder and modifier words
    unsigned ibit;
    for(ibit=0;ibit<32;++ibit){
      if(modifiermap & (1 << ibit)){
// extend the vector if this bit hasn't been seen before
	status->addHistory(TrkErrCode(),modifier(ibit).c_str());
      }
    }
  }
}
