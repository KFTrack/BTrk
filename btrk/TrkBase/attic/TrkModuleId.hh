//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkModuleId.hh,v 1.3 2007/02/05 22:16:38 brownd Exp $
//
// Description:
//     Trivial class to identify modules used in tracking.
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 2004	Lawrence Berkeley Laboratory
//
// Author(s): David Brown, 11/22/04
//
//------------------------------------------------------------------------

#ifndef TRKMODULEID_HH
#define TRKMODULEID_HH

class TrkHistory;
class TrkRecoTrk;
#include <map>
#include <vector>
#include <string>

class TrkModuleId {
public:
// enums defining module mapping.  existing values should NEVER BE CHANGED,
// new numbers may be added as needed

  enum trkFinders {nofinder=0,dchl3trkconverter=1,dchtrackfinder=2,
		   dchradtrkfinder=3,dcxtrackfinder=4,dcxsparsefinder=5,
		   svttrackfinder=6,svtcircle2helix=7,
		   nfinders};
  enum trkModifiers {nomodifier=0,dcxhitadder=1,dchtrkfitupdater=2,trksettrackt0=3,
		     trksettrackt1=4,dchkal1dfit=5,dchkalrx=6,svtkal1dfit=7,
		     svtkalrx=8,dchkalfinalfit=9,
		     svtkalfinalfit=10,trksvthitadder=11,trackmerge=12,
		     trkdchhitadder=13,trkdchradhitadder=14,defaultkalrx=15,
                     trkhitfix=16,trkloopfix=17,trksvthafix=18,trkfailedfix=19,trkmomfix=20,
		     nmodifiers};
// single accessor
  static const TrkModuleId& instance();
// singleton class, so constructor is protected.  Use 'instance' method
  ~TrkModuleId();
// translate enum values to strings
  const std::string& finder(int ifnd) const;
  const std::string& modifier(int imod) const;
// translate string to enum value
  int finder(const std::string& fndmod) const;
  int modifier(const std::string& modmod) const;
// same for TrkHistory.  If the module specified isn't a finder (modifier)
// the appropriate version of 0 will be returned
  int finder(const TrkHistory& hist) const;
  int modifier(const TrkHistory& hist) const;
// access maps
  const std::map<std::string,int>& finders() const { return _finders;}
  const std::map<std::string,int>& modifiers() const { return _modifiers;}
// allow extending the finder and modifier maps
  void addFinder(const TrkHistory& hist);
  void addModifier(const TrkHistory& hist);
// fill and decode a bitmap of finders/modifiers from a track.
  unsigned finderMap(const TrkRecoTrk*) const;
  void setFinders(TrkRecoTrk* trk,unsigned findermap) const;
  unsigned modifierMap(const TrkRecoTrk*) const;
  void setModifiers(TrkRecoTrk* trk,unsigned findermap) const;
private:
// static instance
  static TrkModuleId* _instance;
  TrkModuleId();
// pre-empt
  TrkModuleId(const TrkModuleId&);
  TrkModuleId& operator =(const TrkModuleId&);
// map of finder module names to numbers
  std::map<std::string,int> _finders;
// map of modifier module names to numbers
  std::map<std::string,int> _modifiers;
// reverse-mapping
  std::vector<std::string> _findernames;
  std::vector<std::string> _modifiernames;
};

#endif
