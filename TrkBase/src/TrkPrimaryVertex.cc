//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkPrimaryVertex.cc,v 1.8 2007/05/18 03:39:52 kelsey Exp $
//
// Description:
//
// Utility to find a 'primary' vertex in a set of tracks. A 'seed'
// point can be provided. Tracks are ordered in pt and added
// consecutively, provided that the chisquare increament is smaller
// than maxDeltaChisqPerDof.
//
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Wouter Hulsbergen
//
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
#include "TrkBase/TrkPrimaryVertex.hh"
#include "TrkBase/TrkPocaVertex.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "PDT/Pdt.hh"
#include <algorithm>


TrkPrimaryVertex::TrkPrimaryVertex(const BbrPointErr& seedvertex,
				   const trkcontainer& seedtrklist,
				   double maxDeltaChisqPerDof,
				   bool sort)
{
  _position = seedvertex;
  _status = TrkErrCode(TrkErrCode::succeed);
  trkcontainer tracks;
  if(sort)
    tracks = sortTrackList(seedtrklist);
  else
    tracks = seedtrklist;
  double chisq(0.0);
  int nDof(0);
  if(addTracks(tracks,maxDeltaChisqPerDof,chisq,nDof))
    setConsistency(chisq,nDof);
}

TrkPrimaryVertex::TrkPrimaryVertex(const TrkPocaVertex& seedvertex,
				   const trkcontainer& seedtrklist,
				   double maxDeltaChisqPerDof,
				   bool sort) :
  TrkVertex(seedvertex)
{
  trkcontainer tracks;
  if(sort)
    tracks = sortTrackList(seedtrklist);
  else
    tracks = seedtrklist;
  double chisq(0.0);
  int nDof(0);
  if(addTracks(tracks,maxDeltaChisqPerDof,chisq,nDof))
    setConsistency(chisq,nDof);
}


TrkPrimaryVertex::TrkPrimaryVertex(const trkcontainer& seedtrklist,
				   double maxSeedChisqPerDof,
				   double maxDeltaChisqPerDof)
{
  bool foundseed(false) ;
  trkcontainer tracks = sortTrackList(seedtrklist) ;
  double chisq(0.0);
  int nDof(0);
  for( trkcontainer::iterator it1 = tracks.begin() ;
       it1 != tracks.end() && !foundseed; ++it1) {
    trkcontainer::iterator it2 = it1 ;
    for( ++it2; it2 != tracks.end() && !foundseed; ++it2) {
      const TrkRecoTrk* trk1 = *it1;
      const TrkRecoTrk* trk2 = *it2;

      TrkPocaVertex pocavertex( trk1, 0, trk1->defaultType(), 
				trk2, 0, trk2->defaultType()) ;
      if(pocavertex.status().success() && 
	 pocavertex.chisq()<maxSeedChisqPerDof) {
	_status   = pocavertex.status() ;
	_position = pocavertex.position() ;
	_p4       += pocavertex.p4() ;
	chisq    += pocavertex.chisq() ;
	nDof     += pocavertex.nDof() ;
	_usedtracks.push_back( *it1 ) ;
	_usedtracks.push_back( *it2) ;
	foundseed = true ;
	*it1 = 0 ;
	*it2 = 0 ;
      } 
    }
  }

  if(!foundseed) {
    // problem
    _status.setFailure(1,"no seed found") ;
  } else {
    // remove the tracks we just added and try to add remaining tracks
    trkcontainer::iterator newend =
      std::remove( tracks.begin(), tracks.end(), trkcontainer::value_type(0)) ;
    tracks.erase(newend,tracks.end()) ;
    if(addTracks( tracks,maxDeltaChisqPerDof,chisq,nDof))
      setConsistency(chisq,nDof);
  }
}

TrkPrimaryVertex::TrkPrimaryVertex()
{}

TrkPrimaryVertex::TrkPrimaryVertex(const TrkPrimaryVertex& other) :
  TrkVertex(other)
{}

TrkPrimaryVertex*
TrkPrimaryVertex::clone() const {
  return new TrkPrimaryVertex(*this);
}


TrkPrimaryVertex&
TrkPrimaryVertex::operator = (const TrkPrimaryVertex& other) {
  if(&other != this)
    TrkVertex::operator =(other);
  return *this;
}

bool
TrkPrimaryVertex::addTracks(const trkcontainer& tracks, double maxdchisq,
			    double& chisq, int& nDof)
{
  bool retval(false);
  for( trkcontainer::const_iterator it = tracks.begin() ;
       it != tracks.end() ; it++) {
    const TrkRecoTrk* trk = *it;
    TrkPocaVertex pocavertex( trk, 0, trk->defaultType(), _position) ;
    if( pocavertex.status().success() && 
	pocavertex.chisq() < pocavertex.nDof()*maxdchisq) {
      _usedtracks.push_back(trk) ;
      _info[trk] = TrkVtxInfo(pocavertex.doca(),pocavertex.flt1(),trk->defaultType());
      nDof  += pocavertex.nDof() ;
      chisq += pocavertex.chisq() ;
      _p4    += pocavertex.p4() ;
      _position = pocavertex.position() ;
      retval = true;
    }
  }
  return retval;
}

inline bool trkptsort(const TrkRecoTrk* lhs, const TrkRecoTrk* rhs) 
{
  return fabs(lhs->fitResult()->pt(0))>fabs(rhs->fitResult()->pt(0)) ;
}

TrkPrimaryVertex::trkcontainer
TrkPrimaryVertex::sortTrackList(const trkcontainer& trklist)
{
  trkcontainer tracks = trklist ;
  std::sort( tracks.begin(), tracks.end(), trkptsort ) ;
  return tracks ;
}
