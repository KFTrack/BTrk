//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkPrimaryVertex.hh,v 1.5 2005/10/06 22:37:40 brownd Exp $
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

#ifndef __TRKPRIMARYVERTEX_HH__
#define __TRKPRIMARYVERTEX_HH__

#include "TrkBase/TrkVertex.hh"

class TrkRecoTrk ;
class TrkPocaVertex;

class TrkPrimaryVertex : public TrkVertex
{
public:
  TrkPrimaryVertex();
  TrkPrimaryVertex(const TrkPrimaryVertex& other);
  TrkPrimaryVertex(const TrkPocaVertex& seedvertex,
		   const trkcontainer& seedtrklist,
		   double maxDeltaChisqPerDof=9,
		   bool sort=true) ;
  TrkPrimaryVertex(const BbrPointErr& seedvertex,
		   const trkcontainer& seedtrklist,
		   double maxDeltaChisqPerDof=9,
		   bool sort=true) ;
  TrkPrimaryVertex(const trkcontainer& seedtrklist,
		   double maxSeedChisqPerDof=9,
		   double maxDeltaChisqPerDof=9);
  virtual ~TrkPrimaryVertex() {}
  TrkPrimaryVertex& operator = (const TrkPrimaryVertex& other);
  virtual TrkPrimaryVertex* clone() const;
  
private:
  bool addTracks(const trkcontainer& tracks, double maxDeltaChisqPerDof,
		 double& chisq, int& nDof) ;
  trkcontainer sortTrackList(const trkcontainer& trklist) ;

} ;

#endif
