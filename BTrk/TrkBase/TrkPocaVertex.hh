//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkPocaVertex.hh,v 1.7 2006/03/25 15:15:56 brownd Exp $
//
// Description:
// Utility to vertex two tracks, or a track and a point, using TrkPoca.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Wouter Hulsbergen
//
//------------------------------------------------------------------------

#ifndef __TRKPOCAVERTEX_HH__
#define __TRKPOCAVERTEX_HH__

#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/TrkBase/TrkVertex.hh"
#include "BTrk/TrkBase/TrkParticle.hh"

class TrkRep ;

class TrkPocaVertex : public TrkVertex {  
public:
  TrkPocaVertex();
  TrkPocaVertex(const TrkRep* trk1, double flt1, TrkParticle const& tpart1,
		const TrkRep* trk2, double flt2, TrkParticle const& tpart2,
		double precision=1.e-5) ;
  
  TrkPocaVertex(const TrkRep* trk, double flt, TrkParticle const& tpart,
 		const BbrPointErr& pt,
		double precision=1.e-5);


  TrkPocaVertex(const TrkPocaVertex& other);
  TrkPocaVertex& operator =(const TrkPocaVertex& other);
  virtual TrkPocaVertex* clone() const;

  const TrkPoca& poca() const { return _poca; }
  double doca() const { return _poca.doca(); }
  double flt1() const { return _poca.flt1() ; }
  double flt2() const { return _poca.flt2() ; }

  static CLHEP::HepVector crossproduct(const CLHEP::HepVector& v1,const CLHEP::HepVector& v2);
  template<class T> CLHEP::HepVector convert(const T& vec);

private:
  TrkPoca _poca ;
} ;

#endif

