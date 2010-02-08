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

#include "TrkBase/TrkPoca.hh"
#include "TrkBase/TrkVertex.hh"
#include "PDT/PdtPid.hh"

class TrkRecoTrk ;

class TrkPocaVertex : public TrkVertex {  
public:
  TrkPocaVertex();
  TrkPocaVertex(const TrkRecoTrk* trk1, double flt1, PdtPid::PidType pid1,
		const TrkRecoTrk* trk2, double flt2, PdtPid::PidType pid2,
		const PdtEntry* vtype=0,
		double precision=1.e-5) ;
  
  TrkPocaVertex(const TrkRecoTrk* trk, double flt, PdtPid::PidType pid,
 		const BbrPointErr& pt,
		const PdtEntry* vtype=0,double precision=1.e-5);


  TrkPocaVertex(const TrkPocaVertex& other);
  TrkPocaVertex& operator =(const TrkPocaVertex& other);
  virtual TrkPocaVertex* clone() const;

  const TrkPoca& poca() const { return _poca; }
  double doca() const { return _poca.doca(); }
  double flt1() const { return _poca.flt1() ; }
  double flt2() const { return _poca.flt2() ; }

  static HepVector crossproduct(const HepVector& v1,const HepVector& v2);
  template<class T> HepVector convert(const T& vec);

private:
  TrkPoca _poca ;
} ;

#endif

