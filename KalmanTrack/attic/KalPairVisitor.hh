//--------------------------------------------------------------------------
// File and Version Information:
//
// Description:  KalPairVisitor is an implementation of the Visitor pattern
//               used to determine the trajectory type of a PairSite
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Doug Roberts
//
//------------------------------------------------------------------------
#ifndef KALPAIRVISITOR_HH
#define KALPAIRVISITOR_HH

#include "TrkBase/TrkVisitor.hh"

class TrkSimpTraj;
class HelixTraj;
class TrkCircleTraj;
class NeutTraj;
class TrkDifLineTraj;
class HepPoint;
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Vector.h"
class BbrVectorErr;
class BbrDoubleErr;
class TrkExchangePar;
class TrkParams;

// Class interface //
class KalPairVisitor : public TrkVisitor {

public:

  KalPairVisitor(const TrkSimpTraj&);

  virtual ~KalPairVisitor();

  // ******************
  // data member access
  // ******************

  const HelixTraj*      helix() const      {return _ht;}
  const TrkCircleTraj*  circle() const     {return _ct;}
  const NeutTraj*       neut() const       {return _nt;}   
  const TrkDifLineTraj* line() const       {return _lt;}

  //********************************
  // The visitor functions:
  //********************************

  virtual void trkVisitHelixTraj(const HelixTraj*);
  virtual void trkVisitCircleTraj(const TrkCircleTraj*);
  virtual void trkVisitNeutTraj(const NeutTraj*);
  virtual void trkVisitLineTraj(const TrkDifLineTraj*);

  // Transform parameters
  bool transformParams(BbrDoubleErr fltlen, BbrVectorErr mom,
		       TrkParams& params) const;

  // Fit for best flightlength to perform the transformation
  BbrDoubleErr fitFlightLength(double fltlen, BbrVectorErr mom,
			       TrkParams otherParams, 
			       int level=0) const;

  enum Coord {pX = 0, pY, pZ, X, Y, Z, nCoord};

  // Calculate derivatives of parameters w.r.t. coordinates {p,x}
  const HepMatrix derivsParWrtCoord(Hep3Vector mom, HepPoint pos,
				    int charge) const;

private:

  const HelixTraj*      _ht;
  const TrkCircleTraj*  _ct;
  const NeutTraj*       _nt;
  const TrkDifLineTraj* _lt;

  // Calculate derivatives of coordinates w.r.t. parameters
  const HepMatrix derivsCoordWrtPar(TrkExchangePar params, int charge,
				    double fltlen) const;
    
  // Calcualte derivatives of coordinates w.r.t. beam constraints
  const HepMatrix derivsCoordWrtBeam() const;

  // Calcualte derivatives of parameters w.r.t. flight length
  const HepVector derivsCoordWrtFlt(TrkExchangePar params, int charge,
				    double fltlen) const;

  double calcSigma2FltLen(const double fltIn, const BbrPointErr x) const;

};

#endif
