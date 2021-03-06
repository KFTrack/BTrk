//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkHelixUtils.hh,v 1.3 2008/04/28 06:02:01 brownd Exp $
//
// Description:  package of utility routines for doing things to helices.
//  No data members.  I'll probably want to put this someplace else 
//  eventually.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner
//
//------------------------------------------------------------------------

#ifndef TRKHELIXUTILS_HH
#define TRKHELIXUTILS_HH

#include "CLHEP/Matrix/Matrix.h"
class HelixParams;
class HepPoint;
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/SymMatrix.h"
class BField;
class BbrPointErr;
class BbrVectorErr;
class NeutParams;
class Trajectory;
// Class interface //
class TrkHelixUtils {

public:
  TrkHelixUtils() { };
// basic function to conver point and momentum to 5-parameters, plus
// the flightlength at the vertex point.
  static void helixFromMom(CLHEP::HepVector& parvec, double& fltlen, const HepPoint& vertex, 
		const CLHEP::Hep3Vector& p, double sign, const BField&);
// same, inputing value of Bz
	static void helixFromMom(CLHEP::HepVector& parvec, double& fltlen, const HepPoint& vertex, 
  	const CLHEP::Hep3Vector& p, double sign, double Bz);

     // Create a helix-set from a position and a momentum.
     // Uses nominal B field to determine curvature. 
  static HelixParams helixFromMom(const HepPoint& vertex, 
		const CLHEP::Hep3Vector& p, double sign, const BField&);

     // Does the same, but gives *real* errors on the parameters
     // (instead of a default error matrix)
     // Uses nominal B field to determine curvature. 
  static HelixParams helixFromMomErr(const BbrPointErr& vertex,
                const BbrVectorErr& p,const CLHEP::HepMatrix& cxp,  double sign, const BField&);
     // Does the same but for neutrals
  static NeutParams lineFromMomErr(const BbrPointErr& vertex,
                const BbrVectorErr& p,const CLHEP::HepMatrix& cxp,  double sign, const BField&); 

     // Jacobian for transforming std helix params to new set defined at fltNew
  static CLHEP::HepMatrix jacobianExtrapolate(const HelixParams&, double fltNew);

     // Actually transform the error matrix, as above
  static CLHEP::HepSymMatrix extrapolateCov(HelixParams &, double fltNew);

     // Path length (3-d) to intersection with cylinder at radius rad.
  static double fltToRad(const HelixParams& hel, double rad);

     // path length (3-d) to intersection with given Z plane, starting from estimate
  static bool findZFltlen(Trajectory const& traj, double zplane, double& zflt, double tol=0.1);

private:	
  // Preempt 
  TrkHelixUtils&   operator= (const TrkHelixUtils&);
  TrkHelixUtils(const TrkHelixUtils &);
	static const double small;
};

#endif

