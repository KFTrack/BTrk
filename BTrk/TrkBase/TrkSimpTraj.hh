//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkSimpTraj.hh,v 1.32 2007/09/24 21:56:27 gapon Exp $
//
// Description:  TrkSimpTraj is an abstract class for trajectories that 
//     contain a simple set of parameters (i.e. it's not a composite 
//     trajectory).  It is a base class for things like HelixTraj and 
//     LineTraj, and inherits from TrkDifTraj.  The reference point (what
//     the parameters implicitly use to relate to BaBar 3-space) is now
//     an explicit member, as opposed to being implicitly set to 0,0,0.
//
// Note: the HP compiler is too stupid to process the default ctor, so instead
// the origin will be hard-coded to 0,0,0 for the HP (only).  This could
// well cause problems for HP users someday, but complaints should be addressed
// to the manufacturer of those boxes, not me.
//
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner, Dave Brown
//
//------------------------------------------------------------------------

#ifndef TRKSIMPTRAJ_HH
#define TRKSIMPTRAJ_HH

#include "BTrk/TrkBase/TrkDifTraj.hh"
#include "BTrk/TrkBase/TrkKalDeriv.hh"
#include "BTrk/TrkBase/TrkParams.hh"
#include "BTrk/BbrGeom/HepPoint.h"
#include <vector>

class TrkVisitor;
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include <iosfwd>

//  Prototype the point translation functions.  This is kind of ugly
typedef void (*TranslateParams)(const HepPoint& oldpoint,
                                const HepPoint& newpoint,
                                const CLHEP::HepVector& oldpar,
                                const CLHEP::HepSymMatrix& oldcov,
                                CLHEP::HepVector& newpar,
                                CLHEP::HepSymMatrix& newcov,
                                double fltlen);

// Class interface //
class TrkSimpTraj : public TrkDifTraj, public TrkKalDeriv {

public:

  static HepPoint _theOrigin; // define the origin as a HepPoint


public:

  //-----------------------
  // Constructors and such
  //-----------------------
  TrkSimpTraj(const CLHEP::HepVector& params, const CLHEP::HepSymMatrix& cov,
              const double startRange = -99999.,const double endRange =99999.,
              const HepPoint& refpoint = _theOrigin);
  TrkSimpTraj(const TrkParams& params,
              const double startRange = -99999.,const double endRange =99999.,
              const HepPoint& refpoint = _theOrigin);
  virtual ~TrkSimpTraj();

  virtual TrkSimpTraj* clone() const = 0;

  //--------------------------------------------------
  //  Access to parameters, errors and reference point
  //--------------------------------------------------
  TrkParams*                 parameters()                   {return &_dtparams;}
  const TrkParams*           parameters() const             {return &_dtparams;}
  virtual const TrkSimpTraj* localTrajectory(double fltLen, double& localFlt)
    const;
  const HepPoint&            referencePoint() const         {return _refpoint;}
  virtual void               print(std::ostream& os) const;
  virtual                    void printAll(std::ostream& os) const;

  //--------------------------------------------------
  // Change contents
  //--------------------------------------------------
  // Change the reference point and the parameters
  void           changePoint(const HepPoint& newpoint,double& fltlen);
  // Set the ref point and don't change the params
  void           setPoint(const HepPoint& newpoint)    {_refpoint = newpoint;}
  // inversion function: this inverts both the flight range and the parameters
  // so that the same traj is described but going in the opposite direction.
  TrkSimpTraj& invert();
  // invert the track parameters passed in newparams.
  // Returns true in flags if the inversion requires a sign change in the
  // covariance matrix as well.
  virtual void invertParams(TrkParams* newparams, std::vector<bool>& flags) const = 0;
  // Provide function which translates the reference point of parameters
  virtual TranslateParams paramFunction() const = 0;

  //--------------------------------------------------
  // Visitor access for momentum functions
  //--------------------------------------------------

  virtual void visitAccept(TrkVisitor* vis) const = 0;

  bool operator==(const TrkSimpTraj&) const; // return equivalence, not identy


protected:
  TrkParams _dtparams;
  HepPoint _refpoint; // reference point for parameters
private:
  // Preempt 
  TrkSimpTraj(const TrkSimpTraj &);
  TrkSimpTraj&   operator= (const TrkSimpTraj&);
};
#endif
