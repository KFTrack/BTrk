//--------------------------------------------------------------------------
// File and Version Information:
//
// Description:
//     Defines output interface for any "fitted" track (recon'ed charged,
//     hypothesized charged, or hypothesized uncharged).
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner, Justin Albert
//------------------------------------------------------------------------

#ifndef TRKABSFIT_HH
#define TRKABSFIT_HH

class HepPoint;
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
class BbrPointErr;
class BbrVectorErr;
class TrkDifTraj;
class TrkSimpTraj;
#include <iosfwd>
class HelixParams;
class ChisqConsistency;


// Class interface //
class TrkAbsFit {

public:
  //********************************
  //Global track quantities:
  //********************************
  virtual int                 charge()                    const = 0;
  virtual double              chisq()                     const = 0;
  virtual int                 nDof()                      const = 0;
  virtual const TrkDifTraj&   traj()                      const = 0;


  //********************************
  //Information about track at a given position (flight length)
  //********************************
  virtual HepPoint          position(double fltL)         const = 0;
  virtual Hep3Vector        direction(double fltL)        const = 0;
  virtual Hep3Vector        momentum(double fltL=0.)      const = 0;
  virtual double            pt(double fltL=0.)            const = 0;        
  virtual BbrPointErr       positionErr(double fltL)      const = 0;
  virtual BbrVectorErr      directionErr(double fltL)     const = 0;
  virtual BbrVectorErr      momentumErr(double fltL)      const = 0;
  // Interface to vertexing algorithms (M.Bondioli 7/17/98)
  // covariance matrices of the track at fixed flight length
  virtual HepMatrix         posmomCov(double fltL)        const = 0;
  virtual void              getAllCovs(double fltL,
                                       HepSymMatrix& xxCov,
                                       HepSymMatrix& ppCov,
                                       HepMatrix& xpCov)  const = 0;

  // accessors to 2nd derivative of chi2 wrt x and p.
  // x0 and p0 are filled with the pos and 3mom around which expansion
  // takes place, whilst Weights are filled with the 2nd deriv of chi2.
  // so that:
  //  dx=x-x0   dp=p-p0
  //
  //  chi2(x,p)=0.5 dx^t*Wxx*dx + dx^t*Wxp*dp + 0.5 dp^t*Wpp*dp
  // where:
  //  pos and mom are 3-dim Vectors,
  //  xxWeight ppWeight and xpWeight are 3 by 3 matrices 
  //
  virtual void              getAllWeights(double fltL,
                                          HepVector& pos,
                                          HepVector& mom,
                                          HepSymMatrix& xxWeight,
                                          HepSymMatrix& ppWeight,
                                          HepMatrix&    xpWeight)
                                            const = 0;


  //********************************
  // Valid flight range
  //********************************
  //  validRange = trajectory's valid range = flight length range 
  //   over which the track can be used. 
  //  foundRange = range over which the track is known to 
  //   exist (i.e. from first hit to last hit).  See derived class
  //   TrkFit, which is the base for all recon'ed track fits.

  virtual double            startValidRange()             const = 0;
  virtual double            endValidRange()               const = 0;

  //******************************************
  // Printing
  //******************************************
  virtual void printAll(std::ostream& ostr) const = 0;
  virtual void print(std::ostream& ostr) const = 0;

protected:
  TrkAbsFit();
  virtual ~TrkAbsFit();

private:
  // Preempt
  TrkAbsFit&   operator= (const TrkAbsFit&);
  TrkAbsFit(const TrkAbsFit &);
};

#endif
