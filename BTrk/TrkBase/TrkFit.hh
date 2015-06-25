//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkFit.hh,v 1.15 2004/08/06 06:31:41 bartoldu Exp $
//
// Description:
//     Implements a few basic TrkAbsFit functions (most are left for TrkRep and
//     subclasses).
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner
//            modified by Justin Albert
//------------------------------------------------------------------------

#ifndef TRKFIT_HH
#define TRKFIT_HH

#include "BTrk/TrkBase/TrkAbsFit.hh"
#include "BTrk/TrkBase/TrkParticle.hh"

// Class interface //
class TrkFit : public TrkAbsFit {

public:
  virtual ChisqConsistency    chisqConsistency() const = 0;
  virtual bool validFlightLength(double fltL,double tolerance=0.0)      const;
  void printType(std::ostream& ostr) const;

  virtual int                 nActive()                   const = 0;
  virtual TrkParticle const&     particleType()              const = 0;
  virtual HelixParams    helix(double fltL)            const = 0;
  virtual double            arrivalTime(double fltL)      const = 0;
  virtual double            startFoundRange()             const = 0;
  virtual double            endFoundRange()               const = 0;

protected:
  TrkFit();
  virtual ~TrkFit();
private:
  // Preempt
  TrkFit&   operator= (const TrkFit&);
  TrkFit(const TrkFit &);
};
#endif
