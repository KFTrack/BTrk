//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkRep.hh,v 1.69 2007/09/24 21:56:27 gapon Exp $
//
// Description: Base class for internal track representation classes -- e.g. 
//   HelixRep, KalRep.  Owns and maintains a TrkHotList; and keeps a 
//   pointer to the track that owns the Rep.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Authors: Steve Schaffner
//
//------------------------------------------------------------------------

#ifndef TRKREP_HH
#define TRKREP_HH
#include <memory>
#include "BTrk/TrkBase/TrkDirection.hh"
#include "BTrk/TrkBase/TrkFitStatus.hh"
#include "BTrk/TrkBase/TrkFit.hh"
#include "BTrk/TrkBase/TrkHotList.hh"
#include "BTrk/TrkBase/TrkHitOnTrkUpdater.hh"
#include "BTrk/TrkBase/TrkT0.hh"
#include "BTrk/TrkBase/TrkId.hh"
// These 3 are needed for the OSF compiler:
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/BbrGeom/BbrPointErr.hh"
#include "CLHEP/Vector/ThreeVector.h"
// the following is needed by the Sun compiler
#include "BTrk/ProbTools/ChisqConsistency.hh"

class TrkHitOnTrk;
#include <iosfwd>
class TrkDifTraj;
class TrkErrCode;
class HelixParams;
#include "CLHEP/Vector/ThreeVector.h"
class HepPoint;
class TrkVolume;
class BField;

// Class interface //
class TrkRep : public TrkFitStatus, public TrkFit, public TrkHitOnTrkUpdater {

  public:
    //******************************************
    // Constructors and such
    //******************************************
    // construct from a hotlist -- The rep will not take ownership of anything passed
    //                          -- it will _clone_ the hots passed on inHots
    TrkRep(const TrkHotList& inHots, 
	TrkParticle const& tpart);
    TrkRep(TrkHotList& inHots, 
	TrkParticle const& tpart, bool stealHots=false);
    // construct from a hotlist -- The rep will always _TAKE OWNERSHIP_ of the hots
    //                             and if takeownership==true, ALSO of the list.
    TrkRep(const TrkHotList* inHots, 
	TrkParticle const& tpart);
    TrkRep(TrkHotList* inHots, 
	TrkParticle const& tpart, bool takeownership=false);
    // rep without explicit hotlist
    TrkRep(TrkParticle const& tpart, bool createHotList=false);
    // copy ctor
    TrkRep(const TrkRep& oldRep, TrkParticle const& tpart);
    virtual ~TrkRep();
    // clone() used to copy tracks; cloneNewHypo() for new hypos within track
    virtual TrkRep* clone() const = 0;
    virtual TrkRep* cloneNewHypo(TrkParticle const& tpart) = 0;
    bool operator== (const TrkRep&);

    //******************************************
    // Global quantities:
    //******************************************
    virtual ChisqConsistency    chisqConsistency() const;
    virtual int               nActive()      const;
    double                    startValidRange() const;
    double                    endValidRange()   const;
    virtual double            startFoundRange() const;
    virtual double            endFoundRange()   const;
    //  t0
    TrkT0 const& t0() const { return _trkt0; }
    // flightlength associated with t0.
    double flt0() const { return _flt0; }
    TrkParticle const& particleType() const { return _tpart; }
    // for now, allow setting t0.  Must be coherent with flt0
    void setT0(TrkT0 const& t0,double flt0) { _trkt0 = t0; _flt0 = flt0;}

    //******************************************
    // Information about track at a given position
    //******************************************
    virtual HepPoint                  position(double fltL)         const;
    virtual CLHEP::Hep3Vector                direction(double fltL)        const;
    virtual double                    arrivalTime(double fltL)      const;
    virtual BbrPointErr               positionErr(double fltL)      const;
    virtual BbrVectorErr              directionErr(double fltL)     const;

    //******************************************
    // Hit (list) handling
    //******************************************
    // Simple implementations of these are present in the base class; 
    //   complicated reps (e.g. Kalman) may wish to override.
    virtual void            addHot(TrkHitOnTrk *theHot);
    virtual void            removeHot(TrkHitOnTrk *theHot);
    virtual void            activateHot(TrkHitOnTrk *theHot);
    virtual void            deactivateHot(TrkHitOnTrk *theHot);
    virtual TrkHotList*             hotList()                 {return _hotList.get();}
    virtual const TrkHotList*       hotList() const           {return _hotList.get();}
    virtual void            updateHots();
    virtual bool            resid(const TrkHitOnTrk *theHot,
	double &residual, double &residErr,
	bool exclude=false) const;

    // Distinguishes hotLists that actually hold hits from those that just hold 
    //   info (e.g. nActive) about the hits used to create the fit
    bool hitCapable() const                     {return hotList()->hitCapable();}


    //****************************************** 
    // Fitting stuff
    //******************************************
    virtual TrkErrCode fit() = 0;

    TrkId const& id() const { return _id; }
  protected:
    TrkRep&   operator= (const TrkRep&);
  private:
    // new data members for mu2e
    TrkT0 _trkt0; // t0 value
    double _flt0; // flightlength associated with t=t0.
    TrkParticle _tpart;
    TrkId _id;  // track identifier
  protected:
    //protected, not private, so derived classes can create in copy ctor
    // note: cloning a hotlist requires parentRep to be set first, so
    // this must be declared _after _parentRep
    std::auto_ptr<TrkHotList>  _hotList;
};

#endif
