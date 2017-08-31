//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalRep.hh,v 1.134 2008/05/18 23:00:33 brownd Exp $
//
// Description:
//      Implementation class for TrkRep for a Kalman fit.  See documentation at
//      http://online04.lbl.gov/~brownd/tracking/chep2k_328.ps
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 1997	Lawrence Berkeley Laboratory
//
// Author List:
//      Dave Brown 3/15/97
//------------------------------------------------------------------------
#ifndef KALREP_HH
#define KALREP_HH
#include "BTrk/BbrGeom/HepPoint.h"
#include "BTrk/KalmanTrack/KalSite.hh"
#include "BTrk/KalmanTrack/KalEndSite.hh"
#include "BTrk/KalmanTrack/KalContext.hh"
#include "BTrk/TrkBase/TrkRep.hh"
#include "BTrk/TrkBase/TrkErrCode.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
#include "BTrk/TrkBase/TrkDirection.hh"
#include "BTrk/TrkBase/TrkDifPieceTraj.hh"
#include "BTrk/DetectorModel/DetIntersection.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include <iostream>

class TrkDifTraj;
class HelixParams;
class TrkHit;
class TrkFundHit;
class TrkVolume;
class TrkSimpTraj;
class KalMaterial;
class KalHit;
class KalEndSite;
class BField;
#include <vector>

class KalRep : public TrkRep {
public:
//
//
//  General interface.  The functions below are of general interest and
//  will not affect the state of the KalRep.
//
//
//******************************************
// Local helix parameters and trajectory access
//******************************************
  HelixParams helix(double fltLen) const;
  const TrkSimpTraj* localTrajectory(double fltlen,double& loclen) const;
  const TrkDifTraj& traj() const { return (TrkDifTraj&)pieceTraj(); }
  const TrkDifPieceTraj& pieceTraj() const;
//******************************************
// momentum information
//******************************************
  virtual int charge() const { return _charge; }
  virtual CLHEP::Hep3Vector        momentum(double fltL=0.)      const;
  virtual double            pt(double fltL=0.)            const;
  virtual BbrVectorErr      momentumErr(double fltL)      const;
// covariance matrices of the track at fixed flight length 
  virtual CLHEP::HepMatrix         posmomCov(double fltL)            const;
  virtual void              getAllCovs(double fltL,
				       CLHEP::HepSymMatrix& xxCov,
				       CLHEP::HepSymMatrix& ppCov,
				       CLHEP::HepMatrix& xpCov)      const;
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
					  CLHEP::HepVector& pos,
					  CLHEP::HepVector& mom,
					  CLHEP::HepSymMatrix& xxWeight,
					  CLHEP::HepSymMatrix& ppWeight,
					  CLHEP::HepMatrix&    xpWeight) const;
// override the definition of the valid flightlength. 
  virtual bool validFlightLength(double fltL,double tolerance=0.0)      const;
//******************************************
// Chisquared and related info
//******************************************
  double chisq() const;
// simply count Degrees of freedom
  virtual int nDof() const;
// nDof for a part of a track, up to the given flight length in the given direction
  virtual int nDof(double fltlen,trkDirection tdir) const;
//******************************************
// Special chi(squared) information, mostly for hits
//******************************************
// hit chisquared functions: these require the rep to have already been fit.
// the hot itself on the fit can be optionally excluded.  The activity state
// of the hot is ignored (ie the user should test on it before calling this function)
  double hitChisq(const TrkHit* hot,bool exclude=true) const;
// same as above, except returning chi (=residual/error).  Since this
// is signed, the return value indicates success separately.  The returned
// error includes the track error and the hit error (squared!!) (in chi space).
  bool hitChi(const TrkHit* hot,
	      double& chival,
	      double& chierr2,
	      bool exclude=true) const;
// 1-directional version of the above; the fit must be valid in the
// direction specified (but not necessarily in the other).  This gives
// only a _partial_ consistency test (the hits beyond the one under
// test in the specified direction are not taken into account).  The hit
// is always excluded in this case
  double hitChisq(const TrkHit* hot,trkDirection tdir) const;
// Compute the chisquared consistency of the forward and backward fits, at the
// specified flight length.  Optionally specify a vector of masks which specify
// which parameters to compare (default=all).
  double matchChisq(double fltlen,bool* tparams = 0) const;
// physical-dimension residual; this uses the above.  Note the different default
// value for 'exclude' compared to hitChi!!!!!
  bool resid(const TrkHit *hot,
             double& residual, double& residErr,
             bool exclude=false) const;
//******************************************
// diagnostic info
//******************************************
// amount of material traversed over a given range;
  double radiationFraction(double* range=0) const;
  double energyLoss(double* range=0) const;
// bending over a given range
  double thetaBend(double* range=0) const;
  double phiBend(double* range=0) const;
//******************************************
// cinematic info
//******************************************
  double transitTime(double flt0, double flt1) const;
//******************************************
// Printing; these can produce a _lot_ of output, so be careful!
//******************************************
  void printAll(std::ostream& ostr = std::cout) const;
  void print(std::ostream& ostr = std::cout) const;
//
//
//   Expert only section; the functions below allow changing the state of
//   a KalRep, and as such should only be used by experts.  Several expert-level
//   accessors are also defined here.
//
//
//******************************************
// Constructors and such
//******************************************
// constructor
 KalRep(TrkSimpTraj const&,TrkHitVector const& hots,
	std::vector<DetIntersection> const& dlist,
 	KalContext const& context,TrkParticle const& tpart,
	TrkT0 const& t0,double flt0);
// destructor
  virtual ~KalRep();
//******************************************
// Fitting and such
//******************************************
  virtual TrkErrCode fit(); // will fit the track to convergence, as
// defined in the KalContext parmeters
  TrkErrCode fit(trkDirection tdir); // fit only in a given direction.  The
// the result is not generally useable 
//  extend the track validity as specified, if possible
  TrkErrCode extendThrough(const TrkVolume&,trkDirection trkDir = trkOut );
// extend the rep up to a given flight length (if possible)
  TrkErrCode extendThrough(double fltlen);
// reset the fit: this makes the fit not current, and allows forcing additional iterations
  void resetFit();
// The same as above, but restoring the fit to its initial construction settings (ie seed traj).
// optionally invert the direction of the fit
  TrkErrCode resetAll(bool invert=false);
//******************************************
// Hit manipulation
//******************************************
  void addHit(TrkHit *);
  void removeHit(TrkHit *);
  void activateHit(TrkHit *);
  void deactivateHit(TrkHit *);
  void updateHit(TrkHit *);
  void updateHits();
  //  Use the current output trajectory to update the KalSites.
  //  This removes fitCurrent, but not fitValid.
  TrkErrCode updateSites();
//******************************************
// Constraint manipulation
//******************************************
  void addConstraint(const TrkSimpTraj* traj,double cfltlen,bool* cparams);
// the following will modify _all_ existing constraint sites on the rep.  It will
// also make the fit notcurrent IFF the constraints have actually changed.
  void setConstraints(bool* cparams);
// find nearest material site with rad thickness about the threshold
  const KalMaterial* findMaterial(double flt,double minrad) const;
// add a material site
  void addInter(DetIntersection const& detinter);
// add a generic site.  This function will take ownership of the site
  void addSite(KalSite* site);
//******************************************
// special-purpose trajectory access
//******************************************
  TrkDifTraj& traj() { return (TrkDifTraj&)pieceTraj(); }
  TrkDifPieceTraj& pieceTraj();
// reference trajectory
  const TrkDifPieceTraj* referenceTraj() const {return _reftraj; }
// filter results at a given flight length.  Valid trajectories must be
// passed in.  Return value is success.  Note that the inner trajectory
// represents the _outward_ filter result at the specified flightlength,
// the outer trajectory the _inward_ filter result.  Yes this makes sense,
// since the outer filter result is based on all the track info inwards of
// the flightlength, etc.
  bool filterTrajs(double fltlen,TrkSimpTraj* intraj,TrkSimpTraj* outtraj) const;
// access filter traj in a single direction
  bool filterTraj(double fltlen,trkDirection dir,TrkSimpTraj* traj) const;
// smoothed trajectory EXCLUDING the measurement content of the given site
  bool smoothedTraj(const KalHit* hit,TrkSimpTraj* traj) const;
  bool smoothedTraj(std::vector<KalSite*>::const_iterator ifirst,
    std::vector<KalSite*>::const_iterator isecond, TrkSimpTraj* traj) const;
// 
//******************************************
// iteration convergence testing
//******************************************
  int iterations() const { return _niter; }
  int intersections() const { return _ninter; }
// spatial trajectory separation between iterations, in terms of spatial distance
// and flightlength
  double maximumDistance() const { return _maxdist; }
  double maximumFlightDiff() const { return _maxfltdif; }
// evaluate the change in parameters between iterations: this computes
// a 'chisquared' using the most recent fit's covariance matrix, and the
// difference in parameter vectors.  Note that this is not a true chisquared
// as the parameters are 100% correlated, it just has the mathematical form of
// a chisquared, and forms a convenient dimensionless measure of convergence.
// Can be computed at either end of the track.
  double parameterDifference(trkDirection) const;
// estimate how much the momentum at the origin would change on the _next_
// iteration, given the difference between the current reference and fit
// momentum value, propogated through a crude dE/dx model
  double estimatedMomDiff() const;
// convergence test
  bool converged() const;
// pre-test before fitting, to catch crazys.  Sets a TrkErrCode appropriately if false
  bool isFitable(TrkErrCode& err) const;
// check DOFs are sufficient for KalContext requirements
  bool enoughDofs() const;
//******************************************
// reference quantities
//******************************************
  const TrkSimpTraj* seed() const {return _seedtraj;}
  const TrkSimpTraj* seedTrajectory() const {return _seedtraj;}
  double refMomentum() const { return _refmom; }
  double refMomFltLen()const {return _refmomfltlen;}
//******************************************
// status information
//******************************************
  bool needsFit(trkDirection idir) const { return !_siteflag[idir];}
  bool hasFit(trkDirection idir) const { return _siteflag[idir];}
  bool isExtendable(trkDirection idir) const { return _extendable[idir];}
// if the track stops, return the material site at which it stops
  const KalMaterial* stoppingSite() const { return _stopsite; }
// allow attempting to recover stopped tracks.  This will activate all materials past the
// stopsite (not including it), and reset the fit.  It does _not_ refit the track.
  TrkErrCode removeStopSite();
//******************************************
// miscelaneous
//******************************************
// diirectional chisquared; useful only as a diagnostic
  double chisquared(trkDirection idir=trkIn) const;
// partial chisquared, using hits only up to the given flight length in the
// given direction.
  double chisquared(double fltlen,trkDirection tdkr) const;
// find the sites which bound a given flightlength.  Can return null pointers if
// fltlen is not in range
  void findBoundingSites(double fltlen,const KalSite*& insite,const KalSite*& outsite) const;
// find the site with flightlength nearest to (but not beyond) the input.  Returns
// -1 if the input fltlen is before any site.
  int findNearestSite(double fltlen) const;
// find the first site of a given type
  const KalSite* findSite(KalSite::siteType stype) const;
// config information
  const KalContext& kalContext() const { return _kalcon; }
//******************************************
// access to KalSites
//******************************************
  const std::vector<KalSite*>& siteList() const { return _sites;}
// access first and last hit site
  const KalSite* firstHit() const {
    return _hitrange[0] >= 0 ? 
      _sites[_hitrange[0]] : 0; }
  const KalSite* lastHit() const {
    return _hitrange[1] < int(_sites.size()) ? 
      _sites[_hitrange[1]] : 0; }
// find hit site give the hit
  const KalHit* findHitSite(const TrkHit*) const;
// find material site(s) given an element
  void findMaterialSites(const DetElem*,std::vector<const KalMaterial*>& matsites) const;
// fill vectors of trajectories corresponding to fit requests and constraints
  void fitTrajectories(std::vector<TrkSimpTraj*>& fits,const char* stream="Default") const;
  void constraintTrajectories(std::vector<TrkSimpTraj*>& fits) const;
  double lowFitRange() const { return _fitrange[0];}
  double hiFitRange() const { return _fitrange[1];}
private:
  double _maxdist; // maximum distance between new and ref. trajectory, used to measure convergence
  double _maxfltdif; // maximum distance between new and ref. trajectory, used to measure convergence
protected:
  unsigned _niter; // number of times fit iterated
  unsigned _ninter; // number of times material was intersected
private:
  TrkDifPieceTraj* _ptraj; // piecewise trajectory
  TrkDifPieceTraj* _reftraj; // reference for derivatives + residuals
  KalContext const& _kalcon; // context for this rep
  TrkSimpTraj* _seedtraj; // seed trajectory for the fit
  KalMaterial* _stopsite;// stopping info

  std::vector<KalSite*> _sites; // sorted list of sites
  KalEndSite _endsites[2]; // sites to seed fit
  bool _siteflag[2]; // flag to keep track of the state of the sites
  double _fitrange[2]; // flight range to the edge of the tracking volume
  int _hitrange[2]; // range of sites from 1st to last hit.
  bool _extendable[2]; // is the track extendable in the different directions?
  double _chisq[2]; // chisquared sum for the different fit directions
  double _refmom; // momentum of reference trajectory
  double _refmomfltlen;  // flight length along reference traj where seed mom is valid
  int _charge; // particle charge
//
// private functions
//
  KalRep& operator = (const KalRep&);
  void setFit(trkDirection idir) {_siteflag[idir] = true;}
// set the fit range according to the specified volumes and trajectory
  void setFitRange(const Trajectory* traj);
// build the reference trajectory from the seed.  This fills the list of intersections
  void buildRefTraj();
// build the material sites
  void buildMaterialSites(double[2],std::vector<DetIntersection>& dinter);
// build the end sites
  void buildEndSites();
// return the site at which the track stops while adding material (or 0, if it doesn't stop)
  KalMaterial* createMaterialSites(std::vector<DetIntersection> tlist,
				   std::vector<KalSite*>& sites,
				   int lowindex,int highindex,
				   double mom,trkDirection dedxdir) const;
// build the hit sites
  void buildHitSites();
// build the bend sites
  void buildBendSites(double[2]);
  void createBendSites(double range[2],std::vector<KalSite*>& sites) const;
//  process sites in a given direction (= filtering).  The return
//  value tells whether any sites were changed.
  bool process(KalSite*,int,int,trkDirection);
// Code to do iteration in fit().  Can be overridded by subclasses
  virtual TrkErrCode iterateFit();
// update intersections if needed.  This is VERY EXPENSIVE
  void reIntersect();
// compute the momentum at the reference point, using the current fit result
  void updateRefMom();
// extend the sites, either for track extension, or for finishing the fit
// after iterating on the core
  TrkErrCode extendSites(int isite,trkDirection);
// similarly extend the piecewiset trajectory.
  TrkErrCode extendTraj(int isite,trkDirection);
// helper function for hot manipulations
  KalHit* findHitSite(const TrkHit*,unsigned& isite) const;
// underlying chisquared function
  double chisquared(int startsite,int nsites,trkDirection tdir) const;
// direct manipulation, for convenience
  void setFitRange(double newrange[2]);
  void setFitRange(double lowrange, double highrange);
  void setExtendable(trkDirection idir); 
// helper allowing similar constructors to share code
  void init(const HelixParams&);
// same, using seed parameters
  void initFromSeed();
	void initRep();
	void initSites();

//
//  Protected functions
//
protected:
  std::vector<KalSite*>& kalSites() { return _sites;}
  TrkDifPieceTraj* refTraj() { return _reftraj; }
  //  Build the trajectory from the sites and the old trajectory
  TrkErrCode buildTraj();
// build a 1-direction traj.  The valid range on this traj will be from the first/last
// hit inwards/outwards
  TrkErrCode buildTraj(trkDirection);
// special function to update end sites
  void updateEndSites(double smear,bool diagonly=false);
// direction-dependent version of the above
  void updateSites(int startindex,int endindex,
		   double initialmom,trkDirection dedxdir);
// cleanup after updating
  void fixupSites();
  // fit existing sites
  TrkErrCode fitSites(); 
  // allow subclasses to set extent
  void setUnextendable(trkDirection idir) {_extendable[idir] = false;}
  // allow subclasses to set fit range
  // find 1st + last hit site
  void findHitSites();
  bool stopsIn(const KalMaterial* mat) const;
  KalSite* nextActive(unsigned index,trkDirection tdir) const;

};
#endif
