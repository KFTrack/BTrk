// $Id: WireHitOnTrack.cc 89 2010-01-14 12:34:14Z stroili $
//
#include "BaBar/BaBar.hh"
#include "DchData/WireHitOnTrack.hh"
#include "DchData/DchHitBase.hh"
#include "DchData/DchHit.hh"
#include "DchData/WireHit.hh"

//#include "AbsEnv/AbsEnv.hh"
//#include "DchEnv/DchEnv.hh"

#include "BaBar/Constants.hh"
#include "BbrGeom/BbrAngle.hh"
#include "CLHEP/Matrix/Vector.h"

// DchGeom needed to verify if hit is inside of chamber...
// and to find the trajectory describing the hit, i.e. wire
#include "DchGeom/DchLayer.hh"
#include "DchGeom/DchSWire.hh"

#include "TrkBase/TrkPoca.hh"
#include "TrkBase/TrkDifTraj.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "TrkBase/TrkRep.hh"
#include "TrkBase/TrkSimpTraj.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "DetectorModel/DetCylinder.hh"
#include "DetectorModel/Intersection.hh"
#include "CLHEP/Geometry/Transformation.h"
#include "ErrLogger/ErrLog.hh"
using std::cout;
using std::endl;


//-------------
// Constructors
//-------------
WireHitOnTrack::WireHitOnTrack(const WireHit& baseHit,const Trajectory* wireTr,
			       int ambig)
  : TrkHitOnTrk(&baseHit,10.e-4),
    _ambig(ambig),
    _hitTraj(wireTr),
    _dHit(&baseHit)
{
        // need to flag somehow that that we haven't computed things yet...
        //   now, we know that _drift[0] is for ambig <0, and _drift[1] is ambig >0
        //   and _drift is signed according to ambig. Thus _drift[0] should be -tive
        //   and _drift[1] should be +tive. To indicate this are not yet initialized,
        //   put 'impossible' values here...
  _drift[0] = +9999;
  _drift[1] = -9999;
  setHitResid(-21212121.0);
  setHitRms( 0.02 );
  setHitLen(0.5 * 1000/*baseHit.layer()->zLength()*/);
  setFltLen(0.);
  _startLen = hitTraj()->lowRange() - 5.;
  _endLen   = hitTraj()->hiRange()  + 5.;
}

WireHitOnTrack::WireHitOnTrack(const WireHitOnTrack &hot,TrkRep *newRep,
                             const TrkDifTraj* trkTraj, const WireHit *hb)
  : TrkHitOnTrk(hot,newRep,trkTraj)
{
  _ambig = hot._ambig;
  _hitTraj  = hot._hitTraj;
  _drift[0] = hot._drift[0];
  _drift[1] = hot._drift[1];
  _startLen = hot._startLen;
  _endLen   = hot._endLen;
  _dHit = (hb==0?hot._dHit:hb);
}

WireHitOnTrack::~WireHitOnTrack()
{ ; }
 
TrkHitOnTrk*
WireHitOnTrack::clone(TrkRep *rep, const TrkDifTraj *trkTraj) const
{
  return new WireHitOnTrack(*this,rep,trkTraj);
}

void WireHitOnTrack::print(std::ostream& o) const{_dHit->print(o);}

double
WireHitOnTrack::dcaToWire() const
{
  double dca = -9999.;
  if ( getParentRep() == 0 ) {
//     cout << "no parent rep" << endl;
    return dca;
  }
  // WARNING: cannot use the internal _poca, as it lags one iteration
  //          behind the fit... therfore use _EXTERNAL_ residual
  if (isActive())  {  // FIXME: currently can only use 'resid()' if isActive..
    dca = resid()+drift();
  } else {
    TrkPoca poca(getParentRep()->traj(), fltLen(), *hitTraj(), hitLen(),
                 _tolerance);
    if (poca.status().success()) dca = poca.doca();
  }
  return dca;
}

bool
WireHitOnTrack::updateAmbiguity(double dca)
{
  if (dca < 0 && ambig() >= 0) {
    setAmbig(-1); return isActive();
  } else if (dca > 0 && ambig() <= 0) {
    setAmbig(1); return isActive();
  } else {
    return false;
  }
}

const WireHitOnTrack*
WireHitOnTrack::dchHitOnTrack() const
{
  return this;
}

double
WireHitOnTrack::entranceAngle() const
{
  static Hep3Vector dir;
  static HepPoint pos;
  if (getParentRep() == 0) return 0.;
  getParentRep()->traj().getInfo(fltLen(), pos, dir);
  return BbrAngle(dir.phi() - pos.phi());
}

unsigned
WireHitOnTrack::layerNumber() const
{
  return _dHit->layernumber();
}

double
WireHitOnTrack::dipAngle() const
{
  return getParentRep()==0?0:Constants::pi/2-getParentRep()->traj().direction(fltLen()).theta();
}

TrkErrCode
WireHitOnTrack::updateMeasurement(const TrkDifTraj* traj, bool maintainAmb)
{
  TrkErrCode status=updatePoca(traj,true);
  if (status.failure()) {
    ErrMsg(warning) << "WireHitOnTrack::updateMeasurement failed " << status << endmsg;
    return status;
  }
  assert (_poca!=0);
  double dca=_poca->doca();
  bool forceIteration = (maintainAmb&&ambig()!=0)?false:updateAmbiguity(dca);
  assert(ambig()!=0);
  // Check for hits beyond end plates.  !!Turn off hit if it is == temp. hack
  if (isBeyondEndflange()) setUsability(false);
  if (forceIteration || !driftCurrent() ) {
    updateCorrections(); // force recomputation of drift for current ambig(), setting of hitRms
    forceIteration=true;
  }
  setHitResid(dca-drift());
  return !forceIteration?status:
    TrkErrCode(TrkErrCode::succeed, 11, "Ambiguity flipped");
}

void
WireHitOnTrack::updateCorrections()
{
  const TrkRep* tkRep = getParentRep();
  assert(tkRep != 0);
  double tof = tkRep->arrivalTime(fltLen());
  // at this point, since dcaToWire is computed, _ambig must be either -1 or +1
  assert( ambig() == -1 || ambig() == 1 );
  static HepPoint pos; static Hep3Vector dir;
  _trkTraj->getInfo(fltLen(), pos, dir);
  double eAngle = BbrAngle(dir.phi() - pos.phi());
  double dAngle =Constants::pi/2 - dir.theta();
  double z = pos.z();
  // note the special case for INCOMING tracks:
  //    wire left of track implies hit on the right side of wire
  int wireAmb= fabs(eAngle)<Constants::pi/2?_ambig:-_ambig;
  // provide the underlying hit with the *external* information
  // needed to compute the drift distance, i.e. those numbers that
  // the hit cannot figure out by itself...
  double dist =  _dHit->driftDist(tof, wireAmb, eAngle, dAngle, z); assert(dist>0);
  _drift[ambig()<0?0:1] = ambig() * dist;
  assert( driftCurrent() );

  double newSig = _dHit->sigma(dist, wireAmb, eAngle, dAngle, z);
  assert(newSig>0);
  setHitRms(newSig);
}

double
WireHitOnTrack::driftVelocity() const
{
  const TrkRep* tkRep = getParentRep();
  assert(tkRep != 0);
  double tof =  tkRep->arrivalTime(fltLen());
  static HepPoint pos; static Hep3Vector dir;
  _trkTraj->getInfo(fltLen(), pos, dir);
  double eAngle = BbrAngle(dir.phi() - pos.phi());
  double dAngle =Constants::pi/2 - dir.theta();
  double z = pos.z();
  // note the special case for INCOMING tracks:
  //    wire left of track implies hit on the right side of wire
  int wireAmb= fabs(eAngle)<Constants::pi/2?_ambig:-_ambig;

  static const double epsilon = 0.5e-9; // tof are in s;
  double dist1 =  _dHit->driftDist(tof+epsilon, wireAmb, eAngle, dAngle, z);
  double dist2 =  _dHit->driftDist(tof-epsilon, wireAmb, eAngle, dAngle, z);

  return (dist2-dist1)/(2*epsilon); // velocity in cm/s
}

bool
WireHitOnTrack::timeResid(double &t, double &tErr) const
{
    double v = driftVelocity();
    if (v <= 0) return false;
    t = (fabs(drift())-fabs(dcaToWire()))/v;
    tErr= hitRms()/v;
    return true;
}

bool
WireHitOnTrack::timeAbsolute(double &t, double &tErr) const
{
  double tresid(-1.0);
  if(timeResid(tresid,tErr)){
// add back the track time
    t = tresid + getParentRep()->parentTrack()->trackT0();
    return true;
  } else
    return false;
}

TrkEnums::TrkViewInfo
WireHitOnTrack::whatView() const
{
  return _dHit->whatView();
}

int
WireHitOnTrack::wire() const
{
  return _dHit->wire();
}

int
WireHitOnTrack::whichView() const
{
  return _dHit->whichView();
}

const Trajectory*
WireHitOnTrack::hitTraj() const
{
    return _hitTraj;
}

unsigned
WireHitOnTrack::status() const
{
  return _dHit->status();
}

