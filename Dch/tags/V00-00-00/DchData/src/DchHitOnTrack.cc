// $Id: DchHitOnTrack.cc 89 2010-01-14 12:34:14Z stroili $
//
#include "BaBar/BaBar.hh"
#include "DchData/DchHitOnTrack.hh"
#include "DchData/DchHitBase.hh"

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
DchHitOnTrack::DchHitOnTrack(const TrkFundHit& fundHit,
                             const DchHitBase& baseHit,
                             int ambig,  double t0,
                             bool isProtoII)
  : TrkHitOnTrk(&fundHit,10.e-4),
    _ambig(ambig),
    _isProtoII(isProtoII),
    _hitTraj(baseHit.hitTraj()),
    _fitTime(baseHit.rawTime()-t0*1e-9),
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
  setHitLen(0.5 * baseHit.layer()->zLength());
  setFltLen(0.);
  _startLen = hitTraj()->lowRange() - 5.;
  _endLen   = hitTraj()->hiRange()  + 5.;
}

DchHitOnTrack::DchHitOnTrack(const DchHitOnTrack &hot,TrkRep *newRep,
                             const TrkDifTraj* trkTraj, const DchHitBase *hb)
  : TrkHitOnTrk(hot,newRep,trkTraj)
{
  _ambig = hot._ambig;
  _isProtoII = hot._isProtoII;
  _hitTraj  = hot._hitTraj;
  _fitTime  = hot._fitTime;
  _drift[0] = hot._drift[0];
  _drift[1] = hot._drift[1];
  _startLen = hot._startLen;
  _endLen   = hot._endLen;
  _dHit = (hb==0?hot._dHit:hb);
}

DchHitOnTrack::~DchHitOnTrack()
{ ; }
 
void DchHitOnTrack::print(std::ostream& o) const{
  //layer()->getWire(0)->print(o);
  _hitTraj->print(o);
  if(_poca) o<<"poca flt1 = "<<_poca->flt1()<<" flt2= "<<_poca->flt2()<<" doca "<<_poca->doca()<<"\n";
  baseHit()->print(o);
}

void
DchHitOnTrack::setT0(double t0)
{ 
  _fitTime= _dHit->rawTime()-t0*1e-9; 
}

double
DchHitOnTrack::dcaToWire() const
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
DchHitOnTrack::updateAmbiguity(double dca)
{
  if (dca < 0 && ambig() >= 0) {
    setAmbig(-1); return isActive();
  } else if (dca > 0 && ambig() <= 0) {
    setAmbig(1); return isActive();
  } else {
    return false;
  }
}

const DchHitOnTrack*
DchHitOnTrack::dchHitOnTrack() const
{
  return this;
}

double
DchHitOnTrack::entranceAngle() const
{
  static Hep3Vector dir;
  static HepPoint pos;
  if (getParentRep() == 0) return 0.;
  getParentRep()->traj().getInfo(fltLen(), pos, dir);
  return BbrAngle(dir.phi() - pos.phi());
}

unsigned
DchHitOnTrack::layerNumber() const
{
  return layernumber();
}

double
DchHitOnTrack::dipAngle() const
{
  return getParentRep()==0?0:Constants::pi/2-getParentRep()->traj().direction(fltLen()).theta();
}

TrkErrCode
DchHitOnTrack::updateMeasurement(const TrkDifTraj* traj, bool maintainAmb)
{
  TrkErrCode status=updatePoca(traj,true);
  if (status.failure()) {
    ErrMsg(warning) << "DchHitOnTrack::updateMeasurement failed " << status << endmsg;
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
DchHitOnTrack::updateCorrections()
{
  const TrkRep* tkRep = getParentRep();
  assert(tkRep != 0);
  double tof = tkRep->arrivalTime(fltLen());
  if ( _isProtoII ) tof = protoIItof(tof);
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
  _fitTime = _dHit->driftTime(tof);
  _drift[ambig()<0?0:1] = ambig() * dist;
  assert( driftCurrent() );

  double newSig = _dHit->sigma(dist, wireAmb, eAngle, dAngle, z);
  assert(newSig>0);
  setHitRms(newSig);
}

double
DchHitOnTrack::driftVelocity() const
{
  const TrkRep* tkRep = getParentRep();
  assert(tkRep != 0);
  double tof =  tkRep->arrivalTime(fltLen());
  if ( _isProtoII ) tof = protoIItof(tof);
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
DchHitOnTrack::timeResid(double &t, double &tErr) const
{
    double v = driftVelocity();
    if (v <= 0) return false;
    t = (fabs(drift())-fabs(dcaToWire()))/v;
    tErr= hitRms()/v;
    return true;
}

bool
DchHitOnTrack::timeAbsolute(double &t, double &tErr) const
{
  double tresid(-1.0);
  if(timeResid(tresid,tErr)){
// add back the track time
    t = tresid + getParentRep()->parentTrack()->trackT0();
    return true;
  } else
    return false;
}

double
DchHitOnTrack::protoIItof(double tof) const
{
    // Note explicit mass hypothesis here: 2 GeV muon!!!!!!!!!!!!
    double mass2   = .01116 ;
    double ptot2   = 4.;
    double betainv = sqrt( (ptot2 +  mass2)/ ptot2);
    // The flight for a protoII track is
    // fltProtoTrk = fltLen(radTriggerCounter) - fltLen()
    double radiusTr=69.0 ;                        // radius  Trigger Counter
    double zPM = 29.5;                            // z of photo multiplier
    HepTransformation unit;
    DetCylinder cyl(unit,radiusTr);
    double fltTrig(0);
    Intersection inters(getParentRep()->traj(),cyl);
    inters.intersect(fltTrig);
    //       Intersection(parentRep()->traj(),cyl).intersect(fltTrig);
    // Then take into account the time needed by the light to fly in the
    // Trigger counter
    // from the track trajectory to the PM - coordinates (0,0,29.5) cm
    // velocity in the scintillator = c/1.6
    //xyz track @ trigger count:
    HepPoint coordAtTrig= getParentRep()->position(fltTrig);
    HepPoint PMposition(0., radiusTr, zPM);
    double flInTrig=(coordAtTrig-PMposition).mag();
    double speedinv=betainv/Constants::c;
    double timeInScint = flInTrig*speedinv*1.6;
    return (fltTrig - fltLen())*speedinv -timeInScint;
}

TrkEnums::TrkViewInfo
DchHitOnTrack::whatView() const
{
  return _dHit->whatView();
}

int
DchHitOnTrack::layernumber() const
{
  return _dHit->layernumber();
}

const DchLayer*
DchHitOnTrack::layer() const
{
  return _dHit->layer();
}

int
DchHitOnTrack::wire() const
{
  return _dHit->wire();
}

int
DchHitOnTrack::whichView() const
{
  return _dHit->whichView();
}

double
DchHitOnTrack::rawTime() const
{
  return _dHit->rawTime();
}

double
DchHitOnTrack::charge() const
{
  return _dHit->charge();
}

const Trajectory*
DchHitOnTrack::hitTraj() const
{
    return _hitTraj;
}

const DchHit*
DchHitOnTrack::dchHit() const
{
  return 0;
}


// Replace underlying hit pointer (base class doesn't own it)

void 
DchHitOnTrack::changeBase(DchHitBase* newBase)
{
  _dHit = newBase;
}
