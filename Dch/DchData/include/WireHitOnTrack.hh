//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: WireHitOnTrack.hh 89 2010-01-14 12:34:14Z stroili $
//
// Description:
//   Contains drift chamber hit info, as hit is used on a particular track
//   Inherits from TrkHitOnTrk.  The drift distance is stored as an 
//   absolute value, but returned as |drift|*ambig.  Ambiguity stored as +/- 1; 
//   ambiguity = 0 => pick better value @ first call to updateMeasurement.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Authors: Steve Schaffner
//
// Revision History:
//	20011018  M. Kelsey -- Make "rawTime()" public.
//	20030923  M. Kelsey -- Add function to replace _dHit pointer
//------------------------------------------------------------------------

#ifndef WIREHITONTRACK_HH
#define WIREHITONTRACK_HH

#include "TrkBase/TrkHitOnTrk.hh"
#include "TrkBase/TrkEnums.hh"
#include "BaBar/Constants.hh"
#include <math.h>

class WireHit;
class DchLayer;
class Trajectory;
class TrkRecoTrk;
#include "CLHEP/Matrix/Vector.h"

class WireHitOnTrack : public TrkHitOnTrk {

public:
  WireHitOnTrack(const WireHit& hit,const Trajectory* wireTr=0,int ambig=0);
  virtual ~WireHitOnTrack();

  virtual TrkHitOnTrk* clone(TrkRep *, const TrkDifTraj *trkTraj=0) const;


  // WireHitOnTrack specific functions
  double  entranceAngle() const;
        // the entrance Angle is the difference in phi between the *direction* 
        // of the track at the hit, and the phi of the *location* of hit
        // Signing convention is such that tracks going to the 'left' 
        // are positive and the ones going to the 'right' are negative
        // The definition was chosen so that radially outgoing tracks have 
        // zero entrance angle (by construction), regardless of 'doca'

  double dipAngle() const;
        // dipangle is just pi/2 - theta of the track at the location of the hit

  // In general ambiguities are: -1 -> right, 0 -> don't know, +1 -> left
  //   note: maybe this should just be an 'enum' or even a seperate class...
  // Note the special case of incoming tracks (i.e. |entranceAngle()|>pi/2) 
  // where a track having a wire on the LEFT ( ambig()=+1) has the hit on the
  // RIGHT ( wireAmbig()=-1 ) of the wire
  int      ambig()  const  { return _ambig; }  // wire wrt track direction
  int     wireAmbig()  const   { // hit wrt the wire location
    return fabs(entranceAngle())<Constants::pi/2?ambig():-ambig();} 

  //   if ambiguity unknown, return the average of the 
  //   absolute value of both ambiguities, i.e. pick a positive number...
  double   drift()  const  { return _ambig!=0 ? _drift[_ambig<0 ? 0:1] 
                                              : (_drift[1]-_drift[0])*0.5; }
  double  drift(double dca)  const  { return _drift[dca<0?0:1]; }
  double  dcaToWire() const;

  //   generic virtual functions (required by TrkHitOnTrk)
  virtual const Trajectory* hitTraj() const;
  virtual const WireHitOnTrack* dchHitOnTrack() const;
  virtual bool timeResid(double& t,double& tErr) const;
  virtual bool timeAbsolute(double& t,double& tErr) const;

  double driftVelocity() const;
  // specific virtual functions (required by WireHitOnTrack)
  unsigned                      status()          const;

  // Forwarded to DchHitBase
  int              wire()     const ;
  unsigned       layerNumber() const ;
  int              whichView() const; // 0 for axial, +/- 1 for stereo

  TrkEnums::TrkViewInfo whatView() const;

  // Set used during persistant -> transient and internally
  void                     setAmbig(int a) { _ambig = a<0?-1:a>0?1:0; }

  void print(std::ostream& ) const;
protected:

  WireHitOnTrack(const WireHitOnTrack &hitToBeCopied, TrkRep *newRep,
                const TrkDifTraj* trkTraj,const WireHit *hb=0);

  bool              isBeyondEndflange() const 
                    { return (hitLen() < _startLen || hitLen() > _endLen); }
  // return forceIteration: ambiguity flipped && hit is active
  bool              updateAmbiguity(double dca); 

  virtual TrkErrCode updateMeasurement(const TrkDifTraj* traj, bool maintainAmbiguity);

private:
  friend class DchHOTData;

  void              updateCorrections();
  bool              driftCurrent() const { return ambig()*_drift[ambig()<0?0:1]>0; }


  //Data members
  int               _ambig;  // this is the LR ambiguity wrt the TRACK 
                             // direction;
  double            _drift[2];  // corrected version of what's in the FundHit,
                                // one for each ambiguity: left, right
  const Trajectory* _hitTraj;

  //   cached information to improve code speed
  double            _startLen;  // start hitlen traj
  double            _endLen;    // end hitlen traj
  const WireHit *_dHit;


  //  hide copy constructor and assignment operator
  WireHitOnTrack(const WireHitOnTrack&);
  WireHitOnTrack& operator=(const WireHitOnTrack&);
};

#endif
