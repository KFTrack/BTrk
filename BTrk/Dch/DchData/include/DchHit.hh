//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchHit.hh 89 2010-01-14 12:34:14Z stroili $
//
// Description:
//	Class Header  for |DchHit| 
//      Fundamental Hit class for use by drift chamber fits and pattern
//      recogition code
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      S. Schaffner - Original
//	A. Snyder - Modifications to use |DchGeom| and to construct
//                  from |DchDigi|s
//
// Copyright Information:
//	Copyright (C) 1996	SLAC
//
//------------------------------------------------------------------------

#ifndef DCHHIT_HH
#define DCHHIT_HH
#include "BaBar/Constants.hh"
#include "DchData/DchHitBase.hh"
#include "TrkBase/TrkFundHit.hh"
#include "TrkBase/TrkEnums.hh"

class TrkDetElemId;
class GTrack;
class DchDigi{};
class DchDetector;
class DchTimeToDistList;
class DchHitData;

class DchHit : public DchHitBase, public TrkFundHit
{
public:
  /*  DchHit(const DchDigi& digi,
         const DchDetector& det,
         const DchTimeToDistList &t2d,
         unsigned iTdc=0);*/
  /* DchHit(int layer,int wire, double tdctime,
         const DchDetector& det,
         const DchTimeToDistList &t2d);*/

  DchHit(int layer,int wire,double dist,double sigma,const DchDetector& det);
  DchHit();

  virtual ~DchHit();
// dnb 12/14/01 I need these for hit filtering (unassigned hit persistence)
// compiler implementation is OK.
  DchHit& operator=(const DchHit&);
  DchHit(const DchHit&);
//
  bool operator==(const DchHit&) const;

  //io
  void print(std::ostream &o) const;

  const DchDigi*     digi()  const             {return _digiPtr;}
  double phi() const                           {return _phi;}	// phi at chamber center
  double x() const                             {return _rmid * _cosphi;} // x at chamber center
  double x(double z) const;                                                // x  at global z
  double y() const                             {return _rmid * _sinphi;} // y at chamber center
  double y(double z) const;                                                // y at global z
  double rMid() const                          {return _rmid;} // R at chamber center
  double zlen() const                          {return _zlen;} // chamber extent in z
  // drift distances are *ALWAYS* positive; (they get their sign in 
  // DchHitOnTrack, not here: what would we do if a '0' ambig is specified?)
  double driftDist(double   bunchTime, int ambig) const  
          {return driftDist(bunchTime+crudeTof(), ambig, 0, 0, 0 );}

  // arrivalTime in seconds
  double driftDist(double arrivalTime, int ambig, double entranceAngle,
                   double dipAngle, double z) const;
  // tof in seconds
  // driftTime returned in _ns_ (yes, NANOseconds)
  double driftTime(double tof) const
  { return rawTime()-1.e9*tof; } // z dependence has moved into the t->d code

  virtual const GTrack* getGTrack() const;

  virtual TrkDetElemId elemId() const;
  TrkEnums::TrkViewInfo whatView() const  { return DchHitBase::whatView();}

  unsigned status() const;

private:
  const DchDigi* _digiPtr;      // pointer to digi (0 when running on mini!)

  double _phi;			// phi of wire
  double _cosphi;		// cos(phi)
  double _sinphi;		// sin(phi)
  double _rmid;			// radius of wire at mid chamber
  double _zlen;			// wire length
  unsigned _status;             // status flag;


  // Calculate crude (no track info) correction of time for flight 
  // delay (assumes tracks from origin); wire propagation assumes z=0;
  // return time in seconds
  double crudeTof() const { return  _rmid/Constants::c; }

  //hide the copy ctor and assignment op, at least until somebody needs them
  //DchHit();
};
extern std::ostream& operator<<(std::ostream &o,const DchHit& aHit);

// Might need this again someday:
  // Spawn a HitOnTrk object (with pointer to FundHit of appropriate type)
  //  TrkHitOnTrk* makeHot(TrkRecoTrk *track);

#endif
