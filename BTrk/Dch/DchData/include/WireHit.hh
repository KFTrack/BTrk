//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: WireHit.hh 89 2010-01-14 12:34:14Z stroili $
//
// Description:
//	Class Header  for |WireHit| 
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

#ifndef WIREHIT_HH
#define WIREHIT_HH
#include "BaBar/Constants.hh"
#include "DchData/DchHitBase.hh"
#include "TrkBase/TrkFundHit.hh"
#include "TrkBase/TrkEnums.hh"

class TrkDetElemId;
class GTrack;
class DchDetector;
class DchTimeToDistList;
class WireHitData;

class WireHit : public TrkFundHit
{
public:
  WireHit(double dist=0,double sigma=0);
  

  virtual ~WireHit();

  double phi() const                           {return _phi;}	// phi at chamber center
  double x() const                             {return _rmid * _cosphi;} // x at chamber center
  double x(double z) const;                                                // x  at global z
  double y() const                             {return _rmid * _sinphi;} // y at chamber center
  double y(double z) const;                                                // y at global z
  double rMid() const                          {return _rmid;} // R at chamber center
  double zlen() const                          {return _zlen;} // chamber extent in z

  virtual const GTrack* getGTrack() const;

  virtual TrkDetElemId elemId() const;

  int layernumber() const                      {return _layer;} // layer number
  int wire()  const                            {return _wire;}  // wire number

  unsigned status() const;

  // drift distances are *ALWAYS* positive; 
  // tof in s, rawtime in ns..
  inline double driftDist(double tof, int ambig, double entranceAngle,
                          double dipAngle, double z) const;
  // "fitting" sigma (multiplies by _fittingSigMult [= 1])
  inline double sigma(double driftdist, int ambig, double entranceAngle,
                      double dipAngle, double z) const;
  // "finding" sigma (also adds _addFindingSig in quadrature)
  inline double sigma(double driftdist, int ambig = 0) const;


  int whichView() const {return 1;}
  TrkEnums::TrkViewInfo whatView() const  { return whichView()==0?TrkEnums::xyView:TrkEnums::bothView; }

  //io
  void print(std::ostream &o) const;

private:
  double _dist;
  double _sigma;

  double _phi;			// phi of wire
  double _cosphi;		// cos(phi)
  double _sinphi;		// sin(phi)
  double _rmid;			// radius of wire at mid chamber
  double _zlen;			// wire length
  unsigned _status;             // status flag;
  int _wire;
  int _layer;

};
extern std::ostream& operator<<(std::ostream &o,const WireHit& aHit);

double 
WireHit::driftDist(double tof, int ambig, double entranceAngle,
                   double dipAngle, double z) const 
{ // tof in s, rawtime in ns..
  return _dist;
}

double
WireHit::sigma(double driftdist, int ambig) const
{
  return _sigma;
}

double
WireHit::sigma( double driftdist, int ambig, double entranceAngle,
                   double dipAngle, double z) const
{
  return _sigma;
}

#endif
