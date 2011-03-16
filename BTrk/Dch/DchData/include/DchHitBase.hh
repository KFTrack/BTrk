//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchHitBase.hh 89 2010-01-14 12:34:14Z stroili $
//
// Description:
//	Class Header  for |DchHitBase| 
//	DchHitBase contains the functionality provided by stuff from
//	the environment which is required from a DchHit by DchHitOnTrack
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

#ifndef DCHHITBASE_HH
#define DCHHITBASE_HH
//---------------------
//    C++ Headers    --
//---------------------
#include <math.h>

//-------------------------------------
//    Collaborating Class Headers    --
//-------------------------------------
#include "DchGeomBase/DchCellAddr.hh"
#include "TrkBase/TrkEnums.hh"
#include "DchCalib/DchTimeToDist.hh"

//----------------------------
//    Class Declarations    --
//----------------------------
class DchDetector;  
class DchLayer;
class DchTimeToDist;
class DchTimeToDistList;
class Trajectory;

class DchHitBase 
{
public:

  DchHitBase(unsigned layer, unsigned wire, double rawTime, unsigned iTdc, double charge,
             const DchDetector&, const DchTimeToDistList&);
  virtual ~DchHitBase();

  bool operator==(const DchHitBase& rhs) const { return (this == &rhs); }

  const DchLayer*    layer() const             {return _layerPtr;}
  const DchDetector* geom()  const             {return _geomPtr;}
  int layernumber() const                      {return _layer;} // layer number
  int wire()  const                            {return _wire;}  // wire number

  double charge() const                        { return _charge;}
  double rawTime() const                       { return _rawTime;}
  unsigned tdcIndex() const                       { return _iTdc;}
  double driftTime(double tof) const           { return _rawTime-1.e9*tof;}
  const Trajectory* hitTraj() const;

  // drift distances are *ALWAYS* positive; 
  // tof in s, rawtime in ns..
  inline double driftDist(double tof, int ambig, double entranceAngle,
                          double dipAngle, double z) const;
  // "fitting" sigma (multiplies by _fittingSigMult [= 1])
  inline double sigma(double driftdist, int ambig, double entranceAngle,
                      double dipAngle, double z) const;
  // "finding" sigma (also adds _addFindingSig in quadrature)
  inline double sigma(double driftdist, int ambig = 0) const;


  int whichView() const;
  TrkEnums::TrkViewInfo whatView() const  { return whichView()==0?TrkEnums::xyView:TrkEnums::bothView; }

protected:
  double _rawTime;
  double _charge;
  const DchDetector* _geomPtr;  // pointer to geometry
  const DchLayer*   _layerPtr;  // pointer to layer
  const DchTimeToDist*   _t2d;  // t->d calibration; needs to be AFTER _layer and _wire 
                                //                   as those are used for initialization 
                                //                   of _t2d...
  unsigned _layer;
  unsigned _wire;
  unsigned _iTdc;

  static double _addFindingSig; //add to nom. sigma in quad => finding sig
  static double _fittingSigMult;//scale factor to manipulate nom. sigma


  DchHitBase(const DchHitBase& rhs):
          _rawTime(rhs._rawTime), _charge(rhs._charge),
          _geomPtr(rhs._geomPtr), _layerPtr(rhs._layerPtr),
          _t2d(rhs._t2d), _layer(rhs._layer),
          _wire(rhs._wire), _iTdc(rhs._iTdc)
  {
          ;
  }
  DchHitBase& operator=(const DchHitBase& rhs){
    if(this != &rhs){
      _layer = rhs._layer;
      _wire = rhs._wire;
      _geomPtr = rhs._geomPtr;
      _layerPtr = rhs._layerPtr;
      _t2d = rhs._t2d;
      _rawTime = rhs._rawTime;
      _iTdc = rhs._iTdc;
      _charge = rhs._charge;
    }
    return *this;
  }

private:
  friend class DchMakeHits;
  static void setSigScale(double addFind, double fitMult) {
    _addFindingSig = addFind;
    _fittingSigMult = fitMult;}

  //hide the copy ctor and assignment op, at least until somebody needs them
  DchHitBase();

};


// FIXME: these fcn's were inlined, but DchHitBase.hh is pulled into some DDL files, 
//        and the DDL compiler cannot deal with DchCalib/DchTimeToDist.hh...
//        As a result, these fcn's are NO longer inline...
// NOTE: there is now a precomiler hack in DchDataP/DchMiniHitListP.ddl that
//       avoids the above ooddlx problem. As a result, restore these inlines
//// inline bunch of often called small functions for speed
double 
DchHitBase::driftDist(double tof, int ambig, double entranceAngle,
                   double dipAngle, double z) const 
{ // tof in s, rawtime in ns..
  return _t2d->timeToDist(driftTime(tof),ambig,entranceAngle,dipAngle,z);
}

double
DchHitBase::sigma(double driftdist, int ambig) const
{
  double sig = _fittingSigMult * _t2d->resolution(driftdist,ambig,0,0,0);
  return sqrt(sig*sig + _addFindingSig * _addFindingSig);
}

double
DchHitBase::sigma( double driftdist, int ambig, double entranceAngle,
                   double dipAngle, double z) const
{
  return _fittingSigMult * _t2d->resolution(driftdist,ambig,entranceAngle,
                                            dipAngle,z);
}

#endif
