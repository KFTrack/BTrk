//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchLayer.hh 123 2010-04-29 14:41:45Z stroili $
//
// Description:
//	Class DchLayer.
//      Dch layer class, it's just a "navigation" class for pattern
//      recognition
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	R. Stroili		originator
//	
//
// Copyright Information:
//	Copyright (C) 1997	INFN & Padova University
//
// Revision History:
//	20020619  M. Kelsey -- Break circular dependence w/DchEnv by adding
//		  a data member for the layer's owner, DchDetector.
//------------------------------------------------------------------------

#ifndef DCHLAYER_HH
#define DCHLAYER_HH
//----------------
// BaBar header --
//----------------
#if defined( HP1022 ) && !defined( BABAR_HH )
#include "BaBar/BaBar.hh"
#endif // HP1022 && !BABAR_HH
//-------------
// C Headers --
//-------------

//---------------
// C++ Headers --
//---------------
#include <vector>

//----------------------
// Base Class Headers --
//----------------------
#include "DetectorModel/DetSurfaceElem.hh"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "BaBar/Constants.hh"
#include "DchGeom/DchCell.hh"
#include "DchGeom/DchSWire.hh"
#include "TrkBase/TrkEnums.hh"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
class DchHyperType;
class DchSuperLayer;
class DchDetector;
class Trajectory;
#include <iosfwd>
#include <memory>

//		---------------------
// 		-- Class Interface --
//		---------------------

class DchLayer : public DetSurfaceElem {

  //--------------------
  // Instance Members --
  //--------------------

public:

  // Constructors
  DchLayer(DchHyperType* itsLayerType, const char* name, int id,
      const DetSurface& surf, int lay, const std::vector<DchSWire*>& wires,
      DchCell& thecell, const DchDetector& theDet, bool deb = false);

  // Copy Constructor
  //    DchLayer( const DchLayer& );

  // Destructor
  virtual
  ~DchLayer(void);

  // Operators
  virtual bool
  operator==(const DchLayer& compare) const;

  // Accessors (const)
  //   const DchHyperType* layerType() const
  //     { return (const DchHyperType*) detectorType(); }
  bool
  exist(void) const
  {
    return _exist;
  }
  bool
  debug(void) const
  {
    return _debug;
  }
  int
  layNum(void) const
  {
    return _layer;
  }
  int
  nWires(void) const
  {
    return _wires.size();
  }
  int
  superLayer(void) const
  {
    return 1 + ((layNum() - 1) / 4);
  }
  int
  subLayer(void) const
  {
    return layNum() - 4 * (superLayer() - 1);
  }
  double
  cellHeight(void) const
  {
    return _cellHeight;
  }
  double
  stDip(void) const
  {
    return _stdip;
  }
  double
  zEnd(void) const
  {
    return _zend;
  }
  double
  rMid(void) const
  {
    return _rmid;
  }
  double
  rEnd(void) const
  {
    return _rend;
  }
  double
  rIn(void) const
  {
    return rMid() - 0.5 * cellHeight();
  }
  double
  rOut(void) const
  {
    return rMid() + 0.5 * cellHeight();
  }
  double
  stereo(void) const
  {
    return _stereo;
  }
  double
  dPhiz(void) const
  {
    return _delphi;
  }
  double
  zLength(void) const
  {
    return getWire(0)->zLength();
  }
  double
  sag(void) const
  {
    return getWire(0)->getSag();
  }
  //  this function returns the phi offset in the middle of the chamber
  double
  phiOffset(void) const
  {
    return phiWire(0);
  }
  //  this function returns the cell offset w.r.t rear end-plate
  double
  phiEPOffset(void) const
  {
    return _phiOffset;
  }
  double
  dPhizDC(double z) const
  {
    return phiWireDC(0, z);
  }
  double
  radiusDC(double z) const
  {
    return getWire(0)->radiusDC(z);
  }
  DchSWire *
  getWire(int wire) const
  {
    return _wires[wire];
  }
  int
  view(void) const
  {
    return _view;
  }
  TrkEnums::TrkViewInfo
  whatView() const
  {
    return _view == 0 ? TrkEnums::xyView : TrkEnums::bothView;
  }
  double
  phiWireDC(int cell, double z) const;
  double
  phiWire(int cell) const;
  double
  xWire(int cell) const;
  double
  yWire(int cell) const;
  double
  dPhi(void) const
  {
    return Constants::twoPi / nWires();
  }
  // return the width of the cell at mid chamber corrected by the 
  // stereo angle
  double
  cellWidth(void) const
  {
    return dPhi() * rMid() * cos(stereo());
  }
  // return the width of the cell at z along the chamber
  double
  cellWidth(double z) const
  {
    return dPhi() * radiusDC(z) * cos(stereo());
  }

  //  const DchSuperLayer SLayer(void) const {return *_sLay;}

  // keep same name for backward compatibility
  //   TrkLineTraj makeHitTrajInGlobalCoords(int wire, double z) const;
  const Trajectory*
  makeHitTrajInGlobalCoords(int wire, double z = 0) const;

  //    DchLayer&       operator= ( const DchLayer& );
  //    virtual int operator==( const DchLayer& ) const;
  //            int operator!=( const DchLayer& ) const;

  // Selectors (const)
  const DchCell&
  getCell(int cell) const
  {
    return *_cells[cell];
  }
  int
  cellNum(const HepPoint& global) const;

  // returns list of cells in this layer traversed by the trajectory
  // the int returned by the three functions below refers to the number
  // of cells passed (in the case of whichClosestCell it is just 0 or 1).

  int
  whichCells(const Trajectory* traj, double range[2], std::vector<int>& cells,
      bool rec = false) const;

  // same as above but return also vector of paths

  int
  whichCells(const Trajectory* traj, double range[2], std::vector<int>& cells,
      std::vector<double>& paths, bool rec = false) const;

  // a somewhat faster function if you just want to know what is almost
  // certainly the closest channel
  int
  whichClosestCell(const Trajectory* traj, double range[2], int& closestCell,
      bool rec = false) const;

  // a function that returns the # of the closest cell to the given (phi,z) 
  int
  whichClosestCell(double phi, double z, int& closestCell) const;

  // a function which returns the list of cells in this layer within a
  // specified azimuthal distance ("road width") of the closest cell
  int
  whichCells(const Trajectory* traj, double range[2], double halfRoad,
      std::vector<int>& cells, bool rec = false) const;

  // returns a point on the layer surface at given phi and z in global
  // reference frame
  //   HepPoint surfacePoint( const double phi, const double z ) const;

  // Modifiers
  void
  setDebug(bool deb)
  {
    _debug = deb;
  }

protected:

  //   void setSlayer(const DchSuperLayer* sl) const { _sLay = sl; }
  void
  setCellHeight(double d)
  {
    _cellHeight = d;
  }
  // Helper functions

private:

  // Friends
  friend class DchSuperLayer;
  friend class DchDetector;

private:

  bool
  checkIntersect(int cellnum, const Trajectory* traj, double range[2],
      double path) const;
  // 		       std::vector<int>& cells ) const;
  // Data members
  const DchDetector& _dchDet; // Reference to layer's owner (for whichCells)

  bool _debug;
  bool _exist; // this layer exist ?
  // whichCells
  int _layer; // layer number
  int _nwires; // number of wires
  //  const DchSuperLayer* _sLay;   // pointer to SuperLayer
  std::vector<DchSWire *> _wires; // array of pointers to hits on wires
  // DchDriftCalib dCalib;         // time-to-distance object;   MOVED to DchEnv
  DchCell _theCell; // base cell for the layer
  std::vector<DchCell*> _cells; // cells in layer

  double _stdip; // change in radius from mid to end - nominal
  // value
  double _rend; // radius on rear end-plate
  double _rmid; // radius in the center of the chamber
  double _zend; // z position of end-plate (in the local
  // reference system
  double _stereo; // stereo angle
  double _delphi; // twist angle between mid and rear chamber
  double _phiOffset; // offset of first cell
  double _cellHeight; // cell height
  int _view;

  std::auto_ptr<DchHyperType> _type; // surface type

  DchCell*
  buildCell(int cellNum);

  DchLayer&
  operator=(const DchLayer&);
};

std::ostream&
operator <<(std::ostream& o, const DchLayer&);

#endif // DCHLAYER_HH
