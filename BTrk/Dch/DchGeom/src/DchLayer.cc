//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchLayer.cc 123 2010-04-29 14:41:45Z stroili $
//
// Description:
//	Class DchLayer
//      Do not use this for DchLayerd class (foo<T>).  use DchLayerDchLayer.hh
//      instead.
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	R. Stroili	        originator
//	
//
// Copyright Information:
//	Copyright (C) 1996	INFN & Padova University
//
// Revision History:
//	20020411 M. Kelsey -- Add sense wire to construction of DchCell
//	20020417 M. Kelsey -- Set DchSWire::_id in ::buildCell
//	20020619  M. Kelsey -- Break circular dependence w/DchEnv by adding
//		  a data member for the layer's owner, DchDetector.  Also,
//		  reorder ctor initializers to match .hh declarations
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "DchGeom/DchLayer.hh"

//-------------------------------
//  C++ headers
//-------------------------------
#include <vector>
#include <algorithm>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "BaBar/Constants.hh"
#include "BbrStdUtils/BbrCollectionUtils.hh"
#include "BbrGeom/Trajectory.hh"
#include "DchGeom/DchDetector.hh"
#include "DchGeom/DchHyperType.hh"
#include "DchGeom/DchHyperb.hh"
#include "DchGeomBase/DchCellAddr.hh"
#include "DetectorModel/DetIntersection.hh"
#include "DetectorModel/DetSurface.hh"
#include "DetectorModel/DetSurfaceElem.hh"
#include "ErrLogger/ErrLog.hh"
using std::ostream;

// NOTE:  DO NOT TERMINE #define STATEMENTS WITH SEMI-COLONS!
#define LOWER_PATH_RANGE 20.
#define UPPER_PATH_RANGE 250.
#define CELL_INTERVAL 5.

DchLayer::DchLayer(DchHyperType* lType, const char* name, int id,
    const DetSurface& surf, int lay, const std::vector<DchSWire*>& wires,
    DchCell& thecell, const DchDetector& theDet, bool deb) :
  DetSurfaceElem(lType, name, id, surf), _dchDet(theDet), _debug(deb), _exist(
      !wires.empty()), _layer(lay), _wires(wires), _theCell(thecell), _type(
      lType)
{
  // set pointer to layer for wires
  for (int wire = 0; wire < nWires(); ++wire) {
    getWire(wire)->setLayerPtr(this);
    _cells.push_back(buildCell(wire));
  }

  // get nominal values from wire # 0
  _rend = wires[0]->rEnd();
  _rmid = wires[0]->rMid();
  _stdip = _rend - _rmid;
  // in the local (DC) system frame
  _zend = wires[0]->zEndDC();
  _stereo = wires[0]->stereo();
  // half twist angle
  _delphi = wires[0]->dPhiz();
  // stereo view
  _view = _stereo == 0 ? 0 : (int) (_stereo / (fabs(_stereo)));
  _phiOffset = wires[0]->phiE();
}

DchLayer::~DchLayer(void)
{
  std::for_each(_cells.begin(), _cells.end(), babar::Collection::DeleteObject());
}

double
DchLayer::phiWireDC(int cell, double z) const
{
  if (cell >= 0 && cell < nWires()) {
    return getWire(cell)->phiDC(z);
  } else {
    ErrMsg(error) << "phiWireDC: wrong cell # " << cell << "\n"
        << "  number of cells in this layer is " << nWires() << endmsg;
    return 0.;
  }
}

double
DchLayer::phiWire(int cell) const
{
  // in the middle of the chamber
  if (cell >= 0 && cell < nWires()) {
    return getWire(cell)->phi();
  } else {
    ErrMsg(error) << "phiWire: wrong cell # " << cell << "\n"
        << "  number of cells in this layer is " << nWires() << endmsg;
    return 0.;
  }
}

double
DchLayer::xWire(int cell) const
{
  // in the middle of the chamber
  if (cell >= 0 && cell < nWires()) {
    return getWire(cell)->xMid();
  } else {
    ErrMsg(error) << "xWire: wrong cell # " << cell << "\n"
        << "  number of cells in this layer is " << nWires() << endmsg;
    return 0.;
  }
}

double
DchLayer::yWire(int cell) const
{
  // in the middle of the chamber 
  if (cell >= 0 && cell < nWires()) {
    return getWire(cell)->yMid();
  } else {
    ErrMsg(error) << "YWire: wrong cell # " << cell << "\n"
        << "  number of cells in this layer is " << nWires() << endmsg;
    return 0.;
  }
}

// TrkLineTraj 
const Trajectory*
DchLayer::makeHitTrajInGlobalCoords(int wire, double z) const
{
  if (wire >= 0 && wire < nWires()) {
    return getWire(wire)->getTraj();
  } else {
    ErrMsg(error) << "makeHitTrajInGlobalCoords: wrong cell # " << wire << "\n"
        << "  number of cells in this layer is " << nWires() << endmsg;
    return getWire(0)->getTraj();
  }
}

DchCell*
DchLayer::buildCell(int cellNum)
{
  double phi = cellNum * dPhi();

  if (cellNum < 0 || cellNum >= nWires()) {
    ErrMsg(error) << "getCell: wrong cell # " << cellNum << "\n"
        << "  number of cells in this layer is " << nWires() << endmsg;
    phi = 0; // default this cell to cell # 0
    cellNum = 0;
  }
  int cellID = DchCellAddr::cellIs(cellNum, _layer);

  getWire(cellNum)->_id = cellID; // Set sense wire ID number

  return new DchCell(_theCell, phi, cellID, getWire(cellNum));
}

int
DchLayer::whichCells(const Trajectory* traj, double range[2],
    std::vector<int>& cells, std::vector<double>& paths, bool recursion) const
{

  // first of all clear list of cells
  cells.clear();
  paths.clear();

  // some local variables
  int nCells = 0;
  double lastCellRange = range[0];
  double newrange[2];

  // set some control flags
  bool increment = true;
  bool decrement = true;

  int cellnum = 0;
  int flag = whichClosestCell(traj, range, cellnum);

  double path = 0.;

  if (flag) {

    // check if trajectory intersects found cell
    //     DchCell cell = getCell(cellnum);
    newrange[0] = range[0] - CELL_INTERVAL;
    newrange[1] = range[0] + CELL_INTERVAL;
    if (checkIntersect(cellnum, traj, newrange, path)) {
      if (newrange[0] > lastCellRange) lastCellRange = newrange[0];
      //     std::vector<HepPoint> vec;
      //     double path = cell.intersect( traj, newrange, vec );
      //     if ( path != 0 ) {
      //       if ( newrange[0] > lastCellRange ) lastCellRange =
      // 					   newrange[0];
      cells.push_back(cellnum);
      paths.push_back(path);
    } else {
      if (ErrLogging(debugging)) {
        ErrMsg(debugging) << "not a valid intersection, " << "cell/layer "
            << cellnum << "/" << layNum() << "  "
        // 			 <<vec.size()
            // 			 <<" intersection points found"
            << endmsg;
      }
    }
  } else if (!recursion) {
    //     cout<<"\t\ttry again"<<endl;
    const DchLayer* prevLayer = 0;
    const DchLayer* nextLayer = 0;
    if (layNum() > 1) {
      prevLayer = _dchDet.getDchLayer(layNum() - 1);
    }
    if (layNum() < 40) {
      nextLayer = _dchDet.getDchLayer(layNum() + 1);
    }

    if (prevLayer != 0 && prevLayer->whichClosestCell(traj, range, cellnum)) {
      HepPoint intPoint = traj->position(range[0]);
      whichClosestCell(intPoint.phi(), intPoint.z(), cellnum);
      // 	cout << "\tfound in previous layer"<<endl;
    } else if (nextLayer != 0 && nextLayer->whichClosestCell(traj, range,
        cellnum)) {
      HepPoint intPoint = traj->position(range[0]);
      whichClosestCell(intPoint.phi(), intPoint.z(), cellnum);
    } else {
      ErrMsg(debugging) << "no intersection found on nearby layers!" << endmsg;
      return 0;
    }

    if (checkIntersect(cellnum, traj, range, path)) {
      if (range[0] > lastCellRange) lastCellRange = range[0];
      //     std::vector<HepPoint> vec;
      //     double path = cell.intersect( traj, newrange, vec );
      //     if ( path != 0 ) {
      //       if ( newrange[0] > lastCellRange ) lastCellRange =
      // 					   newrange[0];
      cells.push_back(cellnum);
      paths.push_back(path);
    } else {
      if (ErrLogging(debugging)) {
        ErrMsg(debugging) << "not a valid intersection, " << "cell/layer "
            << cellnum << "/" << layNum() << "  "
        // 			 <<vec.size()
            // 			 <<" intersection points found"
            << endmsg;
      }
    }
    //     if ( ErrLogging(debugging) ) {
    //       ErrMsg(debugging)<<"trajectory does not intersect layer "<<layNum()
    // 		       <<" "<<flag<<" "<<range[0]<<" "<<range[1]<<endmsg;
    //     }
    //     return 0;


  }

  // look neighboring cells
  int cellindex = cellnum;
  // first decrement index
  while (decrement) {
    cellindex--;
    if (cellindex < 0) cellindex += nWires();

    DchCell nearcell = getCell(cellindex);
    double newrange[2];
    newrange[0] = range[0] - CELL_INTERVAL;
    newrange[1] = range[0] + CELL_INTERVAL;
    double path = nearcell.intersect(traj, newrange);
    if (path != 0) {
      if (newrange[0] > lastCellRange) lastCellRange = newrange[0];
      cells.push_back(cellindex);
      paths.push_back(path);
    } else {
      decrement = false;
    }
  }

  // then increment index
  cellindex = cellnum;
  while (increment) {
    cellindex++;
    if (cellindex >= nWires()) cellindex -= nWires();

    DchCell nearcell = getCell(cellindex);
    double newrange[2];
    newrange[0] = range[0] - CELL_INTERVAL;
    newrange[1] = range[0] + CELL_INTERVAL;
    double path = nearcell.intersect(traj, newrange);
    if (path != 0) {
      if (newrange[0] > lastCellRange) lastCellRange = newrange[0];
      cells.push_back(cellindex);
      paths.push_back(path);
    } else {
      increment = false;
    }
  }
  nCells = cells.size();

  //  // add some interval to avoid loops
  //  newrange[0] = lastCellRange + CELL_INTERVAL;
  //  newrange[1] = UPPER_PATH_RANGE;
  //
  //  if ( nCells > 0 && ( newrange[1] > newrange[0] ) ) {
  //    std::vector<int> moreCells;
  //    std::vector<double> morePaths;
  //    int more = whichCells(traj,newrange,moreCells,morePaths,true);
  //    for ( int index=0; index<more; index++ ) {
  //      cells.push_back( moreCells[index] );
  //      paths.push_back( morePaths[index] );
  //    }
  //  }
  return cells.size();
}

int
DchLayer::whichCells(const Trajectory* traj, double range[2],
    std::vector<int>& cells, bool recursion) const
{
  static std::vector<double> paths;

  //     int retval = whichCells(traj, range, cells, paths, recursion);
  return whichCells(traj, range, cells, paths, recursion);
}

int
DchLayer::whichClosestCell(const Trajectory* traj, double range[2],
    int& closestCell, bool recursion) const
{
  closestCell = 0;
  // intersect with layer surface
  DetIntersection inters;
  inters.pathrange[0] = range[0];
  inters.pathrange[1] = range[1];
  inters.pathlen = range[0];

  int flag = intersect(traj, inters);
  if (!flag) return 0; // Exit if no intersection

  range[0] = inters.pathrange[0];
  range[1] = inters.pathrange[1];

  HepPoint point = traj->position(inters.pathrange[0]);
  closestCell = cellNum(point);
  if (closestCell < 0 || closestCell >= nWires()) {
    if (ErrLogging(debugging)) {
      ErrMsg(debugging) << "trajectory not traversing Dch layer " << layNum()
          << " cell number " << closestCell << endmsg;
    }
    return 0;
  }

  return 1;
}

int
DchLayer::whichClosestCell(double phi, double z, int& closestCell) const
{
  //  WARNING: this function works in the global BaBar reference system,
  //  phi and z are in that frame, before using them we need to convert 
  //  them to the local reference frame (just shift z)
  double locZ = _dchDet.gotoLoc(z);
  //  now project phi on the rear end-plate
  double rPhi = phi - dPhiz() * (1. + 2. * locZ / zLength());
  //  and now it's easy to get the cell number ((rPhi-phiOffset())/dPhi())
  //  ... just remember that cells are staggered
  closestCell = (int) (0.5 + (rPhi - phiEPOffset()) / dPhi());
  if (closestCell < 0) {
    closestCell += nWires() - 1;
  }
  return 1;
}

int
DchLayer::whichCells(const Trajectory* traj, double range[2], double halfRoad,
    std::vector<int>& cells, bool recursion) const
{
  // Initialize list if non-recursive call
  if (!recursion) cells.clear();

  // Find cell nearest trajectory as seed
  int middleCell;
  if (0 == whichClosestCell(traj, range, middleCell, recursion)) {
    if (ErrLogging(debugging)) {
      ErrMsg(debugging) << " whichCells: trajectory doesn't intersect layer"
          << endmsg;
    }
    return 0;
  }

  cells.push_back(middleCell); // Intersection point

  // Get radius of layer at intersection (for stereo)
  HepPoint point = traj->position(range[0]);
  double rLyr = sqrt(point.x() * point.x() + point.y() * point.y());

  // Collect cells adjacent to trajectory out to road width
  int ncells = (int) (halfRoad / (dPhi() * rLyr));
  ncells = (ncells < nWires() / 2) ? ncells : nWires() / 2 - 1;
  int i, cLeft, cRight;
  for (i = 1; i <= ncells; i++) { // Expressions wrap cell# at 2-pi
    cLeft = (middleCell + i) % nWires();
    cRight = (middleCell - i + nWires()) % nWires();

    cells.push_back(cLeft);
    cells.push_back(cRight);
  }

  return 1 + 2 * ncells; // Number of cells added to list
}

int
DchLayer::cellNum(const HepPoint& point) const
{
  double projected = ((DchHyperb*) surface())->projectPhiToRear(point);

  if (projected < -1000.) {
    if (ErrLogging(debugging)) {
      ErrMsg(debugging) << " given point does not belong to layer surface "
          << endmsg;
    }
    return -1;
  }
  double angle = projected - _phiOffset;

  double d = angle / dPhi() + 0.5;

  if (d < 0.) d += nWires();
  if (d >= nWires()) d -= nWires();

  return (int) d;
}

bool
DchLayer::checkIntersect(int cellnum, const Trajectory* traj, double range[2],
    double path) const
// 			  std::vector<int>& cells ) const 
{
  bool retval = false;
  // check if trajectory intersects found cell
  DchCell cell = getCell(cellnum);
  std::vector<HepPoint> vec;
  path = cell.intersect(traj, range, vec);
  if (path != 0) {
    retval = true;
    //       if ( range[0] > lastCellRange ) lastCellRange =
    // 					newrange[0];
    //     cells.push_back(cellnum);
  } else {
    if (ErrLogging(debugging)) {
      ErrMsg(debugging) << "not a valid intersection, " << "cell/layer "
          << cellnum << "/" << layNum() << "  " << vec.size()
          << " intersection points found" << endmsg;
    }
  }
  return retval;
}
//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

bool
DchLayer::operator==(const DchLayer&) const
{
  bool answer = false;
  return answer;
}

ostream&
operator<<(ostream& o, const DchLayer& l)
{
  o << " Layer #: " << l.layNum() << " # of wires: " << l.nWires() << "\n"
      << " SuperLayer #: " << l.superLayer() << " view: " << l.view() << "\n"
      << " rEnd: " << l.rEnd() << " rMid: " << l.rMid() << "\n" << " stereo: "
      << l.stereo() << " Delta phi: " << l.dPhiz() << "\n" << " Offset: "
      << l.phiOffset() << " length: " << l.zLength() << endmsg;
  l.getWire(0)->print(o);
  return o;
}
