//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchCellPlaneType.cc 123 2010-04-29 14:41:45Z stroili $
//
// Description:
//	Class DchCellPlaneType
//      Do not use this for DchCellPlaneTyped class (foo<T>).  use DchCellPlaneTypeDchCellPlaneType.hh
//      instead.
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	<Author1>		<originator/contributor etc.>
//	<Author2>		<originator/contributor etc.>
//
// Copyright Information:
//	Copyright (C) 1996	<Institution>
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "DchGeom/DchCellPlaneType.hh"

//-------------
// C Headers --
//-------------
extern "C" {
}

//---------------
// C++ Headers --
//---------------
#include <iostream>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "AbsEnv/AbsEnv.hh"
#include "GenEnv/GenEnv.hh"

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------
#define PLANE_TOLERANCE 0.01

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------

DchCellPlaneType::DchCellPlaneType(const char * name, int id, std::vector<
    TwoDCoord>& corners) :
  DetSurfaceType(name, id, corners), _thick(0.)
{
  _material = gblEnv->getGen()->findDetMaterial("GasWire");

  // find surface boundaries

  _incorner[0] = 99999.;
  _incorner[1] = 99999.;
  _outcorner[0] = -99999.;
  _outcorner[1] = -99999.;
  for (int i = 0; i < corners.size(); i++) {
    if (_incorner[0] > (corners[i])[0]) _incorner[0] = (corners[i])[0];
    if (_incorner[1] > (corners[i])[1]) _incorner[1] = (corners[i])[1];
    if (_outcorner[0] < (corners[i])[0]) _outcorner[0] = (corners[i])[0];
    if (_outcorner[1] < (corners[i])[1]) _outcorner[1] = (corners[i])[1];
  }
}

//--------------
// Destructor --
//--------------
DchCellPlaneType::~DchCellPlaneType()
{
  ;
}

//-------------
// Operators --
//-------------
bool
DchCellPlaneType::physicalMaterial(const TypeCoord* thispoint) const
{
  //
  //  First check if the point is in the rectangle
  //

  double uvalue = (*thispoint)[0];
  double vvalue = (*thispoint)[1];

  //   cout <<  _incorner[0]<< " - " << uvalue << " - " <<_outcorner[0]<<endl;
  //   cout <<  _incorner[1]<< " - " << vvalue << " - " <<_outcorner[1]<<endl;

  if (((uvalue - _incorner[0]) > -PLANE_TOLERANCE && (uvalue - _outcorner[0])
      < PLANE_TOLERANCE) && ((vvalue - _incorner[1]) > -PLANE_TOLERANCE
      && (vvalue - _outcorner[1]) < PLANE_TOLERANCE)) {
    return true;
  } else {
    return false;
  }
}

double
DchCellPlaneType::thickness(const TwoDCoord*) const
{
  return _thick;
}

const DetMaterial&
DchCellPlaneType::material(const TypeCoord*) const
{
  return *_material;
}
