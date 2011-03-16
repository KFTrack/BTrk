//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchHyperType.cc 123 2010-04-29 14:41:45Z stroili $
//
// Description:
//	Class DchHyperType
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
#include "DchGeom/DchHyperType.hh"

//-------------
// C Headers --
//-------------
extern "C" {
}

//---------------
// C++ Headers --
//---------------

//-------------------------------
// Collaborating class Headers --
//-------------------------------
#include "BaBar/Constants.hh"

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

static const char rscid[] =
    "$Id: DchHyperType.cc 123 2010-04-29 14:41:45Z stroili $";

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------
const int nperimeter = 107;

//----------------
// Constructors --
//----------------
DchHyperType::DchHyperType(const char* name, const DetMaterial* material,
    int id, double twist, double length) :
  DetSurfaceType(name, id), _twist(twist), _length(length), _thematerial(
      material)
{
  std::vector<TypeCoord*>* PointStore = myPointStore();
  std::vector<TypeCoord*>* Outline = myOutline();

  _perimeter = new TwoDCoord[nperimeter];
  int ipoint;
  double phistep = 2 * M_PI / 26;
  double phi;
  int jpoint = 0;
  for (ipoint = 0, phi = 0.0; ipoint < 27; ipoint++, phi += phistep) {
    _perimeter[jpoint][0] = phi;
    _perimeter[jpoint][1] = -_length / 2;
    jpoint++;
  }
  for (ipoint = 0, phi = 0.0; ipoint < 52; ipoint += 4, phi += 2 * phistep) {
    _perimeter[jpoint][0] = phi;
    _perimeter[jpoint][1] = -_length / 2;

    _perimeter[jpoint + 1][0] = phi + _twist;
    _perimeter[jpoint + 1][1] = +_length / 2;

    _perimeter[jpoint + 2][0] = phi + phistep + _twist;
    _perimeter[jpoint + 2][1] = +_length / 2;

    _perimeter[jpoint + 3][0] = phi + phistep;
    _perimeter[jpoint + 3][1] = -_length / 2;
    jpoint += 4;
  }
  _perimeter[jpoint][0] = 0.0;
  _perimeter[jpoint][1] = -_length / 2;
  jpoint++;
  for (ipoint = 0, phi = 0.0; ipoint < 27; ipoint++, phi += phistep) {
    _perimeter[jpoint][0] = phi + _twist;
    _perimeter[jpoint][1] = _length / 2;

    jpoint++;
  }
  TwoDCoord* point;
  for (int i = 0; i < nperimeter; i++) {
    point = new TwoDCoord(_perimeter[i][0], _perimeter[i][1]);
    PointStore->push_back(point);
    Outline->push_back(point);
  }
  // Add first point as last point
  Outline->push_back((*PointStore)[0]);
}
//--------------
// Destructor --
//--------------
DchHyperType::~DchHyperType()
{
  delete[] _perimeter;
}

bool
DchHyperType::physicalMaterial(const TypeCoord* point) const
{
  double z = (*point)[1];
  if (z > -_length * .5 && z < _length * .5) {
    return true;
  } else {
    return false;
  }
}

double
DchHyperType::thickness(const TwoDCoord* point) const
{
  return Constants::epsilon;
}

double
DchHyperType::effectiveThickness(const TwoDCoord* point, double dircos1,
    double dircos2) const
{
  return 0;
}

const DetMaterial&
DchHyperType::material(const TypeCoord*) const
{
  return *_thematerial;
}

