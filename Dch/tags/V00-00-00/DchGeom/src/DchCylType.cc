//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchCylType.cc 123 2010-04-29 14:41:45Z stroili $
//
// Description:
//	Class DchCylType
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	R. Stroili		originator
//	
//
// Copyright Information:
//	Copyright (C) 1997	INFN-Pd
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "DchGeom/DchCylType.hh"

//-------------
// C Headers --
//-------------

//---------------
// C++ Headers --
//---------------
#include <math.h>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "ErrLogger/ErrLog.hh"
#include "DchGeomBase/DchSimpleCyl.hh"

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

static const char rscid[] = "$Id: DchCylType.cc 123 2010-04-29 14:41:45Z stroili $";

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------
const int nperimeter=107;
//----------------
// Constructors --
//----------------
DchCylType::DchCylType(const char* name, double radii[2], double length,
		       const DetMaterial* themat, int idnum) :
  DetSurfaceType(name,idnum), _thick(radii[1]-radii[0]), _length(length), _radius(radii[0]),
      _perimeter(0), _thematerial(themat)
{
  std::vector< TypeCoord* >* PointStore = myPointStore();
  std::vector< TypeCoord* >* Outline = myOutline();

  _perimeter = new TwoDCoord[nperimeter];
  int ipoint;
  double phistep = 2*M_PI/26;
  double phi;
  int jpoint=0;
  for(ipoint=0,phi=0.0;ipoint<27;ipoint++,phi+= phistep){
    _perimeter[jpoint][0] = phi;
    //_perimeter[jpoint][1] = _radius;
    _perimeter[jpoint][1] = -_length/2 ;
    jpoint++;
  }
  for(ipoint=0,phi=0.0;ipoint<52;ipoint += 4,phi+= 2*phistep){
    _perimeter[jpoint][0] = phi;
    //_perimeter[jpoint][1] = _radius;
    _perimeter[jpoint][1] = -_length/2 ;
    _perimeter[jpoint+1][0] = phi;
    //_perimeter[jpoint+1][1] = _radius+_thick;
    _perimeter[jpoint+1][1] = +_length/2 ;
    _perimeter[jpoint+2][0] = phi+phistep;
    //_perimeter[jpoint+2][1] = _radius+_thick;
    _perimeter[jpoint+2][1] = +_length/2;
    _perimeter[jpoint+3][0] = phi+phistep;
    //_perimeter[jpoint+3][1] = _radius;
    _perimeter[jpoint+3][1] = -_length/2 ;
    jpoint += 4;
  }
  _perimeter[jpoint][0] = 0.0;
  _perimeter[jpoint][1] = -_length/2;
  jpoint++;
  for(ipoint=0,phi=0.0;ipoint<27;ipoint++,phi+= phistep){
    _perimeter[jpoint][0] = phi;
    //    _perimeter[jpoint][1] = _radius+_thick;
    _perimeter[jpoint][1] = _length/2;
    
    jpoint++;
  }

 
  TwoDCoord* point ;
  for (int i = 0; i < nperimeter; i++){
    point = new TwoDCoord(_perimeter[i][0],_perimeter[i][1]);
    PointStore->push_back((TypeCoord*)point);
    Outline->push_back((TypeCoord*)point);
  }
  point = new TwoDCoord(_perimeter[0][0],_perimeter[0][1]);
  Outline->push_back((TypeCoord*)point);
  PointStore->push_back((TypeCoord*)point);
  
}

DchCylType::DchCylType(const char* name, const DchSimpleCyl* cyl,
		       const DetMaterial* themat, int idnum) :
  DetSurfaceType(name,idnum), _perimeter(0), _thematerial(themat)
{
  _thick = cyl->getOuterRadius()-cyl->getInnerRadius();
  _radius = cyl->getInnerRadius();
  _length = cyl->getLength();

  std::vector< TypeCoord* >* PointStore = myPointStore();
  std::vector< TypeCoord* >* Outline = myOutline();

  
  _perimeter = new TwoDCoord[nperimeter];
  int ipoint;
  double phistep = 2*M_PI/26;
  double phi;
  int jpoint=0;
  for(ipoint=0,phi=0.0;ipoint<27;ipoint++,phi+= phistep){
    _perimeter[jpoint][0] = phi;
    //_perimeter[jpoint][1] = _radius;
    _perimeter[jpoint][1] = -_length/2 ;
    jpoint++;
  }
  for(ipoint=0,phi=0.0;ipoint<52;ipoint += 4,phi+= 2*phistep){
    _perimeter[jpoint][0] = phi;
    //_perimeter[jpoint][1] = _radius;
    _perimeter[jpoint][1] = -_length/2 ;
    _perimeter[jpoint+1][0] = phi;
    //_perimeter[jpoint+1][1] = _radius+_thick;
    _perimeter[jpoint+1][1] = +_length/2 ;
    _perimeter[jpoint+2][0] = phi+phistep;
    //_perimeter[jpoint+2][1] = _radius+_thick;
    _perimeter[jpoint+2][1] = +_length/2;
    _perimeter[jpoint+3][0] = phi+phistep;
    //_perimeter[jpoint+3][1] = _radius;
    _perimeter[jpoint+3][1] = -_length/2 ;
    jpoint += 4;
  }
  _perimeter[jpoint][0] = 0.0;
  _perimeter[jpoint][1] = -_length/2;
  jpoint++;
  for(ipoint=0,phi=0.0;ipoint<27;ipoint++,phi+= phistep){
    _perimeter[jpoint][0] = phi;
    //    _perimeter[jpoint][1] = _radius+_thick;
    _perimeter[jpoint][1] = _length/2;
    
    jpoint++;
  }

 
  TwoDCoord* point ;
  for (int i = 0; i < nperimeter; i++){
    point = new TwoDCoord(_perimeter[i][0],_perimeter[i][1]);
    PointStore->push_back((TypeCoord*)point);
    Outline->push_back((TypeCoord*)point);
  }
  // Put first point in as last point
  Outline->push_back((*PointStore)[0]);
  PointStore->push_back((TypeCoord*)point);
  
}

DchCylType::DchCylType( const DchCylType& copy ) 
{
  ErrMsg(fatal)<<" constructor not allowed "<<endmsg;
}

//--------------
// Destructor --
//--------------
DchCylType::~DchCylType()
{
  if (_perimeter != 0) delete[] _perimeter;
}

//-------------
// Operators --
//-------------
bool 
DchCylType::physicalMaterial(const TypeCoord* thispoint) const 
{
//
//  Check if the point is on the strut
//
  double zval = (*thispoint)[1];
  return zval >= -_length/2 && zval <= _length/2; 
}

const DetMaterial& 
DchCylType::material(const TypeCoord*) const 
{ 
  return *_thematerial;
}

double 
DchCylType::thickness(const TwoDCoord*) const 
{ 
  return _thick; 
}
