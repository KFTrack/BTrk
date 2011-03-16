//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchSuperLayer.cc 123 2010-04-29 14:41:45Z stroili $
//
// Description:
//	Class DchSuperLayer
//      Do not use this for DchSuperLayerd class (foo<T>).  use DchSuperLayerDchSuperLayer.hh
//      instead.
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	R. Stroili		originator
//	
//
// Copyright Information:
//	Copyright (C) 1997	INFN - Pd
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "DchGeom/DchSuperLayer.hh"

//-------------
// C Headers --
//-------------
extern "C" {
}

//---------------
// C++ Headers --
//---------------
#include <assert.h>
using std::endl;
using std::ostream;

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

static const char rscid[] =
    "$Id: DchSuperLayer.cc 123 2010-04-29 14:41:45Z stroili $";
static const int _layInSuper = 4;

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------
DchSuperLayer::DchSuperLayer(const char* detname, int number) :
  DetSet(detname, number)
{
  _radius = _delphi = _delphiinv = 0;
  layers[0] = layers[1] = layers[2] = layers[3] = 0;
  _view = 0;
  _next = _prev = _nextInView = _prevInView = 0;
  _exist = false;
  _slayer = number - 1000;
}

//--------------
// Destructor --
//--------------
DchSuperLayer::~DchSuperLayer()
{
  //  delete [] layers;
}

//-------------
// Modifiers --
//-------------
void
DchSuperLayer::addLayer(int index, const DchLayer* lay)
{
  //---------------------------------------------------------- 
  // here |index| is the index of array of pointers to layers
  // belonging to the superlayer, so this ramges from 0 to 3
  //----------------------------------------------------------
  // check on index number
  assert ( index>=0 && index <4);
  // check that it was not already set
  assert ( layer(index) == 0 );
  // chack that layer stays in this superlayer
  assert ( (int)((lay->layNum()-1)/_layInSuper+1) == slayNum() );

  //   lay->setSlayer(this);
  layers[index] = lay;
}

void
DchSuperLayer::updateInfo(const DchSuperLayer* prev, const DchSuperLayer* next)
{
  //
  // function to set the data-members of this class
  //
  _exist = true;
  _radius = 0.5 * (firstLayer()->rEnd() + lastLayer()->rEnd());
  _view = firstLayer()->view();
  _delphi = firstLayer()->dPhiz();
  _delphiinv = 0.0;
  if (_delphi != 0.) _delphiinv = 1. / _delphi;
  // now the pointers
  _next = next;
  _prev = prev;
}

void
DchSuperLayer::print(ostream& o) const
{
  o << " SuperLayer #: " << slayNum() << "\n" << " rEnd: " << rEnd()
      << " rad0: " << rad0() << "\n" << " view: " << whichView() << " stDip: "
      << stDip() << "\n" << " delphi: " << delPhi() << " zEnd: " << zEnd()
      << "\n" << layers[0]->layNum() << " " << layers[1]->layNum() << " "
      << layers[2]->layNum() << " " << layers[3]->layNum() << "\n" << endl;
}

ostream&
operator<<(ostream& o, DchSuperLayer& sl)
{
  sl.print(o);
  return o;
}
