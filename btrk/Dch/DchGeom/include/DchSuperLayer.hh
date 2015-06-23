//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchSuperLayer.hh 123 2010-04-29 14:41:45Z stroili $
//
// Description:
//	Class DchSuperLayer.
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

#ifndef DCHSUPERLAYER_HH
#define DCHSUPERLAYER_HH

//-------------
// C Headers --
//-------------

//---------------
// C++ Headers --
//---------------
#include <assert.h>

//----------------------
// Base Class Headers --
//----------------------
#include "DetectorModel/DetSet.hh"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "DchGeom/DchLayer.hh"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
class DchLayer;
//		---------------------
// 		-- Class Interface --
//		---------------------
static const int _layInSuper = 1;

class DchSuperLayer : public DetSet {

//--------------------
// Declarations     --
//--------------------

    // Typedefs, consts, and enums

//--------------------
// Instance Members --
//--------------------

public:

    // Constructors
  DchSuperLayer(const char* detname, int number);
  DchSuperLayer() {}

    // Copy Constructor
//    DchSuperLayer( const DchSuperLayer& );

    // Destructor
  virtual ~DchSuperLayer( );

    // Operators
  int whichView(void) const                                  {return _view;}
  bool exist(void) const                                    {return _exist;}
  double rEnd(void) const                                  {return _radius;}
  double rad0(void) const { return 0.5 * (firstLayer()->rMid() + 
		  lastLayer()->rMid()); }
  double zEnd(void) const { return 0.5 * (firstLayer()->zEnd() + 
		  lastLayer()->zEnd()); }
  double stDip(void) const { return 0.5 * (firstLayer()->stDip() + 
		  lastLayer()->stDip()); }
  double delPhi(void) const                                {return _delphi;}
  double delPhiinv(void) const                          {return _delphiinv;}
  //---------------------------------------------------------- 
  // here |index| is the index of array of pointers to layers
  // belonging to the superlayer, so this ramges from 0 to 3
  //----------------------------------------------------------
  const DchLayer* layer(int i) const { assert ( i>=0 && i <_layInSuper ) ; 
                                       return layers[i]; }
  const DchLayer* firstLayer(void) const                 {return layers[0];}
  const DchLayer* lastLayer(void) const                  {return layers[_layInSuper-1];}
  const DchSuperLayer* next(void) const                      {return _next;}
  const DchSuperLayer* prev(void) const                      {return _prev;}
  const DchSuperLayer* nextInView(void) const          {return _nextInView;}
  const DchSuperLayer* prevInView(void) const          {return _prevInView;}
  void setNextInView(DchSuperLayer* sl)                  {_nextInView = sl;}
  void setPrevInView(DchSuperLayer* sl)                  {_prevInView = sl;}
  int slayNum(void) const                                  {return _slayer;}
  // to maintain bacward compatibility
  int index(void) const                                  {return _slayer-1;}
// One line printout
  void print(std::ostream& o=std::cout) const ;
    
//    DchSuperLayer&       operator= ( const DchSuperLayer& );
//    virtual int operator==( const DchSuperLayer& ) const;
//            int operator!=( const DchSuperLayer& ) const;

    // Selectors (const)

    // Modifiers

private:

  // Friends
  friend class DchDetector;
  
  void addLayer(int index, const DchLayer* lay);
  void updateInfo(const DchSuperLayer* prev, const DchSuperLayer* next);

  // Data members
  bool _exist;
  double _radius;  // mean rad.
  double _delphi; // diff in phi between z=0 and zend (=0 for axial
  double _delphiinv; 
  const DchLayer* layers[_layInSuper];
  int _view;  // +1, 0, -1 = U, axial, V
  int _slayer;  // superlayer number
  const DchSuperLayer* _next;
  const DchSuperLayer* _prev;
  const DchSuperLayer* _nextInView;
  const DchSuperLayer* _prevInView;

};

std::ostream&  operator << (std::ostream& o, DchSuperLayer&);

#endif // DCHSUPERLAYER_HH
