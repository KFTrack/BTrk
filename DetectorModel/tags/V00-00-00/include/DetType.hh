// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetType.hh,v 1.25 2004/08/06 05:58:33 bartoldu Exp $
//
//  Description:
//  Base class to describe detector types.  These are generic
//  and are described only in a local (arbitrarily defined) coordinate
//  system.  Physical objects in space are defined by the DetElem
//  class (which includes a pointer to a DetType).
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 8/27/96
// ------------------------------------------------------------------------------
#ifndef DETECTORTYPE_HH
#define DETECTORTYPE_HH
//----------------
// BaBar header --
//----------------
#if defined( HP1022 ) && !defined( BABAR_HH )
#include "BaBar/BaBar.hh"
#include "ErrLogger/ErrLog.hh"
#endif // HP1022 && !BABAR_HH
//
//  Include files
//
#include "DetectorModel/DetTypeCoord.hh"
#include <string>
#include <iostream>
#include <math.h>
#include "BaBar/BaBar.hh"
#include "ErrLogger/ErrLog.hh"
#include <vector>

class DetMaterial;
//
//  Define the class
//
class DetType {
public:
//
//  Constructors; only simple ones, as this is an abstract base class
//
  DetType();
  DetType(const char*,int);
//  Destructor
  virtual ~DetType();
//
// Access functions; these are mostly virtual
//
  void print(std::ostream& os) const;
  int typeNumber() const { return _itype; } // type number
  const std::string& typeName() const {return _tname; } // type name
//
//  Determine whether a given (local) coordinate lies 'within' the object.
//
  virtual bool physicalMaterial(const TypeCoord*) const = 0; // physical size
//
//  Describe the material at a given coordinate
//
  virtual const DetMaterial& material(const TypeCoord*) const = 0;
//
//  Describe the 'outline' of this object
//
  virtual const std::vector< TypeCoord* >* outline() const;
  virtual const std::vector< TypeCoord* >* pointStore() const;
//
//  Functions needed for RogueWave
//
  virtual bool operator == ( const DetType& other ) const;

protected:
  virtual std::vector< TypeCoord* >* myOutline();
  virtual std::vector< TypeCoord* >* myPointStore();
private:
  int _itype; // unique identifier for each type (defined during creation)
  std::string _tname; // helpful name
  // This contains the information to draw elements, with null pointers
  // to define end of line segments. This one is for drawings.
  std::vector< TypeCoord* >* _outline;
  // This list owns the points which the previous list only uses
  std::vector< TypeCoord* >* _pointStore;
// prohibit
  DetType& operator = (const DetType&);
  DetType(const DetType&);
};
#endif
