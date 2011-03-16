//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchPlateDefl.hh 91 2010-01-14 12:37:23Z stroili $
//
// Description:
//	Class DchPlateDefl.
//      This class describes the geometry variations due to the
//      Dch end-plate deflections. It assumes a cylindrical symmetry
//      for these corrections
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	R. Stroili		originator
//	
//
// Copyright Information:
//	Copyright (C) 1998	INFN & Padova University
//
//------------------------------------------------------------------------

#ifndef DCHPLATEDEFL_HH
#define DCHPLATEDEFL_HH

#include <vector>
class DchWirePar;

class DchPlateDefl {
public:

  // Constructors
  DchPlateDefl( );
  DchPlateDefl( const char* ); // construct from file
  DchPlateDefl( const DchPlateDefl& );       // Copy Constructor
  DchPlateDefl( std::vector<DchWirePar>& );

  // Destructor
  virtual ~DchPlateDefl( );

  // Selectors (const)
  const DchWirePar& wireCorr( int ) const;
  const std::vector<DchWirePar>& corrList() const { return _deflCorr; }

  DchPlateDefl& operator=(const DchPlateDefl&);

private:
  // Data members
  std::vector<DchWirePar> _deflCorr;   // correction for plate deflections
                              // the correction is global for one layer
};

#endif // DCHPLATEDEFL_HH
