//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchBBCorr.hh 88 2010-01-14 12:32:57Z stroili $
//
// Description:
//	Class DchBBCorr
//	Dch Transient class containing DchPID Bethe-Bloch parameters
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	F. Galeazzi		Originator
//
// Copyright Information:
//	Copyright (C) 1999	INFN & University of Padova
//
//------------------------------------------------------------------------
#ifndef DCHBBCORR_HH
#define DCHBBCORR_HH

//---------------
// C++ Headers --
//---------------
#include <vector>
#include <string>
#include <iostream>

#include "PDT/PdtPid.hh"


class DchBBCorr 
{
public:
  // Constructors
  DchBBCorr( const std::vector<double>& parVector,
                 const char* functionName );
  // Destructor
  virtual ~DchBBCorr( );

  const char* channelName( ) const { return _functionName.c_str(); }
  unsigned numberOfParameters( ) const { return _par.size(); }
  double fieldValue(unsigned i) const { return _par[i]; }
  virtual double ionization(double momentum, PdtPid::PidType hypo) const = 0;
  void print( std::ostream& o) const;
private:
  // Data members
  std::vector<double> _par;
  std::string _functionName;
  DchBBCorr();

};

#endif // DCHBBCORR_HH
