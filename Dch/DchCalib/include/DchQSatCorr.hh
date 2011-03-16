//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchQSatCorr.hh 88 2010-01-14 12:32:57Z stroili $
//
// Description:
//	Class DchQSatCorr
//	Dch Transient class for charge saturation correction
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	R. Stroili		Originator
//
// Copyright Information:
//	Copyright (C) 1999	INFN & University of Padova
//
//------------------------------------------------------------------------
#ifndef DCHQSATCORR_HH
#define DCHQSATCORR_HH

#include <iostream>
#include <vector>
#include <string>

class DchQSatCorr 
{
public:
  DchQSatCorr( const std::vector<double>& parVector,
	       const char* functionName );
  virtual ~DchQSatCorr( );
  virtual double correction( const std::vector<double>& arg ) const = 0;
  virtual DchQSatCorr* clone() const =0;
  
  unsigned numberOfParameters( ) const { return _par.size(); }
  const char* channelName( ) const { return _functionName.c_str(); }
  double fieldValue(unsigned i) const { return _par[i]; }
  void print( std::ostream& o) const;

protected:
  // Data members
  std::vector<double> _par;
  std::string _functionName;
};

#endif // DCHQSATCORR_HH
