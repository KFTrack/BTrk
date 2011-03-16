//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchQSatCorr.cc 88 2010-01-14 12:32:57Z stroili $
//
// Description:
//	Class DchQSatCorr
//	Dch Transient class for charge saturation corrections
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
#include "BaBar/BaBar.hh"
#include "DchCalib/DchQSatCorr.hh"
#include "ErrLogger/ErrLog.hh"
using std::ostream;

DchQSatCorr::DchQSatCorr( const std::vector<double>& parVector,
			  const char* functionName ) 
  : _par( parVector ), _functionName(functionName)
{;}

DchQSatCorr::~DchQSatCorr( ) 
{}

void
DchQSatCorr::print( ostream& o ) const
{
  o <<" Dch correction for Q saturation:\t"<<_functionName <<endmsg;
  for ( unsigned i=0; i<_par.size(); i++ ) o << "\t" << _par[i];
  o << endmsg;
}
