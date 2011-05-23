//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchBetheBloch.cc 132 2010-09-21 11:15:12Z stroili $
//
// Description:
//	Class DchBetheBloch
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
// Revision History:
//	20070328  Move ::getDchMomentum(TrkFit*) here from DchPidInfo
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "DchCalib/DchBetheBloch.hh"

#include <iostream>
using std::endl;
using std::ostream;


DchBetheBloch::DchBetheBloch( const std::vector<double>& parVector,
                              const char* functionName )
  : _par( parVector ), _functionName(functionName)
{;}

DchBetheBloch::~DchBetheBloch()
{}

void
DchBetheBloch::print( ostream& o ) const
{
  o <<" Dch PID Bethe-Bloch:\t"<<_functionName << endl;
  for ( unsigned i=0; i<_par.size(); ++i ) o << "\t" << _par[i];
  o << endl;
}
