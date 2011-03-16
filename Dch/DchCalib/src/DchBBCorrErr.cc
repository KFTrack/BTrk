//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchBBCorrErr.cc 88 2010-01-14 12:32:57Z stroili $
//
// Description:
//	Class DchBBCorrErr
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
#include "BaBar/BaBar.hh"
#include "DchCalib/DchBBCorrErr.hh"
#include "PDT/Pdt.hh"
using std::endl;
using std::ostream;

DchBBCorrErr::DchBBCorrErr( const std::vector<double>& parVector, 
                                    const char* functionName ) 
  : _par( parVector ), _functionName(functionName)
{;}

DchBBCorrErr::~DchBBCorrErr( ) 
{}

void
DchBBCorrErr::print( ostream& o ) const
{
  o <<" Dch Bethe-Bloch error correction parametrization:\t"<<_functionName <<endl;
  for ( unsigned i=0; i<_par.size(); ++i ) o << "\t" << _par[i];
  o << endl;
}
