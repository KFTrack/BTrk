//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchBetheBlochErr.cc 88 2010-01-14 12:32:57Z stroili $
//
// Description:
//	Class DchBetheBlochErr
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
#include "DchCalib/DchBetheBlochErr.hh"
#include "PDT/Pdt.hh"
using std::endl;
using std::ostream;

DchBetheBlochErr::DchBetheBlochErr( const std::vector<double>& parVector, 
                                    const char* functionName ) 
  : _par( parVector ), _functionName(functionName)
{;}

DchBetheBlochErr::~DchBetheBlochErr( ) 
{}

void
DchBetheBlochErr::print( ostream& o ) const
{
  o <<" Dch Bethe-Bloch error parametrization:\t"<<_functionName <<endl;
  for ( unsigned i=0; i<_par.size(); ++i ) o << "\t" << _par[i];
  o << endl;
}
