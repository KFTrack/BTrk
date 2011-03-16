//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchGFileParser.cc 92 2010-01-14 12:38:30Z stroili $
//
// Description:
//	Class DchGFileParser
//      Do not use this for DchGFileParserd class (foo<T>).  use DchGFileParserDchGFileParser.hh 
//      instead.
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//
// Copyright Information:
//	Copyright (C) 1997	<Institution>
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "DchGeomBase/DchGFileParser.hh"

//-------------
// C Headers --
//-------------
extern "C" {
}

//---------------
// C++ Headers --
//---------------
#include <assert.h>
#include <iomanip>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "CLHEP/String/Strings.h"
#include "ErrLogger/ErrLog.hh"

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------
DchGFileParser::DchGFileParser( const char* filename )
  : _inFile(filename), _filename(filename)
{
  assert( _inFile.good() );
  _begin = _inFile.tellg();
  ErrMsg(trace) << "DchGFileParser opened " << filename
		<< ": _begin = " << _begin << endmsg;
}

//--------------
// Destructor --
//--------------
DchGFileParser::~DchGFileParser()
{}

