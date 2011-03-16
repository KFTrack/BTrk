//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchGeomFileReader.cc 92 2010-01-14 12:38:30Z stroili $
//
// Description:
//	Class DchGeomFileReader
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	<R. Stroili>		<originator>
//	
//
// Copyright Information:
//	Copyright (C) 19997	<INFN>
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "DchGeomBase/DchGeomFileReader.hh"

//-------------
// C Headers --
//-------------
extern "C" {
}

//---------------
// C++ Headers --
//---------------
#include <assert.h>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "DchGeomBase/DchGFileParser.hh"
#include "DchGeomBase/DchDbioParser.hh"
#include "ErrLogger/ErrLog.hh"

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//-----------------
// Constructor   --
//-----------------

DchGeomFileReader::DchGeomFileReader()
{
}

//-------------
// Methods   --
//-------------
DchGFileParser*
DchGeomFileReader::fileParser( const std::string& file ) const
{

//    if ( file.find(".dat") != std::string::npos ) {
//      return 0;
//    } 

  if ( file.find(".db") != std::string::npos ) {
    return new DchDbioParser(file.c_str());
  } 
  
  ErrMsg(fatal) << "No valid filename was given:\t" << file << endmsg;
  return 0;
}




