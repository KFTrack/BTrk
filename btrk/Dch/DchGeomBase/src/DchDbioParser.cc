//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchDbioParser.cc 92 2010-01-14 12:38:30Z stroili $
//
// Description:
//	Class DchDbioParser
//      Do not use this for DchDbioParserd class (foo<T>).  use DchDbioParserDchDbioParser.hh
//      instead.
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Roberto Stroili
//
// Copyright Information:
//	Copyright (C) 1997	<Institution>
//
// Revision History:
//	Converted to STL in 2002 by Roberto Stroili
//	20040409  M. Kelsey -- add missing EOF check to seek() loop
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "DchGeomBase/DchDbioParser.hh"

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
#include <vector>
#include <string>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
//#include "CLHEP/String/Strings.h"
#include "ErrLogger/ErrLog.hh"
using std::streampos;

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------
DchDbioParser::DchDbioParser( const char* filename )
  : DchGFileParser(filename)
{
}

//--------------
// Destructor --
//--------------
DchDbioParser::~DchDbioParser()
{}

//-------------
// Methods   --
//-------------
bool
DchDbioParser::seek( const std::string& theString ) 
{
  std::string line("xxx");
  streampos mark = _inFile.tellg();
//   while ( line.length() > 0 ) {
  while ( !_inFile.eof() ) {
    mark = _inFile.tellg();
    do {
      getline(_inFile, line);
    } while (line == "" && !_inFile.eof());
    if ( line.find(theString) != std::string::npos ) {
      //  found the line, step back
      _inFile.seekg(mark);
      return true;
    }
    line = "xxx";
  }
  ErrMsg(debugging)
    << " DchDbioParser::seek(" << theString << ") failed" << endmsg;
  return false;
}    

bool
DchDbioParser::seekTemplate( const std::string& theString ) 
{
  std::string line("xxx");
  streampos mark;
  while ( !_inFile.eof() ) {
    mark = _inFile.tellg();
    do {
      getline(_inFile, line);
    } while (line == "" && !_inFile.eof());
    if ( line.find(theString) != std::string::npos ) {
      //  found the line, step back
      _inFile.seekg(mark);
      return true;
    } else {
      line = "xxx";
      mark = _inFile.tellg();
    }
  } 
  ErrMsg(debugging)
    << " DchDbioParser::seekTemplate(" << theString << ") failed" << endmsg;
  return false;
}    

bool
DchDbioParser::locate( const std::string& templ, const std::string& name ) 
{
  ErrMsg(debugging) << " DchDbioParser looking for Template "
		    << templ << " name " << name << endmsg;

  if ( seek(templ) && seek(name) && check(templ, name) ) {
    return true;
  }
  //  if not found rewind the stream and search from beginning
  ErrMsg(debugging) << " rewind input stream" << endmsg;
  //    in reality close and reopen file
  fileRewind();
  if ( seek(templ) && seek(name) && check(templ,name) ) {
    return true;
  }
  ErrMsg(debugging) << " Template " << templ << " name " << name
		    << " not found" << endmsg;
  return false;
}    

bool
DchDbioParser::check( const std::string& templ, const std::string& name ) 
{
  //  mark current position
  streampos mark = _inFile.tellg();
  
  std::string str1, str2, str3;
  _inFile >> str1;
  _inFile >> str2;
  _inFile >> str3;
  //  step back to original position
  _inFile.seekg(mark);
  //  check current line
  if ( str1 != "make" || str2 != templ || str3 != name  ) {
    //  rewind stream ... bad entry
    ErrMsg(debugging) << " Bad data line. It does not begin with MAKE " 
		      << "or the second token (" << str2 << ") is not equal "
		      << "to " << templ << " " << name << endmsg;
    return false;
  }
  return true;
}    

bool
DchDbioParser::getDchVol( const std::string& volName, 
			  std::vector<double>& par,  
			  std::vector<double>& position )
{
  std::string templ("DchVol");

  if ( locate(templ, volName) ) {
    streampos mark = _inFile.tellg();
    par.clear();
    position.clear();
    
    std::string str1,str2,tag,name,shape,medium,mother,rotm,only;
    int npar, nin;
    float f1, f2, f3, f4, x, y, z;

    _inFile >> str1 >> str2 >> tag >> name >> shape >> medium >> npar >> f1 
	    >> f2 >> f3 >> f4 >> mother >> nin >> x >> y >> z >> rotm >> only;
    _inFile >> str1;
    par.push_back(f1); par.push_back(f2); par.push_back(f3); 
    par.push_back(f4); position.push_back(x); position.push_back(y); 
    position.push_back(z); 

    _inFile.seekg(mark);
    return true;
  } else {
    return false;
  }
}

bool
DchDbioParser::getDchSWire( int lay, std::vector<double>& par )  
{
  assert( lay > 0 );
  std::string layName;
  if ( lay <10 ) {
    //    layName = HepString("layer0")+HepString(lay);
  } else {
    // layName = HepString("layer")+HepString(lay);
  }

  if ( locate(std::string("DchSwir"), layName) ) {
    par.clear();

    std::string str1,str2,layer;
    int layerNum, nCells;
    float rEnd, phiEnd, twist, stereo, dphi, sag, diam, volt;

    _inFile >> str1 >> str2 >> layer >> layerNum >> nCells >> rEnd >> phiEnd
	    >> twist >> stereo >> dphi >> sag >> diam >> volt;

    par.push_back(nCells); par.push_back(rEnd); par.push_back(phiEnd); 
    par.push_back(twist); par.push_back(stereo); par.push_back(dphi); 
    par.push_back(sag); par.push_back(diam); par.push_back(volt); 
    return true;
  } else {
    return false;
  }
}

bool
DchDbioParser::getDchFWire( const std::string& type, int lay, 
			    std::vector<double>& par )  
{
  assert( lay > 0 );
  std::string layName;
  if ( lay <10 ) {
    //layName = HepString("layer0")+HepString(lay);
  } else {
    //  layName = HepString("layer")+HepString(lay);
  }

  if ( locate(type, layName) ) {
    par.clear();
    std::string str1,str2,layer;
    
    if ( type != "DchLfw" ) {
      int layerNum, nwirs;
      float rEnd, phiEnd, twist, stereo, dphi, sag, diam, volt;

      _inFile >> str1 >> str2 >> layer >> layerNum >> nwirs >> rEnd >> phiEnd
	      >> twist >> stereo >> dphi >> sag >> diam >> volt;

      par.push_back(nwirs); par.push_back(rEnd); par.push_back(phiEnd); 
      par.push_back(twist); par.push_back(stereo); par.push_back(dphi); 
      par.push_back(sag); par.push_back(diam); par.push_back(volt); 
    } else {
      int layerNum;
      float rEnd1, rEnd2, phiEnd, twist, stereo1, stereo2, dphi, sag, diam, 
	volt;

      _inFile >> str1 >> str2 >> layer >> layerNum >> rEnd1 >> rEnd2 >> phiEnd
	      >> twist >> stereo1 >> stereo2 >> dphi >> sag >> diam >> volt;

      par.push_back(rEnd1); par.push_back(rEnd2); par.push_back(phiEnd); 
      par.push_back(twist); par.push_back(stereo1); par.push_back(stereo2); 
      par.push_back(dphi); par.push_back(sag); par.push_back(diam); 
      par.push_back(volt); 
    }
    return true;
  } else {
    return false;
  }
}

