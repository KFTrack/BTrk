//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchGDch.cc 125 2010-05-19 07:31:02Z stroili $
//
// Description:
//	Class DchGDch
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Roberto Stroili    -    INFN
//
// History (add to end):
//      Stroili   May 11, 1997  - creation
//
// Copyright Information:
//	Copyright (C) 1997		INFN
//	
// Revision History:
//	20040408  M. Kelsey -- Fix bad code!  Don't put actions in assert();
//		  Check for  pointers before use; add trace messages during
//		  file processing.
//------------------------------------------------------------------------

//----------------
// BaBar Header --
//----------------
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "DchGeomBase/DchGDch.hh"

//---------------
// C++ Headers --
//---------------
#include <math.h>
#include <iostream>
#include <assert.h>
#include <vector>

//-------------------------------
// collaborating class headers --
//-------------------------------
#include "BaBar/Constants.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "DchGeomBase/DchGFileParser.hh"
#include "DchGeomBase/DchGeomFileReader.hh"
#include "ErrLogger/ErrLog.hh"  
using std::ostream;

//-------------------
// BaBar C Headers --
//-------------------
// #include "gnbase/dch_dbio.h"

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

//              -----------------------------------
//              -- External Function Declaration --
//              -----------------------------------

// extern "C" void bbgeomDbioInit(); // geometry constants initializations

//----------------
// Constructor  --
//----------------

DchGDch::DchGDch( const std::string& fileName ) 
{
  ErrMsg(trace) << "DchGDch: using" << fileName << endmsg;
  DchGeomFileReader filereader;

  initWires();

  DchGFileParser* parser = filereader.fileParser(fileName);
  assert(0 != parser);

  bool parseOk = true;

  double z_forward, z_rear;

  std::vector<double> par;
  std::vector<double> mainPosition;

  //  get main volume position
  parseOk = parser->getDchVol("DCHA", par, mainPosition);
  assert(parseOk);

  //  get gas volume
  std::vector<double> position;
  parseOk = parser->getDchVol("DCHG", par, position);
  assert(parseOk);

  _GasCyl = DchSimpleCyl(par[0],par[1],par[2]*2.);
  for (int i=0; i<3; i++ ) {
    _GasPos[i] = position[i]+mainPosition[i];
  }
  z_forward = GasLength()/2.+_GasPos[2];
  z_rear = -GasLength()/2.+_GasPos[2];

  //  get inner cylinder
  parseOk = parser->getDchVol("DCHI", par, position);
  assert(parseOk);
  _InnerCyl = DchSimpleCyl(par[0],par[1],par[2]*2.);
  for (int i=0; i<3; i++ ) {
    _ICPos[i] = position[i]+mainPosition[i];
  }

  //  get outer cylinder
  parseOk = parser->getDchVol("DCHO", par, position);
  assert(parseOk);
  _OuterCyl = DchSimpleCyl(par[0],par[1],par[2]*2.);
  for (int i=0; i<3; i++ ) {
    _OCPos[i] = position[i]+mainPosition[i];
  }

  std::vector<double> plateOffset;

  //  FORWARD END PLATE
  //      first some calculations: offset
  parseOk = parser->getDchVol("DCFP", par, plateOffset);
  assert(parseOk);

  parseOk = parser->getDchVol("DCF1", par, position);
  assert(parseOk);

  _ForwCyl = DchSimpleCyl(par[0],par[1],par[2]*2.);
  for (int i=0; i<3; i++ ) {
    _FEPPos[i] = position[i]+mainPosition[i]+plateOffset[i];
  }

  //  REAR END PLATE
  //      first some calculations: offset
  parseOk = parser->getDchVol("DCRP", par, plateOffset);
  assert(parseOk);

  parseOk = parser->getDchVol("DCR1", par, position);
  assert(parseOk);

  _RearCyl = DchSimpleCyl(par[0],par[1],par[2]*2.);
  for (int i=0; i<3; i++ ) {
    _REPPos[i] = position[i]+mainPosition[i]+plateOffset[i];
  }


  std::vector< double > ipar;
  std::vector< double > lpar;
  std::vector< double > opar;
  std::string iwir("DchIfw");
  std::string lwir("DchLfw");
  std::string owir("DchOfw");

  int i(0);
  while ( parser->getDchSWire(i+1,par) &&
	  parser->getDchFWire(lwir,i+1,lpar) &&
	  parser->getDchFWire(iwir,i+1,ipar) &&
	  parser->getDchFWire(owir,i+1,opar) ) {
    
    _sWires[i].nCells   = (int)par[DchGFileParser::cells];
    _sWires[i].rEnd     = par[DchGFileParser::Srad];
    _sWires[i].rForw    = par[DchGFileParser::Srad];
    _sWires[i].phiEnd   = par[DchGFileParser::Sphi];
    _sWires[i].phiForw  = par[DchGFileParser::Sphi]
      +par[DchGFileParser::Stwist];
    _sWires[i].twist    = par[DchGFileParser::Stwist];
    _sWires[i].stereo   = par[DchGFileParser::Sstereo];
    _sWires[i].deltaPhi = Constants::twoPi/(int)par[DchGFileParser::cells];
    _sWires[i].zFw      = z_forward;
    _sWires[i].zR       = z_rear;
    _sWires[i].sag      = par[DchGFileParser::Ssag];
    _sWires[i].diameter = par[DchGFileParser::Sdiam];
    _sWires[i].volts    = par[DchGFileParser::Svolt];
    
    _cellStruc[i].nofFWires = 6;
    int wOrder = 0;
    // ordering goes from 1 to 7
    _cellStruc[i].FWires[wOrder].order = wOrder+1;
    _cellStruc[i].FWires[wOrder].rEnd = lpar[DchGFileParser::LFrad1];
    _cellStruc[i].FWires[wOrder].rForw = lpar[DchGFileParser::LFrad1];
    _cellStruc[i].FWires[wOrder].phiEnd = lpar[DchGFileParser::LFphi];
    _cellStruc[i].FWires[wOrder].phiForw = lpar[DchGFileParser::LFphi]+
      lpar[DchGFileParser::LFtwist];
    _cellStruc[i].FWires[wOrder].twist = lpar[DchGFileParser::LFtwist];
    _cellStruc[i].FWires[wOrder].stereo = 
      lpar[DchGFileParser::LFstereo1];
    _cellStruc[i].FWires[wOrder].deltaPhi = Constants::twoPi/
      par[DchGFileParser::cells];
    _cellStruc[i].FWires[wOrder].zFw = z_forward;
    _cellStruc[i].FWires[wOrder].zR = z_rear;
    _cellStruc[i].FWires[wOrder].sag = lpar[DchGFileParser::LFsag];
    _cellStruc[i].FWires[wOrder].diameter = lpar[DchGFileParser::LFdiam];
    _cellStruc[i].FWires[wOrder].volts = lpar[DchGFileParser::LFvolt];
    
    wOrder++;
    
    _cellStruc[i].FWires[wOrder].order = wOrder+1;
    _cellStruc[i].FWires[wOrder].rEnd = lpar[DchGFileParser::LFrad2];
    _cellStruc[i].FWires[wOrder].rForw = lpar[DchGFileParser::LFrad2];
    _cellStruc[i].FWires[wOrder].phiEnd = lpar[DchGFileParser::LFphi];
    _cellStruc[i].FWires[wOrder].phiForw = lpar[DchGFileParser::LFphi]+
      lpar[DchGFileParser::LFtwist];
    _cellStruc[i].FWires[wOrder].twist = lpar[DchGFileParser::LFtwist];
    _cellStruc[i].FWires[wOrder].stereo = 
      lpar[DchGFileParser::LFstereo2];
    //     _cellStruc[i].FWires[wOrder].deltaPhi = DchLfw[i].dphi;
    _cellStruc[i].FWires[wOrder].deltaPhi = Constants::twoPi/
      par[DchGFileParser::cells];
    _cellStruc[i].FWires[wOrder].zFw = z_forward;
    _cellStruc[i].FWires[wOrder].zR = z_rear;
    _cellStruc[i].FWires[wOrder].sag = lpar[DchGFileParser::LFsag];
    _cellStruc[i].FWires[wOrder].diameter = lpar[DchGFileParser::LFdiam];
    _cellStruc[i].FWires[wOrder].volts = lpar[DchGFileParser::LFvolt];
    
    wOrder++;
    
    _cellStruc[i].FWires[wOrder].order = wOrder+1;
    _cellStruc[i].FWires[wOrder].rEnd = opar[DchGFileParser::Fradius];
    _cellStruc[i].FWires[wOrder].rForw = opar[DchGFileParser::Fradius];
    _cellStruc[i].FWires[wOrder].phiEnd = opar[DchGFileParser::Fphi];
    _cellStruc[i].FWires[wOrder].phiForw = opar[DchGFileParser::Fphi]+
      opar[DchGFileParser::Ftwist];
    _cellStruc[i].FWires[wOrder].twist = opar[DchGFileParser::Ftwist];
    _cellStruc[i].FWires[wOrder].stereo = opar[DchGFileParser::Fstereo];
    //     _cellStruc[i].FWires[wOrder].deltaPhi = DchOfw[i].dphi;
    _cellStruc[i].FWires[wOrder].deltaPhi = Constants::twoPi/
      par[DchGFileParser::cells];
    _cellStruc[i].FWires[wOrder].zFw = z_forward;
    _cellStruc[i].FWires[wOrder].zR = z_rear;
    _cellStruc[i].FWires[wOrder].sag = opar[DchGFileParser::Fsag];
    _cellStruc[i].FWires[wOrder].diameter = opar[DchGFileParser::Fdiam];
    _cellStruc[i].FWires[wOrder].volts = opar[DchGFileParser::Fvolt];
    
    wOrder++;
    
    if ( opar[DchGFileParser::nwires] == 2) {
      _cellStruc[i].nofFWires = 7;
      _cellStruc[i].FWires[wOrder].order = wOrder+1;
      _cellStruc[i].FWires[wOrder].rEnd = opar[DchGFileParser::Fradius];
      _cellStruc[i].FWires[wOrder].rForw = opar[DchGFileParser::Fradius];
      _cellStruc[i].FWires[wOrder].deltaPhi = opar[DchGFileParser::Fdphi];
      _cellStruc[i].FWires[wOrder].phiEnd = opar[DchGFileParser::Fphi] +
	_cellStruc[i].FWires[wOrder].deltaPhi;
      _cellStruc[i].FWires[wOrder].phiForw = opar[DchGFileParser::Fphi] +
	_cellStruc[i].FWires[wOrder].deltaPhi + 
	opar[DchGFileParser::Ftwist];
      _cellStruc[i].FWires[wOrder].twist = opar[DchGFileParser::Ftwist];
      _cellStruc[i].FWires[wOrder].stereo = 
	opar[DchGFileParser::Fstereo];
      //       _cellStruc[i].FWires[wOrder].deltaPhi = DchOfw[i].dphi;
      _cellStruc[i].FWires[wOrder].zFw = z_forward;
      _cellStruc[i].FWires[wOrder].zR = z_rear;
      _cellStruc[i].FWires[wOrder].sag = opar[DchGFileParser::Fsag];
      _cellStruc[i].FWires[wOrder].diameter = opar[DchGFileParser::Fdiam];
      _cellStruc[i].FWires[wOrder].volts = opar[DchGFileParser::Fvolt];
      
      wOrder++;
      
    }
    
    _cellStruc[i].FWires[wOrder].order = wOrder+1;
    _cellStruc[i].FWires[wOrder].rEnd = lpar[DchGFileParser::LFrad2];
    _cellStruc[i].FWires[wOrder].rForw = lpar[DchGFileParser::LFrad2];
    _cellStruc[i].FWires[wOrder].deltaPhi = Constants::twoPi/
      par[DchGFileParser::cells];
    _cellStruc[i].FWires[wOrder].phiEnd = lpar[DchGFileParser::LFphi] +
      _cellStruc[i].FWires[wOrder].deltaPhi;
    _cellStruc[i].FWires[wOrder].phiForw = lpar[DchGFileParser::LFphi] +
      _cellStruc[i].FWires[wOrder].deltaPhi + 
      lpar[DchGFileParser::LFtwist];
    _cellStruc[i].FWires[wOrder].twist = lpar[DchGFileParser::LFtwist];
    _cellStruc[i].FWires[wOrder].stereo = 
      lpar[DchGFileParser::LFstereo2];
    //     _cellStruc[i].FWires[wOrder].deltaPhi = DchLfw[i].dphi;
    _cellStruc[i].FWires[wOrder].zFw = z_forward;
    _cellStruc[i].FWires[wOrder].zR = z_rear;
    _cellStruc[i].FWires[wOrder].sag = lpar[DchGFileParser::LFsag];
    _cellStruc[i].FWires[wOrder].diameter = lpar[DchGFileParser::LFdiam];
    _cellStruc[i].FWires[wOrder].volts = lpar[DchGFileParser::LFvolt];
    
    wOrder++;
    
    _cellStruc[i].FWires[wOrder].order = wOrder+1;
    _cellStruc[i].FWires[wOrder].rEnd = lpar[DchGFileParser::LFrad1];
    _cellStruc[i].FWires[wOrder].rForw = lpar[DchGFileParser::LFrad1];
    _cellStruc[i].FWires[wOrder].deltaPhi = Constants::twoPi/
      par[DchGFileParser::cells];
    _cellStruc[i].FWires[wOrder].phiEnd = lpar[DchGFileParser::LFphi] +
      _cellStruc[i].FWires[wOrder].deltaPhi;
    _cellStruc[i].FWires[wOrder].phiForw = lpar[DchGFileParser::LFphi] +
      _cellStruc[i].FWires[wOrder].deltaPhi + 
      lpar[DchGFileParser::LFtwist];
    _cellStruc[i].FWires[wOrder].twist = lpar[DchGFileParser::LFtwist];
    _cellStruc[i].FWires[wOrder].stereo = 
      lpar[DchGFileParser::LFstereo1];
    //     _cellStruc[i].FWires[wOrder].deltaPhi = DchLfw[i].dphi;
    _cellStruc[i].FWires[wOrder].zFw = z_forward;
    _cellStruc[i].FWires[wOrder].zR = z_rear;
    _cellStruc[i].FWires[wOrder].sag = lpar[DchGFileParser::LFsag];
    _cellStruc[i].FWires[wOrder].diameter = lpar[DchGFileParser::LFdiam];
    _cellStruc[i].FWires[wOrder].volts = lpar[DchGFileParser::LFvolt];
    
    wOrder++;
    
    _cellStruc[i].FWires[wOrder].order = wOrder+1;
    _cellStruc[i].FWires[wOrder].rEnd = ipar[DchGFileParser::Fradius];
    _cellStruc[i].FWires[wOrder].rForw = ipar[DchGFileParser::Fradius];
    _cellStruc[i].FWires[wOrder].phiEnd = ipar[DchGFileParser::Fphi] + 
      ipar[DchGFileParser::Fdphi];
    _cellStruc[i].FWires[wOrder].phiForw = ipar[DchGFileParser::Fphi] +
      ipar[DchGFileParser::Fdphi] + ipar[DchGFileParser::Ftwist];
    _cellStruc[i].FWires[wOrder].twist = ipar[DchGFileParser::Ftwist];
    _cellStruc[i].FWires[wOrder].stereo = ipar[DchGFileParser::Fstereo];
    //     _cellStruc[i].FWires[wOrder].deltaPhi = ipar[DchGFileParser::Fdphi];
    _cellStruc[i].FWires[wOrder].deltaPhi = Constants::twoPi/
      par[DchGFileParser::cells];
    _cellStruc[i].FWires[wOrder].zFw = z_forward;
    _cellStruc[i].FWires[wOrder].zR = z_rear;
    _cellStruc[i].FWires[wOrder].sag = ipar[DchGFileParser::Fsag];
    _cellStruc[i].FWires[wOrder].diameter = ipar[DchGFileParser::Fdiam];
    _cellStruc[i].FWires[wOrder].volts = ipar[DchGFileParser::Fvolt];
    
    wOrder++;
    
    if ( ipar[DchGFileParser::nwires] == 2) {
      _cellStruc[i].nofFWires = 7;
      _cellStruc[i].FWires[wOrder].order = wOrder+1;
      _cellStruc[i].FWires[wOrder].rEnd = ipar[DchGFileParser::Fradius];
      _cellStruc[i].FWires[wOrder].rForw = ipar[DchGFileParser::Fradius];
      _cellStruc[i].FWires[wOrder].phiEnd = ipar[DchGFileParser::Fphi];
      _cellStruc[i].FWires[wOrder].phiForw = ipar[DchGFileParser::Fphi] +
	ipar[DchGFileParser::Ftwist];
      _cellStruc[i].FWires[wOrder].twist = ipar[DchGFileParser::Ftwist];
      _cellStruc[i].FWires[wOrder].stereo = ipar[DchGFileParser::Fstereo];
      //       _cellStruc[i].FWires[wOrder].deltaPhi = ipar[DchGFileParser::Fdphi];
      _cellStruc[i].FWires[wOrder].deltaPhi = Constants::twoPi/
	par[DchGFileParser::cells];
      _cellStruc[i].FWires[wOrder].zFw = z_forward;
      _cellStruc[i].FWires[wOrder].zR = z_rear;
      _cellStruc[i].FWires[wOrder].sag = ipar[DchGFileParser::Fsag];
      _cellStruc[i].FWires[wOrder].diameter = ipar[DchGFileParser::Fdiam];
      _cellStruc[i].FWires[wOrder].volts = ipar[DchGFileParser::Fvolt];
    }
    i++;
  }
  _nLayers = i;


}

// copy constructor
DchGDch::DchGDch( const DchGDch& other ) 
  : _version(999)
{
  // volumes for Kalman filter
  // gas volume 
  _GasCyl = other._GasCyl;
  _InnerCyl = other._InnerCyl;
  _OuterCyl = other._OuterCyl;
  _RearCyl = other._RearCyl;
  _ForwCyl = other._ForwCyl;

  // elements for tracking
  _nLayers = other._nLayers;
  _zoffset = other._zoffset;
  _zlen = other._zlen;
  _cellHeight = other._cellHeight;

  for ( int i=0; i<3; i++ ) {
    _GasPos[i] = other._GasPos[i];
    _ICPos[i] = other._ICPos[i];
    _OCPos[i] = other._OCPos[i];
    _REPPos[i] = other._REPPos[i];
    _FEPPos[i] = other._FEPPos[i];
  }

  
  // copy wires
  for ( int layer=0; layer<_nLayers; layer++ ) {
    _rOffset[layer] = other._rOffset[layer];

    // first sense wires ...
    _sWires[layer].nCells   = other._sWires[layer].nCells;
    _sWires[layer].rEnd     = other._sWires[layer].rEnd;
    _sWires[layer].rForw    = other._sWires[layer].rForw; 
    _sWires[layer].phiEnd   = other._sWires[layer].phiEnd;
    _sWires[layer].phiForw  = other._sWires[layer].phiForw;
    _sWires[layer].twist    = other._sWires[layer].twist;
    _sWires[layer].stereo   = other._sWires[layer].stereo;
    _sWires[layer].deltaPhi = other._sWires[layer].deltaPhi;
    _sWires[layer].zR       = other._sWires[layer].zR;
    _sWires[layer].zFw      = other._sWires[layer].zFw;
    _sWires[layer].sag      = other._sWires[layer].sag;
    _sWires[layer].diameter = other._sWires[layer].diameter;
    _sWires[layer].volts    = other._sWires[layer].volts;

    // ... then field wires to define cell structure
    int nofwires = other._cellStruc[layer].nofFWires;
    _cellStruc[layer].nofFWires = nofwires;
    for ( int fwire=0; fwire<nofwires; fwire++ ) {
      _cellStruc[layer].FWires[fwire].order = 
	other._cellStruc[layer].FWires[fwire].order;
      _cellStruc[layer].FWires[fwire].rEnd = 
	other._cellStruc[layer].FWires[fwire].rEnd;
      _cellStruc[layer].FWires[fwire].rForw = 
	other._cellStruc[layer].FWires[fwire].rForw;
      _cellStruc[layer].FWires[fwire].phiEnd = 
	other._cellStruc[layer].FWires[fwire].phiEnd;
      _cellStruc[layer].FWires[fwire].phiForw = 
	other._cellStruc[layer].FWires[fwire].phiForw;
      _cellStruc[layer].FWires[fwire].twist = 
	other._cellStruc[layer].FWires[fwire].twist;
      _cellStruc[layer].FWires[fwire].stereo = 
	other._cellStruc[layer].FWires[fwire].stereo;
      _cellStruc[layer].FWires[fwire].deltaPhi = 
	other._cellStruc[layer].FWires[fwire].deltaPhi;
      _cellStruc[layer].FWires[fwire].zR = 
	other._cellStruc[layer].FWires[fwire].zR;
      _cellStruc[layer].FWires[fwire].zFw = 
	other._cellStruc[layer].FWires[fwire].zFw;
      _cellStruc[layer].FWires[fwire].sag = 
	other._cellStruc[layer].FWires[fwire].sag;
      _cellStruc[layer].FWires[fwire].diameter = 
	other._cellStruc[layer].FWires[fwire].diameter;
      _cellStruc[layer].FWires[fwire].volts = 
	other._cellStruc[layer].FWires[fwire].volts;
    }
  }
}

DchGDch::DchGDch( int lay ) 
  : _version(999), _nLayers(lay)
{
  DchSimpleCyl cyl;
  _name         = "DCH Detector";
  _GasCyl = cyl;
  _InnerCyl = cyl;
  _OuterCyl = cyl;
  _RearCyl = cyl;
  _ForwCyl = cyl;

  for ( int j=0; j<3; j++ ) {
    _GasPos[j] = 0.;
    _ICPos[j] = 0.;
    _OCPos[j] = 0.;
    _REPPos[j] = 0.;
    _FEPPos[j] = 0.;
  }

  for ( int i=0; i<_nLayers; i++ ) {
    _rOffset[i] = 0.;
  }
  _zoffset = 0.;
  _zlen = 0.;
  _cellHeight = 0.;
}


//--------------
// Destructor --
//--------------

DchGDch::~DchGDch() {
}
  
//              -----------------------------------------
//              -- Public  Function Member Definitions --
//              -----------------------------------------

void 
DchGDch::initWires(void) {
  for (int layer=0; layer < 40; ++layer  ) {
    _sWires[layer].nCells   = 0;
    _sWires[layer].rEnd     = 0;
    _sWires[layer].rForw    = 0;
    _sWires[layer].phiEnd   = 0;
    _sWires[layer].phiForw  = 0;
    _sWires[layer].twist    = 0;
    _sWires[layer].stereo   = 0;
    _sWires[layer].deltaPhi = 0;
    _sWires[layer].zFw      = 0;
    _sWires[layer].zR       = 0;
    _sWires[layer].sag      = 0;
    _sWires[layer].diameter = 0;
    _sWires[layer].volts    = 0;
  }
}

void 
DchGDch::print(ostream& o) const {
  o << _name << endmsg;
}

void DchGDch::printAll(ostream& o) const {
  o << _name << " version # " <<_version<<"\n"
    << " Gas volume:\n"
    << "     inner radius: "<<GasInRad()<<"\n"
    << "     outer radius: "<<GasOutRad()<<"\n"
    << "     length:       "<<GasLength()<<"\n"
    << "     position:     "<<GasPosX()<<" "<<GasPosY()<<" "<<GasPosZ()
    <<endmsg;
  if ( _version >= 10 ) {
    o << " Inner cylinder:\n"
      << "     inner radius: "<<ICInRad()<<"\n"
      << "     outer radius: "<<ICOutRad()<<"\n"
      << "     length:       "<<ICLength()<<"\n"
      << "     position:     "<<ICPosX()<<" "<<ICPosY()<<" "<<ICPosZ()
      <<endmsg;
    o << " Outer cylinder:\n"
      << "     inner radius: "<<OCInRad()<<"\n"
      << "     outer radius: "<<OCOutRad()<<"\n"
      << "     length:       "<<OCLength()<<"\n"
      << "     position:     "<<OCPosX()<<" "<<OCPosY()<<" "<<OCPosZ()
      <<endmsg;
    o << " Rear End Plate:\n"
      << "     inner radius: "<<REPInRad()<<"\n"
      << "     outer radius: "<<REPOutRad()<<"\n"
      << "     length:       "<<REPLength()<<"\n"
      << "     position:     "<<REPPosX()<<" "<<REPPosY()<<" "<<REPPosZ()
      <<endmsg;
    o << " Forward End Plate:\n"
      << "     inner radius: "<<FEPInRad()<<"\n"
      << "     outer radius: "<<FEPOutRad()<<"\n"
      << "     length:       "<<FEPLength()<<"\n"
      << "     position:     "<<FEPPosX()<<" "<<FEPPosY()<<" "<<FEPPosZ()
      <<endmsg;
  }
  o << "print sense wire radius, phi at rear end-plate, twist, phi at "
    << "forward end-plate, stereo and sag" 
    << endmsg;
  for ( int i=0; i<_nLayers; ++i ) {
    o << i << "\t" << _sWires[i].rEnd << "\t" <<  _sWires[i].phiEnd << "\t"
      << _sWires[i].twist << "\t" << _sWires[i].phiForw << "\t" 
      << _sWires[i].stereo << "\t" << _sWires[i].sag << endmsg;
  }
  o << "print cell structure sense wire radius, phi at rear end-plate, twist, "
    << "phi at forward end-plate, stereo and sag" 
    << endmsg;
  for ( int j=0; j<_nLayers; ++j ) {
    for ( int k=0; k<_cellStruc[j].nofFWires; ++k ) {
      o << j << " " << k <<"\t" << _cellStruc[j].FWires[k].rEnd << "\t" 
	<< _cellStruc[j].FWires[k].phiEnd << "\t"
	<< _cellStruc[j].FWires[k].twist << "\t" 
	<< _cellStruc[j].FWires[k].phiForw << "\t"
	<< _cellStruc[j].FWires[k].stereo << "\t" 
	<< _cellStruc[j].FWires[k].sag << endmsg;
    }
  }
}


ostream&  operator << (ostream& o, const DchGDch& a) {
  a.printAll(o); 
  return o;
}

