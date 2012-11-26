//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchGDch.hh 92 2010-01-14 12:38:30Z stroili $
//
// Description:
//	DrcGDch Class - 
//        Simple data class containing nominal parameters needed
//        to set-up Dch geometry. The structure is very similar to
//        the structure of the dch.db file.
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Roberto Stroili     -    INFN
//
// History (add to end):
//      Stroili   May 11, 1997  - creation
//
// Copyright Information:
//	Copyright (C) 1997		INFN
//	
//
//------------------------------------------------------------------------

#ifndef DCHGDCH_HH
#define DCHGDCH_HH

#include <iostream>
#include <fstream>
#include <assert.h>
#include <string>
#include <vector>

// collaborating class headers
#include "DchGeomBase/DchSimpleCyl.hh"
#include "DchGeomBase/DchPhiSegmCyl.hh"

#include <iosfwd>
#include "CLHEP/Vector/ThreeVector.h"
class DchGDchP;
class DchSympleCyl;

struct senseWire {    // sense wire structure
  int nCells;         // # of cells in this layer
  double rEnd;        // radius at rear end plate
  double rForw;       // radius at forward end plate
  double phiEnd;      // phi of sense wire #0 at rear end plate
  double phiForw;     // phi of sense wire #0 at forward end plate
  double twist;       // twist of layer
  double stereo;      // tangent of stero angle
  double deltaPhi;    // phi aperture of cell
  double zFw;         // forward zeta coordinate
  double zR;          // rear zeta coordinate
  double sag;         // sagitta
  double diameter;    // wire diameter
  double volts;       // wire voltage
};

struct fieldWire {    // field wire structure
  int order;          // numbering order of wire in cell
  double rEnd;        // radius at rear end plate
  double rForw;       // radius at forward end plate
  double phiEnd;      // phi of sense wire #0 at rear end plate
  double phiForw;     // phi of sense wire #0 at forward end plate
  double twist;       // twist of layer
  double stereo;      // tangent of stero angle
  double deltaPhi;    // phi aperture of cell
  double zFw;         // forward zeta coordinate
  double zR;          // rear zeta coordinate
  double sag;         // sagitta
  double diameter;    // wire diameter
  double volts;       // wire voltage
};

struct cellStruct{     // basic cell structure
  int nofFWires;      // number of field wires defining the cell (6 or 7)
  fieldWire FWires[7];  // pointer to the field wires
};

static const int nLayer = 40;

class DchGDch {

public:

// Constructor
//  DchGDch( const std::string& fileName ); // construct from ascii file db
  DchGDch( int lay );  // "dummy" constructor, used to build the geometry
                       // from the persistent geometry 
  DchGDch( const DchGDch& ); // copy construct

// Destructor
  virtual ~DchGDch( );

// accessor functions
  //  gas volume (Gas)
  double GasInRad(void) const {return _GasCyl.getInnerRadius();}
  double GasOutRad(void) const {return _GasCyl.getOuterRadius();}
  double GasLength(void) const {return _GasCyl.getLength();}
  const DchSimpleCyl& getGasCyl(void) const {return _GasCyl;}
  const double* GasPos(void) const {return _GasPos;};
  double GasPosX(void) const {return _GasPos[0];}
  double GasPosY(void) const {return _GasPos[1];}
  double GasPosZ(void) const {return _GasPos[2];}
  //  inner cylinder (IC)
  double ICInRad(void) const {return _InnerCyl.getInnerRadius();}
  double ICOutRad(void) const {return _InnerCyl.getOuterRadius();}
  double ICLength(void) const {return _InnerCyl.getLength();}
  const DchSimpleCyl& getInnerCyl(void) const {return _InnerCyl;}
  const double* ICPos(void) const {return _ICPos;}
  double ICPosX(void) const {return _ICPos[0];}
  double ICPosY(void) const {return _ICPos[1];}
  double ICPosZ(void) const {return _ICPos[2];}
  //  cylinder layer of inner guard wires (IWGrd)
  double IWGrdInRad(void) const {return _InnerWGrdCyl.getInnerRadius();}
  double IWGrdOutRad(void) const {return _InnerWGrdCyl.getOuterRadius();}
  double IWGrdLength(void) const {return _InnerWGrdCyl.getLength();}
  const DchSimpleCyl& getInnerWGrdCyl(void) const {return _InnerWGrdCyl;}
  const double* IWGrdPos(void) const {return _IWGrdPos;}
  double IWGrdPosX(void) const {return _IWGrdPos[0];}
  double IWGrdPosY(void) const {return _IWGrdPos[1];}
  double IWGrdPosZ(void) const {return _IWGrdPos[2];}
  //  outer cylinder (OC)
  double OCInRad(void) const {return _OuterCyl.getInnerRadius();}
  double OCOutRad(void) const {return _OuterCyl.getOuterRadius();}
  double OCLength(void) const {return _OuterCyl.getLength();}
  const DchSimpleCyl& getOuterCyl(void) const {return _OuterCyl;}
  const double* OCPos(void) const {return _OCPos;}
  double OCPosX(void) const {return _OCPos[0];}
  double OCPosY(void) const {return _OCPos[1];}
  double OCPosZ(void) const {return _OCPos[2];}
  //  rear end plate (REP)
  double REPInRad(void) const {return _RearCyl.getInnerRadius();}
  double REPOutRad(void) const {return _RearCyl.getOuterRadius();}
  double REPLength(void) const {return _RearCyl.getLength();}
  const DchSimpleCyl& getRearCyl(void) const {return _RearCyl;}
  const double* REPPos(void) const {return _REPPos;}
  double REPPosX(void) const {return _REPPos[0];}
  double REPPosY(void) const {return _REPPos[1];}
  double REPPosZ(void) const {return _REPPos[2];}
  //  forward end plate (FEP)
  double FEPInRad(void) const {return _ForwCyl.getInnerRadius();}
  double FEPOutRad(void) const {return _ForwCyl.getOuterRadius();}
  double FEPLength(void) const {return _ForwCyl.getLength();}
  const DchSimpleCyl& getForwCyl(void) const {return _ForwCyl;}
  const double* FEPPos(void) const {return _FEPPos;};
  double FEPPosX(void) const {return _FEPPos[0];};
  double FEPPosY(void) const {return _FEPPos[1];};
  double FEPPosZ(void) const {return _FEPPos[2];};
  //  layer numbering goes from 1 to 40...
  const senseWire& sWire(int lay) const { assert ( lay>0 && lay<=40 );
                                  return _sWires[lay-1];}
  const cellStruct& cell(int lay) const { assert ( lay>0 && lay<=40 );
                                  return _cellStruc[lay-1];}
  const fieldWire& fWire(int lay, int wire) const { assert ( lay>0 && lay<=40 );  
                                  assert (wire>=0 && wire<7);  
                                  return _cellStruc[lay-1].FWires[wire];}
  int nLayers(void) const { return _nLayers; }

  int version(void) const { return _version; }


  void initWires(void);
// One line printing
  void print(std::ostream& o) const;
// more complete printing
  void printAll(std::ostream& o=std::cout) const;

  std::string GasMat() const { return _GasMaterial; }
  std::string ICMat()  const { return _ICMaterial; }
  std::string IWGrdMat()  const { return _IWGrdMaterial; }
  std::string OCMat()  const { return _OCMaterial; }
  std::string REPMat(int iSub=0) const { return _REPMaterialSubs.at(iSub); }
  std::string FEPMat(int iSub=0) const { return _FEPMaterialSubs.at(iSub); }
  //  rear end plate (REP)
  size_t REPnSubs() const { return _RearCylSubs.size(); }
  double REPInRad(int iSub) const {return _RearCylSubs.at(iSub).getInnerRadius();}
  double REPOutRad(int iSub) const {return _RearCylSubs.at(iSub).getOuterRadius();}
  double REPLength(int iSub) const {return _RearCylSubs.at(iSub).getLength();}
  bool   REPisPhiSegmented(int iSub) const {return _RearCylSubs.at(iSub).isPhiSegmented();}
  double REPPhi0(int iSub) const {
          if (_RearCylSubs.at(iSub).getInnerRadius()) {
                  //return ((DchPhiSegmCyl*) &_RearCylSubs.at(iSub))->getPhi0();
                  return _RearCylSubs.at(iSub).getPhi0();
          } else {
                  return 0.0;
          }
  }
  double REPsolidDPhi(int iSub) const {
          if (_RearCylSubs.at(iSub).getInnerRadius()) {
                  //return ((DchPhiSegmCyl*) &_RearCylSubs.at(iSub))->getSolidDeltaPhi();
                  return _RearCylSubs.at(iSub).getSolidDeltaPhi();
          } else {
                  return 0.0;
          }
  }
  double REPhollowDPhi(int iSub) const {
          if (_RearCylSubs.at(iSub).getInnerRadius()) {
                  //return ((DchPhiSegmCyl*) &_RearCylSubs.at(iSub))->getHollowDeltaPhi();
                  return _RearCylSubs.at(iSub).getHollowDeltaPhi();
          } else {
                  return 0.0;
          }
  }
  const DchSimpleCyl& getRearCyl(int iSub) const {return _RearCylSubs.at(iSub);}
  const double* REPPos(int iSub) const {return _REPPosSubs.at(iSub);}
  double REPPosX(int iSub) const {return _REPPosSubs.at(iSub)[0];}
  double REPPosY(int iSub) const {return _REPPosSubs.at(iSub)[1];}
  double REPPosZ(int iSub) const {return _REPPosSubs.at(iSub)[2];}
  //  forward end plate (FEP)
  size_t FEPnSubs() const { return _ForwCylSubs.size(); }
  double FEPInRad(int iSub) const {return _ForwCylSubs.at(iSub).getInnerRadius();}
  double FEPOutRad(int iSub) const {return _ForwCylSubs.at(iSub).getOuterRadius();}
  double FEPLength(int iSub) const {return _ForwCylSubs.at(iSub).getLength();}
  const DchSimpleCyl& getForwCyl(int iSub) const {return _ForwCylSubs.at(iSub);}
  const double* FEPPos(int iSub) const {return _FEPPosSubs.at(iSub);}
  double FEPPosX(int iSub) const {return _FEPPosSubs.at(iSub)[0];}
  double FEPPosY(int iSub) const {return _FEPPosSubs.at(iSub)[1];}
  double FEPPosZ(int iSub) const {return _FEPPosSubs.at(iSub)[2];}
  bool   FEPisPhiSegmented(int iSub) const {return _ForwCylSubs.at(iSub).isPhiSegmented();}
  double FEPPhi0(int iSub) const {
          if (_ForwCylSubs.at(iSub).getInnerRadius()) {
                  //return ((DchPhiSegmCyl*) &_ForwCylSubs.at(iSub))->getPhi0();
                  return _ForwCylSubs.at(iSub).getPhi0();
          } else {
                  return 0.0;
          }
  }
  double FEPsolidDPhi(int iSub) const {
          if (_ForwCylSubs.at(iSub).getInnerRadius()) {
                  //return ((DchPhiSegmCyl*) &_ForwCylSubs.at(iSub))->getSolidDeltaPhi();
                  return _ForwCylSubs.at(iSub).getSolidDeltaPhi();
          } else {
                  return 0.0;
          }
  }
  double FEPhollowDPhi(int iSub) const {
          if (_ForwCylSubs.at(iSub).getInnerRadius()) {
                  //return ((DchPhiSegmCyl*) &_ForwCylSubs.at(iSub))->getHollowDeltaPhi();
                  return _ForwCylSubs.at(iSub).getHollowDeltaPhi();
          } else {
                  return 0.0;
          }
  }

private:

  friend class DchGDchP;
  friend class DchGDchCdbR;
  friend bool testCdb(const DchGDch*, const DchGDch*);

  DchGDch( void );  // default constructor not implemented

// Data
//  General quantities
  std::string  _name;        // detector name and version number
  int _version;
  int _nLayers;
  DchSimpleCyl _GasCyl;
  double _GasPos[3];
  std::string _GasMaterial;
  DchSimpleCyl _InnerCyl;
  double _ICPos[3];
  std::string _ICMaterial;
  DchSimpleCyl _InnerWGrdCyl;
  double _IWGrdPos[3];
  std::string _IWGrdMaterial;
  DchSimpleCyl _OuterCyl;
  double _OCPos[3];
  std::string _OCMaterial;
  DchSimpleCyl _RearCyl;
  double _REPPos[3];
  //std::vector<DchSimpleCyl> _RearCylSubs;
  std::vector<DchPhiSegmCyl> _RearCylSubs;
  std::vector<double*> _REPPosSubs;
  std::vector<std::string> _REPMaterialSubs;
  DchSimpleCyl _ForwCyl;
  double _FEPPos[3];
  //std::vector<DchSimpleCyl> _ForwCylSubs;
  std::vector<DchPhiSegmCyl> _ForwCylSubs;
  std::vector<double*> _FEPPosSubs;
  std::vector<std::string> _FEPMaterialSubs;
  senseWire _sWires[nLayer];
  cellStruct _cellStruc[nLayer];
  double _rOffset[nLayer];
  double _zoffset;
  double _zlen;
  double _cellHeight;


};

std::ostream&  operator << (std::ostream& o, const DchGDch&);

#endif









