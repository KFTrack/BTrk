//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: DchDetector.hh 91 2010-01-14 12:37:23Z stroili $ 
//
// Description:
//      DchDetector Class - 
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Roberto Stroili
//
// History (add to end):
//      Stroili May, 11 1997  - creation
//
// Copyright Information:
//      Copyright (C) 1997              INFN
//      
//
//------------------------------------------------------------------------

#ifndef DCHDETECTOR_HH
#define DCHDETECTOR_HH

#include <assert.h>

#include "DchGeom/DchLayer.hh"
#include "DchGeom/DchCell.hh"
//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

class DetCylinder;
class DetSet;
class DetSurfaceSet;
class DetElem;
//class DetAlignElem;
class DchGasSet;
class DchCylType;
class DchVolType;
class DchVolElem;
class DchHyperb;
class DchHyperType;
class DchSuperLayer;
class DchSWire;
//class DchGlobalAlign;
//class DchPlateAlign;
//class DchPlateDefl;
//class DchWireAlign;
class TrkVolume;
class TrkFit;
class DchGDch;
#include <iosfwd>
#include <vector>
//class HepString;

class DchDetector {
public:

  //  Constructors
  static DchDetector* GetInstance(){return fgInstance;} 
  DchDetector( const DchGDch&, bool deb=false ); // builds the whole 
                                                 // transient tree from dbio
  //  Destructor
  virtual ~DchDetector( );

  //  accessors
  const DetSet& dchSet() const  { return *_dchset; }
  const DetSet& dchGasSet() const  { return *_dchgasset; }
  const DetSet& dchLayerSet() const  { return *_dchlayerset; }
  DchVolElem* getDchGasVol(void) const;
  const TrkVolume& trkVolume() const { return *_trkVolume; }

  bool debug(void) { return _debug; }
  int version(void) const { return _version; }
  void calcZOffset(void);

  double zlen(void) const { return _dclayer[0]->zLength(); }
  double zOffSet(void) const {return _zOffset;}
  double rIn(void) const;
  double rOut(void) const;
  //  these are in the BaBar refernece frame
  double xWire(int lay, int wire, double z) const;
  double yWire(int lay, int wire, double z) const;
    
  // access superlayer info
  int firstSLayNum(void) const                        {return _firstSlayNum;}
  int lastSLayNum(void) const                          {return _lastSlayNum;}
  int nSuper(void) const {return _nSlay;}
  int nAxialSuper(void) const {return _nAxSlay;}
  int nStereoSuper(int iview) const { return iview == 0 ? 0 : 
    _nSterSlay[(iview+1)/2]; }
  const DchSuperLayer* firstSlay(void) const             {return _firstSlay;}
  const DchSuperLayer* lastSlay(void) const               {return _lastSlay;}
  const DchSuperLayer* firstSlayInView(int iview) const
                                          {return _firstSlayInView[iview+1];}
  const DchSuperLayer* lastSlayInView(int iview) const 
                                           {return _lastSlayInView[iview+1];}
  const DchSuperLayer* slayer(int index) const;

  //  access layer info
  int firstLayNum(void) const                           {return _firstLayer;}
  int lastLayNum(void) const              {return _firstLayer + _nLayer - 1;}
  int nLayer(void) const                                    {return _nLayer;}
  DchLayer* dchLayer(int laynum) const;  // "slow" access
  DchLayer* getDchLayer(int laynum) const { return _dclayer[laynum-
							    firstLayNum()]; }
  const DchCell&  getDchCell(int laynum, int cell) const {
    return getDchLayer(laynum)->getCell(cell); }
  const DchLayer* firstLayer(void) const {return _dclayer[0];}
  const DchLayer* lastLayer(void) const   {return _dclayer[lastLayNum()-
					                     firstLayNum()];}
  const DchLayer* nextLayer(int lay) const 
                                      { return _nextlay[lay-firstLayNum()]; }
  const DchLayer* prevLayer(int lay) const   
                                      { return _prevlay[lay-firstLayNum()]; }
  const DchLayer* nextLayer(const DchLayer* layer) const
                          { return _nextlay[layer->layNum()-firstLayNum()]; } 
  const DchLayer* prevLayer(const DchLayer* layer) const
                          { return _prevlay[layer->layNum()-firstLayNum()]; } 
  const DetCylinder* innerCylPID(void) const            { return _innerCyl; }
  bool exist(int lay) const;
  bool existDet(int lay) const;
  //   int view(int lay) const;
  int nWires(int lay) const            { return getDchLayer(lay)->nWires(); }
  double rMid(int lay) const             { return getDchLayer(lay)->rMid(); }
  double rEnd(int lay) const             { return getDchLayer(lay)->rEnd(); }
  double zLength(int lay) const       { return getDchLayer(lay)->zLength(); }
  //  some functions to get cell widths
  double cellWidNormFact(int lay) const { assert(lay>=firstLayNum() && 
						 lay<=lastLayNum());
                                          return _cellWidthNormFact[lay-1]; }
  double cellWidth(int lay) const { assert(lay>=firstLayNum() && 
					   lay<=lastLayNum());
                                      return getDchLayer(lay)->cellWidth(); }

  double cellWidth(int lay, double z) const { assert(lay>=firstLayNum() && 
						     lay<=lastLayNum());
                                     return getDchLayer(lay)->cellWidth(z); }

  // alignment stuff
  // ---------------
  // apply alignment
  /* void apply(const DchGlobalAlign&); // align the whole DCH
  void apply(const DchPlateAlign&); // align the end plates
  void apply(const DchPlateDefl&); // correct for end-plate deflections
  void apply(const DchWireAlign&); // single wire corrections
  // access current alignment
  const DchGlobalAlign* globalAlignment() const { return _globalalign; }
  const DchPlateDefl* plateDeflection() const { return _platedefl; }
  */
  // modifiers
  void setDebug(bool deb);

  // ad hoc transformation... for backward compatibility
  double zBaBar(double z) const { return z + _zOffset; }
  double gotoLoc(double z) const { return z - _zOffset; }

  // One line printout
  void print(std::ostream& o=std::cout) const;
  // DetectorModel printout
  void printAll(std::ostream& o) const;
  // print coordinate for some wires; for comparison with MC geoemtry
  void printDebug(std::ostream& o) const;
  // dump to file
  void writeToFile(std::string filename) const;
  // get Dch entrance momentum
  Hep3Vector entranceMomentum( const TrkFit* fit ) const;  
  Hep3Vector entranceMomentum( const TrkFit* fit, double& flt ) const;
  bool entranceMomentum( const TrkFit* fit, double& flt, Hep3Vector& mom ) 
    const;

private:
  // private member functions
  // 
  DchDetector( void ) {;}               // dummy constructor
  void buildpointers(void);		// make the Layer & SuperLayer pointers
  void buildSuperLayers(void);          // nuild super-layers
  const DchSuperLayer* getDchSLayer(int slaynum) const;

  //   data members
  // DETECTORS
  DetSet*           _dchset;  // full Dch detector set
  DetSet*        _dchgasset;  // Dch gas volume set
  DetSet*      _dchlayerset;  // Dch gas volume set
  DetSet*     _dchslayerset;  // Dch gas volume set
  TrkVolume*     _trkVolume;  // tracking volume for the DCH
  DetCylinder*    _innerCyl;  // inner cylinder (for PID reconstruction)

  //  main volumes
  DchCylType*     _inCylType;   // inner cylinder
  DchCylType*     _outCylType;  // outer cylinder
  DchVolType*     _rEPType;     // rear End Plate
  DchVolType*     _fEPType;     // forward End Plate
  std::vector<DchVolType*>     _rEPSubsType;     // rear End Plate sub emelents
  std::vector<DchVolType*>     _fEPSubsType;     // forward End Plate sub emelents

  DchVolType*     _gasVolType;
  std::vector<std::vector<DchSWire*> >  _senseWire;  // pointer to all the sense wires
                                                     //  _senseWire[layer][wire]

  double _cellWidthNormFact[40]; // normalization factor for cell width

  //  pointer arrays for Layers and SuperLayers
  DchLayer* _dclayer[40];	// 40 layers for chamber
  DchLayer* _nextlay[40];	// pointer to next layer
  DchLayer* _nextlayinvw[40];	// pointer to next layer in same view
  DchLayer* _prevlay[40];	// pointer to previous layer in view
  DchLayer* _prevlayinvw[40];	// pointer to previous layer in same view
  DchSuperLayer** slayList;    // will be array of pointers to superlayers
  const DchSuperLayer* _firstSlayInView[3];
  const DchSuperLayer* _lastSlayInView[3];
  const DchSuperLayer* _firstSlay;
  const DchSuperLayer* _lastSlay;
  
  // description
  int               _version;
  bool              _debug;
  int               _firstLayer;
  int               _nLayer;
  //  some cached information
  int               _firstSlayNum;
  int               _lastSlayNum;
  int               _nSlay;
  int               _nAxSlay;
  int               _nSterSlay[2];
  double            _zOffset;

  // alignment
  //DchGlobalAlign* _globalalign; // last global alignment applied
  //int _nglobal; // number of times global alignment has been applied
  //DchPlateDefl*   _platedefl;

  static DchDetector *fgInstance;
};

std::ostream&  operator << (std::ostream& o, const DchDetector&);

#endif
