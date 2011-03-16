//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: DchDetector.cc 123 2010-04-29 14:41:45Z stroili $ 
//
// Description:
//      Class DchDetector
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Roberto Stroili  -  INFN
//
// History (add to end):
//      Stroili   May 11, 1997  - creation
//
// Copyright Information:
//      Copyright (C) 1997             INFN
//      
// Revision History:
//	20000511  Michael Kelsey -- Make _dchlayerset as DetSurfaceSet;
//		  move inner cylinder reference surface to top of ctor.
//	20020619  M. Kelsey -- Pass self into DchLayer ctor.
//------------------------------------------------------------------------

//----------------
// BaBar Header --
//----------------
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "DchGeom/DchDetector.hh"

//---------------
// C++ Headers --
//---------------
#include <iostream>
#include <fstream>
#include <iomanip>
#include <assert.h>
#include <stdlib.h>
#include <string>
#include <vector>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "BaBar/Constants.hh"
#include "ErrLogger/ErrLog.hh"
#include "DchGeomBase/DchGDch.hh"
#include "DchGeomBase/DchCellAddr.hh"
#include "DchGeomBase/DchWirePar.hh"
#include "DchGeom/DchCylType.hh"
#include "DchGeom/DchVolType.hh"
#include "DchGeom/DchVolElem.hh"
//#include "DchGeom/DchGlobalAlign.hh"
//#include "DchGeom/DchPlateAlign.hh"
#include "DchGeom/DchPlateDefl.hh"
//#include "DchGeom/DchWireAlign.hh"
#include "DchGeom/DchHyperb.hh"
#include "DchGeom/DchHyperType.hh"
#include "DchGeom/DchSelections.hh"
#include "DchGeom/DchSWire.hh"
#include "DchGeom/DchFWire.hh"
#include "DchGeom/DchCell.hh"
#include "DchGeom/DchLayer.hh"
#include "DchGeom/DchSuperLayer.hh"
#include "DchGeom/DchSearchId.hh"
#include "CLHEP/Geometry/Transformation.h"
#include "CLHEP/Geometry/AngleSets.h"
#include "CLHEP/Vector/ThreeVector.h"
//#include "CLHEP/String/Strings.h"
#include "DetectorModel/DetMaterial.hh"
#include "DetectorModel/DetSet.hh"
#include "DetectorModel/DetSurfaceElem.hh"
#include "DetectorModel/DetSurfaceSet.hh"
#include "DetectorModel/DetCylinder.hh"
#include "DetectorModel/DetAlignElem.hh"
#include "DetectorModel/Intersection.hh"
#include "TrkBase/TrkErrCode.hh"
#include "TrkBase/TrkFit.hh"
#include "TrkBase/TrkDifTraj.hh"
//#include "TrkGeom/TrkSimpVolume.hh"
#include "BaBar/BbrCollectionUtils.hh"
#include "MatEnv/MatDBInfo.hh"
using std::endl;
using std::fstream;
using std::ios;
using std::ofstream;
using std::ostream;

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

static const int _layPerSL = 1;

DchDetector *DchDetector::fgInstance=0;
//----------------
// Constructors --
//----------------

DchDetector::DchDetector(const DchGDch& gdch, bool deb) :
  _inCylType(0), _outCylType(0), _rEPType(0), _fEPType(0), _senseWire(
      gdch.nLayers()), _version(10), _debug(deb), _firstLayer(1), _nLayer(gdch.nLayers())
//, _globalalign(0), _nglobal(0), _platedefl(0)
{
  fgInstance=this;
  MatDBInfo* mtdbinfo=new MatDBInfo;
  
  double radii[2];
  double length;

  if (deb) gdch.printAll();

  //  "common" euler angles
  EulerAngles euler(0., 0., 0.);
  double zR = 0., zF = 0.;

  // Inner cylinder reference surface
  Hep3Vector x0(gdch.ICPosX(), gdch.ICPosY(), gdch.ICPosZ());
  HepTransformation tr0(x0, euler);
  DetCylinder incyl(tr0, gdch.ICInRad());

  // Initialize set pointers

  _dchset = new DetSet("Dch Detector", 1);
  _dchgasset = new DetSet("Dch gas volume Detector", 2);
  _dchlayerset = new DetSurfaceSet("Dch layer set", 3, incyl, 1, 1);
  _dchslayerset = new DetSet("Dch super layer set", 4);
  *_dchset += *_dchgasset;

  // physical model of chamber
  // =========================

  // some physical objects are available only for dbio-generated geometry 
  // (version >= 10)
  if (gdch.version() >= 10) {
    // inner cylinder
    radii[0] = gdch.ICInRad();
    radii[1] = gdch.ICOutRad();
    length = gdch.ICLength();
    const DetMaterial* material = mtdbinfo->findDetMaterial("IT-InRad");// gblEnv->getGen()->findDetMaterial("Beryllium");
    if(material)material->printAll(std::cout); 
    //assert(0 != material);
    _inCylType = new DchCylType("Dch inner cylinder", radii, length, material,
        1);

    DetSurfaceElem* inelem = new DetSurfaceElem(_inCylType,
        "Dch Inner Cylinder", 1, incyl);

    //    *_dchset += *inelem;
    *_dchgasset += *inelem;

    // outer cylinder
    radii[0] = gdch.OCInRad();
    radii[1] = gdch.OCOutRad();
    length = gdch.OCLength();
    material = 0;//gblEnv->getGen()->findDetMaterial("Carbon fiber");
    //assert(0 != material);
    _outCylType = new DchCylType("Dch outer cylinder", radii, length, material,
        2);
    Hep3Vector x1(gdch.OCPosX(), gdch.OCPosY(), gdch.OCPosZ());
    HepTransformation tr1(x1, euler);
    DetCylinder outcyl(tr1, radii[0]);
    DetSurfaceElem* outelem = new DetSurfaceElem(_outCylType,
        "Dch Outer Cylinder", 2, outcyl);
    //    *_dchset += *outelem;
    *_dchgasset += *outelem;

    // rear end plate
    radii[0] = gdch.REPInRad();
    radii[1] = gdch.REPOutRad();
    length = gdch.REPLength();
    material = 0;//gblEnv->getGen()->findDetMaterial("Aluminum");
    //assert(0 != material);

    Hep3Vector rEPpos(gdch.REPPosX(), gdch.REPPosY(), gdch.REPPosZ());
    HepTransformation rEP_tr(rEPpos, euler);

    _rEPType = new DchVolType("Dch rear end-plate", radii[0], radii[1], -length
        / 2., length / 2., material, 31);
    DetVolumeElem* repelem = new DetVolumeElem(_rEPType, "Dch rear endplate",
        31, rEP_tr);

    zR = gdch.REPPosZ() - length * .5;

    *_dchgasset += *repelem;

    // forward end plate
    radii[0] = gdch.FEPInRad();
    radii[1] = gdch.FEPOutRad();
    length = gdch.FEPLength();
    material = 0;//gblEnv->getGen()->findDetMaterial("Aluminum");
    //assert(0 != material);

    Hep3Vector fEPpos(gdch.FEPPosX(), gdch.FEPPosY(), gdch.FEPPosZ());
    HepTransformation fEP_tr(fEPpos, euler);

    _fEPType = new DchVolType("Dch forward end-plate", radii[0], radii[1],
        -length / 2., length / 2., material, 32);
    DetVolumeElem* fepelem = new DetVolumeElem(_fEPType,
        "Dch forward endplate", 32, fEP_tr);

    zF = gdch.FEPPosZ() + length * .5;

    *_dchgasset += *fepelem;

  }
  // gas volume
  const DetMaterial* material = mtdbinfo->findDetMaterial("IT-gas1");//gblEnv->getGen()->findDetMaterial("GasWire");
  material->printAll(std::cout);
  //assert(0 != material);

  radii[0] = gdch.GasInRad();
  radii[1] = gdch.GasOutRad();
  length = gdch.GasLength();

  Hep3Vector vol_pos(gdch.GasPosX(), gdch.GasPosY(), gdch.GasPosZ());
  HepTransformation vol_tr(vol_pos, euler);

  _gasVolType = new DchVolType("Dch gas volume", radii[0], radii[1], -length
      / 2., length / 2., material, 33);
  DchVolElem* gasvol =
      new DchVolElem(_gasVolType, "Dch volume gas", 34, vol_tr);
  _gasVolType->setDebug(_debug);
  gasvol->setDebug(_debug);

  *_dchgasset += *gasvol;

  //  inner fake cylinder, needed by Dch PID reconstruction to get a reasonable
  //  momentum. The radius is 2 cm inside the gas volume
  double inrad = gdch.GasInRad() + 2.;
  _innerCyl = new DetCylinder(vol_tr, inrad);

  // calculate z offset, it needs GasVol
  calcZOffset();

  //  set volume radii == end-plate radii as they're larger than the gas 
  //  volume ones
  radii[0] = gdch.REPInRad();
  radii[1] = gdch.REPOutRad();

  double tol(.2);
  //_trkVolume = new TrkSimpVolume("Dch Tracking Volume", radii[0] - tol,
  //    radii[1] + tol, zR - tol, zF + tol);

  // tracking chamber
  // ================

  if (_debug && ErrLogging(debugging)) {
    ostream& outstream = ErrMsg(debugging);
    for (int lay = 1; lay < nLayer() + 1; ++lay) {
      outstream << " LAYER # " << lay << "\n" << "            # of cells "
          << gdch.sWire(lay).nCells << "\n" << "            rEnd       "
          << gdch.sWire(lay).rEnd << "\n" << "            rForw      "
          << gdch.sWire(lay).rForw << "\n" << "            phiEnd     "
          << gdch.sWire(lay).phiEnd << "\n" << "            phiForw    "
          << gdch.sWire(lay).phiForw << "\n" << "            twist      "
          << gdch.sWire(lay).twist << "\n" << "            stereo     "
          << gdch.sWire(lay).stereo << "\n" << "            deltaPhi   "
          << gdch.sWire(lay).deltaPhi << "\n" << "            sag        "
          << gdch.sWire(lay).sag << "\n" << "            diameter   "
          << gdch.sWire(lay).diameter << "\n" << "            volts      "
          << gdch.sWire(lay).volts << endmsg;
    }
  }
  //  build the "layers", envelop of the sense wires

  // build wire type array

  // here the layer indexing goes between 1 and 40 as the layer numbers !!!
  int ilay; // to be used in another loop
  for (ilay = firstLayNum(); ilay <= lastLayNum(); ilay++) {
    const senseWire& sWir = gdch.sWire(ilay);
    // build layer surface
    DchHyperb LayerSurf(vol_tr, sWir.rEnd, fabs(gdch.GasLength()), sWir.twist);
    //   VOID NOT YET AVAILABLE use GasWire for the moment
    material = 0;//gblEnv->getGen()->findDetMaterial("GasWire");
    //assert(0 != material);
    // build the basic cell
    std::vector<DchFWire*> fieldWires;
    for (int fwires = 0; fwires < gdch.cell(ilay).nofFWires; fwires++) {
      // wire coordinate on rear end-plate
      fieldWire fwir = gdch.fWire(ilay, fwires);
      HepPoint point1(fwir.rEnd * cos(fwir.phiEnd), fwir.rEnd
          * sin(fwir.phiEnd), fwir.zR);
      // wire coordinate on forward end-plate
      HepPoint point2(fwir.rForw * cos(fwir.phiForw), fwir.rForw * sin(
          fwir.phiForw), fwir.zFw);
      fieldWires.push_back(new DchFWire(point1, point2));
    }
    DchCell dchCell(fieldWires);

    HepPoint rearPoint(sWir.rEnd * cos(sWir.phiEnd), sWir.rEnd * sin(
        sWir.phiEnd), sWir.zR);
    HepPoint frntPoint(sWir.rForw * cos(sWir.phiForw), sWir.rForw * sin(
        sWir.phiForw), sWir.zFw);
    int nCells = sWir.nCells;
    double sag = sWir.sag;
    double deltaCell = sWir.deltaPhi;
    _senseWire[ilay - 1].clear();
    _senseWire[ilay - 1].reserve(nCells);
    for (int cell = 0; cell < nCells; cell++) {
      HepPoint r(rearPoint);
      r.rotateZ(deltaCell * cell);
      HepPoint f(frntPoint);
      f.rotateZ(deltaCell * cell);
      _senseWire[ilay - 1].push_back(new DchSWire(r, f, sag));
    }

    DchHyperType* hyperType = new DchHyperType("stereo type", material, 44,
        sWir.twist, gdch.GasLength());
    DchLayer* nlayer = new DchLayer(hyperType, "Dch Layer", 100 + ilay,
        LayerSurf, ilay, _senseWire[ilay - 1], dchCell, *this, _debug);
    *_dchlayerset += *nlayer;
  }

  // add LayerSet to DCHset
  *_dchset += *_dchlayerset;

  // build pointers to make navigation faster
  buildpointers();

  // do superlayers
  // --------------
  buildSuperLayers();

  // add SuperLayerSet to DCHset (_dchslayerset is "built" in buildSuperLayers)
  *_dchset += *_dchslayerset;

  // set "nominal" cell height for each layer
  // =========================================================
  double rOther;
  double height;
  for (ilay = firstLayNum(); ilay <= lastLayNum(); ilay++) {
    if (ilay == firstLayNum()) {
      rOther = getDchLayer(ilay + 1)->rMid();
    } else {
      rOther = getDchLayer(ilay - 1)->rMid();
    }
    height = fabs(getDchLayer(ilay)->rMid() - rOther);
    getDchLayer(ilay)->setCellHeight(height);
    //
    //  compute cell width normalization factors (wrt layer # 16)
    //
    _cellWidthNormFact[ilay - 1] = getDchLayer(ilay)->cellWidth()
        / getDchLayer(16)->cellWidth();
  }

  if (_debug && ErrLogging(debugging)) {
    ostream& outstream = ErrMsg(debugging);
    outstream << " SUPER-LAYERS \n";
    for (int islay = firstSLayNum(); islay <= lastSLayNum(); islay++) {
      getDchSLayer(islay)->print(outstream);
    }
    outstream << " LAYERS \n";
    for (ilay = firstLayNum(); ilay <= lastLayNum(); ilay++) {
      outstream << *getDchLayer(ilay) << "\n";
    }
    outstream << endmsg;
  }
  
}

//--------------
// Destructor --
//--------------
DchDetector::~DchDetector()
{

  // delete sense wires
  typedef std::vector<std::vector<DchSWire*> >::iterator iter;
  for (iter i = _senseWire.begin(); i != _senseWire.end(); ++i) {
    std::for_each(i->begin(), i->end(), babar::Collection::DeleteObject());
  }
  // delete some pointers
  //delete _trkVolume;
  //delete _globalalign;
  //delete _platedefl;
  //delete _innerCyl;
  //
  //  First delete all elements
  //
  DetElemList list;
  _dchset->select(list, dsMatchAll, 0, 0, 0);
  std::for_each(list.begin(), list.end(), babar::Collection::DeleteObject());

  int sectlim[2] = { 0, 1040 };
  std::string setStr("Dch");
  DchSearchId selSet(sectlim, &setStr);
  DetSetList slist;
  _dchset->setSelect(slist, &selectDchSet, (void*) &selSet);
  if (ErrLogging(trace)) {
    ErrMsg(trace) << "DchDetector: \tnumber of sets to destroy: "
        << slist.size() << endmsg;
  }
  if (_debug && ErrLogging(debugging)) {
    ostream& outstream = ErrMsg(debugging);
    DetSetList::iterator siter = slist.begin();
    for (siter = slist.begin(); siter != slist.end(); siter++) {
      outstream << "\t" << (*siter)->setName() << "\n";
    }
    outstream << endmsg;
  }
  std::for_each(slist.begin(), slist.end(), babar::Collection::DeleteObject());

  // delete pointer arrays
  //   delete [] _dclayer;
  //   delete [] _nextlay;
  //   delete [] _nextlayinvw;
  //   delete [] _prevlay;
  //   delete [] _prevlayinvw;
  //   delete [] _firstSlayInView;
  //   delete [] _lastSlayInView;
  delete[] slayList;
  delete _inCylType;
  delete _outCylType;
  delete _rEPType;
  delete _fEPType;
  delete _gasVolType;
  //  delete _dchslayerset;
  //  delete _dchlayerset;
  //  delete _dchgasset;
  //   delete _dchset;

}

DchLayer*
DchDetector::dchLayer(int laynum) const
{
  DetElemList elist;
  _dchset->select(elist, &findDchLayer, 0);
  int nlayers = elist.size();

  // check if number of layers is the expected one...
  assert ( nlayers == nLayer() );
  assert ( laynum>=firstLayNum() && laynum<=lastLayNum() );

  // do not suppose layers are in order...
  DetElemList::iterator diter = elist.begin();
  for (diter = elist.begin(); diter != elist.end(); diter++) {
    if (((DchLayer*) *diter)->layNum() == laynum) {
      return (DchLayer*) *diter;
    }
  }
  return 0;
}

DchVolElem*
DchDetector::getDchGasVol(void) const
{
  DetElemList elist;
  _dchset->select(elist, &selectDchGasVol, 0);
  int nvol = elist.size();

  // check if there's only one volume gas
  assert ( nvol == 1 );

  return (DchVolElem*) *(elist.begin());
}

const DchSuperLayer*
DchDetector::getDchSLayer(int slaynum) const
{
  DetSetList slist;
  _dchset->setSelect(slist, &findDchSuperLayer, 0);
  int nslayers = slist.size();

  // check if number of layers is the expected one...
  assert ( nslayers == nSuper() );
  assert ( slaynum>=firstSLayNum() && slaynum<=lastSLayNum() );

  // do not suppose layers are in order...
  DetSetList::iterator siter = slist.begin();
  for (siter = slist.begin(); siter != slist.end(); siter++) {
    if (((DchSuperLayer*) *siter)->slayNum() == slaynum) {
      return (DchSuperLayer*) *siter;
    }
  }
  return 0;
}

const DchSuperLayer*
DchDetector::slayer(int slaynum) const
{
  assert ( slaynum>=firstSLayNum() && slaynum<=lastSLayNum() );
  return slayList[slaynum - 1];
}

void
DchDetector::buildpointers(void)
{
  // first layers
  for (int index = 0; index < lastLayNum() - firstLayNum() + 1; index++) {
    // initialize
    _nextlay[index] = 0;
    _prevlay[index] = 0;
    _nextlayinvw[index] = 0;
    _prevlayinvw[index] = 0;
    int layi = index + firstLayNum();

    if (!existDet(layi)) {
      ErrMsg(fatal) << " layer # " << layi << " does not exist!" << endmsg;
    }
    _dclayer[index] = dchLayer(layi);

    //next and previous pointers
    if (existDet(layi + 1)) {
      _nextlay[index] = dchLayer(layi + 1);
    }
    if (existDet(layi - 1)) {
      _prevlay[index] = _dclayer[index - 1];
    }

    //next in view pointer
    int iview = dchLayer(layi)->view();
    int jndex;
    for (jndex = index + 1; jndex < lastLayNum() - firstLayNum() + 1; jndex++) {
      int layj = jndex + firstLayNum();
      if (!existDet(layj)) {
        ErrMsg(fatal) << " layer # " << layj << " does not exist!" << endmsg;
      }
      if (iview != dchLayer(layj)->view()) continue;
      _nextlayinvw[index] = dchLayer(jndex + 1);
      break;
    } //(int jndex=index+1; _dclayer[jndex].Exist(); jndex++) 

    //prev in view pointer
    for (jndex = index - 1; jndex >= 0; jndex--) {
      int layj = jndex + firstLayNum();
      if (!existDet(layj)) {
        ErrMsg(fatal) << " layer # " << layj << " does not exist!" << endmsg;
      }
      if (iview != dchLayer(layj)->view()) continue;
      _prevlayinvw[index] = _dclayer[jndex];
      break;
    } //(int jndex=index+1; _dclayer[jndex].exist(); jndex++) 

  } //(int index=0; _dclayer[index].exist(); index++)
}

void
DchDetector::buildSuperLayers(void)
{
  // some initializations
  _nAxSlay = 0;
  _nSterSlay[0] = _nSterSlay[1] = 0;

  _nSlay = nLayer() / _layPerSL;

  _firstSlayNum = 1;
  _lastSlayNum = _nSlay;
  slayList = new DchSuperLayer*[_nSlay];

  // initialize pointers 
  _firstSlayInView[0] = _firstSlayInView[1] = _firstSlayInView[2] = 0;
  _lastSlayInView[0] = _lastSlayInView[1] = _lastSlayInView[2] = 0;

  int islay;

  // build the SuperLayers
  for (islay = 0; islay < _nSlay; islay++) {
    DchSuperLayer* superlay = new DchSuperLayer("Dch SuperLayer", 1001 + islay);
    *_dchslayerset += *superlay;
    slayList[islay] = superlay;
  }

  _firstSlay = slayList[0];
  _lastSlay = slayList[_nSlay - 1];

  // set pointers to Layers in SuperLayers
  for (int lay = firstLayNum(); lay <= lastLayNum(); lay++) {
    int superlayer = (lay - 1) / _layPerSL;
    int index = (lay - 1) % _layPerSL;
    slayList[superlayer]->addLayer(index, getDchLayer(lay));
  }

  // update SuperLayer data members
  DchSuperLayer* oldSlayByView[3] = { 0, 0, 0 };

  for (islay = 0; islay < _nSlay; islay++) {
    const DchSuperLayer* prev = 0;
    const DchSuperLayer* next = 0;

    // SuperLayer view
    int iview = slayList[islay]->layer(0)->view();
    int viewIndex = iview + 1;

    // count SuperLayer types
    if (iview == 0) _nAxSlay++;
    else if (iview == -1) _nSterSlay[0]++;
    else if (iview == 1) _nSterSlay[1]++;

    // build pointer links
    if (islay > 0) prev = slayList[islay - 1];
    if (islay < _nSlay - 1) next = slayList[islay + 1];

    // fill first and last SuperLayer pointers
    if (firstSlayInView(iview) == 0) _firstSlayInView[viewIndex]
        = slayList[islay];
    _lastSlayInView[viewIndex] = slayList[islay];

    slayList[islay]->updateInfo(prev, next);

    // now the poiters to SuperLayers of the same view
    if (oldSlayByView[viewIndex] != 0) {
      oldSlayByView[viewIndex]->setNextInView(slayList[islay]);
      slayList[islay]->setPrevInView(oldSlayByView[viewIndex]);
    }

    oldSlayByView[viewIndex] = slayList[islay];
  }
}

bool
DchDetector::exist(int layer) const
{
  if ((layer < firstLayNum()) || (layer > lastLayNum())) {
    return false;
  } else {
    return getDchLayer(layer)->exist();
  }
}

bool
DchDetector::existDet(int layer) const
{
  if ((layer < firstLayNum()) || (layer > lastLayNum())) {
    return false;
  } else {
    return dchLayer(layer)->exist();
  }
}

double
DchDetector::rIn(void) const
{
  return ((DchVolType*) getDchGasVol()->detectorType())->rmin();
}

double
DchDetector::rOut(void) const
{
  return ((DchVolType*) getDchGasVol()->detectorType())->rmax();
}

double
DchDetector::xWire(int lay, int wire, double z) const
{
  assert (exist(lay));
  return getDchLayer(lay)->getWire(wire)->xWireDC(gotoLoc(z));
}

double
DchDetector::yWire(int lay, int wire, double z) const
{
  assert (exist(lay));
  return getDchLayer(lay)->getWire(wire)->yWireDC(gotoLoc(z));
}

void
DchDetector::calcZOffset(void)
{
  const HepPoint locPoint(0., 0., 0.);
  DchVolElem* gasVol = getDchGasVol();
  HepPoint globPoint = gasVol->transf().transFrom(locPoint);
  _zOffset = globPoint.z();
}

void
DchDetector::setDebug(bool deb)
{
  _debug = deb;
  if (getDchGasVol() != 0) getDchGasVol()->setDebug(deb);
}

//
//  Alignment functions.  These remove the previous alignment,
//  apply the new alignment, and keep a copy of the new alignment.
//  This also counts the number of times the alignment has been applied,
//  in order to be able to rebuild the whole tree when numerical precision
//  is lost (but that's not implemented yet).
//
/*
void
DchDetector::apply(const DchGlobalAlign& glob)
{
  if (ErrLogging(debugging)) {
    ErrMsg(debugging) << " DchDetector: apply global alignment " << endmsg;
  }
  //
  //  First, remove the old alignment (if any)
  //
  if (_globalalign != 0) {
    ErrMsg(error) << " Old global alignment found !!! \n"
        << "          it shouldn't be there..." << endmsg;
    delete _globalalign;
  }
  //
  //  Apply the new
  //
  _dchset->applyGlobal(glob.alignment());
  //
  //  Keep a copy the alignment and Increment the counter
  //
  _globalalign = new DchGlobalAlign(glob);
  _nglobal++;
  //
  calcZOffset();
  //
  // now align the sense wires...
  //
  for (int lay = 0; lay < nLayer(); lay++) { //loop over layers
    // now loop over wires
    for (int wire = 0; wire < getDchLayer(lay + 1)->nWires(); wire++) {
      _senseWire[lay][wire]->wireAlign(_globalalign->alignment());
    }
  }
}

void
DchDetector::apply(const DchPlateAlign& palign)
{
  if (ErrLogging(debugging)) {
    ErrMsg(debugging) << " DchDetector: apply endplate alignment " << endmsg;
  }
  // 
  //  First, remove the old alignment (if any) 
  // 

  // 
  //  Apply the new 
  // 

  // 
  //  Keep a reference to the alignment and Increment the counter 
  // 
  const DetAlignElem& rearAlign = palign.alignElements(0);
  const DetAlignElem& forwAlign = palign.alignElements(1);

  // 
  // now align the sense wires... 
  // 
  for (unsigned lay = 0; lay < nLayer(); ++lay) { //loop over layers
    // now loop over wires 
    for (unsigned wire = 0; wire < getDchLayer(lay + 1)->nWires(); ++wire) {
      _senseWire[lay][wire]->wireAlign(rearAlign.transform(),
          forwAlign.transform());
    }
  }
}

void
DchDetector::apply(const DchPlateDefl& defl)
{
  if (ErrLogging(trace)) {
    ErrMsg(trace) << " DchDetector: apply endplate deflection " << endmsg;
  }
  //
  //  First, remove the old alignment (if any)
  //
  if (_platedefl != 0) {
    ErrMsg(error) << " Old plate deflection found !!! \n"
        << "          it shouldn't be there..." << endmsg;
    delete _platedefl;
  }

  _platedefl = new DchPlateDefl(defl);
  // apply corrections
  const std::vector<DchWirePar>& v = _platedefl->corrList();
  for (unsigned i = 0; i < v.size(); ++i) {
    int id = v[i].getID();
    if (id > 0) { // check if there's a correction
      // default (null) correction
      // has ID=-1
      int layer = DchCellAddr::layerIs(id);
      assert ( layer>0 && layer<41); // check if the layer is correct
      int wire = DchCellAddr::wireIs(id);
      assert ( wire==0 ); // wire ID must be 0 for this correction to be valid

      // loop over wires of a layer
      unsigned nw = getDchLayer(layer)->nWires();
      for (wire = 0; wire < nw; ++wire)
        _senseWire[layer - 1][wire]->wireAlign(v[i]);
    }
  } // loop over pseudo-layers (they don't need to be ordered)
}

void
DchDetector::apply(const DchWireAlign& walign)
{
  if (ErrLogging(debugging)) {
    ErrMsg(debugging) << " DchDetector: apply single wire alignment " << endmsg;
  }

  // apply correction to single wire

  // loop over corrections
  const std::vector<DchWirePar>& v = walign.corrList();
  for (unsigned i = 0; i < v.size(); ++i) {
    const DchWirePar& wirecorr = v[i];
    int id = wirecorr.getID();

    // first some checks
    int layer = DchCellAddr::layerIs(id);
    int wire = DchCellAddr::wireIs(id);
    assert ( layer>0 && layer<41 );
    assert ( wire>=0 && wire<getDchLayer(layer)->nWires() );

    // apply correction
    _senseWire[layer - 1][wire]->wireAlign(wirecorr);
  }
}
*/
Hep3Vector
DchDetector::entranceMomentum(const TrkFit* fit, double& fltDch) const
{
  fltDch = 0.;

  if (fit == 0) {
    ErrMsg (error) << "null TrkFit pointer, return 0 momentum" << endmsg;
    return Hep3Vector(0., 0., 0.);
  }

  //  check if the track has hit's in the DCH

  //  if (fit->nDch() > 0) {

    //  get radius of low-range point
    HepPoint x = fit->traj().position(fit->traj().lowRange());
    double lowRangeRad = x.mag() * sin(x.theta());

    //--- Extrapolate the tracks to 'radius' (2 cm. inside the gas volume) ---
    if (lowRangeRad > innerCylPID()->radius()) {
      fltDch = fit->traj().lowRange();
    } else {
      TrkErrCode error = Intersection(fit->traj(), *innerCylPID()).intersect(
          fltDch);
      if (!error.success()) {
        ErrMsg(debugging) << "Failed to intersect: " << error << endmsg;
        fltDch = 0;
      }
    }
    //  }
  return fit->momentum(fltDch);
}

Hep3Vector
DchDetector::entranceMomentum(const TrkFit* fit) const
{
  double fltDch(0);
  return entranceMomentum(fit, fltDch);
}

bool
DchDetector::entranceMomentum(const TrkFit* fit, double& flt, Hep3Vector& mom) const
{
  mom = entranceMomentum(fit, flt);
  return (flt == 0.) ? false : true;
}

void
DchDetector::printAll(ostream& o) const
{
  o << "The Dch Detector contains the following set : " << "\n\t";
  _dchset->printAll(o);
  o << "\n";

  o << " SUPER-LAYERS \n";
  for (int islay = firstSLayNum(); islay <= lastSLayNum(); islay++) {
    getDchSLayer(islay)->print(o);
  }
  o << " LAYERS \n";
  for (int ilay = firstLayNum(); ilay <= lastLayNum(); ilay++) {
    o << *getDchLayer(ilay) << "\n";
  }

  int prec = o.precision();
  o.precision(6);
  // output global quantities: zOffset, chamber length, inner & outer radius
  o << zOffSet() << "   " << zlen() << "   " << rIn() << "   " << rOut()
      << endl;
  o.setf(ios::fixed, ios::floatfield);
  // now the layer quantities
  for (int lay = firstLayNum(); lay <= lastLayNum(); lay++) {
    // layer number, # of wires, radius on rear end-plate
    o.precision(4);
    o << lay << "\t" << nWires(lay) << "\t" << rEnd(lay) << "   ";
    o.precision(7);
    // phi offset of first wire on rear end-plate
    o << getDchLayer(lay)->phiOffset() - getDchLayer(lay)->dPhiz()
    // tan(stereo)
        << "   " << getDchLayer(lay)->stereo()
    // twist angle between mid and rear chamber
        << "   " << getDchLayer(lay)->dPhiz() << "   "
    // sag (for future use)
        << getDchLayer(lay)->sag() << endl;
  }
  // geometry version
  o << version() << endl;
  // reset precision
  o.precision(prec);
}

void
DchDetector::writeToFile(std::string filename) const
{
  ofstream outfile(filename.c_str());
  if (outfile) {
    print(outfile);
  } else {
    ErrMsg(error) << " dump(file): error opening output file" << endmsg;
  }
}

void
DchDetector::print(ostream& o) const
{
  o << "The Dch Detector contains the following set : ";
  _dchset->print(o);
  o << endl;
}

ostream&
operator<<(ostream& o, const DchDetector& a)
{
  a.print(o);
  return o;
}

void
DchDetector::printDebug(ostream& o) const
{
  o << "\tLayer #1\t Cell # 0\n";
  getDchLayer(1)->getWire(0)->printInfo(o);
  o << "\tLayer #1\t Cell # 75\n";
  getDchLayer(1)->getWire(75)->printInfo(o);
  o << "\tLayer #5\t Cell # 0\n";
  getDchLayer(5)->getWire(0)->printInfo(o);
  o << "\tLayer #5\t Cell # 75\n";
  getDchLayer(5)->getWire(75)->printInfo(o);
  o << "\tLayer #20\t Cell # 0\n";
  getDchLayer(20)->getWire(0)->printInfo(o);
  o << "\tLayer #20\t Cell # 75\n";
  getDchLayer(20)->getWire(75)->printInfo(o);
  o << "\tLayer #36\t Cell # 0\n";
  getDchLayer(36)->getWire(0)->printInfo(o);
  o << "\tLayer #36\t Cell # 75\n";
  getDchLayer(36)->getWire(75)->printInfo(o);
  o << "\tLayer #40\t Cell # 0\n";
  getDchLayer(40)->getWire(0)->printInfo(o);
  o << "\tLayer #40\t Cell # 75\n";
  getDchLayer(40)->getWire(75)->printInfo(o);
}
