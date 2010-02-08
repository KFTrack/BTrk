//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkHotListEmpty.cc,v 1.22 2004/07/05 05:43:27 hulsberg Exp $
//
// Description:
//     
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner
//
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
#include "TrkBase/TrkHotListEmpty.hh"
#include "TrkBase/TrkHitOnTrk.hh"
#include <assert.h>
#include <algorithm>

TrkHotListEmpty::TrkHotListEmpty(int nact, int nsv, int ndc,
                                 double sfr, double  efr)
{
  _nPhi = nsv/2; _nZ = nsv-_nPhi;
  _nAxial = ndc/3; _nStereo = ndc-_nAxial;
  _stFndRng = sfr;
  _endFndRng = efr;
  _firstdch = _lastdch = 0;
  for(unsigned isvt=0;isvt<5;isvt++)
    _svtpat[isvt] = TrkView(TrkEnums::noView);
}

TrkHotListEmpty::TrkHotListEmpty(unsigned nPhi, unsigned nZ,
                                 unsigned nAxial, unsigned nStereo,
                                 double sfr, double efr,
                                 unsigned firstdch,unsigned lastdch, 
                                 TrkView svtpattern[5],
				 const std::vector<unsigned>& inactive) :
  _nAxial(nAxial),_nStereo(nStereo),_nPhi(nPhi),_nZ(nZ),
  _stFndRng(sfr),_endFndRng(efr),_firstdch(firstdch),_lastdch(lastdch),
  _inactive(inactive)
{
  for(unsigned j=0;j<5;j++)
    _svtpat[j] = svtpattern[j];
}

TrkHotListEmpty::TrkHotListEmpty(const TrkHotListEmpty& rhs)
{
  _nAxial    = rhs._nAxial;
  _nStereo   = rhs._nStereo;
  _nPhi      = rhs._nPhi;
  _nZ        = rhs._nZ;
  _stFndRng  = rhs._stFndRng;
  _endFndRng = rhs._endFndRng;
  _firstdch  = rhs._firstdch;
  _lastdch   = rhs._lastdch;
  _inactive  = rhs._inactive;
  for (unsigned i=0;i<5;++i) _svtpat[i] = rhs._svtpat[i];

}


TrkHotListEmpty&
TrkHotListEmpty::operator = (const TrkHotList& other) {
  if(this != &other){
    _nPhi = other.nSvt(TrkEnums::xyView);
    _nZ = other.nSvt(TrkEnums::zView);
    _nAxial = other.nDch(TrkEnums::xyView);
    _nStereo = other.nDch(TrkEnums::zView);
    _stFndRng = other.startFoundRange();
    _endFndRng = other.endFoundRange();
    _firstdch = other.firstDchLayer();
    _lastdch = other.lastDchLayer();
    for(unsigned ilay=0;ilay<5;ilay++)
      _svtpat[ilay] = other.svtView(ilay+1);// layer numbering starts at 1
    unsigned nhits = other.nHit();
    std::vector<unsigned> inactive;
    for(unsigned ihit=0;ihit<nhits;ihit++)
      if(!other.isActive(ihit))
	inactive.push_back(ihit);
    _inactive = inactive;
  }
  return *this;
}
  

TrkHotListEmpty::~TrkHotListEmpty()
{
}

TrkHotList*
TrkHotListEmpty::clone(TrkBase::Functors::cloneHot) const
{
  return new TrkHotListEmpty(*this);
}

TrkHotListEmpty::TrkHotListEmpty(const TrkHotList& other) :
  _nAxial(other.nDch(TrkEnums::xyView)),
  _nStereo(other.nDch(TrkEnums::zView)),
  _nPhi(other.nSvt(TrkEnums::xyView)),
  _nZ(other.nSvt(TrkEnums::zView)),
  _stFndRng(other.startFoundRange()),
  _endFndRng(other.endFoundRange()),
  _firstdch(other.firstDchLayer()),
  _lastdch(other.lastDchLayer())
{
  for(unsigned ilay=0;ilay<5;ilay++)
    _svtpat[ilay] = other.svtView(ilay+1);// layer numbering starts at 1
  unsigned nhit = other.nHit();
  for(unsigned ihit=0;ihit<nhit;ihit++)
    if(!isActive(ihit))
      _inactive.push_back(ihit);
}

int
TrkHotListEmpty::nActive(TrkEnums::TrkViewInfo view) const
{
  return nSvt(view)+nDch(view);
}

int
TrkHotListEmpty::nSvt(TrkEnums::TrkViewInfo view) const
{
  switch (view) {
  case TrkEnums::bothView:
    return _nPhi  + _nZ;
  case TrkEnums::xyView:
    return _nPhi;
  case TrkEnums::zView:
    return _nZ;
  default:
    return -1;
  }
}

int
TrkHotListEmpty::nDch(TrkEnums::TrkViewInfo view) const
{
  switch (view) {
  case TrkEnums::zView:
    return _nStereo;
  case TrkEnums::xyView:
    return _nAxial;
  case TrkEnums::bothView:
    return _nStereo + _nAxial;
  default:
    return -1;
  }
}

int
TrkHotListEmpty::nHit(TrkEnums::TrkViewInfo view) const
{
  unsigned nhit = nActive(view);
  if(view == TrkEnums::bothView)
    nhit += _inactive.size();
  else
    nhit += _inactive.size()/2;
  return nhit;
}

double
TrkHotListEmpty::startFoundRange() const
{
  return _stFndRng;
}

double
TrkHotListEmpty::endFoundRange() const
{
  return _endFndRng;
}


const std::vector<TrkHitOnTrk*>&
TrkHotListEmpty::hotlist() const
{
  static const std::vector<TrkHitOnTrk*> dummy;
  assert(dummy.empty());
  return dummy;

}

std::vector<TrkHitOnTrk*>&
TrkHotListEmpty::hotlist()
{
  static std::vector<TrkHitOnTrk*> dummy;
  assert(dummy.empty());
  return dummy;
}

void
TrkHotListEmpty::append(TrkHitOnTrk* )
{
}

void
TrkHotListEmpty::remove(TrkHitOnTrk* )
{
}

TrkHitOnTrk*
TrkHotListEmpty::findHot(const TrkFundHit*) const
{
  return 0;
}

bool
TrkHotListEmpty::hitCapable() const
{
  return false;
}

void
TrkHotListEmpty::updateHots()
{
  return;
}

TrkView
TrkHotListEmpty::svtView(int ilay) const
{
  return (ilay>=1 && ilay<=5)?_svtpat[ilay-1]:TrkView(TrkEnums::noView);
}

unsigned
TrkHotListEmpty::firstDchLayer() const
{
  return _firstdch;
}

unsigned
TrkHotListEmpty::lastDchLayer() const
{
  return _lastdch;
}

bool
TrkHotListEmpty::isActive(unsigned ihot) const {
  std::vector<unsigned>::const_iterator ifound =
    std::find(_inactive.begin(),_inactive.end(),ihot);
  return ifound == _inactive.end();
}
