//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkHotListData.cc,v 1.10 2004/07/05 05:43:27 hulsberg Exp $
//
//  Description:
//  Class TrkHotListData
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Infomation;
//	Copyright (C) 2002	Lawrence Berkeley Laboratory
//
// Author(s): Dave Brown 10/202/02
//
//------------------------------------------------------------------------


#include "BaBar/BaBar.hh"
#include "TrkBase/TrkHotListData.hh"
#include "CommonUtils/ComPackFlatFloat.hh"
#include "CommonUtils/ComPackInt.hh"
#include "CommonUtils/ComPackSignedExpFloat.hh"
#include "CommonUtils/ComPackBase.hh"
#include "TrkBase/TrkHotListEmpty.hh"
#include "ErrLogger/ErrLog.hh"
#include <assert.h>
// count the # of layers
const unsigned TrkHotListData::_nsvt(5);
const unsigned TrkHotListData::_ndch(40);
// svtdata shifts and masks
const unsigned TrkHotListData::_vsft(2);
const unsigned TrkHotListData::_vmsk(0x3);
const unsigned TrkHotListData::_nphisft(_nsvt*_vsft);
const unsigned TrkHotListData::_nphimsk(0x7);
const unsigned TrkHotListData::_nzsft(3+_nphisft);
const unsigned TrkHotListData::_nzmsk(0x7);
// dchdata shifts and masks
const unsigned TrkHotListData::_naxsft(0);
const unsigned TrkHotListData::_naxmsk(0x3f);
const unsigned TrkHotListData::_nstsft(6+_naxsft);
const unsigned TrkHotListData::_nstmsk(0x7f);
// 3 bits at the end of the dch word are currently unused
// flight range shifts and masks
const unsigned TrkHotListData::_startsft(0);
const unsigned TrkHotListData::_startmsk(0x3ff);
const unsigned TrkHotListData::_rangesft(10+_startsft);
const unsigned TrkHotListData::_rangemsk(0x3ff);
const unsigned TrkHotListData::_flaysft(10+_rangesft);
const unsigned TrkHotListData::_flaymsk(0x3f);
const unsigned TrkHotListData::_llaysft(6+_flaysft);
const unsigned TrkHotListData::_llaymsk(0x3f);
// packers
const ComPackSignedExpFloat 
TrkHotListData::_packfndstart(10,3,256.0);
const ComPackFlatFloat 
TrkHotListData::_packfndrange(0.0,256.0,10);

TrkHotListData::TrkHotListData() : _svtdata(0),_dchdata(0),_fndrng(0)
{}

TrkHotListData::TrkHotListData(d_UShort svtdata,d_UShort dchdata,d_ULong rangedata) :
  _svtdata(svtdata),_dchdata(dchdata),_fndrng(rangedata)
{}

TrkHotListData::TrkHotListData(const TrkHotList& hotlist) :
  _svtdata(0),_dchdata(0),_fndrng(0)
{
// count the # of hits by view
  int nphi = hotlist.nSvt(TrkEnums::xyView);
  int nz = hotlist.nSvt(TrkEnums::zView);
// construct the pattern, and count the minimum # of hits implied by that
  int mphi(0);
  int mz(0);
  for(unsigned ilay=0;ilay<_nsvt;ilay++){
    TrkView svtv = hotlist.svtView(ilay+1); // layer numbering starts at 1
    _svtdata |= (svtv.viewData() << (ilay*_vsft));
    if(svtv.contains(TrkEnums::xyView))mphi++;
    if(svtv.contains(TrkEnums::zView))mz++;
  }
// store the _difference_ between the # of hits/view and the minimum # implied
// by the pattern.
  nphi -= mphi;
  nz -= mz;
  if(nphi<0 || nz<0){
    ErrMsg(error) << "Inconsistent TrkHotList Svt information" << endmsg;
    if(nphi<0)nphi=0;
    if(nz<0)nz=0;
  }
// truncate if necessary
  if(nphi>_nphimsk)nphi=_nphimsk;
  if(nz>_nzmsk)nz=_nzmsk;
  _svtdata |= (nphi << _nphisft) | (nz << _nzsft);
// pack # of dch hits
  int nax = hotlist.nDch(TrkEnums::xyView);
  int nst = hotlist.nDch(TrkEnums::zView);
  if(nax<0)nax=0;
  if(nst<0)nst=0;
// truncate if necessary
  if(nax > _naxmsk) nax = _naxmsk;
  if(nst > _nstmsk) nst = _nstmsk;
  _dchdata = (nax << _naxsft) | (nst << _nstsft);
// pack the found range as start + range
  unsigned fndstartpack,fndrangepack;
  double fndstart = hotlist.startFoundRange();
  double fndrange = hotlist.endFoundRange()-fndstart;
  ComPackBaseBase::StatusCode stacode = _packfndstart.pack(fndstart,fndstartpack);
  ComPackBaseBase::StatusCode rngcode = _packfndrange.pack(fndrange,fndrangepack);
  assert(stacode != ComPackBaseBase::TAG_BAD &&
	 rngcode != ComPackBaseBase::TAG_BAD);
// put these into the word
  _fndrng = (fndstartpack << _startsft) | (fndrangepack << _rangesft);
// add the inner and outer Dch layer
  unsigned fdch = hotlist.firstDchLayer();
  unsigned ldch = hotlist.lastDchLayer();
  unsigned firstdch = (fdch>0&&fdch<=_ndch) ? 
    ((unsigned)fdch & _flaymsk) : 0;
  unsigned lastdch = (ldch>0&&ldch<=_ndch) ? 
    ((unsigned)ldch & _llaymsk) : 0;
  _fndrng |=  (firstdch << _flaysft) | (lastdch << _llaysft);
}


TrkHotListData::TrkHotListData(const TrkHotListData& other) :
  _svtdata(other._svtdata),_dchdata(other._dchdata),_fndrng(other._fndrng)
{}


TrkHotListData&
TrkHotListData::operator =(const TrkHotListData& other ) {
  if(this != &other){
    _svtdata = other._svtdata;
    _dchdata = other._dchdata;
    _fndrng = other._fndrng;
  }
  return *this;
}

TrkHotListData::~TrkHotListData()
{}


// allow separate unpacking of all parameters;
// the hotlist function _RETURNS OWNERSHIP_
TrkHotListEmpty*
TrkHotListData::hotList() const
{
// unpack the svt information
  TrkView viewpat[5];
  for(unsigned ilay=0;ilay<_nsvt;ilay++){
    viewpat[ilay] = svtView(ilay);
  }
// unpack the found range as start + range
  double fndstart = foundStart();
  double fndrange = foundRange();
  double fndend = fndstart+fndrange;
// empty unassigned vector
  std::vector<unsigned> unassigned;
  TrkHotListEmpty* hotlist = 
    new TrkHotListEmpty(nPhi(),nZ(),nAxial(),nStereo(),
			fndstart,fndend,
			firstDch(),lastDch(),
			viewpat,unassigned);
  return hotlist ;
}

unsigned
TrkHotListData::nPhi() const {
  int mphi(0);
  for(unsigned ilay=0;ilay<_nsvt;ilay++)
    if(svtView(ilay).contains(TrkEnums::xyView))mphi++;
  return ((_svtdata >> _nphisft) & _nphimsk) + mphi;
}

unsigned 
TrkHotListData::nZ() const {
  int mz(0);
  for(unsigned ilay=0;ilay<_nsvt;ilay++)
    if(svtView(ilay).contains(TrkEnums::zView))mz++;
  return ((_svtdata >> _nzsft) & _nzmsk) + mz;
}

unsigned
TrkHotListData::nAxial() const {
  return ((_dchdata >> _naxsft) & _naxmsk);
}

unsigned
TrkHotListData::nStereo() const {
  return ((_dchdata >> _nstsft) & _nstmsk);
}

unsigned
TrkHotListData::firstDch() const {
  return ((_fndrng >> _flaysft) & _flaymsk);
}


unsigned
TrkHotListData::lastDch() const {
  return ((_fndrng >> _llaysft) & _llaymsk);
}

TrkView
TrkHotListData::svtView(int ilay) const {
  unsigned short svtd =  ( (_svtdata >> (ilay*_vsft)) & _vmsk);
  return TrkView(svtd);
}

double
TrkHotListData::foundStart() const {
  unsigned fndstartpack = ((_fndrng >> _startsft) & _startmsk);
  double fndstart;
  ComPackBaseBase::StatusCode stacode = _packfndstart.unpack(fndstartpack,fndstart);
  assert(stacode != ComPackBaseBase::TAG_BAD);
  return fndstart;
}

double
TrkHotListData::foundRange() const {
  unsigned fndrangepack = ((_fndrng >> _rangesft) & _rangemsk);
  double fndrange;
  ComPackBaseBase::StatusCode rngcode = _packfndrange.unpack(fndrangepack,fndrange);
  assert(rngcode != ComPackBaseBase::TAG_BAD);
  return fndrange;
}
