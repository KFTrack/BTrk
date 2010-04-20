// File and Version Information:
//      $Id: TrkRecoTrk.cc,v 1.109 2004/10/12 15:06:06 raven Exp $
//
// Description:
//      implementation for TrkRecoTrk
//
//
// Author List: Stephen Schaffner
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "BaBar/Constants.hh"
#include <assert.h>
#include <functional>
#include <algorithm>
#include "CLHEP/Alist/AIterator.h"
#include "TrkBase/TrkRecoTrk.hh"
#include "TrkBase/TrkHitOnTrk.hh"
#include "TrkBase/TrkFundHit.hh"
#include "TrkBase/TrkDifTraj.hh"
#include "TrkBase/TrkExchangePar.hh"
#include "TrkBase/TrkErrCode.hh"
#include "TrkBase/TrkExtInterface.hh"
#include "TrkBase/TrkFitStatus.hh"
#include "TrkBase/TrkRepIter.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "PDT/Pdt.hh"
#include "BbrGeom/BbrPointErr.hh"
#include "BbrGeom/BbrVectorErr.hh"
#include "difAlgebra/DifPoint.hh"
#include "difAlgebra/DifVector.hh"
#include "difAlgebra/DifIndepPar.hh"
#include "TrkBase/TrkRep.hh"
#include "TrkBase/TrkContext.hh"
#include "ErrLogger/ErrLog.hh"
#include "BaBar/PdtPid.hh"
#include "boost/shared_ptr.hpp"
#include <vector>
using std::endl;
using std::ostream;


class TrkRecoTrkImpl {
public:
  friend class TrkRecoTrk;
private:
  TrkRecoTrkImpl()
     : _reps(PdtPid::nPidType,boost::shared_ptr<TrkRep>((TrkRep *)0)),
       _hitInterfaces(PdtPid::nPidType,boost::shared_ptr<TrkHitList>((TrkHitList *)0))
  { }

  // the reps live here; owned by trk;
  // Pointers to the Reps, indexed by particle type; always has nPart entries,
  // which may point to as few as one single distinct Rep...
  typedef std::vector<boost::shared_ptr<TrkRep> > repList_t;
  typedef repList_t::const_iterator repConstIter;
  typedef repList_t::iterator repIter;
  repList_t  _reps;
  typedef std::vector<boost::shared_ptr<TrkHitList> > hitList_t;
  typedef hitList_t::iterator hitListIter;
  hitList_t _hitInterfaces;

};

TrkRecoTrk::TrkRecoTrk(PdtPid::PidType defaultPart, const TrkContext& ctext,
                       double t0) :
  _impl(new TrkRecoTrkImpl), 
  _id(ctext.getId()),
  _fitNumber(PdtPid::nPidType,(int)0),
  _defaultType(defaultPart),
  _trackT0(t0),
  _bField( ctext.bField() )
{
  // No TrkRep is defined here; must be created in appropriate FitMaker.
  unsigned i=0;
  for (TrkRecoTrkImpl::hitListIter iface = _impl->_hitInterfaces.begin();
       iface!= _impl->_hitInterfaces.end(); ++iface) {
    iface->reset( new TrkHitList(this, (PdtPid::PidType)i++) );  //cast
  }
}

// persistence reconstitution.  This sets a nul value for BField and IdManager
TrkRecoTrk::TrkRecoTrk(PdtPid::PidType defaultPart,long idnum,double t0) :
  _impl(new TrkRecoTrkImpl),
  _id(idnum,0),
  _fitNumber(PdtPid::nPidType,(int)0),
  _defaultType(defaultPart),
  _trackT0(t0),
  _bField(0)
{
  // No TrkRep is defined here; must be created in appropriate FitMaker.
  unsigned i=0;
  for (TrkRecoTrkImpl::hitListIter iface = _impl->_hitInterfaces.begin();
       iface!= _impl->_hitInterfaces.end(); ++iface) {
    iface->reset( new TrkHitList(this, (PdtPid::PidType)i++) );  //cast
  }
}

//-- Copy constructor
TrkRecoTrk::TrkRecoTrk(const TrkRecoTrk& rhs) :
  _impl(new TrkRecoTrkImpl),
  _id(rhs._id.idManager()),
  _fitNumber(PdtPid::nPidType,(int)0),
  _storage(rhs._storage),
  _trackT0( rhs._trackT0),
  _bField(rhs._bField)
{
  _defaultType = rhs.defaultType();
  unsigned i=0;
  for (TrkRecoTrkImpl::hitListIter iface = _impl->_hitInterfaces.begin();
       iface!= _impl->_hitInterfaces.end(); ++iface) {
    iface->reset( new TrkHitList(this, (PdtPid::PidType)i++) );  //cast
  }
  copyReps(rhs);
}

TrkRecoTrk::~TrkRecoTrk()
{
        delete _impl;
}

const TrkRecoTrk&
TrkRecoTrk::operator=(const TrkRecoTrk &right)
{
  if (&right == this) return *this;
  _trackT0 = right._trackT0;
  _defaultType = right.defaultType();
  copyReps(right);
  _bField = right._bField;
  AbsEvtObj::operator=(right);
  _id.setNewValue(right._id);
  _storage = right._storage;
  return *this;
}

const TrkId&
TrkRecoTrk::id() const
{
  return _id;
}

double
TrkRecoTrk::trackT0() const
{
  return _trackT0;
}

PdtPid::PidType
TrkRecoTrk::whichFit(PdtPid::PidType hypo) const
{
  if (_impl->_reps[hypo].get() == 0) {
    return hypo;
  }
//  return _repPtr[hypo]->fitValid() ? _repPtr[hypo]->particleType()
//                                   : PdtPid::null;
  return _impl->_reps[hypo]->particleType();
}

int
TrkRecoTrk::fitNumber(PdtPid::PidType hypo) const
{
  PdtPid::PidType used = whichFit(hypo);
  if (used == PdtPid::null) return -1;
  int index = used;
  return _fitNumber[index];
}

void
TrkRecoTrk::print(ostream& ostr) const
{
  ostr << "Trk: " << id() << " def: "
       << Pdt::lookup(defaultType())->name()
       << " fitNumber:" << fitNumber(defaultType());
  const TrkFit* daFit = fitResult();
  const TrkHitList* daList = hits();
  const TrkFitStatus* daStatus = status();
  if (daList != 0) {
    ostr << " nhit: " << daList->nHit();
  }
  if (daFit != 0) {
    ostr << " nactive: " << daFit->nActive() << " chisq: " << daFit->chisq();
  }
  if (daStatus != 0) {
    ostr << " 3-d: " << (daStatus->is2d() == 0);
  }
  ostr << " t0: " << _trackT0 << "\n";
  if (daFit != 0) {
    TrkExchangePar par = daFit->helix(0.);
    ostr << "phi0: " << par.phi0()
       << " om: "<< par.omega()
       << " d0: " << par.d0() << " z0: " << par.z0()
       << " ct: " << par.tanDip();
  }
  ostr << endl;
}

void
TrkRecoTrk::printAll(ostream& ostr) const
{
  //  This should be expanded to print other hypotheses as well
  ostr << "Trk: " << id() << " def: "
       << Pdt::lookup(defaultType())->name()
       << " fitNumber:" << fitNumber(defaultType());
  const TrkFit* daFit = fitResult();
  const TrkHitList* daList = hits();
  const TrkFitStatus* daStatus = status();
  if (daList != 0) {
    ostr << " nhit: " << daList->nHit();
  }
  if (daFit != 0) {
    ostr << " nactive: " << daFit->nActive() << " chisq: " << daFit->chisq();
  }
  if (daStatus != 0) {
    ostr << " 3-d: " << (daStatus->is2d() == 0);
  }
  ostr << " t0: " << _trackT0 << "\n";
  getRep(defaultType())->printAll(ostr);
  ostr << endl;
}

TrkErrCode
TrkRecoTrk::addFit(PdtPid::PidType hypo,bool fit)
{
  // If there is no fit, create one.  If hypo points to a fit for a different
  //   particle type, create a fit of type hypo, and point at that.  Carry
  //   out the fit if needed.
  if (hits() == 0) {
    // Unfittable rep
    return TrkErrCode(TrkErrCode::fail, 11,
                      "TrkRecoTrk::addFit(): cannot add a fit to this track.");
  }
  if (whichFit(hypo) == hypo) {
    return TrkErrCode(TrkErrCode::succeed, 11,
                      "TrkRecoTrk::addFit(): requested fit already exists.");
  }
  _impl->_reps[hypo].reset( _impl->_reps[defaultType()]->cloneNewHypo(hypo) );
  TrkErrCode fitErr(TrkErrCode::succeed, 1);
  if (fit && !_impl->_reps[hypo]->fitCurrent()) {
    fitErr = _impl->_reps[hypo]->fit();
  }
  ++_fitNumber[hypo];
  return fitErr;
}

void
TrkRecoTrk::resetT0(double t)
{
  _trackT0 = t;
  updateReps();
}

void
TrkRecoTrk::updateReps()
{
  std::pair<TrkRepIter,TrkRepIter> x = uniqueReps();
  std::for_each(x.first,x.second,std::mem_fun_ref(&TrkRep::updateHots));
  std::transform(_fitNumber.begin(),_fitNumber.end(),
                 _fitNumber.begin(),
                 std::bind2nd(std::plus<int>(),1));
}

bool
TrkRecoTrk::operator==(const TrkRecoTrk &other) const
{
  return _id == other._id;
}

bool
TrkRecoTrk::operator<(const TrkRecoTrk &other) const
{
  return _id < other._id;
}


//***********************
// Protected functions
//***********************

TrkRep*
TrkRecoTrk::getRep(PdtPid::PidType hypo)
{
  assert(hypo>=PdtPid::electron && hypo<= PdtPid::proton);
  TrkRep* theRep = _impl->_reps[hypo].get();
// insist the default rep exist
  if(hypo == defaultType())assert(0 != theRep);
  return theRep;
}

const TrkRep*
TrkRecoTrk::getRep(PdtPid::PidType hypo) const
{
  assert(hypo>=PdtPid::electron && hypo<= PdtPid::proton);
  const TrkRep* theRep = _impl->_reps[hypo].get();
  if(hypo == defaultType())assert(0 != theRep);
  return theRep;
}

void
TrkRecoTrk::changeDefault(PdtPid::PidType newHypo)
{
  if (newHypo == defaultType()) return;
  assert(whichFit(newHypo) != PdtPid::null);

  TrkHotList *oldList = getRep(defaultType())->hotList();
  std::for_each(oldList->begin(),
                oldList->end(),
                std::mem_fun_ref(&TrkHitOnTrk::setUnusedHit));
  assert(getRep(newHypo) != 0);
  TrkHotList *newList= getRep(newHypo)->hotList();
  std::for_each(newList->begin(),
                newList->end(),
                std::mem_fun_ref(&TrkHitOnTrk::setUsedHit));
  _defaultType = newHypo;
}

void
TrkRecoTrk::copyReps(const TrkRecoTrk& rhs)
{
  TrkRecoTrkImpl::repIter lhs = _impl->_reps.begin();
  for (TrkRecoTrkImpl::repIter i=rhs._impl->_reps.begin();i!=rhs._impl->_reps.end();++i,++lhs) {
      TrkRecoTrkImpl::repIter j=std::find(rhs._impl->_reps.begin(),i,*i);
      if (j==i) { // first time this one is seen
          lhs->reset( (*i)->clone(this) );
          (*lhs)->setValid((*i)->fitValid());
          (*lhs)->setCurrent((*i)->fitCurrent());
      } else {
          *lhs = *(_impl->_reps.begin()+(j-rhs._impl->_reps.begin()));
      }
  }
  assert(_fitNumber.size()==rhs._fitNumber.size());
  std::copy(rhs._fitNumber.begin(),rhs._fitNumber.end(),_fitNumber.begin());
}

void
TrkRecoTrk::setRep(TrkRep* r)
{
  // Sets the default rep to be r, clears out other reps., and sets
  //  non-default rep ptrs to point at default.  Increments all fit numbers.
  std::fill(_impl->_reps.begin(),_impl->_reps.end(),boost::shared_ptr<TrkRep>(r));
  std::transform(_fitNumber.begin(),_fitNumber.end(),
                 _fitNumber.begin(),
                 std::bind2nd(std::plus<int>(),1));
}

void
TrkRecoTrk::repointHypo(PdtPid::PidType hypo, PdtPid::PidType fit)
{
  // Do we have to do anything?
  if (fit == hypo || getRep(fit) == getRep(hypo)) return;

  if (hypo == defaultType()) {
    ErrMsg(error) <<
      "TrkRecoTrk: can't make default hypothesis point at different fit"
                  << endmsg;
    return;
  }
  _impl->_reps[hypo] = _impl->_reps[fit];
}

bool
TrkRecoTrk::attach(TrkExtInterface& interface, PdtPid::PidType hypo) const
{
  const TrkRep* rp = getRep(hypo);
  return rp!=0?interface.attach(rp):0;
}


bool
TrkRecoTrk::attach(TrkExtInterface& interface, PdtPid::PidType hypo)
{
  TrkRep* rp = getRep(hypo);
  return rp!=0?interface.attach(rp):0;
}

TrkHitList*
TrkRecoTrk::hits(PdtPid::PidType hypo)
{
  const TrkRep* rp = getRep(hypo);
  return rp==0?0:_impl->_hitInterfaces[hypo].get();
}

const TrkHitList*
TrkRecoTrk::hits(PdtPid::PidType hypo) const
{
  const TrkRep* rp = getRep(hypo);
  return rp==0?0:_impl->_hitInterfaces[hypo].get();
}

const TrkFit*
TrkRecoTrk::fitResult() const
{
  return fitResult(defaultType());
}

const TrkFit*
TrkRecoTrk::fitResult(PdtPid::PidType hypo) const
{
  const TrkRep* rp = getRep(hypo);
  return rp==0?0:(rp->fitValid()?rp:0);
}

const TrkFitStatus*
TrkRecoTrk::status() const
{
  return status(defaultType());
}

const TrkFitStatus*
TrkRecoTrk::status(PdtPid::PidType hypo) const
{
  return getRep(hypo);
}

TrkFitStatus*
TrkRecoTrk::status()
{
  return status(defaultType());
}

TrkFitStatus*
TrkRecoTrk::status(PdtPid::PidType hypo)
{
  return getRep(hypo);
}

void
TrkRecoTrk::setFitNumber(PdtPid::PidType hypo, int newNumber)
{
  _fitNumber[hypo] = newNumber;
}

void
TrkRecoTrk::addHypoTo(TrkRep* newRep, PdtPid::PidType hypo)
{
  _impl->_reps[hypo].reset( newRep );
}

void
TrkRecoTrk::setIdManager(TrkIdManager* idMan)
{
  _id.setIdManager(idMan);
}

ostream&
operator<<(ostream& os, const TrkRecoTrk& tk)
{
  tk.print(os);
  return os;
}

void
TrkRecoTrk::setBField(const BField* field)
{
  _bField = field;
}

void
TrkRecoTrk::markForStore(PdtPid::PidType hypo,double fltlen,
			 const char* listname) {
// first, translate to the real hypo
  PdtPid::PidType realhypo = whichFit(hypo);
// Then, make sure this hypo has a valid fit
  if(getRep(realhypo)!= 0 && getRep(realhypo)->fitValid())
// add an entry in the toStore list (if it's unique)
    _storage[std::string(listname)].insert(TrkStoreHypo(realhypo,fltlen));
  else
// It's an error to try to store invalid fits
    ErrMsg(error) << "Invalid fits cannot be marked for storage" << endmsg;
}

const std::set<TrkStoreHypo>&
TrkRecoTrk::storageRequests(const char* listname) const {
  static std::set<TrkStoreHypo> empty; // empty set to return if list doesn't exist
  std::map<std::string,std::set<TrkStoreHypo> >::const_iterator foundit =
    _storage.find(std::string(listname));
  if(foundit != _storage.end())
    return foundit->second;
  else
    return empty;
}

void
TrkRecoTrk::clearStorageRequests(const char* listname) {
  _storage[std::string(listname)].clear();
}

void
TrkRecoTrk::storageLists(std::set<std::string>& storage) const {
// clear the output set
  storage.clear();
// iterate over all the storage requests
  std::map<std::string,std::set<TrkStoreHypo> >::const_iterator miter = _storage.begin();
  while(miter != _storage.end()){
    storage.insert(miter->first);
    miter++;
  }
}

TrkHotList*
TrkRecoTrk::hots(PdtPid::PidType hypo)
{
  TrkRep* rp = getRep(hypo);
  return rp==0?0:rp->hotList();
}

const TrkHotList*
TrkRecoTrk::hots(PdtPid::PidType hypo) const
{
  const TrkRep* rp = getRep(hypo);
  return rp==0?0:rp->hotList();
}

std::pair<TrkRepIter,TrkRepIter>
TrkRecoTrk::allReps() const
{
  typedef std::vector<TrkRep*> RPL;
  boost::shared_ptr<RPL> x( new RPL );
  for (TrkRecoTrkImpl::repConstIter i=_impl->_reps.begin();i!=_impl->_reps.end();++i) {
          x->push_back(i->get());
  }
  return std::make_pair(TrkRepIter(x,0),TrkRepIter(x,x->size()));
}

std::pair<TrkRepIter,TrkRepIter>
TrkRecoTrk::uniqueReps() const
{
  typedef std::vector<TrkRep*> RPL;
  boost::shared_ptr<RPL> x( new RPL );
  for (TrkRecoTrkImpl::repConstIter i=_impl->_reps.begin();i!=_impl->_reps.end();++i) {
     if (std::find(x->begin(),x->end(),i->get())==x->end()) x->push_back(i->get());
  }
  return std::make_pair(TrkRepIter(x,0),TrkRepIter(x,x->size()));
}
