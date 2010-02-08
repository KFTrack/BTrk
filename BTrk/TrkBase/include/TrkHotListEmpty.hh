//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkHotListEmpty.hh,v 1.16 2003/10/25 21:58:48 brownd Exp $
//
// Description: List of hits (as HitOnTrk objects) associated with a 
//  reconstructed track.  Abstract base class.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner
//
//------------------------------------------------------------------------
#ifndef TRKHOTLISTEMPTY_HH
#define TRKHOTLISTEMPTY_HH

#include "TrkBase/TrkHotList.hh"
#include "TrkBase/TrkView.hh"

class TrkFundHit;
class TrkHitOnTrk;
class TrkRep;

// Class interface //
class TrkHotListEmpty : public TrkHotList {
public:
// constructors and such.  The following assumes a uniform distribution
// of hots between views, and default values for start/end layer.  This
// constructor is DEPRECATED, please don't use it in new code
  TrkHotListEmpty(int nActive, int nSvt, int nDch,
                  double startFoundRange, double endFoundRange);
// full constructor; specify hots by type, and include in
  TrkHotListEmpty(unsigned nPhi,unsigned nZ, unsigned nAxial, unsigned nStereo,
                  double startFoundRange, double endFoundRange,
                  unsigned firstdchlay,unsigned lastdchlay,
                  TrkView svtpattern[5],
		  const std::vector<unsigned>& inactive);
// copy constructor; this is the prefered way to create an empty Hot list
  TrkHotListEmpty(const TrkHotList& other);
// equivalence is OK
  TrkHotListEmpty&   operator= (const TrkHotList&);  

  virtual TrkHotList* clone(TrkBase::Functors::cloneHot) const;
  virtual ~TrkHotListEmpty();

  virtual bool      hitCapable()      const;
  virtual int        nActive(TrkEnums::TrkViewInfo view=TrkEnums::bothView) const;
  virtual int        nDch(TrkEnums::TrkViewInfo view=TrkEnums::bothView) const;
  virtual int        nSvt(TrkEnums::TrkViewInfo view=TrkEnums::bothView) const;
  virtual int        nHit(TrkEnums::TrkViewInfo view=TrkEnums::bothView) const;
  virtual TrkView svtView(int layer) const;
  virtual unsigned  firstDchLayer() const;
  virtual unsigned  lastDchLayer() const;
  virtual double    startFoundRange() const;
  virtual double    endFoundRange()   const;
  virtual bool      isActive(unsigned ihot) const;
  virtual void append(TrkHitOnTrk* );
  virtual void remove(TrkHitOnTrk* );
  TrkHitOnTrk* findHot(const TrkFundHit*) const;
  virtual void updateHots();

protected:
  virtual const std::vector<TrkHitOnTrk*>&   hotlist()   const;
  virtual       std::vector<TrkHitOnTrk*>&   hotlist();
private:
  unsigned _nAxial,_nStereo;
  unsigned _nPhi,_nZ;
  double _stFndRng;
  double _endFndRng;
  unsigned _firstdch;
  unsigned _lastdch;
  TrkView _svtpat[5];
  std::vector<unsigned> _inactive; // indices of inactive hots
  TrkHotListEmpty(const TrkHotListEmpty& rhs);  //copy ctor
// allow persistent to set inactive
  void setInactive(std::vector<unsigned>& inactive ) {
    _inactive = inactive; }
  friend class KalMiniTrkK;
};

#endif
