//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkRecoTrk.hh,v 1.83 2004/08/06 06:31:43 bartoldu Exp $
//
// Description:
//      This is the standard reconstructed charged track class.  Only a few
//   functions, describing the track as a whole, are in this interface.
//   The remainder of the information can be obtained through one of the
//   four interfaces available here (TrkFit, TrkFitStatus, TrkHitList,
//   TrkExtInterface -- see comments below for details).
//   The interface you get will be for a particular mass hypothesis (if
//   leave out the particle-type argument, you will get the interface for
//   the default type for that track).  Some interfaces may not be available
//   for some tracks: some kinds of tracks don't have hit lists, for example,
//   and sometimes fits fail for some hypotheses.  In such cases, you get
//   a null pointer returned to you -- so test it.
//   All hypotheses should have the same kind of internal representation
//   A single fit may represent more than one particle
//   type; whichFit(hypo) tells you which fit is actually being used when
//   you ask about "hypo".
//
// Track creation:
//   Only FitMaker objects are allowed to create tracks, and only they can
//   change the TrkRep inside an existing track.  The RecoTrk ctor permits
//   tracks to be created without Reps, but FitMakers are required (by
//   fiat, not by syntax) to install a Rep in a track before finishing
//   with it.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Steve Schaffner
//------------------------------------------------------------------------
#ifndef TRKRECOTRK_HH
#define TRKRECOTRK_HH

#include <functional>
#include <map>
#include <set>
#include "TrkBase/TrkId.hh"
#include "BaBar/PdtPid.hh"
#include "TrkBase/TrkDirection.hh"
#include "TrkBase/TrkHitList.hh"
#include "TrkBase/TrkStoreHypo.hh"
#include "TrkBase/TrkFit.hh"
#include "TrkBase/TrkFitStatus.hh"

class TrkRep;
class TrkRepIter;
class TrkContext;
class TrkFundHit;
#include <iosfwd>
class TrkExchangePar;
class TrkVolume;
class TrkFit;
class TrkFitStatus;
class TrkErrCode;
class BField;
class TrkRecoTrkImpl;

class TrkRecoTrk  {
public:
  typedef std::unary_function<TrkRecoTrk,bool> predicate_type;
  //*********************************
  //Global track quantities:
  //*********************************
  const TrkId&      id()                            const;
  PdtPid::PidType   defaultType()                   const {return _defaultType;}
  PdtPid::PidType   whichFit(PdtPid::PidType hypo)  const;
  int               fitNumber(PdtPid::PidType hypo) const;
  double            trackT0()                       const { return _trackT0; }
  double            trackT0err()		    const { return _trackT0err; }
  const BField&     bField()                        const {return *_bField;}

  TrkErrCode        addFit(PdtPid::PidType hypo,bool fit=true);  // also fits if requested
  // Note: resetT0() requires refit to make fit current
  void              resetT0(double t0, double t0err);         // also updates hits

  //**********************************************************
  // To get information about the track as fitted to a particular mass
  //   hypothesis, use one of the following interfaces.  In each case,
  //   you can either specify a hypothesis, or omit the argument and get the
  //   default hypothesis for this track.
  //**********************************************************

  //**********************************************************
  // (1) Standard information about the fitted track (momentum, position,
  //   charge, chisq etc).
  //**********************************************************
  const TrkFit* fitResult() const;
  const TrkFit* fitResult(PdtPid::PidType hypo) const;


  //**********************************************************
  // (2) Interface for accessing and manipulating the track's list of hits,
  //  and for fitting the track.
  //**********************************************************
        TrkHitList*   hits()                      {return hits(defaultType());}
  const TrkHitList*   hits() const                {return hits(defaultType());}
        TrkHitList*   hits(PdtPid::PidType hypo);
  const TrkHitList*   hits(PdtPid::PidType hypo) const;

// same for hots.  This is more direct
        TrkHotList*   hots()                      {return hots(defaultType());}
  const TrkHotList*   hots() const                {return hots(defaultType());}
        TrkHotList*   hots(PdtPid::PidType hypo);
  const TrkHotList*   hots(PdtPid::PidType hypo) const;



  //**********************************************************
  // (3) Specialized information about the fit; of interest mostly to experts
  //**********************************************************
  const TrkFitStatus* status() const;
  const TrkFitStatus* status(PdtPid::PidType hypo) const;
  TrkFitStatus* status();
  TrkFitStatus* status(PdtPid::PidType hypo);


  //**************************************************
  // Constructors and such (normal ctor is protected)
  //**************************************************
  // Copy constructor (leaves original unchanged):
  TrkRecoTrk(const TrkRecoTrk& right);
  // Destructor
  virtual ~TrkRecoTrk();
  const TrkRecoTrk& operator=(const TrkRecoTrk& right);
  bool operator==(const TrkRecoTrk &other) const;
  bool operator<(const TrkRecoTrk &other) const;

  //**********************************************************
  // Printing
  //**********************************************************
  virtual void        print(std::ostream& ) const;
  virtual void        printAll(std::ostream& ) const;
  //**********************************************************
  // Persistence
  //**********************************************************
// Mark a particular hypo at a particular flight length for storage.
// This only works for the mini.  Several 'lists' of storage requests
// may be associated with each track
  void markForStore(PdtPid::PidType hypo,double fltlen,const char* listname="Default");
  const std::set<TrkStoreHypo>& storageRequests(const char* listname="Default") const;
// clear out all marked stores
  void clearStorageRequests(const char* listname="Default");
// return the set of fit storage lists known to this track
  void storageLists(std::set<std::string>& storage) const;
private:
  //*** Data members ***
  TrkRecoTrkImpl* _impl; // the reps live here; owned by trk; 
                         // the reason they're stashed away in this class
                         // is because (for #($*)( reasons) ooddlx must
                         // parse code which needs to know that TrkRecoTrk
                         // inherits from AbsEvtObj -- but ooddlx cannot
                         // deal with somewhat complex ANSI C++ constructions
                         // (e.g. namespaces). Stashing the reps into 
                         // TrkRecoTrkImpl, which can be fwd declared as we
                         // only use a pointer to it insures that ooddlx 
                         // doesn't get 'confused' by real ANSI C++...
  TrkId _id;                                // unique id # in event
  std::vector<int> _fitNumber;   //number of times fit has been altered
// keep track of the storage requests, by list.  This is sorted first by hypo, then
// by flightlength.
  std::map<std::string,std::set<TrkStoreHypo> > _storage;
  PdtPid::PidType _defaultType;
  double _trackT0;
  double _trackT0err;
  const BField* _bField;
protected:
public:
  TrkRep* getRep(PdtPid::PidType hypo);
  const TrkRep* getRep(PdtPid::PidType hypo) const;
  // protected functio
  void copyReps(const TrkRecoTrk& rhs);
  // The following takes ownership of the argument; it replaces the default Rep
  //    and zeroes others.
  void setRep(TrkRep*);
  // Make hypothesis <hypo> use fit currently used for hypo <fit>:
  void repointHypo(PdtPid::PidType hypo, PdtPid::PidType fit);
  void changeDefault(PdtPid::PidType newHypo);

  // return the list of unique distinct Reps attached to this track
  std::pair<TrkRepIter,TrkRepIter> uniqueReps() const;
  // return the list of the 5 reps this track is pointing at
  std::pair<TrkRepIter,TrkRepIter> allReps() const;
  void setFitNumber(PdtPid::PidType hypo, int newNumber);
  void updateReps();
  // addHypoTo takes ownership of newRep!
  void addHypoTo(TrkRep* newRep, PdtPid::PidType hypo);
// a couple of lame functions to limp past inconsistencies between the persistent
// design and the tracking design.  Ugh
  void setBField(const BField* field);
// Constructors are protected (construct through FitMaker)
  TrkRecoTrk(PdtPid::PidType defaultPart, const TrkContext&, double t0);
// persistence constructor.  BField must be set later
  TrkRecoTrk(PdtPid::PidType defaultPart,long idnum,double t0);

  // Access to TrkRep, for testing only; use it at your peril
  const TrkRep* testRep( PdtPid::PidType hypo ) const { return getRep(hypo);}
  friend class TrkFitMaker;
  friend class TrkHitOnTrk;
  friend class TrkHitList;
};

std::ostream& operator<<(std::ostream& os, const TrkRecoTrk& tk);

#endif
