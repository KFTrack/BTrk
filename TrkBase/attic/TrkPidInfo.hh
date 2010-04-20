//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkPidInfo.hh 62 2009-03-27 08:58:45Z stroili $
//
// Description:
//     
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner
//
//------------------------------------------------------------------------
#ifndef TRKPIDINFO_HH
#define TRKPIDINFO_HH

#include "AbsPid/AbsPidInfo.hh"

class TrkRecoTrk;
#include <iosfwd>

// Class interface //
class TrkPidInfo : public AbsPidInfo {

public:
  TrkPidInfo(const TrkRecoTrk*);  // keeps ptr to trk
  virtual ~TrkPidInfo();

  virtual TrkPidInfo* clone()     const;
  const   TrkRecoTrk* recoTrack() const;
  void    print(std::ostream&)         const;
  
protected:	
  
private:	
  const TrkRecoTrk* _trk;

  TrkPidInfo&   operator= (const TrkPidInfo&);
  TrkPidInfo(const TrkPidInfo &);

  friend class TrkMakePid;
};

#endif







