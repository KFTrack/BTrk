//-----------------------------------------------------------------------------
// File and Version Information
//      $Id: MatElmDictionary.hh 516 2010-01-15 08:22:00Z stroili $
//
// Description:
//      Class MatElmDictionary (Element Dictionary)
//      Header file
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Mossadek Talby  (SLAC - CPPM/IN2P3 University of Aix-Marseille II)
//
// Modification History:
//   October 15, 1998 - Talby : created
//-----------------------------------------------------------------------------

#ifndef MATELMDICTIONARY_HH
#define MATELMDICTIONARY_HH

//----------------------
// Base Class Headers --
//----------------------
#include <string>
#include <map>
#include "BTrk/BaBar/BbrCollectionUtils.hh"

#include "BTrk/MatEnv/MatElementObj.hh"
#include "BTrk/MatEnv/MatElementList.hh"
//-------------------------------
// Collaborating Class Headers --
//-------------------------------

class MatElmDictionary : public std::map<std::string*, MatElementObj*, babar::Collection::PtrLess>
{

public:

// Constructor 
  MatElmDictionary();

// Destructor
  virtual ~MatElmDictionary();

  void FillElmDict(MatElementList* elmlist);

};

#endif /* MATELMDICTIONARY_HH */

