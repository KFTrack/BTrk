#ifndef DCHCELLADDR_HH
#define DCHCELLADDR_HH

//--------------------------------------------------------------------------
//
// Environment:
//      This software was developed for the BaBar collaboration.  If you
//      use all or part of it, please give an appropriate acknowledgement.
//
// Copyright Information:
//      Copyright (C) 1999      <INFN>
//
//------------------------------------------------------------------------

//		---------------------
// 		-- Class Interface --
//		---------------------


/**
 *  C++ source file code DchCellAddr. 
 *  This class is only a container for three static functions to map a Dch 
 *  cell address to its layer/wire numbers. Needed to break some circular
 *  dependency within the Dch code
 *
 *  This software was developed for the BaBar collaboration.  If you
 *  use all or part of it, please give an appropriate acknowledgement.
 *
 *  Copyright (C) 1999 [INFN & Padova University]
 *
 *  @see DchCellAddrDchCellAddr
 *
 *  @version $Id: DchCellAddr.hh 92 2010-01-14 12:38:30Z stroili $ 
 *
 *  @author (R. Stroili)	(originator);
 *  
 */

class DchCellAddr {

//--------------------
// Instance Members --
//--------------------

public:

  // Constructors
  DchCellAddr( void );

  // Destructor
  virtual ~DchCellAddr( ) {;}

  // Operators
    
//------------------
// Static Members --
//------------------

public:

  // Selectors (const)
  static int wireIs(const int &cell)  { return cell%1000; }
  static int layerIs(const int &cell) { return cell/1000; }
  static int cellIs(const int &wire, const int &layer) { return 
							   layer*1000+wire; }

};

#endif // DCHCELLADDR_HH
