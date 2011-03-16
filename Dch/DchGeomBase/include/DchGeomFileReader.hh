#ifndef DCHGEOMFILEREADER_HH
#define DCHGEOMFILEREADER_HH

//--------------------------------------------------------------------------
//
// Environment:
//      This software was developed for the BaBar collaboration.  If you
//      use all or part of it, please give an appropriate acknowledgement.
//
// Copyright Information:
//      Copyright (C) 2000      <INFN Padova>
//
//------------------------------------------------------------------------

//-------------
// C Headers --
//-------------

//---------------
// C++ Headers --
//---------------
#include <string>

//----------------------
// Base Class Headers --
//----------------------

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
class DchGFileParser;


//		---------------------
// 		-- Class Interface --
//		---------------------


/**
 *  C++ source file code DchGeomFileReader. Simple factory for
 *  providing the appropriate parser for the ASCII geometry files.
 *  At the moment the only files accepted are the "old" style,
 *  i.e. the ProtoII geometry file, and the simulation dbio file.
 *  The first kind of files is defined by the extension _dat_ while
 *  the second one is defined by the estension _db_
 *
 *  Copyright (C) 2000 [INFN]
 *
 *  @see DchGeomFileReaderDchGeomFileReader
 *
 *  @version $Id: DchGeomFileReader.hh 92 2010-01-14 12:38:30Z stroili $ 
 *
 *  @author (R. Stroili)		(originator);
 *
 */

class DchGeomFileReader {

//--------------------
// Declarations     --
//--------------------

    // Typedefs, consts, and enums

//--------------------
// Instance Members --
//--------------------

public:

  // Constructors
  DchGeomFileReader( );

  // Destructor
  virtual ~DchGeomFileReader( ) {;}


  /**
   *  This class has only one function that gets as input parameter a file
   *  name and returns a pointer to the appropriate parser for that kind of 
   *  file
   *
   *  @param is           input file name
   *  @return             the pointer to a file parser
   *
   *  @see DchGeomFileReaderDchGeomFileReader#fileParser
   *
   *  Note: YOU are responsible for deleting the DchGFileParser 
   *        returned by the factory
   */

  DchGFileParser* fileParser( const std::string& file ) const;

};

#endif // DCHGEOMFILEREADER_HH
