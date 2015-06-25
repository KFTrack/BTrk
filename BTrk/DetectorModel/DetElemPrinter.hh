
class DetSet;
class DetType;
#include <iosfwd>
class HepPoint;
#include <string>
#include <list>
#include <vector>

class DetElemPrinter
{

public:
  // constructor
  DetElemPrinter( const DetSet& set ) : _set( set ) {}

  // destructor
  ~DetElemPrinter() {}

  // print
  void print( std::ostream& ) const;

private:
  
  // the set
  const DetSet& _set;


  // print one face
  void printFace( const std::vector< HepPoint > & vertexList, std::ostream& )  const;

  // replace blank by underscores
  void removeBlanks( std::string& c ) const;

};

// the type hash function
unsigned detElemPrinterTypeHash( const DetType* const & p );















