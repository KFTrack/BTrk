#ifndef BFMAP_HH
#define BFMAP_HH
//
// Class to hold one magnetic field map. The map is defined on a regular cartesian grid.
// All field maps are given in the standard Mu2e coordinate system.
// Units are: space point in mm, field values in tesla.
//
// Original Rob Kutschke, based on work by Julie Managan and Bob Bernstein.
//

#include <iosfwd>
#include <string>
#include <vector>
#include "mu2eFast/Container3D.hh"
#include "CLHEP/Vector/ThreeVector.h"


class BFMap {

public:


  BFMap();
  BFMap( const std::string& key, bool warnIfOutside=false);
  BFMap(const std::string& filename, 
    const Hep3Vector& origin,
    int const nx, 
    int const ny, 
    int const nz,
    bool warnIfOutside=false);
  ~BFMap();

    // Accessors
  Hep3Vector getBField(Hep3Vector const& point) const;

  int nx() const { return _nx; }
  int ny() const { return _ny; }
  int nz() const { return _nz; }

  double xmin() const {return _xmin;}; double xmax() const {return _xmax;};
  double ymin() const {return _ymin;}; double ymax() const {return _ymax;};
  double zmin() const {return _zmin;}; double zmax() const {return _zmax;};

  double dx() const {return _dx;}; 
  double dy() const {return _dy;};
  double dz() const {return _dz;};

  const std::string& getKey() const { return _key; };


  void print( std::ostream& os) const;
private:

    // Filename, database key or other id information that describes
    // where this map came from.
  std::string _key;

    // If true, then print a warning message when a point is outside the region 
    // in which the map is defined; else return a field with a value of (0.,0.,0.);
    // This does happen under normal operation of G4 so we should not warn be default.
  bool _warnIfOutside;

    // Grid dimensions
  unsigned int _nx, _ny, _nz;

    // Min and Max values.
  double _xmin, _xmax, _ymin, _ymax, _zmin, _zmax;

    // Distance between points.
  double _dx, _dy, _dz;

    // Origin from external setup
  Hep3Vector _origin;

    // Vector arrays for gridpoints and field values
  Container3D<Hep3Vector> _grid;
  Container3D<Hep3Vector> _field;

    // Functions used internally and by the code that populates the maps.
public:
  const Hep3Vector& origin() { return _origin;}
  unsigned int nx() { return _nx;}
  unsigned int ny() { return _ny;}
  unsigned int nz() { return _nz;}
  
  Container3D<Hep3Vector>& grid() { return _grid; }
  Container3D<Hep3Vector>& field() { return _field; }
  
    // Validity checker
  bool isValid(Hep3Vector const& point) const;

    // method to store the neighbors
  void getNeighbors(int ix, int iy, int iz, Container3D<Hep3Vector>& neighborsBF) const;

    // Interpolater
  Hep3Vector interpolate(Container3D<Hep3Vector> const vec,
    Hep3Vector const frac) const;

    // Polynomial fit function used by interpolater
  double gmcpoly2(std::vector<double> const& f1d, double const& x) const;

    // Define the limits and step sizes for the maps.
  void setLimits ( double xmin, double xmax, double ymin, double ymax,
    double zmin, double zmax);

    // Compute grid indices for a given point.
  std::size_t iX( double x){
    return static_cast<int>((x - _xmin)/_dx + 0.5);
  }

  std::size_t iY( double y){
    return static_cast<int>((y - _ymin)/_dy + 0.5);
  }

  std::size_t iZ( double z){
    return static_cast<int>((z - _zmin)/_dz + 0.5);
  }

};


#endif
