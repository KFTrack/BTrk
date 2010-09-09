#ifndef BFMAP_HH
#define BFMAP_HH
//
// Class to hold one magnetic field map. The map is defined on a regular cartesian grid.
// All field maps are given in the standard Mu2e coordinate system.
// Units are: space point in mm, field values in tesla.
//
// $Id: BFMap.hh,v 1.4 2010/09/01 20:29:02 genser Exp $
// $Author: genser $
// $Date: 2010/09/01 20:29:02 $
//
// Original Rob Kutschke, based on work by Julie Managan and Bob Bernstein.
// Rewritten in part by Krzysztof Genser to save execution time
//

#include <iosfwd>
#include <string>
#include <vector>
#include "mu2eFast/Container3D.hh"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {
  class BFMap {

  public:

    friend class BFieldManagerMaker;

    BFMap():
      _key(),
      _warnIfOutside(false),
      _nx(),
      _ny(),
      _nz(),
      _origin(),
      _grid(),
      _field(),
      _isDefined(){
    }

    BFMap( const std::string& key, bool warnIfOutside=false):
      _key(key),
      _warnIfOutside(warnIfOutside),
      _nx(),
      _ny(),
      _nz(),
      _origin(),
      _grid(),
      _field(),
      _isDefined(){
    }

    BFMap(std::string filename, 
          CLHEP::Hep3Vector const& origin,
          int const nx, 
          int const ny, 
          int const nz,
          bool warnIfOutside=false):
      _key(filename),
      _warnIfOutside(warnIfOutside),
      _nx(nx),
      _ny(ny),
      _nz(nz),
      _origin(origin),
      _grid(_nx,_ny,_nz),
      _field(_nx,_ny,_nz),
      _isDefined(_nx,_ny,_nz,false){
    };
    
    ~BFMap(){};

    // Accessors
    CLHEP::Hep3Vector getBField(CLHEP::Hep3Vector const& point) const;

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

    const CLHEP::Hep3Vector& origin() const { return _origin; }
    const mu2e::Container3D<CLHEP::Hep3Vector>& grid() const { return _grid; }
    const mu2e::Container3D<CLHEP::Hep3Vector>& field() const { return _field; }
    
    mu2e::Container3D<CLHEP::Hep3Vector>& grid() { return _grid; }
    mu2e::Container3D<CLHEP::Hep3Vector>& field() { return _field; }
    
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
    CLHEP::Hep3Vector _origin;

    // Vector arrays for gridpoints and field values
    mu2e::Container3D<CLHEP::Hep3Vector> _grid;
    mu2e::Container3D<CLHEP::Hep3Vector> _field;
    mu2e::Container3D<bool> _isDefined;

    // Functions used internally and by the code that populates the maps.

    // Validity checker
    bool isValid(CLHEP::Hep3Vector const& point) const;

    // method to store the neighbors
    void getNeighbors(int ix, int iy, int iz, CLHEP::Hep3Vector neighborsBF[3][3][3]) const;

    // Interpolator
    CLHEP::Hep3Vector interpolate(CLHEP::Hep3Vector const vec[3][3][3],
                                  double const frac[3]) const;

    // Polynomial fit function used by interpolator
    double gmcpoly2(double const f1d[3], double const& x) const;
public:
    // Define the limits and step sizes for the maps.
    void setLimits ( double xmin, double xmax, double ymin, double ymax,
                     double zmin, double zmax);
    
    void setDefined(unsigned ix, unsigned iy, unsigned iz, bool value=true) {
      _isDefined.set(ix,iy,iz,value);
    }
    
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
private:

  };

} // end namespace mu2e

#endif
