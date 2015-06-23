// $id$
#ifndef DCHWIREARRAYBASE_HH
#define DCHWIREARRAYBASE_HH

#include "DchGeom/DchLayout.hh"
#include <vector>
class DchDetector;

class DchWireArrayBase : public DchLayout
{
protected:
        DchWireArrayBase(const DchDetector& det);
        DchWireArrayBase(const DchLayout& layout);
public:
        virtual ~DchWireArrayBase() {;}
        bool isValid(unsigned layer) const { return layer>=firstLayer() && layer<=lastLayer();}
        bool isValid(unsigned layer,unsigned wire) const { return isValid(layer)
                                                               && wire>=firstWire(layer) && wire<=lastWire(layer); };
        size_t nLayers() const { return _offset.size()-1; }
        size_t nWires() const { return _offset.back(); }
        size_t nWires(unsigned layer) const { return lo(layer+1) - lo(layer); }

        size_t firstLayer() const { return 1;} //FIXME
        size_t firstWire(unsigned layer) const { return _first[le(layer)];}

        size_t lastLayer() const { return firstLayer()+nLayers()-1;}
        size_t lastWire(unsigned layer) const { return firstWire(layer)+nWires(layer)-1;}

        bool hasSameShape(const DchWireArrayBase& i) { return _first==i._first && _offset==i._offset;}
protected:
        unsigned index(unsigned layer,unsigned wire) const { return lo(layer)+(wire-firstWire(layer));}
        unsigned le(unsigned layer) const { return layer-firstLayer();}
        unsigned lo(unsigned layer) const { return _offset[le(layer)];}

private:
        std::vector<unsigned>  _first;   // wire # of the first wire in each layer
        std::vector<ptrdiff_t> _offset;  // offset in _t for the start of each layer

};



class DchWireArrayIteratorBase
{
protected:
        DchWireArrayIteratorBase(unsigned l,unsigned w)
                :_l(l),_w(w) {}
public:
        unsigned layer() const { return _l;}
        unsigned wire() const { return _w;}
protected:
        bool eq(const DchWireArrayIteratorBase& i) const { return _l==i._l && _w==i._w ;}
        bool gt(const DchWireArrayIteratorBase& i) const { return _l>i._l || (_l==i._l&&_w>i._w);}

        bool atEnd(const DchWireArrayBase& p) const { return _l==p.lastLayer()&&_w==p.lastWire(_l)+1;}
        bool atBegin(const DchWireArrayBase& p) const { return _l==p.firstLayer()&&_w==p.firstWire(_l);}
        void inc(const DchWireArrayBase& p) { if(!atEnd(p)   && ++_w>p.lastWire(_l)  &&!atEnd(p))   {++_l;_w=p.firstWire(_l);}}
        void dec(const DchWireArrayBase& p) { if(!atBegin(p) && --_w<p.firstWire(_l) &&!atBegin(p)) {--_l;_w=p.lastWire(_l); }}
        void inc(const DchWireArrayBase& p,int i);
        void dec(const DchWireArrayBase& p,int i);
private:
        size_t _l,_w; // current location
};


#endif
