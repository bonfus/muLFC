#include <iostream>
#include "types.h"


void Crys2Cart(const Mat3& trmat, const MatX& pos, RefMatX out, bool reverseDirection);
void recips (const Mat3& lat, RefMatX rec);

class DistanceCalc {
  public:
    DistanceCalc(const Mat3& lattice, const MatX& atomicPositions);
    void GetMinDistancesFromAtoms(const IVecX& atoms_type, const MatX& intPositionFrac, RefVecX distances);
    T GetMinDistanceFromAtoms(const Vec3& intPosition);
  private:
    Mat3 scLattice;
    MatX scAtomsPosCart;
    VecX scAtomNum;
};
