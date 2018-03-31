#include <iostream>
#include "types.h"

T GetMinDistanceFromAtoms(const Mat3& lattice, const MatX& atomicPositions, const Vec3& intPosition);
void Crys2Cart(const Mat3& trmat, const MatX& pos, RefMatX out, bool reverseDirection);
void recips (const Mat3& lat, RefMatX rec);
T GetMinDistanceFromAtoms(const Mat3& lattice, const MatX& atomicPositions, const Vec3& intPosition);
