#pragma once
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

class Lattice {
  public:
    Lattice(const Mat3& unitCell,
            const MatX& atomicPositions,
            const VecX& atomicOccupations,
            const IVecX& atomicOccupationsGroups,
            const MatX& sitesCorrelation,
            const VecX& _phi,
            const CMatX& _FC,
            const Vec3& _K);
    Lattice(const Mat3& unitCell, 
                  const MatX& atomicPositions, 
                  const VecX& _Phi,
                  const CMatX& _FC,
                  const Vec3& _K);
    void GetCell(RefMat3 cell);
    void GetReciprocalCell(RefMat3 cell);
    void MaterializeOccupationsInCell(RefIVecX occupations);
    void SetOccupations(VecX& atomicOccupations, IVecX& atomicOccupationsGroups, MatX& sitesCorrelation);
    T GetCorrelationTemperature(const IVecX& occupations);
  public:
    Mat3 directCell;
    Mat3 recirpcalCell;
    MatX atomPosFrac;
    MatX atomPosCart;
    VecX atomFracOcc;
    IVecX atomFracOccGroups;
    bool anyPartialOccupation;
    IMatX atomFracTable;
    MatX correlation;
  public:
    unsigned int nAtoms;
    VecX Phi;
    CMatX FC;
    Vec3 K;
};

class LatticeException {
  public:
    explicit LatticeException( const std::string & mesg ) : mesg_(mesg) {}
    std::string mesg_;
};

extern "C" {
#include "clattice.h"
}
