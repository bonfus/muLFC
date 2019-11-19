#pragma once
#include <iostream>
#include "types.h"


void Crys2Cart(const Mat3& trmat, const MatX& pos, RefMatX out);
void Cart2Crys(const Mat3& trmat, const MatX& pos, RefMatX out);
void recips (const Mat3& lat, RefMatX rec);

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
