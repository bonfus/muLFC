#include <iostream>
#include "lattice.h"

void Crys2Cart(const Mat3& trmat, const MatX& pos, RefMatX out, bool reverseDirection)
{
  Vec3 aux;
  if (reverseDirection) {
    out = trmat.inverse() * pos;
  } else {
    out = trmat * pos;
  }
}


void recips (const Mat3& lat, RefMatX rec)
{
  int ii,lattype;
  T bzvol, ucvol;
  Vec3 a1,a2,a3,b1,b2,b3;
  T pi = 3.1415926535897932384626433;

  a1 = lat.col(0);
  a2 = lat.col(1);
  a3 = lat.col(2);

  ucvol=a1(0)*(a2(1)*a3(2)-a2(2)*a3(1))+a1(1)*(a2(2)*a3(0)-a2(0)*a3(2))+a1(2)*(a2(0)*a3(1)-a2(1)*a3(0));
  bzvol=8*pi*pi*pi/ucvol;

  //calculate reciprocal-space lattice vectors b1-b3
  b1(0)=2.0*pi*(a2(1)*a3(2)-a2(2)*a3(1))/ucvol;
  b1(1)=2.0*pi*(a2(2)*a3(0)-a2(0)*a3(2))/ucvol;
  b1(2)=2.0*pi*(a2(0)*a3(1)-a2(1)*a3(0))/ucvol;
  b2(0)=2.0*pi*(a3(1)*a1(2)-a3(2)*a1(1))/ucvol;
  b2(1)=2.0*pi*(a3(2)*a1(0)-a3(0)*a1(2))/ucvol;
  b2(2)=2.0*pi*(a3(0)*a1(1)-a3(1)*a1(0))/ucvol;
  b3(0)=2.0*pi*(a1(1)*a2(2)-a1(2)*a2(1))/ucvol;
  b3(1)=2.0*pi*(a1(2)*a2(0)-a1(0)*a2(2))/ucvol;
  b3(2)=2.0*pi*(a1(0)*a2(1)-a1(1)*a2(0))/ucvol;

  rec.col(0) = b1;
  rec.col(1) = b2;
  rec.col(2) = b3;
}

DistanceCalc::DistanceCalc(const Mat3& lattice, const MatX& atomicPositions) {
    Mat3 aux;
    Vec3 vaux;
    Mat3 sctmp;
    MatX atomicPosSCFrac;
    double nrm;
    int natoms = atomicPositions.cols();


    sctmp << (T) 3,     0,     0, 
                 0, (T) 3,     0,
                 0,     0, (T) 3;
    
    atomicPosSCFrac.resize(3, 3*3*3*natoms);
    scAtomsPosCart.resize(3, 3*3*3*natoms);
    scAtomNum.resize(3*3*3*natoms);

    // calculate scaled lattice parameters
    scLattice = lattice*sctmp;

    // scale lattice coordinates by 3
    int idx = 0;
    for (int l=0; l < natoms; l++)
    {
        for (int i=0; i < 3; i ++) {
            for (int j=0; j < 3; j ++) {
                for (int k=0; k < 3; k ++) {
                    vaux << ((T) i)/3. ,((T) j)/3. ,((T) k)/3.;
                    
                    atomicPosSCFrac.col(idx) = atomicPositions.col(l)/3. + vaux;
                    scAtomNum(idx) = l;
                    idx++;
                }
            }
        }
    }

    // apos in Cart coordinates
    Crys2Cart(scLattice, atomicPosSCFrac, scAtomsPosCart, false);
}

void DistanceCalc::GetMinDistancesFromAtoms(const IVecX& atoms_type, const MatX& intPositionFrac, RefVecX distances)
{
    MatX intPositionSCFrac;
    MatX intPositionSCCart;
    Vec3 vauxCart;
    double nrm;
    int npositions = intPositionFrac.cols();

    intPositionSCFrac.resize(3, npositions);
    intPositionSCCart.resize(3, npositions);

    // lattice scaled by 3 in all directions
    intPositionSCFrac = intPositionFrac/3. + 0.3333333333*MatX::Ones(3, npositions);
    
    // pos in cartesian coordinates
    Crys2Cart(scLattice, intPositionSCFrac, intPositionSCCart, false);

    
    for (int i=0; i < intPositionSCCart.cols(); i ++) {
        nrm = std::numeric_limits<double>::max();
        for (int l=0; l < scAtomsPosCart.cols(); l++) {
            if (atoms_type(scAtomNum(l)) == 0) continue;
            vauxCart = scAtomsPosCart.col(l) - intPositionSCCart.col(i);
            nrm = std::min(vauxCart.norm(), nrm);
        }
        distances(i) = nrm;
    }
}

T DistanceCalc::GetMinDistanceFromAtoms(const Vec3& intPosition)
{
    Vec3 vauxCart;
    Vec3 intPositionSCFrac;
    Vec3 intPositionSCCart;
    T nrm;

    // lattice scaled by 3 in all directions
    intPositionSCFrac = intPosition/3. + 0.3333333333*Vec3::Ones();
    
    // pos in cartesian coordinates
    Crys2Cart(scLattice, intPositionSCFrac, intPositionSCCart, false);

    nrm = std::numeric_limits<T>::max();

    // atomic positions in 3x3x3 lattice
    for (int l=0; l < scAtomsPosCart.cols(); l++) {
        vauxCart = scAtomsPosCart.col(l) - intPositionSCCart;
        nrm = std::min(vauxCart.norm(), nrm);
    }

    return nrm;
}
