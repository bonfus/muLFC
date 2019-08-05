#include <iostream>
#include <random>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include "lattice.h"
#include "config.h"

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

/* Lattice Class
 *
 */

Lattice::Lattice(const Mat3& unitCell, 
                  const MatX& atomicPositions, 
                  const VecX& atomicOccupations, 
                  const IVecX& atomicOccupationsGroups,
                  const MatX& sitesCorrelation, 
                  const VecX& _Phi,
                  const CMatX& _FC,
                  const Vec3& _K) : Phi(_Phi), FC(_FC), K(_K)
{
    int label;

    directCell = unitCell;
    recips(unitCell, recirpcalCell);

    /* === Atoms data === */
    nAtoms = atomicPositions.cols();
    atomPosCart.resize(3, nAtoms);
    atomPosFrac.resize(3, nAtoms);
    atomFracOcc.resize(nAtoms);
    atomFracOccGroups.resize(nAtoms);
    atomFracTable.resize(nAtoms, nAtoms);
    correlation.resize(nAtoms, nAtoms);

    /* Atomic positions: reduced coordinates */
    atomPosFrac = atomicPositions;

    /* Atomic positions: reduced coordinates -> Cartesian Coordinates */
    Crys2Cart(unitCell, atomicPositions, atomPosCart, false);
    

    atomFracOcc       = atomicOccupations;
    atomFracOccGroups = atomicOccupationsGroups;
    correlation       = sitesCorrelation;

    anyPartialOccupation = (atomFracOcc.array() < 1.0-EPS).any();

    if ( anyPartialOccupation ) {
        atomFracTable.setOnes();
        atomFracTable *= -1;

        for (unsigned int i = 0; i < nAtoms; i++)
        {
            label = atomFracOccGroups(i);
            if ( label > 0 ){
                /* Add atom i in the first available position */
                for (unsigned int j = 0; j < nAtoms; j++) {
                    if (atomFracTable(j, label-1) < 0) {
                        atomFracTable(j, label-1) = i;
                        break;
                    }
                }
            }
        }
    }
}
Lattice::Lattice(const Mat3& unitCell, 
                  const MatX& atomicPositions, 
                  const VecX& _Phi,
                  const CMatX& _FC,
                  const Vec3& _K) : Phi(_Phi), FC(_FC), K(_K)
{

    directCell = unitCell;
    recips(unitCell, recirpcalCell);

    /* === Atoms data === */
    nAtoms = atomicPositions.cols();
    atomPosCart.resize(3, nAtoms);
    atomPosFrac.resize(3, nAtoms);
    atomFracOcc.resize(nAtoms);
    atomFracOccGroups.resize(nAtoms);
    atomFracTable.resize(nAtoms, nAtoms);
    correlation.resize(nAtoms, nAtoms);

    /* Atomic positions: reduced coordinates */
    atomPosFrac = atomicPositions;

    /* Atomic positions: reduced coordinates -> Cartesian Coordinates */
    Crys2Cart(unitCell, atomicPositions, atomPosCart, false);
    

    atomFracOcc.setOnes();
    atomFracOccGroups.setZero();
    correlation.setZero();

    anyPartialOccupation = false;
    atomFracTable.setOnes();
    atomFracTable *= -1;
}

void Lattice::SetOccupations(VecX& atomicOccupations, IVecX& atomicOccupationsGroups, MatX& sitesCorrelation) {

    int label;
    atomFracOcc       = atomicOccupations;
    atomFracOccGroups = atomicOccupationsGroups;
    correlation       = sitesCorrelation;

    anyPartialOccupation = (atomFracOcc.array() < 1.0-EPS).any();

    if ( anyPartialOccupation ) {
        atomFracTable.setOnes();
        atomFracTable *= -1;

        for (unsigned int i = 0; i < nAtoms; i++)
        {
            label = atomFracOccGroups(i);
            if ( label > 0 ){
                /* Add atom i in the first available position */
                for (unsigned int j = 0; j < nAtoms; j++) {
                    if (atomFracTable(j, label-1) < 0) {
                        atomFracTable(j, label-1) = i;
                        break;
                    }
                }
            }
        }
    }
}


void Lattice::GetCell(RefMat3 cell) { cell = directCell; }
void Lattice::GetReciprocalCell(RefMat3 cell) { cell = recirpcalCell; }
void Lattice::MaterializeOccupationsInCell(RefIVecX occupations){

    T rndVal;
    unsigned int atmIdx;
    int label;
    occupations.setOnes();
    if (not anyPartialOccupation)
        return;

    std::mt19937 rng;
    rng.seed(std::random_device()());
    /* produces random numbers in a range [a,b) */
    std::uniform_real_distribution<T> distr(0, 1);

    do {
        for (unsigned int i = 0; i < nAtoms; i++)
        {
            /* search within the table allocations that should be assigned */
            label = atomFracTable(0, i);
        
            if ( label >= 0 ) {
                /* generate a PRNG in the range [0, 1) */
                rndVal = distr(rng);
        
                /* This loop is for all the possible entries in having this label.
                 * The maximum value is nAtoms (occupiation split by all atoms in the cell) */
                for (unsigned int j = 0; j < nAtoms; j++)
                {
                    /* all elements except those part of a group are set -1.
                     * When there are no more elements to check, exti the loop
                     * for this label. */
        
                    if (atomFracTable(j, i) < 0) {
                        break;
                    } else {
                        atmIdx = atomFracTable(j, i);
                    }
        
                    /* If the current value of the random number is negative
                     *  we already assigned the occupation to some other atom */
                    if (rndVal < 0.0) {
                        occupations(atmIdx) = 0;
                        continue;
                    }
        
                    /* if the current value of rndVal is in the interval
                     * [0, atomFracOcc(atmIdx) )
                     * site is assigned.*/
                    if ( rndVal < atomFracOcc(atmIdx) ) {
                        occupations(atmIdx) = 1;
                    } else {
                        occupations(atmIdx) = 0;
                    }
                    /* remove the interval that we just checked to see if the
                     * random value is in the next interval (if any is present)
                     */
                    rndVal -= atomFracOcc(atmIdx);
                }
            }
        }
        if ( GetCorrelationTemperature(occupations) < distr(rng) ) break;
    } while ( true );
}

T Lattice::GetCorrelationTemperature(const IVecX& occupations){

    T temp;
    int label;
    //  std::cout << "corr max " << correlation.array().abs() << std::endl;
    if ((correlation.array().abs() <= EPS).all())
        return 0.0;

    temp = 0.0;
    for (unsigned int i = 0; i < nAtoms; i++)
    {
        /* search within the table allocations that should be assigned */
        if ( occupations(i) == 1 ) {
            for (unsigned int j = 0; j < nAtoms; j++)
                temp += correlation(j, i)*occupations(j);
        }
    }
    return temp;
}


// ----------------
// Python interface
// ----------------


namespace py = pybind11;

void init_lattice(py::module &m) {
    py::class_<Lattice>(m, "Lattice")
    .def(py::init<Mat3, MatX, VecX, CMatX, Vec3>())
    ;
}
