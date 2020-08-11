#include <iostream>
#include <random>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include "lattice.h"
#include "config.h"

void Crys2Cart(const Mat3& trmat, const MatX& pos, RefMatX out)
{
    out = trmat * pos;
}

void Cart2Crys(const Mat3& trmat, const MatX& pos, RefMatX out)
{
    out = trmat.inverse() * pos;
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

    /*calculate reciprocal-space lattice vectors b1-b3 */
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
                 const Vec3& _K)
{
    int label;

    directCell = unitCell;
    recips(unitCell, recirpcalCell);

    /* === Magnetic data === */
    Phi = _Phi;
    FC  = _FC;
    K   = _K;
    Crys2Cart(recirpcalCell, K, KCart);

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
    Crys2Cart(unitCell, atomicPositions, atomPosCart);


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
            if ( label > 0 ) {
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
                 const Vec3& _K)
{

    directCell = unitCell;
    recips(unitCell, recirpcalCell);

    /* === Magnetic data === */
    Phi = _Phi;
    FC  = _FC;
    K   = _K;
    Crys2Cart(recirpcalCell, K, KCart);

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
    Crys2Cart(unitCell, atomicPositions, atomPosCart);


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
            if ( label > 0 ) {
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

void Lattice::GetMagneticMoments(const Vec3& r, RefMatX atomMoms) {

    unsigned int a;     /* counter for atoms */
    T c, s, KdotR;      /*cosine and sine of K.R */
    Vec3 R;

    Crys2Cart(directCell, r, R);

    KdotR =  KCart.dot(R); /* R \dot K (done in cartesian coordinates) */
    for (a = 0; a < nAtoms; ++a)
    {
        c = cos (KdotR + 2.0*M_PI*Phi(a) );
        s = sin (KdotR + 2.0*M_PI*Phi(a) );

        atomMoms.col(a) = c*FC.col(a).real() + s*FC.col(a).imag();
    }
}

void Lattice::GetCell(RefMat3 cell) {
    cell = directCell;
}
void Lattice::GetReciprocalCell(RefMat3 cell) {
    cell = recirpcalCell;
}
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
    /*  std::cout << "corr max " << correlation.array().abs() << std::endl; */
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


/* ---------------- */
/* Python interface */
/* ---------------- */


namespace py = pybind11;

void init_lattice(py::module &m) {
    py::class_<Lattice>(m, "Lattice")
    .def(py::init<Mat3, MatX, VecX, CMatX, Vec3>())
    ;
}
