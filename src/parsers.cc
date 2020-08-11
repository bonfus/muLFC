#define _USE_MATH_DEFINES

#include <math.h>
#include <iostream>
#ifndef NAN
#error Nan support is required
#endif

#include "config.h"
#include "types.h"

int ParseAndFilterMagneticAtoms(const T *in_positions,
          const T *in_fc, const T *in_phi,
          const unsigned int in_natoms,
          RefMatX AtomPos, RefMatX MagAtomPos, RefCMatX FC, RefVecX phi)
{

    CVec3 AtomFC;
    unsigned int a;
    int mag_atoms=0;

    // std::cout << "in_natoms" << in_natoms << AtomPos.cols() <<  MagAtomPos.cols() << FC.cols() << std::endl;

    /* Validate input*/
    if (AtomPos.cols() != in_natoms) return -1;
    if (MagAtomPos.cols() != in_natoms) return -2;
    if (FC.cols() != in_natoms) return -3;
    if (phi.size() != in_natoms) return -4;


    for (a = 0; a < in_natoms; a++)
    {
        AtomFC.real() << in_fc[6*a] , in_fc[6*a+2] , in_fc[6*a+4];
        AtomFC.imag() << in_fc[6*a+1] , in_fc[6*a+3] , in_fc[6*a+5];

        // std::cout << "AtomFC" << AtomFC.transpose() << std::endl;

        AtomPos.col(a) << in_positions[3*a] , in_positions[3*a+1] , in_positions[3*a+2];
        if ( AtomFC.norm() >  EPS ){
            FC.col(mag_atoms) = AtomFC;
            MagAtomPos.col(mag_atoms) = AtomPos.col(a);
            phi(mag_atoms) = in_phi[a];
            mag_atoms++;
        }
    }

    return mag_atoms;
}

