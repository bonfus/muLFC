/**
 * @file simplesum.cc
 * @author Pietro Bonfa
 * @date 18 March 2018
 * @brief Dipolar field calculator
 *
 */

#define _USE_MATH_DEFINES
#ifdef _OPENMP
	#include <omp.h>
#endif
#include <stdio.h>
#include <math.h>
#ifndef NAN
#error Nan support is required
#endif

#include "pile.h"
#include "config.h"
#include "dipolesum.h"
#include "parsers.h"

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif


extern "C" {

/**
 * This function calculates the dipolar field for a muon site.
 *
 * @param in_positions positions of the magnetic atoms in fractional
 *         coordinates. Each position is specified by the three
 *         coordinates and the 1D array must be 3*in_natoms long.
 * @param in_fc Fourier components. For each atom in in_positions 6 numbers must be specified.
 *              The standard input is Re(FC_x) Im(FC_x) Re(FC_y) Im(FC_y) Re(FC_z) Im(FC_z)
 *              (the input format can be changed by defining the
 *              _ALTERNATE_FC_INPUT at compile time, but this is highly discouraged).
 *              These values must be provided in the Cartesian coordinate system
 *              defined by in_cell.
 * @param in_K the propagation vector in *reciprocal lattice units*.
 * @param in_phi the phase for each of the atoms given in in_positions.
 * @param in_muonpos position of the muon in fractional coordinates
 * @param in_supercell extension of the supercell along the lattice vectors.
 * @param in_cell lattice cell. The three lattice vectors should be entered
 *         with the following order: a_x, a_y, a_z, b_z, b_y, b_z, c_x, c_y, c_z.
 * @param radius Lorentz sphere radius
 * @param nnn_for_cont number of nearest neighboring atoms to be included
 *                      for the evaluation of the contact field.
 * @param cont_radius only atoms within this radius are eligible to contribute to
 *                      the contact field. This option is redundant but speeds
 *                      up the evaluation significantly
 * @param min_radius_from_atoms: minimal distance from atoms reported in in_positions.
 *                      When the distance from an atom of the system is smaller than required,
 *                      NaN is returned.
 * @param in_natoms: number of atoms in the lattice.
 * @param in_nmounpos: input, number of muon positions to be evaluated.
 * @param out_field_cont Contact filed in Tesla in the Cartesian coordinates system defined by in_cell. A coupling of 1 \f$ \mathrm{Ang} ^{-1} \sim 13.912~\mathrm{mol/emu} \f$ is assumed.
 * @param out_field_dip  Dipolar field in Tesla in the Cartesian coordinates system defined by in_cell.
 * @param out_field_lor  Lorentz field in Tesla in the Cartesian coordinates system defined by in_cell.
 */


void  SimpleSum(const T *in_positions,
          const T *in_fc, const T *in_K, const T *in_phi,
          const T *in_muonpos, const int * in_supercell, const T *in_cell,
          const T radius, const unsigned int nnn_for_cont, const T cont_radius,
          const T min_radius_from_atoms,
          const unsigned int in_natoms, unsigned int in_nmounpos,
          T *out_field_cont, T *out_field_dip, T *out_field_lor)
{

    unsigned int scx, scy, scz; /*supercell sizes */

    MatX atomicPos(3,in_natoms) ;
    MatX MagAtomicPos(3,in_natoms) ;
    Vec3 muonPos;
    Mat3 lattice;


    /* description of the magnetic structure. */
    /* data provided in cartesian coordinates */
    VecX phi(in_natoms) ;
    CMatX FC(3,in_natoms);
    Vec3 K;
    Vec3 B, BLor, BCont;
    unsigned int a, imu;     /* counter for atoms */
    int mag_atoms;
    T r;

    /* define dupercell size */
    scx = in_supercell[0];
    scy = in_supercell[1];
    scz = in_supercell[2];


    lattice.col(0) << in_cell[0], in_cell[1], in_cell[2];
    lattice.col(1) << in_cell[3], in_cell[4], in_cell[5];
    lattice.col(2) << in_cell[6], in_cell[7], in_cell[8];

#ifdef _DEBUG
    for (i=0;i<3;i++)
        printf("Cell is: %i %e %e %e\n",i,in_cell[i*3],in_cell[i*3+1],in_cell[i*3+2]);

    /*printf("a %e %e %e\n", sc_lat.a.x, sc_lat.a.y, sc_lat.a.z); */
#endif

    K.x() = in_K[0]; K.y() = in_K[1]; K.z() = in_K[2];

#ifdef _DEBUG
    std::cout << K.transpose() << std::endl;
    printf("Radius is: %e\n",radius);
#endif

    /* Filter out non magnetic atoms, according to EPS value in config.h */
    mag_atoms = ParseAndFilterMagneticAtoms(in_positions, in_fc, in_phi, in_natoms,
                                  atomicPos, MagAtomicPos, FC, phi);
    if (mag_atoms > 0 ) {
        phi.conservativeResize(mag_atoms);
        FC.conservativeResize(3, mag_atoms);
        MagAtomicPos.conservativeResize(3, mag_atoms);
    }
    /* End of atom filtering */

    DistanceCalc DC(lattice, atomicPos);

    for (imu = 0; imu < in_nmounpos; imu++) {

        /* muon position in reduced coordinates */
        muonPos.x() =  in_muonpos[imu*3 + 0];
        muonPos.y() =  in_muonpos[imu*3 + 1];
        muonPos.z() =  in_muonpos[imu*3 + 2];

        /* check distance from atoms is fine */
        if (min_radius_from_atoms > 0.) {
            r = DC.GetMinDistanceFromAtoms(muonPos);
            if (r < min_radius_from_atoms) {
                out_field_lor[imu*3+0] = NAN;
                out_field_lor[imu*3+1] = NAN;
                out_field_lor[imu*3+2] = NAN;
                out_field_cont[imu*3+0] = NAN;
                out_field_cont[imu*3+1] = NAN;
                out_field_cont[imu*3+2] = NAN;
                out_field_dip[imu*3+0] = NAN;
                out_field_dip[imu*3+1] = NAN;
                out_field_dip[imu*3+2] = NAN;
                continue;
            }
        }


        BCont.setZero(); B.setZero(); BLor.setZero();

        if (mag_atoms > 0) {
            DipoleSum(MagAtomicPos,
              FC, K, phi,
              muonPos, scx, scy, scz,
              lattice, radius, nnn_for_cont, cont_radius,
              BCont, B, BLor);
        }


        out_field_lor[imu*3+0] = BLor.x();
        out_field_lor[imu*3+1] = BLor.y();
        out_field_lor[imu*3+2] = BLor.z();


        out_field_cont[imu*3+0] = BCont.x();
        out_field_cont[imu*3+1] = BCont.y();
        out_field_cont[imu*3+2] = BCont.z();


        out_field_dip[imu*3+0] = B.x();
        out_field_dip[imu*3+1] = B.y();
        out_field_dip[imu*3+2] = B.z();
    }

}

void  SimpleSum2(CLattice * latt,
          const T *in_muonpos, const int * in_supercell,
          const T radius, const unsigned int nnn_for_cont, const T cont_radius,
          const T min_radius_from_atoms,
          unsigned int in_nmounpos,
          T *out_field_cont, T *out_field_dip, T *out_field_lor)
{

    unsigned int scx, scy, scz; /*supercell sizes */

    Vec3 muonPos;
    Lattice *l;

    Vec3 B, BLor, BCont;
    unsigned int imu;     /* counter for muons */
    T r;

    /* define dupercell size */
    scx = in_supercell[0];
    scy = in_supercell[1];
    scz = in_supercell[2];

    l = reinterpret_cast<Lattice*>(latt);
    /* End of atom filtering */

    for (imu = 0; imu < in_nmounpos; imu++) {

        /* muon position in reduced coordinates */
        muonPos.x() =  in_muonpos[imu*3 + 0];
        muonPos.y() =  in_muonpos[imu*3 + 1];
        muonPos.z() =  in_muonpos[imu*3 + 2];


        BCont.setZero(); B.setZero(); BLor.setZero();

        DipoleSum2(*l,
              muonPos, scx, scy, scz,
              radius, nnn_for_cont, cont_radius,
              BCont, B, BLor);


        out_field_lor[imu*3+0] = BLor.x();
        out_field_lor[imu*3+1] = BLor.y();
        out_field_lor[imu*3+2] = BLor.z();


        out_field_cont[imu*3+0] = BCont.x();
        out_field_cont[imu*3+1] = BCont.y();
        out_field_cont[imu*3+2] = BCont.z();


        out_field_dip[imu*3+0] = B.x();
        out_field_dip[imu*3+1] = B.y();
        out_field_dip[imu*3+2] = B.z();
    }

}

} /* extern "C" */

