/**
 * @file rotatesum.c
 * @author Pietro Bonfa
 * @date 9 Sep 2016
 * @brief Dipolar field calculator
 *
 */


#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include "parsers.h"
#include "dipolesum.h"
#include "config.h"

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

void mat3_aangle(const Vec3& v, T r, RefMat3 m);

extern "C" {

/**
 * This function calculates the dipolar field for a set of rotations of the
 * local moments around a specified axis.
 *
 * @param in_positions positions of the magnetic atoms in fractional
 *         coordinates. Each position is specified by the three
 *         coordinates and the 1D array must be 3*in_natoms long.
 * @param in_fc Fourier components. For each atom in in_positions 6 numbers must be specified.
 *              The standard input is Re(FC_x) Im(FC_x) Re(FC_y) Im(FC_y) Re(FC_z) Im(FC_z)
 *              (the input format can be changed by defining the
 *              _ALTERNATE_FC_INPUT at compile time, but this is highly discouraged.)
 *              These values must be provided in the Cartesian coordinate system
 *              defined by in_cell.
 * @param in_K the propagation vector in reciprocal lattice units.
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
 * @param in_natoms: number of atoms in the lattice.
 * @param in_natoms: number of muon sites.
 * @param in_axis: axis for the rotation
 * @param in_nangles: the code will perform in_nangles rotations of 360 deg/in_nangles
 * @param out_field_cont Contact filed in Cartesian coordinates defined by in_cell. A coupling of 1 \f$ \mathrm{Ang} ^{-1} \sim 13.912~\mathrm{mol/emu} \f$ is assumed.
 * @param out_field_dip  Dipolar field in Cartesian coordinates defined by in_cell.
 * @param out_field_lor  Lorentz field in Cartesian coordinates defined by in_cell.
 */
void RotataSum(const double *in_positions,
          const double *in_fc, const double *in_K, const double *in_phi,
          const double *in_muonpos, const int * in_supercell, const double *in_cell,
          const double radius, const unsigned int nnn_for_cont, const double cont_radius,
          const unsigned int in_natoms, const unsigned int in_nmuonpos,
          const double *in_axis, unsigned int in_nangles,
          double *out_field_cont, double *out_field_dip, double *out_field_lor)
{

    unsigned int scx, scy, scz; /*supercell sizes */

    MatX atomicPos(3,in_natoms) ;
    MatX MagAtomicPos(3,in_natoms) ;
    Vec3 muonPos;
    Mat3 lattice;
    Mat3 rmat;
    Vec3 axis;
    T angle;

    /* description of the magnetic structure. */
    /* data provided in cartesian coordinates */
    VecX phi(in_natoms) ;
    CMatX FC(3,in_natoms);
    Vec3 K;
    MatX B, BLor, BCont;
    unsigned int mu, a, angn;     /* counters */
    int mag_atoms;

    MatX tr_mat(3, 3*in_nangles);

    /* define dupercell size */
    scx = in_supercell[0];
    scy = in_supercell[1];
    scz = in_supercell[2];

    /* defines axis */
    axis.x() = in_axis[0];
    axis.y() = in_axis[1];
    axis.z() = in_axis[2];

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


    mag_atoms = ParseAndFilterMagneticAtoms(in_positions, in_fc, in_phi, in_natoms,
                                  atomicPos, MagAtomicPos, FC, phi);


    if (mag_atoms > 0 ) {
        phi.conservativeResize(mag_atoms);
        FC.conservativeResize(3, mag_atoms);
        MagAtomicPos.conservativeResize(3, mag_atoms);
    }

    B.resize(3, in_nangles);
    BLor.resize(3, in_nangles);
    BCont.resize(3, in_nangles);

    B.setZero();
    BLor.setZero();
    BCont.setZero();

    for (angn = 0; angn < in_nangles; ++angn)
    {

        angle = -2*M_PI*((T) angn/ (T) in_nangles);
        mat3_aangle(axis, angle, rmat);

        tr_mat.block(0, angn*3, 3, 3) = rmat;
    }

    for (mu = 0; mu < in_nmuonpos; mu++)
    {

        /* muon position in reduced coordinates */
        muonPos.x() =  in_muonpos[3*mu+0];
        muonPos.y() =  in_muonpos[3*mu+1];
        muonPos.z() =  in_muonpos[3*mu+2];

        BCont.setZero(); B.setZero(); BLor.setZero();

        TransformAndSum(MagAtomicPos,
              FC, K, phi,
              muonPos, scx, scy, scz,
              lattice, radius, nnn_for_cont, cont_radius, tr_mat,
              BCont, B, BLor);

        a = 3*mu*in_nangles;
        for (angn = 0; angn < in_nangles; ++angn)
        {

            out_field_lor[a+0] = BLor.col(angn).x();
            out_field_lor[a+1] = BLor.col(angn).y();
            out_field_lor[a+2] = BLor.col(angn).z();
            //
            out_field_dip[a+0] = B.col(angn).x();
            out_field_dip[a+1] = B.col(angn).y();
            out_field_dip[a+2] = B.col(angn).z();
            //
            out_field_cont[a+0] = BCont.col(angn).x();
            out_field_cont[a+1] = BCont.col(angn).y();
            out_field_cont[a+2] = BCont.col(angn).z();
            a = a + 3;
        }
    }
}

}

/* Get the rotation matrix representation a rotation 'r' around axis 'v' */
/* make sure 'v' is a unit vector */
void mat3_aangle(const Vec3& v, T r, RefMat3 m)
{
	T c, s, C;

	c = cos(r);
	s = -sin(r);
	C = 1 - c;
	m.setZero();

	m.col(0).x() = v.x() * v.x() * C + c;
	m.col(0).y() = v.y() * v.x() * C + v.z() * s;
	m.col(0).z() = v.z() * v.x() * C - v.y() * s;
	m.col(1).x() = v.x() * v.y() * C - v.z() * s;
	m.col(1).y() = v.y() * v.y() * C + c;
	m.col(1).z() = v.z() * v.y() * C + v.x() * s;
	m.col(2).x() = v.x() * v.z() * C + v.y() * s;
	m.col(2).y() = v.y() * v.z() * C - v.x() * s;
	m.col(2).z() = v.z() * v.z() * C + c;
}
