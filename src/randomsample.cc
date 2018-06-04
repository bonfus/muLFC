/**
 * @file randomsample.cc
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
#include "lattice.h"
#include "dipolesum.h"
#include "randomness.h"
#include "parsers.h"

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

#define BUFFER_SIZE 1000

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
 * @param in_supercell extension of the supercell along the lattice vectors.
 * @param in_cell lattice cell. The three lattice vectors should be entered 
 *         with the following order: a_x, a_y, a_z, b_z, b_y, b_z, c_x, c_y, c_z.
 * @param radius Lorentz sphere radius
 * @param nnn_for_cont number of nearest neighboring atoms to be included
 *                      for the evaluation of the contact field.
 * @param cont_radius only atoms within this radius are eligible to contribute to
 *                      the contact field. This option is redundant but speeds
 *                      up the evaluation significantly
 * @param in_constraints: minimum and maximum distance from aech atom reported in in_positions.
 * @param in_natoms: number of atoms in the lattice.
 * @param in_nmounpos: input, number of muon pos to be evaluated.
 * @param out_muonpos: positions of the muon in fractional coordinates.
 * @param out_field_cont Contact filed in Tesla in the Cartesian coordinates system defined by in_cell. A coupling of 1 \f$ \mathrm{Ang} ^{-1} \sim 13.912~\mathrm{mol/emu} \f$ is assumed.
 * @param out_field_dip  Dipolar field in Tesla in the Cartesian coordinates system defined by in_cell.
 * @param out_field_lor  Lorentz field in Tesla in the Cartesian coordinates system defined by in_cell.
 */


void  RandomSample(const T *in_positions, 
          const T *in_fc, const T *in_K, const T *in_phi,
          const int * in_supercell, const T *in_cell, 
          const T radius, const unsigned int nnn_for_cont, const T cont_radius,
          const unsigned int in_natoms, const unsigned int in_nmounpos,
          const T *in_constraints, const int * in_constraint_active, 
          const unsigned int in_nconstraints,
          T* out_muonpos, T *out_field_cont, T *out_field_dip, T *out_field_lor)
{

    unsigned int scx, scy, scz; /*supercell sizes */
    MatX atomicPos(3,in_natoms);
    MatX constraints(2,in_nconstraints);
    IMatX constraint_group(in_natoms, in_nconstraints);
    MatX MagAtomicPos(3,in_natoms) ;
    MatX muonPositions(3, in_nmounpos);
    
    Vec3 muonPos;
    Mat3 lattice;
    Vec3 boxOShift;
    
    
    /* description of the magnetic structure. */
    /* data provided in cartesian coordinates */
    VecX phi(in_natoms) ;
    CVec3 AtomFC;
    CMatX FC(3,in_natoms);
    Vec3 K;
    MatX B, BLor, BCont;       
    unsigned int a, imu;     /* counter for atoms */
    int mag_atoms;
    T r, boxEdge;
    
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

    // parse constraints
    for (int i=0; i < in_nconstraints; i++) {
        for (int j=0; j< in_natoms; j++) {
            constraint_group(j,i) = in_constraint_active[i*in_natoms + j];
        }
        constraints.col(i).x() = in_constraints[i*2];
        constraints.col(i).y() = in_constraints[i*2+1];
    }

    mag_atoms = ParseAndFilterMagneticAtoms(in_positions, in_fc, in_phi, in_natoms, 
                                  atomicPos, MagAtomicPos, FC, phi);

    
    if (mag_atoms > 0 ) {
        phi.conservativeResize(mag_atoms);
        FC.conservativeResize(3, mag_atoms);
        MagAtomicPos.conservativeResize(3, mag_atoms);
    }
    DistanceCalc DC(lattice, atomicPos);
    UniformRandomInsideUnitCell RndGenerator (lattice, BUFFER_SIZE);

    //std::cout << "Box info " << boxEdge << " - " << boxOShift.transpose() << std::endl;
    /* Generate random positions */
    MatX PositionBuffer(3, BUFFER_SIZE);
    VecX distances(BUFFER_SIZE);
    IVecX valid(BUFFER_SIZE);
    int filled=0;
    
    while(filled < in_nmounpos) {
        for (int i=0; i < BUFFER_SIZE; i++) {
            RndGenerator.GetRandomPos(muonPos);
            PositionBuffer.col(i) = muonPos;
        }
        valid.setOnes();
        for (int j=0; j < in_nconstraints; j++) {
            DC.GetMinDistancesFromAtoms(constraint_group.col(j), PositionBuffer, distances);
            for (int i=0; (i < BUFFER_SIZE); i++) {
                if (distances(i) < constraints.col(j).x() ||
                    distances(i) >= constraints.col(j).y()) valid(i) = 0;
            }
        }
        for (int i=0; (i < BUFFER_SIZE) && (filled < in_nmounpos); i++) {
            if (valid(i) == 1) {
                out_muonpos[filled*3+0] = PositionBuffer.col(i).x();
                out_muonpos[filled*3+1] = PositionBuffer.col(i).y();
                out_muonpos[filled*3+2] = PositionBuffer.col(i).z();
                muonPositions.col(filled) = PositionBuffer.col(i);
                filled++;
            }
        }
    }


    B.resize(3, in_nmounpos); 
    BCont.resize(3, in_nmounpos); 
    BLor.resize(3, in_nmounpos); 

    BCont.setZero(); B.setZero(); BLor.setZero();

    if (mag_atoms > 0) {
        DipoleSumMany(MagAtomicPos, 
          FC, K, phi,
          muonPositions, scx, scy, scz,
          lattice, radius, nnn_for_cont, cont_radius,
          BCont, B, BLor);
    }
    
    
    for (imu = 0; imu < in_nmounpos; imu++) {
    
        out_field_lor[imu*3+0] = BLor.col(imu).x();
        out_field_lor[imu*3+1] = BLor.col(imu).y();
        out_field_lor[imu*3+2] = BLor.col(imu).z();
        
        
        out_field_cont[imu*3+0] = BCont.col(imu).x();
        out_field_cont[imu*3+1] = BCont.col(imu).y();
        out_field_cont[imu*3+2] = BCont.col(imu).z();
        
    
        out_field_dip[imu*3+0] = B.col(imu).x();
        out_field_dip[imu*3+1] = B.col(imu).y();
        out_field_dip[imu*3+2] = B.col(imu).z();
    }

}

} /* extern "C" */

