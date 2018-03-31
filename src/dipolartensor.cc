/**
 * @file dipolartensor.c
 * @author Ifeanyi 
 * @date   27 Jan 2017
 * @author Pietro Bonfa
 * @date 9 Sep 2016
 * @brief Dipolar tensor calculator.
 * 
 * 
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
#include "types.h"
#include "lattice.h"

#ifdef _OPENMP
#include <omp.h>
#endif 

extern "C" {
  
/**
 * This function calculates dipolar tensors in fractional coordinates.
 * 
 * @param in_positions positions of the magnetic atoms in fractional 
 *         coordinates. Each position is specified by the three
 *         coordinates and the 1D array must be 3*size long.
 * @param in_muonpos position of the muon in fractional coordinates
 * @param in_supercell extension of the supercell along the lattice vectors.
 * @param in_cell lattice cell. The three lattice vectors should be entered 
 *         with the following order: a_x, a_y, a_z, b_z, b_y, b_z, c_x, c_y, c_z.
 * @param in_radius Lorentz sphere radius
 * @param in_natoms: number of atoms in the lattice.
 * @param out_field the dipolar tensor as 1D array with 9 entries: 
 *          T_11, T_12, T_13, T_21, T_22, T_23, T_31, T_32, T_33
 *          \f{equation}{ T= 
 *               \begin{matrix}
 *               1 & 2 & 3 \\
 *               4 & 5 & 6 \\
 *               7 & 8 & 9
 *               \end{matrix}
 *           \f}
 */
void DipolarTensor(const double *in_positions, 
          const double *in_muonpos, const int * in_supercell, const double *in_cell, 
          const double in_radius, unsigned int in_natoms,
          double *out_field) 
{

    unsigned int scx, scy, scz = 10; /*supercell sizes */
    unsigned int i,j,k; /* counters for supercells */
    
    Vec3 atmpos;
    Vec3 atmposCart;
    Vec3 muonpos;
    Vec3 muonposCart;
    Vec3 r, aux;

        
    Mat3 sc_lat;
    Mat3 lattice;
    
    double n;
    double onebrcube; /* 1/r^3 */
    double onebrfive; /* 1/r^5 */
    
    
    Mat3 A, D;
    
    unsigned int atom;     /* counter for atoms */

    /* define dupercell size */
    scx = in_supercell[0];
    scy = in_supercell[1];
    scz = in_supercell[2];

#ifdef _DEBUG    
    printf("I use: %i %i %i\n",scx, scy, scz);    
    printf("Size is: %i\n",in_natoms);
#endif 
    
    lattice.col(0) << in_cell[0], in_cell[1], in_cell[2];
    lattice.col(1) << in_cell[3], in_cell[4], in_cell[5];
    lattice.col(2) << in_cell[6], in_cell[7], in_cell[8];
    
    aux << (double) scx, (double) scy, (double) scz;
    sc_lat = lattice * aux.asDiagonal();
    
    
    /* muon position in reduced coordinates */
    muonpos.x() =  (in_muonpos[0] + (scx/2) ) / (float) scx;
    muonpos.y() =  (in_muonpos[1] + (scy/2) ) / (float) scy;
    muonpos.z() =  (in_muonpos[2] + (scz/2) ) / (float) scz;

    Crys2Cart(sc_lat, muonpos, muonposCart, false);



    A.setZero();
    D.setZero();


    for (i = 0; i < scx; ++i)
    {
        for (j = 0; j < scy; ++j)
        {
            for (k = 0; k < scz; ++k)
            {
                /* loop over atoms */
                for (atom = 0; atom < in_natoms; ++atom)
                {
                    
                    /* atom position in reduced coordinates */
                    atmpos.x() = ( in_positions[3*atom] + (float) i) / (float) scx;
                    atmpos.y() = ( in_positions[3*atom+1] + (float) j) / (float) scy;
                    atmpos.z() = ( in_positions[3*atom+2] + (float) k) / (float) scz;
                    

                    
                    /* go to cartesian coordinates (in Angstrom!) */
                    Crys2Cart(sc_lat, atmpos, atmposCart, false);
                    
                    /*printf("atompos: %e %e %e\n", atmpos.x, atmpos.y, atmpos.z); */
                    /* difference between atom pos and muon pos (cart coordinates) */
                    
                    r = atmposCart-muonposCart;
                    
                    n = r.norm();
                    if (n < in_radius)
                    {


                        /* vector */
                        onebrcube = 1.0/pow(n,3);
                        onebrfive = 1.0/pow(n,5);
                        
                        
                        /* See uSR bible (Yaouanc Dalmas De Reotier, page 81) */
                        /* alpha = x */
                        D.col(0).x() = -onebrcube+3.0*r.x()*r.x()*onebrfive;
                        D.col(0).y() = 3.0*r.x()*r.y()*onebrfive;
                        D.col(0).z() = 3.0*r.x()*r.z()*onebrfive;
                        
                        /* alpha = y */
                        D.col(1).x() = D.col(0).y();
                        D.col(1).y() = -onebrcube+3.0*r.y()*r.y()*onebrfive;
                        D.col(1).z() = 3.0*r.y()*r.z()*onebrfive;
                        
                        /* alpha = z */
                        D.col(2).x() = D.col(0).z();
                        D.col(2).y() = D.col(1).z(); /* this is pretty stupid */
                        D.col(2).z() = -onebrcube+3.0*r.z()*r.z()*onebrfive;
                        
                        A += D ;

                        
                    }
                }
            }
        }
    }


    /* B = vec3_muls(0.9274009, B); // we should multiply for a volume */
    out_field[0] = A.col(0).x(); out_field[1] = A.col(0).y(); out_field[2] = A.col(0).z();
    out_field[3] = A.col(1).x(); out_field[4] = A.col(1).y(); out_field[5] = A.col(1).z();
    out_field[6] = A.col(2).x(); out_field[7] = A.col(2).y(); out_field[8] = A.col(2).z();

}

} // extern C
