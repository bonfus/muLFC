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
#include "mat3.h"

#ifdef _OPENMP
#include <omp.h>
#endif 


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

    unsigned int scx, scy, scz = 10; //supercell sizes
    unsigned int i,j,k; // counters for supercells
    
    vec3 * atmpos;
    vec3 * muonpos;
    vec3 * r;

        
    mat3 * sc_lat;
    
    double n, rx, ry, rz;
    double onebrcube; // 1/r^3
    double onebrfive; // 1/r^5
    
    
    mat3 * A;

    double Axx=0.0; double Axy=0.0; double Axz=0.0;
    double Ayy=0.0; double Ayz=0.0;
    double Azz=0.0;

    
    unsigned int atom;     // counter for atoms
    
    
    
    
    
    
    // define dupercell size
    scx = in_supercell[0];
    scy = in_supercell[1];
    scz = in_supercell[2];

#ifdef _DEBUG    
    printf("I use: %i %i %i\n",scx, scy, scz);    
    printf("Size is: %i\n",in_natoms);
#endif 
    
    sc_lat = new_mat3(in_cell[0], in_cell[1], in_cell[2],
                      in_cell[3], in_cell[4], in_cell[5],
                      in_cell[6], in_cell[7], in_cell[8]);
 
#ifdef _DEBUG      
    for (i=0;i<3;i++)
        printf("Cell is: %i %e %e %e\n",i,in_cell[i*3],in_cell[i*3+1],in_cell[i*3+2]);
        
    //printf("a %e %e %e\n", sc_lat.a.x, sc_lat.a.y, sc_lat.a.z);
#endif     

    mat3* sctmp = new_mat3_diag((scalar) scx, (scalar) scy, (scalar) scz);
    mat3_mul(sctmp,sc_lat, sc_lat);
    mat3_free(sctmp);
    
    
    // muon position in reduced coordinates
    muonpos = new_vec3( (in_muonpos[0] + (scx/2) ) / (scalar) scx,
                        (in_muonpos[1] + (scy/2) ) / (scalar) scy,
                        (in_muonpos[2] + (scz/2) ) / (scalar) scz);

#ifdef _DEBUG
    printf("Muon pos (frac): %e %e %e\n",muonpos.x,muonpos.y,muonpos.z);
#endif

    mat3_vmul(muonpos,sc_lat,muonpos);

#ifdef _DEBUG
    printf("Muon pos (cart): %e %e %e\n",muonpos.x,muonpos.y,muonpos.z);
#endif


#ifdef _DEBUG
    for (atom = 0; atom < in_natoms; ++atom)
    {
                    
        // atom position in reduced coordinates
        atmpos->x =  in_positions[3*atom] ;
        atmpos->y =  in_positions[3*atom+1] ;
        atmpos->z =  in_positions[3*atom+2] ;
        
        printf("Atom pos (crys): %e %e %e\n",atmpos.x,atmpos.y,atmpos.z);
        
        // go to cartesian coordinates (in Angstrom!)
        atmpos = mat3_vmul(atmpos,sc_lat);  
        
        printf("Atom pos (cart): %e %e %e\n",atmpos.x,atmpos.y,atmpos.z);
    }
        
#endif


atmpos = new_vec3_zero();


#pragma omp parallel
{    
#pragma omp for collapse(3) private(i,j,k,atom,r,n,atmpos,D,onebrcube,onebrfive) reduction(+:Axx,Axy,Axz,Ayy,Ayz,Azz)
    for (i = 0; i < scx; ++i)
    {
        for (j = 0; j < scy; ++j)
        {
            for (k = 0; k < scz; ++k)
            {
                // loop over atoms
                for (atom = 0; atom < in_natoms; ++atom)
                {
                    
                    // atom position in reduced coordinates
                    vec3_set(atmpos, ( in_positions[3*atom] + (double) i) / (double) scx,
                                     ( in_positions[3*atom+1] + (double) j) / (double) scy,
                                     ( in_positions[3*atom+2] + (double) k) / (double) scz);
                    
                    // go to cartesian coordinates (in Angstrom!)
                    mat3_vmul(atmpos, sc_lat, atmpos);
                    
                    //printf("atompos: %e %e %e\n", atmpos.x, atmpos.y, atmpos.z);
                    // difference between atom pos and muon pos (cart coordinates)
                    
                    vec3_sub(atmpos,muonpos);
                    r = atmpos;
                    vec3_muls(-1.,r);
                    n = vec3_norm(r);
                    if (n < in_radius)
                    {


                        // vector
                        onebrcube = 1.0/pow(n,3);
                        onebrfive = 1.0/pow(n,5);
                        
                        vec3_get(r, &rx, &ry, &rz);
                        // See uSR bible (Yaouanc Dalmas De Reotier, page 81)
                        // alpha = x
                        Axx += -onebrcube+3.0*rx*rx*onebrfive;
                        Axy += 3.0*rx*ry*onebrfive;
                        Axz += 3.0*rx*rz*onebrfive;
                        
                        // alpha = y
                        Ayy += -onebrcube+3.0*ry*ry*onebrfive;
                        Ayz += 3.0*ry*rz*onebrfive;
                        
                        // alpha = z
                        Azz += -onebrcube+3.0*rz*rz*onebrfive;
                                           
                        
                    }                    

                }
                
            }
        }
    }
}
    
    A = new_mat3_zero();
    mat3_set(A, Axx, Axy, Axz, Axy, Ayy, Ayz, Axz, Ayz, Azz);
    mat3_getp(A, out_field);
    mat3_free(A);
    vec3_free(muonpos);
    mat3_free(sc_lat);
    vec3_free(atmpos);
}


