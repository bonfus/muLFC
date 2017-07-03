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
#include <math.h>
#include "mat3.h"
#include "pile.h"
#include "config.h"

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif


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
          unsigned int in_natoms, 
          const double *in_axis, unsigned int in_nangles,
          double *out_field_cont, double *out_field_dip, double *out_field_lor)
{

    unsigned int scx, scy, scz = 10; /*supercell sizes */
    unsigned int i,j,k; /* counters for supercells */
    
    struct vec3 atmpos;
    struct vec3 muonpos;
    struct vec3 r;
    struct vec3 m;   /*magnetic moment of atom */
    struct vec3 rm;   /*rotated magnetic moment of atom */
    struct vec3 u;   /* unit vector */
        
    struct mat3 sc_lat;
    
    double n; /* contains norm of vectors */
    double c,s; /*cosine and sine of K.R */
    double onebrcube; /* 1/r^3 */
    double angle;
    
    struct vec3 sk  ;
    struct vec3 isk ;
    double  phi ;
    
    struct vec3 R, K;
    
    unsigned int a, angn;     /* counter for atoms */
    
    /* for rotation */
    struct vec3 axis;
    struct mat3 rmat;
	struct vec3 * B = malloc(in_nangles * sizeof(struct vec3));
	struct vec3 * BLor = malloc(in_nangles * sizeof(struct vec3));
	pile * MCont = malloc(in_nangles * sizeof(pile));
	struct vec3 BCont;
	int NofM = 0;
	double SumOfWeights = 0;

    /* defines axis */
    axis.x = in_axis[0];
    axis.y = in_axis[1];
    axis.z = in_axis[2];
    
    
    
    
    /* define dupercell size */
    scx = in_supercell[0];
    scy = in_supercell[1];
    scz = in_supercell[2];
    
#ifdef _DEBUG    
    printf("I use: %i %i %i\n",scx, scy, scz);    
    printf("Total atoms: %i\n",in_natoms);
    printf("N Angles: %u\n",in_nangles);
#endif
    
    sc_lat.a.x = in_cell[0];
    sc_lat.a.y = in_cell[1];
    sc_lat.a.z = in_cell[2];
    sc_lat.b.x = in_cell[3];
    sc_lat.b.y = in_cell[4];
    sc_lat.b.z = in_cell[5];
    sc_lat.c.x = in_cell[6];
    sc_lat.c.y = in_cell[7];
    sc_lat.c.z = in_cell[8];

#ifdef _DEBUG    
    for (i=0;i<9;i++)
        printf("Cell is: %i %e\n",i,in_cell[i]);
#endif
    
    K.x = in_K[0];
    K.y = in_K[1];
    K.z = in_K[2];

#ifdef _DEBUG
    printf("K is: %e %e %e \n",K.x,K.y,K.z);
    printf("Radius is: %e\n",radius);
#endif

    sc_lat = mat3_mul(
                        mat3_diag((float) scx, (float) scy, (float) scz),
                        sc_lat);
    
    
    /* muon position in reduced coordinates */
    muonpos.x =  (in_muonpos[0] + (scx/2) ) / (float) scx;
    muonpos.y =  (in_muonpos[1] + (scy/2) ) / (float) scy;
    muonpos.z =  (in_muonpos[2] + (scz/2) ) / (float) scz;
    
    muonpos = mat3_vmul(muonpos,sc_lat);    

#ifdef _DEBUG
    printf("Muon pos is: %e %e %e\n",muonpos.x,muonpos.y,muonpos.z);
#endif
   
    for (angn = 0; angn < in_nangles; ++angn)
    {
        B[angn] = vec3_zero();
        BLor[angn] = vec3_zero();
        pile_init(&MCont[angn],nnn_for_cont);
    }
    
    for (i = 0; i < scx; ++i)
    {
        for (j = 0; j < scy; ++j)
        {
            for (k = 0; k < scz; ++k)
            {
                /* loop over atoms */
                for (a = 0; a < in_natoms; ++a)
                {
                    
                    /* atom position in reduced coordinates */
                    atmpos.x = ( in_positions[3*a] + (float) i) / (float) scx;
                    atmpos.y = ( in_positions[3*a+1] + (float) j) / (float) scy;
                    atmpos.z = ( in_positions[3*a+2] + (float) k) / (float) scz;
                    
                    /* go to cartesian coordinates (in Angstrom!) */
                    atmpos = mat3_vmul(atmpos,sc_lat);
                    
                    /*printf("atompos: %e %e %e\n", atmpos.x, atmpos.y, atmpos.z); */
                    /* difference between atom pos and muon pos (cart coordinates) */
                    
                    r = vec3_sub(atmpos,muonpos);
                    
                    n = vec3_norm(r);
                    if (n < radius)
                    {
                        /* calculate magnetic moment */
#ifdef _ALTERNATE_FC_INPUT
                        printf("ERROR!!! If you see this in the Python extension something went wrong!\n");
                         sk.x = in_fc[6*a];   sk.y = in_fc[6*a+1]; sk.z = in_fc[6*a+2];
                        isk.x = in_fc[6*a+3];isk.y = in_fc[6*a+4];isk.z = in_fc[6*a+5];
#else
                         sk.x = in_fc[6*a];   sk.y = in_fc[6*a+2]; sk.z = in_fc[6*a+4];
                        isk.x = in_fc[6*a+1];isk.y = in_fc[6*a+3];isk.z = in_fc[6*a+5];
#endif
                          phi = in_phi[a];
                        /*printf("sk = %e %e %e\n", sk.x, sk.y, sk.z); */
                        /*printf("isk = %e %e %e\n", isk.x, isk.y, isk.z); */
                        
                        
                        R.x = (float) i; R.y = (float) j; R.z = (float) k; 
                        
                        c = cos ( 2.0*M_PI * (vec3_dot(K,R) + phi));
                        s = sin ( 2.0*M_PI * (vec3_dot(K,R) + phi));
                        
                        m = vec3_zero();
                        m = vec3_add ( vec3_muls(c, sk), m);
                        m = vec3_add ( vec3_muls(s, isk), m);
                        
                        
                        
                        /*printf("I sum: r = %e, p = %e %e %e\n",n, r.x, r.y, r.z); */
                        /*printf("I sum: m = %e %e %e\n", m.x, m.y, m.z); */
                        /* sum it */
                        /* B += (( 3.0 * np.dot(nm,atom[1]) * atom[1] - nm ) / atom[0]**3)*0.9274009 */
                        
                        /* unit vector */
                        u = vec3_muls(1.0/n,r);
                        onebrcube = 1.0/pow(n,3);
                        
                        /* do the rotation */
                        for (angn = 0; angn < in_nangles; ++angn)
                        {
                            angle = 2*M_PI*((float) angn/ (float) in_nangles);
                            rmat = mat3_aangle(axis, angle);

#ifdef _DEBUG
                            printf("Rotation matrix is: %e %e %e\n",rmat.a.x,rmat.a.y,rmat.a.z);
                            printf("Rotation matrix is: %e %e %e\n",rmat.b.x,rmat.b.y,rmat.b.z);
                            printf("Rotation matrix is: %e %e %e\n",rmat.c.x,rmat.c.y,rmat.c.z);
#endif

                            /* rotate moment */
                            rm = mat3_mulv(rmat, m);
                            
                            BLor[angn] = vec3_add(BLor[angn],rm);
                            
                            B[angn] = vec3_add(
                                        B[angn],
                                        vec3_muls( onebrcube ,vec3_sub(vec3_muls(3.0*vec3_dot(rm,u),u), rm))
                                    );

							/* Calculate Contact Field */
							if (n < cont_radius) {
#ifdef _DEBUG                      
								printf("Adding moment to Cont: n: %e, m: %e %e %e! (Total: %d)\n", n, rm.x,rm.y,rm.z,nnn_for_cont);
#endif							
								pile_add_element(&MCont[angn], pow(n,CONT_SCALING_POWER), vec3_muls(1./pow(n,CONT_SCALING_POWER),rm));  /* see ass.c for this line */
							}
                                    
                        }
#ifdef _DEBUG               
                        for (angn = 0; angn < in_nangles; ++angn)
                            printf("B %d is now : %e %e %e\n", angn, B[angn].x, B[angn].y, B[angn].z);
#endif                        
                    }                    

                }
                
            }
        }
    }
    
    /* Lorentz Field (explanation of the numbers in ass.c) */
    

    for (angn = 0; angn < in_nangles; ++angn)
    {
		BLor[angn] = vec3_muls(0.33333333333*11.654064, vec3_muls(3./(4.*M_PI*pow(radius,3)),BLor[angn]));/* to tesla units */
		/*printf("The Lorents field contribution is: %e %e %e Tesla and it will NOT be added!!\n",BLor[angn].x,BLor[angn].y,BLor[angn].z);     */

        out_field_lor[3*angn+0] = BLor[angn].x;
        out_field_lor[3*angn+1] = BLor[angn].y;
        out_field_lor[3*angn+2] = BLor[angn].z;
    }
    free(BLor);    
    
    
    /* Dipolar Field */
    
    for (angn = 0; angn < in_nangles; ++angn)
    {
        B[angn] = vec3_muls(0.9274009, B[angn]); /* to tesla units */
        /*B[angn] = vec3_add(B[angn], BLor); */
        out_field_dip[3*angn+0] = B[angn].x;
        out_field_dip[3*angn+1] = B[angn].y;
        out_field_dip[3*angn+2] = B[angn].z;
    }
    free(B);
    
    /* Contact Filed */
    
    /* Contact Field */
    NofM = 0;
    SumOfWeights = 0;
    
    for (angn = 0; angn < in_nangles; ++angn)
    {
		/* (re) initialize */
		BCont = vec3_zero();
		NofM = 0; /* Number of moments considered */
		SumOfWeights = 0;
		
		for (i=0; i < nnn_for_cont; i++) {
			if (MCont[angn].ranks[i] >= 0.0) {
				BCont = vec3_add(BCont, MCont[angn].elements[i]);
				SumOfWeights += (1./MCont[angn].ranks[i]);
				NofM++;
			}
		}
		
		pile_free(&MCont[angn]);
		
		
		if (NofM >0) {
			BCont = vec3_muls((1./SumOfWeights) * 7.769376 , BCont);
		} /* otherwise is zero anyway! */
		
		out_field_cont[3*angn+0] = BCont.x;
		out_field_cont[3*angn+1] = BCont.y;
		out_field_cont[3*angn+2] = BCont.z; 
	}   
}					 


