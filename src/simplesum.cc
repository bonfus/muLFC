/**
 * @file simplesum.cc
 * @author Pietro Bonfa
 * @date 18 March 2018
 * @brief Dipolar field calculator
 *     
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
#include "pile.h"
#include "config.h"
#include "lattice.h"
#include "dipolesum.h"

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

#ifdef _OPENMP
	#include <omp.h>
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
 * @param in_natoms: number of atoms in the lattice.
 * @param inout_nmounpos: input, number of muon pos to be evaluated, output number pos evaluated. If negative, random positions will be generated.
 * @param out_field_cont Contact filed in Tesla in the Cartesian coordinates system defined by in_cell. A coupling of 1 \f$ \mathrm{Ang} ^{-1} \sim 13.912~\mathrm{mol/emu} \f$ is assumed.
 * @param out_field_dip  Dipolar field in Tesla in the Cartesian coordinates system defined by in_cell.
 * @param out_field_lor  Lorentz field in Tesla in the Cartesian coordinates system defined by in_cell.
 */


void  SimpleSum(const T *in_positions, 
          const T *in_fc, const T *in_K, const T *in_phi,
          const T *in_muonpos, const int * in_supercell, const T *in_cell, 
          const T radius, const unsigned int nnn_for_cont, const T cont_radius,
          const T min_radius_from_atoms, 
          unsigned int in_natoms, int inout_nmounpos,
          T *out_field_cont, T *out_field_dip, T *out_field_lor)
{

    unsigned int scx, scy, scz; /*supercell sizes */
    
    MatX atomicPos(3,in_natoms) ;
    Vec3 muonPos;
    Mat3 lattice;
    
    
    
    /* description of the magnetic structure. */
    /* data provided in cartesian coordinates */
    VecX phi(in_natoms) ;
    CMatX FC(3,in_natoms);
    Vec3 K;
    Vec3 B, BLor, BCont;       
    unsigned int a;     /* counter for atoms */    
        
    
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
    
    
    /* muon position in reduced coordinates */
    muonPos.x() =  in_muonpos[0];
    muonPos.y() =  in_muonpos[1];
    muonPos.z() =  in_muonpos[2];

    for (a = 0; a < in_natoms; a++)
    {
        FC.col(a).real() << in_fc[6*a] , in_fc[6*a+2] , in_fc[6*a+4];
        FC.col(a).imag() << in_fc[6*a+1] , in_fc[6*a+3] , in_fc[6*a+5];
        atomicPos.col(a) << in_positions[3*a] , in_positions[3*a+1] , in_positions[3*a+2] ;
        phi(a) = in_phi[a];
    }

    DipoleSum(atomicPos, 
          FC, K, phi,
          muonPos, scx, scy, scz,
          lattice, radius, nnn_for_cont, cont_radius,
          BCont, B, BLor);


    out_field_lor[0] = BLor.x();
    out_field_lor[1] = BLor.y();
    out_field_lor[2] = BLor.z();
    
    
    out_field_cont[0] = BCont.x();
    out_field_cont[1] = BCont.y();
    out_field_cont[2] = BCont.z();
    

    out_field_dip[0] = B.x();
    out_field_dip[1] = B.y();
    out_field_dip[2] = B.z();

}



#if 0

void  SimpleSumOld(const double *in_positions, 
          const double *in_fc, const double *in_K, const double *in_phi,
          const double *in_muonpos, const int * in_supercell, const double *in_cell, 
          const double radius, const unsigned int nnn_for_cont, const double cont_radius,
          const double min_radius_from_atoms, 
          unsigned int in_natoms, int inout_nmounpos,
          double *out_field_cont, double *out_field_dip, double *out_field_lor)
{

    unsigned int scx, scy, scz; /*supercell sizes */
    unsigned int i,j,k; /* counters for supercells */
    
    Vec3 atmpos;
    Vec3 muonpos;
    Vec3 r;
    Vec3 m;   /*magnetic moment of atom */
    Vec3 u;   /* unit vector */
    Vec3 aux;
        
    Mat3 sc_lat;
    
    T n;   /* contains norm of vectors */
    T c,s; /*cosine and sine of K.R */
    T onebrcube; /* 1/r^3 */
    
    
    /* description of the magnetic structure. */
    /* data provided in cartesian coordinates */
    Vec3 sk  ;
    Vec3 isk ;
    double  phi ;
    
    Vec3 R, K, B, BLor;
    pile MCont;
       
    
    unsigned int a;     /* counter for atoms */
	  Vec3 BCont;
	  int NofM = 0; /* Number of moments considered */
	  double SumOfWeights = 0;
    
        
    
    /* define dupercell size */
    scx = in_supercell[0];
    scy = in_supercell[1];
    scz = in_supercell[2];

    
    sc_lat.col(0) << in_cell[0], in_cell[1], in_cell[2];
    sc_lat.col(1) << in_cell[3], in_cell[4], in_cell[5];
    sc_lat.col(2) << in_cell[6], in_cell[7], in_cell[8];
 
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
    aux << (double) scx, (double) scy, (double) scz;
    sc_lat = aux.asDiagonal() * sc_lat;
    
    
    /* muon position in reduced coordinates */
    aux(0) =  (in_muonpos[0] + (scx/2) ) / (double) scx;
    aux(1) =  (in_muonpos[1] + (scy/2) ) / (double) scy;
    aux(2) =  (in_muonpos[2] + (scz/2) ) / (double) scz;

#ifdef _DEBUG
    std::cout << "Muon pos (frac): " << aux.transpose() << std::endl;
#endif
    
    Crys2Cart(sc_lat, aux, muonpos, false);

#ifdef _DEBUG
    std::cout << "Muon pos (cart): " << muonpos.transpose() << std::endl;
#endif


#ifdef _DEBUG
    for (a = 0; a < in_natoms; ++a)
    {
                    
        /* atom position in reduced coordinates */
        atmpos.x =  in_positions[3*a] ;
        atmpos.y =  in_positions[3*a+1] ;
        atmpos.z =  in_positions[3*a+2] ;
        
        printf("Atom pos (crys): %e %e %e\n",atmpos.x,atmpos.y,atmpos.z);
        
        /* go to cartesian coordinates (in Angstrom!) */
        atmpos = mat3_vmul(atmpos,sc_lat);  
        
        printf("Atom pos (cart): %e %e %e\n",atmpos.x,atmpos.y,atmpos.z);
#ifdef _ALTERNATE_FC_INPUT
        printf("ERROR!!! If you see this the extension compilation went wrong!\n");
        printf("FC (real, imag): %e %e %e %e %e %e\n",in_fc[6*a],in_fc[6*a+1],in_fc[6*a+2],in_fc[6*a+3],in_fc[6*a+4],in_fc[6*a+5]);
#else
        printf("FC (real, imag): %e %e %e %e %e %e\n",in_fc[6*a],in_fc[6*a+2],in_fc[6*a+4],in_fc[6*a+1],in_fc[6*a+3],in_fc[6*a+5]);
#endif        
        printf("phi: %e\n",in_phi[a]);
    }
        
#endif


    B = Vec3::Zero();
    BLor = Vec3::Zero();
    pile_init(&MCont, nnn_for_cont);

    

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
                    aux.x() = ( in_positions[3*a] + (double) i) / (double) scx;
                    aux.y()= ( in_positions[3*a+1] + (double) j) / (double) scy;
                    aux.z() = ( in_positions[3*a+2] + (double) k) / (double) scz;
                    
                
                    /* go to cartesian coordinates (in Angstrom!) */
                    Crys2Cart(sc_lat, aux, atmpos, false);
                    
                    
                    /*printf("atompos: %e %e %e\n", atmpos.x, atmpos.y, atmpos.z); */
                    /* difference between atom pos and muon pos (cart coordinates) */
                    
                    r = atmpos - muonpos;
                    
                    n = r.norm();
                    if (n < radius)
                    {
                        /* calculate magnetic moment */
                         sk.x() = in_fc[6*a];   sk.y() = in_fc[6*a+2]; sk.z() = in_fc[6*a+4];
                        isk.x() = in_fc[6*a+1];isk.y() = in_fc[6*a+3];isk.z() = in_fc[6*a+5];

                        phi = in_phi[a];
                        /*printf("sk = %e %e %e\n", sk.x, sk.y, sk.z); */
                        /*printf("isk = %e %e %e\n", isk.x, isk.y, isk.z); */
                        
                        
                        R.x() = (double) i; R.y() = (double) j; R.z() = (double) k; 
                        
                        c = cos ( 2.0*M_PI * (K.dot(R) + phi ));
                        s = sin ( 2.0*M_PI * (K.dot(R) + phi ));
                        
                        m = c*sk + s*isk;
						

                        BLor += m;
                        
                        /* Calculate Contact Field */
                        if (n < cont_radius) {
                            pile_add_element(&MCont, pow(n,CONT_SCALING_POWER), (1./pow(n,CONT_SCALING_POWER))*m);
                        }

                        /* unit vector */
                        u = r.normalized();
                        onebrcube = 1.0/pow(n,3);
                        

                        B += onebrcube * ( u*3.0*m.dot(u) -  m );

                    }                    

                }
                
            }
        }
    }

    /*  1 bohr_magneton/(1angstrom^3) = 9274009.5(amperes ∕ meter)
     *   mu_0 = 0.0000012566371((meter tesla) ∕ ampere)
     *   BLor = (mu_0/3)*M_Lor
     *   Note that befor this line BLor is just the sum of the magnetic moments!
     *      magnetic_constant * 1 bohr_magneton = 11.654064 T⋅Å^3
     */
    
    BLor = 0.33333333333*11.654064 * (3./(4.*M_PI*pow(radius,3))) * BLor;
    /*printf("The Lorents field contribution is: %e %e %e Tesla!!\n",BLor.x,BLor.y,BLor.z); */

    out_field_lor[0] = BLor.x();
    out_field_lor[1] = BLor.y();
    out_field_lor[2] = BLor.z();
    
    /* Contact Field */
    BCont = Vec3::Zero();
    NofM = 0; /* Number of moments considered */
    SumOfWeights = 0;
    for (i=0; i < nnn_for_cont; i++) {
		  if (MCont.ranks[i] > 0.0) {
		  	BCont = BCont + MCont.elements[i];
		  	SumOfWeights += (1./MCont.ranks[i]); /* Add the contribution weighted as norm^3 to the total */
		  	NofM++;
		  }
	  }
	
	  pile_free(&MCont);
	
    /* (2 magnetic_constant/3)⋅1bohr_magneton   = ((2 ⋅ magnetic_constant) ∕ 3) ⋅ (1 ⋅ bohr_magneton)
     *   ≈ 7.769376E-27((g⋅m^3) ∕ (A⋅s^2))
     *   ≈ 7.769376 T⋅Å^3
     */
	
	  if (NofM >0) {
	  	BCont = (1./SumOfWeights) * 7.769376 * BCont;
	  } /* otherwise is zero anyway! */
    
    out_field_cont[0] = BCont.x();
    out_field_cont[1] = BCont.y();
    out_field_cont[2] = BCont.z();
    
    
    /* Dipolar Field */
    /* mu_0/4pi = 0.1E-6(newton ∕ ampere^2) = 0.1E-6((meter tesla) ∕ ampere) */
    B = 0.92740098* B; /* to tesla units */

    out_field_dip[0] = B.x();
    out_field_dip[1] = B.y();
    out_field_dip[2] = B.z();

}
#endif

} /* extern "C" */

