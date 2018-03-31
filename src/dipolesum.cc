/**
 * @file simplesum.cc
 * @author Pietro Bonfa
 * @date 18 March 2018
 * @brief Dipolar field calculator
 *     
 */
#include <stdio.h>
#include <vector>
#include "types.h"
#include "pile.h"
#include "config.h"
#include "lattice.h"

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

void  DipoleSum(const MatX& atomicPositions, 
          const CMatX& FC, const Vec3 K, const VecX& phi,
          const Vec3& muonpos, unsigned int scx, unsigned int scy, unsigned int scz,
          const Mat3& lattice, 
          const double radius, const unsigned int nnn_for_cont, const double cont_radius,
          RefVec3 BC, RefVec3 BD, RefVec3 BL)
{

     /*supercell sizes */
    unsigned int i,j,k; /* counters for supercells */
    
    MatX atomposRed;
    Vec3 muonposCart;
    Vec3 atomposCart;
    Vec3 R;
    Vec3 r;
    Vec3 r_mu;
    Vec3 m;   /*magnetic moment of atom */
    Vec3 u;   /* unit vector */
    Vec3 aux;
        
    Mat3 sc_lat;
    
    T n;   /* contains norm of vectors */
    T c,s; /*cosine and sine of K.R */
    T onebrcube; /* 1/r^3 */
    
    pile MCont;
       
    
    unsigned int a;     /* counter for atoms */
	  int NofM = 0; /* Number of moments considered */
	  double SumOfWeights = 0;
    
    
    aux << (double) scx, (double) scy, (double) scz;
    sc_lat = lattice * aux.asDiagonal();
        
    /* muon position in reduced coordinates */
    aux.x() =  (muonpos.x() + (scx/2) ) / (double) scx;
    aux.y() =  (muonpos.y() + (scy/2) ) / (double) scy;
    aux.z() =  (muonpos.z() + (scz/2) ) / (double) scz;

//    std::cout << "Muon pos (frac): " << aux.transpose() << std::endl;

    
    Crys2Cart(sc_lat, aux, muonposCart, false);

//    std::cout << "Muon pos (cart): " << muonposCart.transpose() << std::endl;


    BC = Vec3::Zero();
    BD = Vec3::Zero();
    BL = Vec3::Zero();
    pile_init(MCont, nnn_for_cont);
    
    // Calculate scaling of atomic position in fractional description
    atomposRed.resize(3, atomicPositions.cols());
    for (a = 0; a < atomicPositions.cols(); ++a) {
        atomposRed(0,a) = atomicPositions(0,a)/ (T) scx;
        atomposRed(1,a) = atomicPositions(1,a)/ (T) scy;
        atomposRed(2,a) = atomicPositions(2,a)/ (T) scz;
    }
    
    for (i = 0; i < scx; ++i)
    {
        r.x() = ( (T) i) / (T) scx;
        for (j = 0; j < scy; ++j)
        {
            r.y() = ( (T) j) / (T) scy;
            for (k = 0; k < scz; ++k)
            {
                r.z() = ( (T) k) / (T) scz;
                
                /* loop over atoms */
                for (a = 0; a < atomposRed.cols(); ++a)
                {
                    
                    /* atom position in reduced coordinates */
                    aux = atomposRed.col(a) + r;
                
                    /* go to cartesian coordinates (in Angstrom!) */
                    Crys2Cart(sc_lat, aux, atomposCart, false);
                    
                    
                    /*printf("atompos: %e %e %e\n", atmpos.x, atmpos.y, atmpos.z); */
                    /* difference between atom pos and muon pos (cart coordinates) */
                    
                    r_mu = atomposCart - muonposCart;
                    
                    n = r_mu.norm();
                    if (n < radius)
                    {
                        /* calculate magnetic moment */
                        
                        R.x() = (T) i; R.y() = (T) j; R.z() = (T) k; 
                        
                        c = cos ( 2.0*M_PI * (K.dot(R) + phi(a) ));
                        s = sin ( 2.0*M_PI * (K.dot(R) + phi(a) ));
                        
                        m = c*FC.col(a).real() + s*FC.col(a).imag();
						
                        /* Calculate Lorentz Field */
                        BL += m;

                        /* Calculate Contact Field */
                        if (n < cont_radius) {
                            pile_add_element(MCont, pow(n,CONT_SCALING_POWER), (1./pow(n,CONT_SCALING_POWER))*m);
                        }

                        /* unit vector */
                        u = r_mu.normalized();
                        onebrcube = 1.0/pow(n,3);
                        

                        BD += onebrcube * ( u*3.0*m.dot(u) -  m );

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
    
    BL = 0.33333333333*11.654064 * (3./(4.*M_PI*pow(radius,3))) * BL;
    /*printf("The Lorents field contribution is: %e %e %e Tesla!!\n",BLor.x,BLor.y,BLor.z); */

    
    /* Contact Field */
    BC.setZero();
    NofM = 0; /* Number of moments considered */
    SumOfWeights = 0;

    for (i=0; i < nnn_for_cont; i++) {
		  if (MCont.ranks(i) > 0.0) {
		  	BC += MCont.elements.col(i);
		  	SumOfWeights += (1./MCont.ranks(i)); /* Add the contribution weighted as norm^3 to the total */
		  	NofM++;
		  }
	  }
	
	  pile_free(MCont);
	
    /* (2 magnetic_constant/3)⋅1bohr_magneton   = ((2 ⋅ magnetic_constant) ∕ 3) ⋅ (1 ⋅ bohr_magneton)
     *   ≈ 7.769376E-27((g⋅m^3) ∕ (A⋅s^2))
     *   ≈ 7.769376 T⋅Å^3
     */
	
	  if (NofM >0) {
	  	BC = (1./SumOfWeights) * 7.769376 * BC;
	  } /* otherwise is zero anyway! */
    
    /* Dipolar Field */
    /* mu_0/4pi = 0.1E-6(newton ∕ ampere^2) = 0.1E-6((meter tesla) ∕ ampere) */
    BD = 0.92740098* BD; /* to tesla units */

}

void  TransformAndSum(const MatX& atomicPositions, 
          const CMatX& FC, const Vec3 K, const VecX& phi,
          const Vec3& muonpos, unsigned int scx, unsigned int scy, unsigned int scz,
          const Mat3& lattice, 
          const double radius, const unsigned int nnn_for_cont, const double cont_radius,
          const MatX& tr_mat,
          RefMatX BC, RefMatX BD, RefMatX BL)
{

     /*supercell sizes */
    unsigned int i,j,k; /* counters for supercells */
    
    MatX atomposRed;
    Vec3 muonposCart;
    Vec3 atomposCart;
    Vec3 R;
    Vec3 r;
    Vec3 r_mu;
    Vec3 m;   /*magnetic moment of atom */
    Vec3 mt;   /*magnetic moment of atom, transformed */
    Vec3 u;   /* unit vector */
    Vec3 aux;
        
    Mat3 sc_lat, tr_mat_loc;
    
    T n;   /* contains norm of vectors */
    T c,s; /*cosine and sine of K.R */
    T onebrcube; /* 1/r^3 */
    
    pile MCont;
       
    
    unsigned int a;     /* counter for atoms */
	  int NofM = 0; /* Number of moments considered */
	  double SumOfWeights = 0;
        
    aux << (double) scx, (double) scy, (double) scz;
    sc_lat = lattice * aux.asDiagonal();
        
    /* muon position in reduced coordinates */
    aux.x() =  (muonpos.x() + (scx/2) ) / (double) scx;
    aux.y() =  (muonpos.y() + (scy/2) ) / (double) scy;
    aux.z() =  (muonpos.z() + (scz/2) ) / (double) scz;

//    std::cout << "Muon pos (frac): " << aux.transpose() << std::endl;

    
    Crys2Cart(sc_lat, aux, muonposCart, false);

//    std::cout << "Muon pos (cart): " << muonposCart.transpose() << std::endl;


    BC = Vec3::Zero();
    BD = Vec3::Zero();
    BL = Vec3::Zero();
    pile_init(MCont, nnn_for_cont);
    
    // Calculate scaling of atomic position in fractional description
    atomposRed.resize(3, atomicPositions.cols());
    for (a = 0; a < atomicPositions.cols(); ++a) {
        atomposRed(0,a) = atomicPositions(0,a)/ (T) scx;
        atomposRed(1,a) = atomicPositions(1,a)/ (T) scy;
        atomposRed(2,a) = atomicPositions(2,a)/ (T) scz;
    }
    
    for (i = 0; i < scx; ++i)
    {
        for (j = 0; j < scy; ++j)
        {
            for (k = 0; k < scz; ++k)
            {
                r.x() = ( (T) i) / (T) scx;
                r.y() = ( (T) j) / (T) scy;
                r.z() = ( (T) k) / (T) scz;
                
                /* loop over atoms */
                for (a = 0; a < atomposRed.cols(); ++a)
                {
                    
                    /* atom position in reduced coordinates */
                    aux = atomposRed.col(a) + r;
                
                    /* go to cartesian coordinates (in Angstrom!) */
                    Crys2Cart(sc_lat, aux, atomposCart, false);
                    
                    
                    /*printf("atompos: %e %e %e\n", atmpos.x, atmpos.y, atmpos.z); */
                    /* difference between atom pos and muon pos (cart coordinates) */
                    
                    r_mu = atomposCart - muonposCart;
                    
                    n = r_mu.norm();
                    if (n < radius)
                    {
                        /* calculate magnetic moment */
                        
                        R.x() = (T) i; R.y() = (T) j; R.z() = (T) k; 
                        
                        c = cos ( 2.0*M_PI * (K.dot(R) + phi(a) ));
                        s = sin ( 2.0*M_PI * (K.dot(R) + phi(a) ));
                        
                        m = c*FC.col(a).real() + s*FC.col(a).imag();

                        /* Collect moments for Contact Field. Will be rotated later */
                        if (n < cont_radius) {
                            pile_add_element(MCont, pow(n,CONT_SCALING_POWER), (1./pow(n,CONT_SCALING_POWER))*m);
                        }

                        /* unit vector */
                        u = r_mu.normalized();
                        onebrcube = 1.0/pow(n,3);

                        for(int ti=0; ti < tr_mat.cols()/3; ti++){
                                                        
                            mt = tr_mat.block(0, ti*3, 3, 3)*m;
                            
                            /* Calculate Lorentz Field */
                            BL.col(ti) += mt;
                            /* Calculate Dipole Field */
                            BD.col(ti) += onebrcube * ( u*3.0*mt.dot(u) -  mt );
                        }
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
    
    BL = 0.33333333333*11.654064 * (3./(4.*M_PI*pow(radius,3))) * BL;
    /*printf("The Lorents field contribution is: %e %e %e Tesla!!\n",BLor.x,BLor.y,BLor.z); */

    
    /* Contact Field */
    BC.setZero();
    NofM = 0; /* Number of moments considered */
    SumOfWeights = 0;

    for (i=0; i < nnn_for_cont; i++) {
		  if (MCont.ranks(i) > 0.0) {
        for(int ti=0; ti < tr_mat.cols()/3; ti++){
          tr_mat_loc = tr_mat.block(0, ti*3, 3, 3);
          mt = tr_mat_loc * MCont.elements.col(i);
          BC.col(ti) += mt;
        }
        SumOfWeights += (1./MCont.ranks(i)); /* Add the contribution weighted as norm^3 to the total */
        NofM++;
		  }
	  }
	
	  pile_free(MCont);
	
    /* (2 magnetic_constant/3)⋅1bohr_magneton   = ((2 ⋅ magnetic_constant) ∕ 3) ⋅ (1 ⋅ bohr_magneton)
     *   ≈ 7.769376E-27((g⋅m^3) ∕ (A⋅s^2))
     *   ≈ 7.769376 T⋅Å^3
     */
	
	  if (NofM >0) {
	  	BC = (1./SumOfWeights) * 7.769376 * BC;
	  } /* otherwise is zero anyway! */
    
    /* Dipolar Field */
    /* mu_0/4pi = 0.1E-6(newton ∕ ampere^2) = 0.1E-6((meter tesla) ∕ ampere) */
    BD = 0.92740098* BD; /* to tesla units */
    
}

