/**
 * @file dipolesum.cc
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
    
    pile MCont(nnn_for_cont);
       
    
    unsigned int a;     /* counter for atoms */
	  int NofM = 0; /* Number of moments considered */
	  double SumOfWeights = 0;
    
    
    aux << (T) scx, (T) scy, (T) scz;
    sc_lat = lattice * aux.asDiagonal();
        
    /* muon position in reduced coordinates */
    aux.x() =  (muonpos.x() + (scx/2) ) / (T) scx;
    aux.y() =  (muonpos.y() + (scy/2) ) / (T) scy;
    aux.z() =  (muonpos.z() + (scz/2) ) / (T) scz;

//    std::cout << "Muon pos (frac): " << aux.transpose() << std::endl;

    
    Crys2Cart(sc_lat, aux, muonposCart, false);

//    std::cout << "Muon pos (cart): " << muonposCart.transpose() << std::endl;


    BC = Vec3::Zero();
    BD = Vec3::Zero();
    BL = Vec3::Zero();
    
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
                            MCont.add_element(pow(n,CONT_SCALING_POWER), (1./pow(n,CONT_SCALING_POWER))*m);
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

void  DipoleSumMany(const MatX& atomicPositions, 
          const CMatX& FC, const Vec3 K, const VecX& phi,
          const MatX& muonpos, unsigned int scx, unsigned int scy, unsigned int scz,
          const Mat3& lattice, 
          const double radius, const unsigned int nnn_for_cont, const double cont_radius,
          RefMatX BC, RefMatX BD, RefMatX BL)
{

     /*supercell sizes */
    unsigned int i,j,k; /* counters for supercells */
    
    MatX atomposRed;
    MatX muonposCart;
    Vec3 atomposCart;
    Vec3 R;
    Vec3 r;
    Vec3 r_mu;
    Vec3 m;   /*magnetic moment of atom */
    Vec3 u;   /* unit vector */
    MatX aux;
    Vec3 vaux;
        
    Mat3 sc_lat;
    
    T n;   /* contains norm of vectors */
    T c,s; /*cosine and sine of K.R */
    T onebrcube; /* 1/r^3 */
       
    
    unsigned int a;     /* counter for atoms */
	  int NofM = 0; /* Number of moments considered */
	  double SumOfWeights = 0;
    int nMuons = muonpos.cols();
    pile ** MCont = new pile*[nMuons];
    for (int i=0; i < nMuons; i++) {
        MCont[i] = new pile(nnn_for_cont);
    }
    
    std::cout << nMuons << std::endl;
    aux.resize(3, nMuons);
    muonposCart.resize(3, nMuons);
    
    vaux << (T) scx, (T) scy, (T) scz;
    sc_lat = lattice * vaux.asDiagonal();
        
    /* muon position in reduced coordinates */
    for (int i=0; i < nMuons; i++) {
        aux.col(i).x() =  (muonpos.col(i).x() + (scx/2) ) / (T) scx;
        aux.col(i).y() =  (muonpos.col(i).y() + (scy/2) ) / (T) scy;
        aux.col(i).z() =  (muonpos.col(i).z() + (scz/2) ) / (T) scz;
    }

//    std::cout << "Muon pos (frac): " << aux.transpose() << std::endl;

    
    Crys2Cart(sc_lat, aux, muonposCart, false);

//    std::cout << "Muon pos (cart): " << muonposCart.transpose() << std::endl;


    BC.setZero();
    BD.setZero();
    BL.setZero();
    
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
                
                R.x() = (T) i; R.y() = (T) j; R.z() = (T) k;
                /* loop over atoms */
                for (a = 0; a < atomposRed.cols(); ++a)
                {
                    
                    /* atom position in reduced coordinates */
                    vaux = atomposRed.col(a) + r;
                
                    /* go to cartesian coordinates (in Angstrom!) */
                    Crys2Cart(sc_lat, vaux, atomposCart, false);

                    /* calculate magnetic moment */
                    
                    
                    
                    c = cos ( 2.0*M_PI * (K.dot(R) + phi(a) ));
                    s = sin ( 2.0*M_PI * (K.dot(R) + phi(a) ));
                    
                    m = c*FC.col(a).real() + s*FC.col(a).imag();
                    
                    /*printf("atompos: %e %e %e\n", atmpos.x, atmpos.y, atmpos.z); */
                    /* difference between atom pos and muon pos (cart coordinates) */
                    for (int imu = 0; imu < nMuons; imu++)
                    {
                        r_mu = atomposCart - muonposCart.col(imu);
                        
                        n = r_mu.norm();
                        if (n < radius)
                        {

                        
                            /* Calculate Lorentz Field */
                            BL.col(imu) += m;
                            //std::cout << "BL " << BL.col(imu).transpose() << std::endl;
                            /* Calculate Contact Field */
                            if (n < cont_radius) {
                                MCont[imu]->add_element(pow(n,CONT_SCALING_POWER), (1./pow(n,CONT_SCALING_POWER))*m);
                            }
                        
                            /* unit vector */
                            u = r_mu.normalized();
                            onebrcube = 1.0/pow(n,3);

                            BD.col(imu) += onebrcube * ( u*3.0*m.dot(u) -  m );
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
    
    for (int imu = 0; imu < nMuons; imu++)
    {
        NofM = 0; /* Number of moments considered */
        SumOfWeights = 0;
        for (i=0; i < nnn_for_cont; i++) {
		        if (MCont[imu]->ranks(i) > 0.0) {
		        	  BC.col(imu) += MCont[imu]->elements.col(i);
		        	  SumOfWeights += (1./MCont[imu]->ranks(i)); /* Add the contribution weighted as norm^3 to the total */
		        	  NofM++;
		        }
	      }
        /* (2 magnetic_constant/3)⋅1bohr_magneton   = ((2 ⋅ magnetic_constant) ∕ 3) ⋅ (1 ⋅ bohr_magneton)
         *   ≈ 7.769376E-27((g⋅m^3) ∕ (A⋅s^2))
         *   ≈ 7.769376 T⋅Å^3
         */
        
	      if (NofM >0) {
	      	BC.col(imu) = (1./SumOfWeights) * 7.769376 * BC.col(imu);
	      } /* otherwise is zero anyway! */
        delete(MCont[imu]);
    }
    delete[] MCont;
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
    
    pile MCont(nnn_for_cont);
       
    
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
                            MCont.add_element(pow(n,CONT_SCALING_POWER), (1./pow(n,CONT_SCALING_POWER))*m);
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

void  FastIncom(const MatX& atomicPositions, 
          const CMatX& FC, const Vec3 K, const VecX& phi,
          const Vec3& muonpos, unsigned int scx, unsigned int scy, unsigned int scz,
          const Mat3& lattice, 
          const double radius, const unsigned int nnn_for_cont, const double cont_radius,
          const unsigned int in_nangles,
          RefMatX BC, RefMatX BD, RefMatX BL)
{
  
    Vec3 atomposCart;
    Vec3 muonposCart;
    
    Vec3 r;
    Vec3 u;   /* unit vector */
    
    Mat3 recLattice;
    Mat3 sc_lat;
    
    double n;
    double c,s; /*cosine and sine of K.R */
    double onebrcube; /* 1/r^3 */
        
    /* tmp value for speed and clearness */
    T KdotR;
    
    VecX stagmom;                    /* this is m_0 */
    MatX refatmpos;             /* reference atom used to produce C and S     */
    MatX refatomposCart;             /* reference atom used to produce C and S     */
    MatX Ahelix;
    MatX Bhelix;/* two unit vectors describing the helix in the m_0 (cos(phi).a +/- sin(phi).b) */
    MatX SDip;
    MatX CDip; /* sums of contribution providing cosine and sine prefactors */
    MatX SLor;  
    MatX CLor;/* sums of contribution providing cosine and sine prefactors */
    
    pile CCont(nnn_for_cont), SCont(nnn_for_cont);
    
    Vec3 KCart;
    
    int natoms;
    unsigned int i,j,k, a, angn;     /* counter for atoms */
    Vec3 aux;
    T angle = 0;
    Vec3 BDip;
    Vec3 BLor;
    Vec3 BCont;

    /* for contact field evaluation */
    Vec3 CBCont = Vec3::Zero();
    Vec3 SBCont = Vec3::Zero();
    int NofM = 0; /* Number of moments considered */
    double SumOfWeights = 0;
    int l;
    /* initialize variables */

    natoms = atomicPositions.cols();



    stagmom.resize(natoms);                    /* this is m_0 */
    refatmpos.resize(3, natoms);             /* reference atom used to produce C and S     */
    refatomposCart.resize(3, natoms);             /* reference atom used to produce C and S     */
    Ahelix.resize(3, natoms ) ;
    Bhelix.resize(3, natoms );/* two unit vectors describing the helix in the m_0 (cos(phi).a +/- sin(phi).b) */
    SDip.resize(3, natoms );
    CDip.resize(3, natoms ); /* sums of contribution providing cosine and sine prefactors */
    SLor.resize(3, natoms );  
    CLor.resize(3, natoms );/* sums of contribution providing cosine and sine prefactors */

    Ahelix.setZero(3, natoms);
    Bhelix.setZero(3, natoms);
    CDip.setZero(3, natoms);
    SDip.setZero(3, natoms);
    CLor.setZero(3, natoms);
    SLor.setZero(3, natoms);


 
    recips(lattice, recLattice);
    
    
    // go to reciprocal cartesia coordinates.
    Crys2Cart(recLattice, K, KCart, false);

    aux << (double) scx, (double) scy, (double) scz;
    sc_lat = lattice * aux.asDiagonal();
    
    
    /* muon position in reduced coordinates */
    aux.x() =  (muonpos.x() + (scx/2) ) / (float) scx;
    aux.y() =  (muonpos.y() + (scy/2) ) / (float) scy;
    aux.z() =  (muonpos.z() + (scz/2) ) / (float) scz;
    
    Crys2Cart(sc_lat, aux, muonposCart, false);


    for (a = 0; a < natoms; ++a)
    {
        /* reference atom in reduced coordinates */
        /*   the first atom is chosen as reference */
        refatmpos.col(a).x() =  (atomicPositions.col(a).x() + (scx/2) ) / (T) scx;
        refatmpos.col(a).y() =  (atomicPositions.col(a).y() + (scy/2) ) / (T) scy;
        refatmpos.col(a).z() =  (atomicPositions.col(a).z() + (scz/2) ) / (T) scz;



        
        
        /* now take care of magntism */

        stagmom(a) = FC.col(a).real().norm();
        Ahelix.col(a) = FC.col(a).real().normalized();

        

        /* check if real and imaginary parts are the same */
        if (fabs(stagmom(a) - FC.col(a).imag().norm())>EPS)
        {
            printf("ERROR!!! Staggered moment is different in real and imag parts of atom %u\n Use another routine!\n",a);
        }

        /* now B */
        Bhelix.col(a) =  FC.col(a).imag().normalized();

        if (fabs(Ahelix.col(a).dot(Bhelix.col(a))) > EPS)
        {
            printf("ERROR!!! Real and imaginary part of atom %u are not orthogonal by %e!\n",a,Ahelix.col(a).dot(Bhelix.col(a)));
        }
   
        
        if (fabs(phi(a)) > EPS)
        {
            printf("WARNING!!! Phi not completely tested! Double check your results.\n");
        }
        
    }
    Crys2Cart(sc_lat, refatmpos, refatomposCart, false);

/* parallel execution starts here */
/* the shared variables are listed just to remember about data races! */
/* other variable shaed by default: refatmpos,atmpos,phi,Ahelix,Bhelix */
//pragma omp parallel shared(SDip,CDip,SLor,CLor,SCont,CCont,scx,scy,scz,in_positions) 
{
//pragma omp for collapse(3) private(i,j,k,a,r,n,c,s,u,crysvec,onebrcube,atmpos)
    for (i = 0; i < scx; ++i)
    {
        for (j = 0; j < scy; ++j)
        {
            for (k = 0; k < scz; ++k)
            {
                /* loop over atoms */
                for (a = 0; a < natoms; ++a)
                {
                    
                    /* atom position in reduced coordinates */
                    aux.x() = ( atomicPositions.col(a).x() + (T) i) / (T) scx;
                    aux.y() = ( atomicPositions.col(a).y() + (T) j) / (T) scy;
                    aux.z() = ( atomicPositions.col(a).z() + (T) k) / (T) scz;
                    
                    /* go to cartesian coordinates (in Angstrom!) */
                    Crys2Cart(sc_lat, aux, atomposCart, false);
                    
                    /* difference between atom pos and muon pos (cart coordinates) */
                    r = atomposCart - muonposCart;
                    
                    n = r.norm();
                    if (n < radius)
                    {
                        
                        /* unit vector */
                        u = r.normalized();
                        onebrcube = 1.0/pow(n,3);
                        
                        // (atomposCart - muonposCart)  - refatomposCart
                        aux = r - refatomposCart.col(a);
                        
                        KdotR =  KCart.dot(aux); // R \dot K (done in cartesian coordinates)
                        
                        /* */
                        c = cos (KdotR + 2.0*M_PI*phi(a) ) ;
                        s = sin (KdotR + 2.0*M_PI*phi(a) ) ;

                        /* sum all data */
                        //pragma omp critical(dipolar)
                        {
                            /* Dipolar */
                            CDip.col(a) += c * onebrcube * (u*3.0*Ahelix.col(a).dot(u) - Ahelix.col(a) ) 
                                         + s * onebrcube * (u*3.0*Bhelix.col(a).dot(u) - Bhelix.col(a) );

                            SDip.col(a) += s * onebrcube * (u*3.0*Ahelix.col(a).dot(u) - Ahelix.col(a) ) 
                                         - c * onebrcube * (u*3.0*Bhelix.col(a).dot(u) - Bhelix.col(a) );

                        }
                        //pragma omp critical(lorentz)
                        {
                            /* Lorentz */
                            CLor.col(a) += c * Ahelix.col(a) + s*Bhelix.col(a);
                            SLor.col(a) += s * Ahelix.col(a) - c*Bhelix.col(a);

                        }
                        /* Contact */
                        if (n < cont_radius) {
                            //pragma omp critical(contact)
                            {
                                aux = stagmom(a) * (c * Ahelix.col(a) + s * Bhelix.col(a));
                                CCont.add_element(pow(n,CONT_SCALING_POWER), aux);
                                
                                aux = stagmom(a) * (s * Ahelix.col(a) - c * Bhelix.col(a));
                                SCont.add_element(pow(n,CONT_SCALING_POWER), aux);
                            }
                        }
                    }                    
                }
            }
        }
    }
}
    
    angle=0;
    /* for contact field evaluation */
    CBCont = Vec3::Zero();
    SBCont = Vec3::Zero();
    NofM = 0; /* Number of moments considered */
    SumOfWeights = 0;    
    
// pragma omp parallel sections private(angn,angle,BDip,BLor,BCont,i) firstprivate(CBCont,SBCont,NofM,SumOfWeights)
{
    /* first portion, dipolar fields and Lorentz */
    // pragma omp section
    {
        for (angn = 0; angn < in_nangles; ++angn)
        {
            angle = 2*M_PI*((T) angn / (T) in_nangles);
            
            /*  === Dipolar Field === */
            BDip = Vec3::Zero();
            /* loop over atoms */
            for (a = 0; a < natoms; ++a)
            {
                BDip +=  stagmom(a) * (cos(angle) * CDip.col(a) - sin(angle) * SDip.col(a));
            }
            
            BD.col(angn) = 0.9274009 * BDip; /* to tesla units */
            
            
            /*  === Lorentz Field === */
            BLor = Vec3::Zero();
            /* loop over atoms */
            for (a = 0; a < natoms; ++a)
            {
                BLor +=  stagmom(a) * (cos(angle) * CLor.col(a) - sin(angle) * SLor.col(a));
            }
            
            BL.col(angn) = 0.33333333333*11.654064*(3./(4.*M_PI*pow(radius,3))) * BLor;            
        }
    }

    BC.setZero();
    /* second portion, contact fields */
    //pragma omp section
    {
        /*  === Contact Field === */
        
        for (i=0; i < nnn_for_cont; i++) {
            if ((CCont.ranks(i) >= 0.0) && (fabs(CCont.ranks(i) - SCont.ranks(i))<EPS)) {
                CBCont += (1./CCont.ranks(i)) * CCont.elements.col(i);
                SBCont += (1./SCont.ranks(i)) * SCont.elements.col(i);
                SumOfWeights += 1./CCont.ranks(i);
                NofM++;
            } else {
                printf("Something VERY odd ! ranks 1: %e ranks 2: %e \n", CCont.ranks(i) , SCont.ranks(i) );
            }
        }
            
        /* (2 magnetic_constant/3)⋅1bohr_magneton   = ((2 ⋅ magnetic_constant) ∕ 3) ⋅ (1 ⋅ bohr_magneton) */
        /*   ≈ 7.769376E-27((g⋅m^3) ∕ (A⋅s^2)) */
        /*   ≈ 7.769376 T⋅Å^3 */
        
        for (angn = 0; angn < in_nangles; ++angn) {
            
            angle = 2*M_PI*((T) angn / (T) in_nangles);
            
            if (NofM >0) {
                BC.col(angn) = (1./SumOfWeights) * 7.769376 * (cos(angle) * CBCont - sin(angle) * SBCont);
            } /* otherwise is zero anyway! */
        }
    }
} /* end of omp parallel sections */


}
