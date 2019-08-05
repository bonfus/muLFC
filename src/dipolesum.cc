/**
 * @file dipolesum.cc
 * @author Pietro Bonfa
 * @date 18 March 2018
 * @brief Dipolar field calculator
 *
 */

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <Eigen/Geometry> 
 
#include <stdio.h>
#include <vector>
#include "types.h"
#include "pile.h"
#include "config.h"
#include "lattice.h"

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

Vec3 EMPTY;

/**
    Returns magnetic moments from Fourier Components given for the cell
    in position R given

    @param digit the single digit to encode.
    @return a bar code of the digit using "|" as the long
    bar and "," as the half bar.
*/
void FCtoMagMom(const CMatX& FC, const Vec3 KCart, const VecX& phi, const Vec3 R, RefMatX m)
{
    unsigned int a;     /* counter for atoms */
    T c, s, KdotR;      /*cosine and sine of K.R */

    KdotR =  KCart.dot(R); // R \dot K (done in cartesian coordinates)
    for (a = 0; a < m.cols(); ++a)
    {
        c = cos (KdotR + 2.0*M_PI*phi(a) );
        s = sin (KdotR + 2.0*M_PI*phi(a) );

        m.col(a) = c*FC.col(a).real() + s*FC.col(a).imag();
    }
}

/* Never really used... */
void  DipoleSum(Lattice latt,
          const Vec3& muonpos, unsigned int scx, unsigned int scy, unsigned int scz,
          const double radius, const unsigned int nnn_for_cont, const double cont_radius,
          RefVec3 BC, RefVec3 BD, RefVec3 BL)
{

     /*supercell sizes */
    unsigned int i,j,k; /* counters for supercells */

    MatX atomposCart;
    IVecX occupations;
    MatX atomMoms;
    Vec3 muonposCart;
    Vec3 R;
    Vec3 KCart;
    Vec3 r;
    Vec3 r_mu;
    Vec3 m;   /*magnetic moment of atom */
    Vec3 u;   /* unit vector */
    Vec3 aux;

    Mat3 lattice;
    Mat3 recLattice;

    T n;   /* contains norm of vectors */
    T onebrcube; /* 1/r^3 */

    pile MCont(nnn_for_cont);

    unsigned int a;     /* counter for atoms */
	  int NofM = 0; /* Number of moments considered */
	  double SumOfWeights = 0;
    bool have_moments_for_cell = false;

    latt.GetCell(lattice);
    
    /* Calculate reciprocal lattice */
    recips(lattice, recLattice);

    /* Propagation vector in Cartesian coordinates */
    Crys2Cart(recLattice, latt.K, KCart, false);

    // Get the position of the central cell
    aux << (T) (scx/2), (T) (scy/2), (T) (scz/2);
    Crys2Cart(lattice, aux, R, false);

    /* muon position: reduced coordinates -> Cartesian Coordinates */
    Crys2Cart(lattice, muonpos, muonposCart, false);

    /* Cartesian coordinates of the muon moved in the central cell of the supercell */
    muonposCart = muonposCart + R;

    /* === Atoms data === */
    atomposCart.resize(3, latt.nAtoms);
    atomMoms.resize(3, latt.nAtoms);
    occupations.resize(latt.nAtoms);

    /* Atomic positions: reduced coordinates -> Cartesian Coordinates */
    Crys2Cart(lattice, latt.atomPosFrac, atomposCart, false);

    /* === Fields data === */
    BC = Vec3::Zero();
    BD = Vec3::Zero();
    BL = Vec3::Zero();

    /* === LOOP over SUPERCELL === */
    for (i = 0; i < scx; ++i)
    {
        r.x() = ( (T) i );
        for (j = 0; j < scy; ++j)
        {
            r.y() = ( (T) j );
            for (k = 0; k < scz; ++k)
            {
                r.z() = ( (T) k );

                /* Radius from origin of cell position under consideration (in Cartesian coordinates) */
                Crys2Cart(lattice, r, R, false);

                have_moments_for_cell = false;
                /* loop over atoms */
                for (a = 0; a < atomposCart.cols(); ++a)
                {

                    /* atom position in cartesian coordinates */
                    aux = R + atomposCart.col(a);
                    r_mu = aux - muonposCart;

                    /* norm of the distance vector between the muon and the atom */
                    n = r_mu.norm();

                    if (n < radius)
                    {
                        /* calculate magnetic moment for cell at R */
                        if ( not have_moments_for_cell ) {
                            latt.MaterializeOccupationsInCell(occupations);
                            //std::cout << "Occupations are: " << occupations.transpose() << std::endl;
                            FCtoMagMom(latt.FC, KCart, latt.Phi, R, atomMoms);
                            have_moments_for_cell = true;
                        }
                        if (occupations(a) == 0) continue;

                        m = atomMoms.col(a);
                        if (m.norm() <= EPS) continue;

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

                    } /* if inside sphere */
                } /* for each atom in unitcell */
            } /* nz */
        } /* ny */
    } /* nx */

    /* === Lorentz Field === */
    /*  1 bohr_magneton/(1angstrom^3) = 9274009.5(amperes ∕ meter)
     *   mu_0 = 0.0000012566371((meter tesla) ∕ ampere)
     *   BLor = (mu_0/3)*M_Lor
     *   Note that befor this line BLor is just the sum of the magnetic moments!
     *      magnetic_constant * 1 bohr_magneton = 11.654064 T⋅Å^3
     */

    BL = 0.33333333333*11.654064 * (3./(4.*M_PI*pow(radius,3))) * BL;
    /*printf("The Lorentz field contribution is: %e %e %e Tesla!!\n",BLor.x,BLor.y,BLor.z); */


    /* === Contact Field === */
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

    /* === Dipolar Field === */
    /* mu_0/4pi = 0.1E-6(newton ∕ ampere^2) = 0.1E-6((meter tesla) ∕ ampere) */
    BD = 0.92740098*BD; /* to tesla units */

}


void  DipoleSumMany(Lattice latt,
          const MatX& muonpos, unsigned int scx, unsigned int scy, unsigned int scz,
          const double radius, const unsigned int nnn_for_cont, const double cont_radius,
          const Vec3& axis, const int nAngles,
          RefMatX BC, RefMatX BD, RefMatX BL)
{

     /*supercell sizes */
    unsigned int i,j,k; /* counters for supercells */

    MatX atomposCart;
    IVecX occupations;
    MatX atomMoms;
    MatX muonposCart;
    Vec3 R;
    Vec3 KCart;
    Vec3 r;
    Vec3 r_mu;
    Vec3 m;   /*magnetic moment of atom */
    Vec3 u;   /* unit vector */
    Vec3 aux;

    Mat3 lattice;
    Mat3 recLattice;
    Mat3 unitRotMat;

    T n;   /* contains norm of vectors */
    T onebrcube; /* 1/r^3 */


    unsigned int a;     /* counter for atoms */
    int NofM = 0; /* Number of moments considered */
	double SumOfWeights = 0;
    int nMuons = muonpos.cols();
    int which_index = 0;

    
    // Set data dimensions
    muonposCart.resizeLike(muonpos);
    
    pile ** MCont = new pile*[nMuons];
    for (int i=0; i < nMuons; i++) {
        MCont[i] = new pile(nnn_for_cont);
    }
    bool have_moments_for_cell = false;

    unitRotMat = AngleAxisd(2.0*M_PI/nAngles, axis);

    latt.GetCell(lattice);
    
    /* Calculate reciprocal lattice */
    recips(lattice, recLattice);

    /* Propagation vector in Cartesian coordinates */
    Crys2Cart(recLattice, latt.K, KCart, false);

    // Cell in the middle
    aux << (T) (scx/2), (T) (scy/2), (T) (scz/2);
    Crys2Cart(lattice, aux, R, false);

    /* muon position: reduced coordinates -> Cartesian Coordinates */
    Crys2Cart(lattice, muonpos, muonposCart, false);

    /* Cartesian coordinates in the center of the supercell */
    muonposCart.colwise() += R;

    /* === Atoms data === */
    atomposCart.resize(3, latt.nAtoms);
    atomMoms.resize(3, latt.nAtoms);
    occupations.resize(latt.nAtoms);

    /* Atomic positions: reduced coordinates -> Cartesian Coordinates */
    Crys2Cart(lattice, latt.atomPosFrac, atomposCart, false);

    /* === Fields data === */
    BC.setZero();
    BD.setZero();
    BL.setZero();

    /* === LOOP over SUPERCELL === */
    for (i = 0; i < scx; ++i)
    {
        r.x() = ( (T) i );
        for (j = 0; j < scy; ++j)
        {
            r.y() = ( (T) j );
            for (k = 0; k < scz; ++k)
            {
                r.z() = ( (T) k );

                /* Cell position in Cartesian coordinates */
                Crys2Cart(lattice, r, R, false);


                have_moments_for_cell = false;
                /* loop over atoms */
                for (a = 0; a < atomposCart.cols(); ++a)
                {
                    /* calculate magnetic moment for cell at R */
                    if ( not have_moments_for_cell ) {
                        //latt.MaterializeOccupationsInCell(occupations);
                        FCtoMagMom(latt.FC, KCart, latt.Phi, R, atomMoms);
                        have_moments_for_cell = true;
                    }
                    //if (occupations(a) == 0) continue;
                    
                    aux = R + atomposCart.col(a);

                    /* loop over muons */
                    for (int imu = 0; imu < nMuons; imu++)
                    {
                        /* atom position in cartesian coordinates */
                        r_mu = aux - muonposCart.col(imu);

                        /* norm of the distance vector between the muon and the atom */
                        n = r_mu.norm();

                        if (n < radius)
                        {
                            m = atomMoms.col(a);
                            if (m.norm() <= EPS) continue;
                            
                            for (int iang = 0; iang < nAngles; iang++) {
    
                                /* Calculate Lorentz Field */
                                BL.col(iang + imu*nAngles) += m;

                                /* unit vector */
                                u = r_mu.normalized();
                                onebrcube = 1.0/pow(n,3);

                                BD.col(iang + imu*nAngles) += onebrcube * ( u*3.0*m.dot(u) -  m );
                                
                                /* increment by rotation */
                                if(nAngles > 1) m = unitRotMat*m;
                                
                            } /* loop on angles */

                            /* Store contact field once, to be rotated later  */
                            if (n < cont_radius) {
                                MCont[imu]->add_element(pow(n,CONT_SCALING_POWER), (1./pow(n,CONT_SCALING_POWER))*m);
                            }

                        } /* if inside sphere */
                    } /* for each muon */
                } /* for each atom in unitcell */
            } /* nz */
        } /* ny */
    } /* nx */


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
        for (int iang = 0; iang < nAngles; iang++) {
            NofM = 0; /* Number of moments considered */
            SumOfWeights = 0;
            for (i=0; i < nnn_for_cont; i++) {
                    if (MCont[imu]->ranks(i) > 0.0) {
                        BC.col(iang + imu*nAngles) += MCont[imu]->elements.col(i);
                        SumOfWeights += (1./MCont[imu]->ranks(i)); /* Add the contribution weighted as norm^3 to the total */
                        NofM++;
                    }
            }
            /* (2 magnetic_constant/3)⋅1bohr_magneton   = ((2 ⋅ magnetic_constant) ∕ 3) ⋅ (1 ⋅ bohr_magneton)
            *   ≈ 7.769376E-27((g⋅m^3) ∕ (A⋅s^2))
            *   ≈ 7.769376 T⋅Å^3
            */
    
            if (NofM >0) {
                BC.col(iang + imu*nAngles) = (1./SumOfWeights) * 7.769376 * BC.col(iang + imu*nAngles);
            } /* otherwise is zero anyway! */
            
            /* now rotate everything */
            for (i=0; i < nnn_for_cont; i++) {
                m = MCont[imu]->elements.col(i);
                MCont[imu]->elements.col(i) = unitRotMat*m;
            }
        }
        delete(MCont[imu]);
    }
    delete[] MCont;
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
void DT(MatX& atomicPositions,  MatX& muonpos,
          int scx, int scy, int scz, 
          Mat3 lattice, const double radius,
          MatX& DT) 
{

    unsigned int i,j,k,a; /* counters for supercells */
    
    MatX atomposCart;
    Vec3 muonposCart;
    Vec3 r, aux, R, r_mu;
    Mat3 sc_lat;
    
    double n;
    double onebrcube; /* 1/r^3 */
    double onebrfive; /* 1/r^5 */
    
    
    Mat3 A, D;
    
    unsigned int atom;     /* counter for atoms */

    
    aux << (double) scx, (double) scy, (double) scz;
    sc_lat = lattice * aux.asDiagonal();

    // Get the position of the central cell
    aux << (T) (scx/2), (T) (scy/2), (T) (scz/2);
    Crys2Cart(lattice, aux, R, false);

    /* muon position: reduced coordinates -> Cartesian Coordinates */
    Crys2Cart(lattice, muonpos, muonposCart, false);

    /* Cartesian coordinates of the muon moved in the central cell of the supercell */
    muonposCart = muonposCart + R;

    /* === Atoms data === */
    atomposCart.resizeLike(atomicPositions);


    /* Atomic positions: reduced coordinates -> Cartesian Coordinates */
    Crys2Cart(lattice, atomicPositions, atomposCart, false);


    A.setZero();
    D.setZero();


    /* === LOOP over SUPERCELL === */
    for (i = 0; i < scx; ++i)
    {
        r.x() = ( (T) i );
        for (j = 0; j < scy; ++j)
        {
            r.y() = ( (T) j );
            for (k = 0; k < scz; ++k)
            {
                r.z() = ( (T) k );
                
                /* Radius from origin of cell position under consideration (in Cartesian coordinates) */
                Crys2Cart(lattice, r, R, false);

                for (a = 0; a < atomposCart.cols(); ++a)
                {


                    /* atom position in cartesian coordinates */
                    aux = R + atomposCart.col(a);
                    r_mu = aux - muonposCart;
                    
                    
                    /* norm of the distance vector between the muon and the atom */
                    n = r_mu.norm();

                    if (n < radius)
                    {


                        /* vector */
                        onebrcube = 1.0/pow(n,3);
                        onebrfive = 1.0/pow(n,5);
                        

                        /* See uSR bible (Yaouanc Dalmas De Reotier, page 81) */
                        /* alpha = x */
                        D.col(0).x() = -onebrcube+3.0*r_mu.x()*r_mu.x()*onebrfive;
                        D.col(0).y() = 3.0*r_mu.x()*r_mu.y()*onebrfive;
                        D.col(0).z() = 3.0*r_mu.x()*r_mu.z()*onebrfive;
                        
                        /* alpha = y */
                        D.col(1).x() = D.col(0).y();
                        D.col(1).y() = -onebrcube+3.0*r_mu.y()*r_mu.y()*onebrfive;
                        D.col(1).z() = 3.0*r_mu.y()*r_mu.z()*onebrfive;
                        
                        /* alpha = z */
                        D.col(2).x() = D.col(0).z();
                        D.col(2).y() = D.col(1).z(); /* this is pretty stupid */
                        D.col(2).z() = -onebrcube+3.0*r_mu.z()*r_mu.z()*onebrfive;
                        
                        A += D ;

                        
                    }
                }
            }
        }
    }

    /* B = vec3_muls(0.9274009, B); // we should multiply for a volume */
    DT=A;

}


// ----------------
// Compatibility interface
// ----------------

namespace py = pybind11;

py::tuple  Fields(std::string s,  MatX& atomicPositions,  CMatX& FC, 
              MatX& K,  VecX& Phi,  
           MatX& muonpos, const IVec3& sc,  Mat3& unitCell,
          const double radius, const unsigned int nnn_for_cont, const double cont_radius,
          const int nangles=1, Vec3& axis=EMPTY)
{
    
    
    atomicPositions.transposeInPlace();
    FC.transposeInPlace();
    K.transposeInPlace();
    //Phi.transposeInPlace();
    //std::cout << "nMuons before"<< muonpos.cols() << std::endl;
    if (muonpos.cols() > 1) muonpos.transposeInPlace();
    unitCell.transposeInPlace();
    
    
    Lattice l = Lattice(unitCell, atomicPositions, Phi, FC, K);
    
    // Set data dimensions
    int nMuons = muonpos.cols();
    //std::cout << "nMuons "<< nMuons << std::endl;
    
    
    //std::cout << "atomicPositions "<< atomicPositions.rows() << " " << atomicPositions.cols() << std::endl;

    MatX BD(3, nMuons*nangles), BL(3, nMuons*nangles), BC(3, nMuons*nangles);
    if (s=="s" || s=="simple")
    {
        DipoleSumMany(l, muonpos, sc(0), sc(1), sc(2), radius,  nnn_for_cont, cont_radius,
                        axis, nangles,
                        BC,  BD,  BL);
    }

    if (s=="r" || s=="rotate")
    {
        DipoleSumMany(l, muonpos, sc(0), sc(1), sc(2), radius,  nnn_for_cont, cont_radius,
                        axis, nangles,
                        BC,  BD,  BL);
    }


    if (s=="i" || s=="inc") {
        for (int i=0; i<nMuons; i++)
        {
            FastIncom(atomicPositions, FC, K, Phi, muonpos.col(i),  sc(0), sc(1), sc(2),
                    unitCell, radius, nnn_for_cont, cont_radius, nangles,
                    BC.block(0, i * nangles, 3,  nangles ),
                    BD.block(0, i * nangles, 3,  nangles ),
                    BL.block(0, i * nangles, 3,  nangles ));            
        }
    }
    
    BC.transposeInPlace();
    BD.transposeInPlace();
    BL.transposeInPlace();

    return py::make_tuple ( BC, BD, BL);
}

py::tuple Simple( MatX& atomicPositions,  CMatX& FC, 
              MatX& K,  VecX& Phi,  
           MatX& muonpos, const IVec3& sc,  Mat3& unitCell,
          const double radius, const unsigned int nnn_for_cont, const double cont_radius)
{

    atomicPositions.transposeInPlace();
    FC.transposeInPlace();
    K.transposeInPlace();
    //Phi.transposeInPlace();

    if (muonpos.cols() > 1) muonpos.transposeInPlace();
    unitCell.transposeInPlace();
    
    Lattice l = Lattice(unitCell, atomicPositions, Phi, FC, K);
    

    
    // Set data dimensions
    int nMuons = muonpos.cols();
    
    MatX BD(3, nMuons), BL(3, nMuons), BC(3, nMuons);
    
    Vec3 axis;
    DipoleSumMany(l, muonpos, sc(0), sc(1), sc(2), radius,  nnn_for_cont, cont_radius,
                    axis, 1, BC,  BD,  BL);

    BC.transposeInPlace();
    BD.transposeInPlace();
    BL.transposeInPlace();
    return py::make_tuple ( BC, BD, BL);
}

Mat3 DipolarTensor( MatX& atomicPositions,  
           MatX& muonpos, const IVec3& sc,  Mat3& unitCell,
          const double radius)
{

    atomicPositions.transposeInPlace();
    unitCell.transposeInPlace();
    if (muonpos.cols() > 1) std::cout << "nMuons> 1 not implemented "<<std::endl;
    
    MatX DTr(3,3);
    
    DT(atomicPositions, muonpos, sc(0), sc(1), sc(2), unitCell, radius, DTr);

    DTr.transposeInPlace();

    return DTr;
}

// ----------------
// Python interface
// ----------------

using namespace pybind11::literals;

void init_dipolesum(py::module &m) {
    m.def("DipoleSum", &DipoleSum);
    m.def("DipoleSumMany", &DipoleSumMany);
    m.def("Fields", &Fields, 
            "Gets you fields!", 
            "s"_a,  "atomicPositions"_a, 
            "FC"_a, "K"_a, "Phi"_a, 
            "muonpos"_a, "sc"_a, "unitCell"_a, 
            "radius"_a, "nnn_for_cont"_a, 
            "cont_radius"_a, 
            "nangles"_a=1, "axis"_a=EMPTY);
    m.def("Simple", &Simple);
    m.def("DipolarTensor", &DipolarTensor,py::return_value_policy::copy);
}

///////// pybind11::gil_scoped_release release;
///////// 
///////// while (true)
///////// {
/////////     // do something and break
///////// }
///////// 
///////// pybind11::gil_scoped_acquire acquire;
