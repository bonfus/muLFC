/**
 * @file simplesum.c
 * @author Pietro Bonfa
 * @date 9 Sep 2016
 * @brief Dipolar field calculator
 *
 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
#include "mat3.h"
#include "pile.h"
#include "config.h"

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

#ifdef _OPENMP
        #include <omp.h>
#endif

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
 * @param in_natoms: number of atoms in the lattice.
 * @param out_field_cont Contact filed in Tesla in the Cartesian coordinates system defined by in_cell. A coupling of 1 \f$ \mathrm{Ang} ^{-1} \sim 13.912~\mathrm{mol/emu} \f$ is assumed.
 * @param out_field_dip  Dipolar field in Tesla in the Cartesian coordinates system defined by in_cell.
 * @param out_field_lor  Lorentz field in Tesla in the Cartesian coordinates system defined by in_cell.
 */
void  SimpleSum(const double *in_positions,
                const double *in_fc, const double *in_K, const double *in_phi,
                const double *in_muonpos, const int * in_supercell, const double *in_cell,
                const double radius, const unsigned int nnn_for_cont, const double cont_radius,
                unsigned int in_natoms,
                double *out_field_cont, double *out_field_dip, double *out_field_lor)
{

    unsigned int scx, scy, scz; /*supercell sizes */
    unsigned int i,j,k; /* counters for supercells */

    vec3 * atmpos;
    vec3 * muonpos;
    vec3 * r;
    vec3 * m;   /*magnetic moment of atom */
    vec3 * u;   /* unit vector */

    mat3 * sc_lat;
    mat3 * sctmp;

    double n;   /* contains norm of vectors */
    double c,s; /*cosine and sine of K.R */
    double onebrcube; /* 1/r^3 */


    /* description of the magnetic structure. */
    /* data provided in cartesian coordinates */
    vec3 * sk;
    vec3 * isk;
    double phi;

    vec3 *R, *K, *B, *BLor, *BCont;
    pile MCont;

#ifdef _OPENMP
    double Bx=0.0;
    double By=0.0;
    double Bz=0.0;

    double BLorx=0.0;
    double BLory=0.0;
    double BLorz=0.0;
#endif

    unsigned int a;     /* counter for atoms */

    int NofM = 0; /* Number of moments considered */
    double SumOfWeights = 0;


    /* */
    atmpos = new_vec3_zero();
    m = new_vec3_zero();   /*magnetic moment of atom */

    sk  = new_vec3_zero();
    isk = new_vec3_zero();

    R = new_vec3_zero();

    /* define dupercell size */
    scx = in_supercell[0];
    scy = in_supercell[1];
    scz = in_supercell[2];

#ifdef _DEBUG
    printf("I use: %i %i %i\n",scx, scy, scz);
    printf("Total atoms: %i\n",in_natoms);
#endif

    sc_lat = new_mat3(in_cell[0], in_cell[1], in_cell[2],
                      in_cell[3], in_cell[4], in_cell[5],
                      in_cell[6], in_cell[7], in_cell[8]);

#ifdef _DEBUG
    for (i=0; i<3; i++)
        printf("Cell is: %i %e %e %e\n",i,in_cell[i*3],in_cell[i*3+1],in_cell[i*3+2]);

    /*printf("a %e %e %e\n", sc_lat.a.x, sc_lat.a.y, sc_lat.a.z); */
#endif

    K = new_vec3( in_K[0], in_K[1], in_K[2]);

#ifdef _DEBUG
    printf("K is: %e %e %e \n",K->x,K->y,K->z);
    printf("Radius is: %e\n",radius);
#endif

    sctmp = new_mat3_diag((double) scx, (double) scy, (double) scz);
    mat3_mul(sctmp,sc_lat, sc_lat);
    mat3_free(sctmp);


    /* muon position in reduced coordinates */
    muonpos = new_vec3(    (in_muonpos[0] + (scx/2) ) / (double) scx,
                           (in_muonpos[1] + (scy/2) ) / (double) scy,
                           (in_muonpos[2] + (scz/2) ) / (double) scz);

#ifdef _DEBUG
    printf("Muon pos (frac): %e %e %e\n",muonpos->x,muonpos->y,muonpos->z);
#endif

    mat3_vmul(muonpos,sc_lat,muonpos);

#ifdef _DEBUG
    printf("Muon pos (cart): %e %e %e\n",muonpos->x,muonpos->y,muonpos->z);
#endif


#ifdef _DEBUG
    for (a = 0; a < in_natoms; ++a)
    {

        /* atom position in reduced coordinates */
        /*atmpos.x =  in_positions[3*a] ; */
        /*atmpos.y =  in_positions[3*a+1] ; */
        /*atmpos.z =  in_positions[3*a+2] ; */

        printf("Atom pos (crys): %e %e %e\n",atmpos->x,atmpos->y,atmpos->z);

        /* go to cartesian coordinates (in Angstrom!) */
        /*atmpos = mat3_vmul(atmpos,sc_lat); */

        printf("Atom pos (cart): %e %e %e\n",atmpos->x,atmpos->y,atmpos->z);
#ifdef _ALTERNATE_FC_INPUT
        printf("ERROR!!! If you see this the extension compilation went wrong!\n");
        printf("FC (real, imag): %e %e %e %e %e %e\n",in_fc[6*a],in_fc[6*a+1],in_fc[6*a+2],in_fc[6*a+3],in_fc[6*a+4],in_fc[6*a+5]);
#else
        printf("FC (real, imag): %e %e %e %e %e %e\n",in_fc[6*a],in_fc[6*a+2],in_fc[6*a+4],in_fc[6*a+1],in_fc[6*a+3],in_fc[6*a+5]);
#endif
        printf("phi: %e\n",in_phi[a]);
    }

#endif


    B = new_vec3_zero();
    BLor = new_vec3_zero();
    pile_init(&MCont, nnn_for_cont);

#pragma omp parallel shared(MCont) /* remember data race! */
    {
#pragma omp for collapse(3) schedule(guided,20)  private(i,j,k,a,r,n,atmpos,sk,isk,phi,R,c,s,m,u,onebrcube) reduction(+:Bx,By,Bz,BLorx,BLory,BLorz)
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
                        vec3_set(atmpos, ( in_positions[3*a] + (double) i) / (double) scx,
                                 ( in_positions[3*a+1] + (double) j) / (double) scy,
                                 ( in_positions[3*a+2] + (double) k) / (double) scz);



                        /* go to cartesian coordinates (in Angstrom!) */
                        mat3_vmul(atmpos, sc_lat, atmpos);

                        /*printf("atompos: %e %e %e\n", atmpos.x, atmpos.y, atmpos.z); */
                        /* difference between atom pos and muon pos (cart coordinates) */

                        vec3_sub(atmpos, muonpos); /* multiplied by -1 below!! */
                        r = atmpos;


                        n = vec3_norm(r);
                        if (n < radius)
                        {
                            /* calculate magnetic moment */
#ifdef _ALTERNATE_FC_INPUT
                            printf("ERROR!!! If you see this in the Python extension something went wrong!\n");
                            vec3_set(sk, in_fc[6*a], in_fc[6*a+1], in_fc[6*a+2]);
                            vec3_set(isk, in_fc[6*a+3], in_fc[6*a+4], = in_fc[6*a+5]);
#else
                            vec3_set(sk, in_fc[6*a], in_fc[6*a+2], in_fc[6*a+4]);
                            vec3_set(isk, in_fc[6*a+1], in_fc[6*a+3], in_fc[6*a+5]);
#endif

                            phi = in_phi[a];
                            /*printf("sk = %e %e %e\n", sk.x, sk.y, sk.z); */
                            /*printf("isk = %e %e %e\n", isk.x, isk.y, isk.z); */


                            vec3_set(R, (double) i, (double) j, (double) k);

                            c = cos ( 2.0*M_PI * (vec3_dot(K,R) + phi ));
                            s = sin ( 2.0*M_PI * (vec3_dot(K,R) + phi ));

                            vec3_set(m, 0., 0., 0.);
                            vec3_daxpy(c, sk, m);
                            vec3_daxpy(s, isk, m);


                            /* calculate Lorentz Field */
#ifdef _OPENMP
                            BLorx += m.x;
                            BLory += m.y;
                            BLorz += m.z;
#else
                            vec3_add(BLor,m);
#endif

                            /* Calculate Contact Field */
                            if (n < cont_radius) {
#ifdef _DEBUG
                                printf("Adding moment to Cont: n: %e, m: %e %e %e! (Total: %d)\n", n, m->x,m->y,m->z,nnn_for_cont);
#endif                                          /* We add the moment multiplied by r^3 and then devide by Sum ^N r^3 */
#pragma omp critical
                                {
                                    vec3_cpy(sk, m); /* sk used as dummy variable */
                                    vec3_muls(1./pow(n,CONT_SCALING_POWER),sk);
                                    pile_add_element(&MCont, pow(n,CONT_SCALING_POWER), sk);
                                }
                            }


                            /*printf("I sum: r = %e, p = %e %e %e\n",n, r.x, r.y, r.z); */
                            /*printf("I sum: m = %e %e %e\n", m.x, m.y, m.z); */
                            /* sum it */
                            /* B += (( 3.0 * np.dot(nm,atom[1]) * atom[1] - nm ) / atom[0]**3)*0.9274009 */

                            /* unit vector */
                            vec3_muls(-1.0/n,r);
                            u = r; /* pointer copy, jusr for readibility */
                            onebrcube = 1.0/pow(n,3);

#ifdef _OPENMP
                            /* m is used as dummy variable for the sum! */
                            m = vec3_muls( onebrcube,vec3_sub(vec3_muls(3.0*vec3_dot(m,u),u), m));
                            Bx += m.x;
                            By += m.y;
                            Bz += m.z;
#else
                            vec3_muls(3.0*vec3_dot(m,u),u);
                            vec3_sub(u, m);
                            vec3_daxpy(onebrcube, u, B);
#endif
#ifdef _DEBUG
                            printf("B is now : %e %e %e %e\n", n, B->x, B->y, B->z);
#endif
                        }

                    }

                }
            }
        }
    }

#ifdef _DEBUG
    printf("Done with iterations!\n");
#endif

#ifdef _OPENMP
    B.x = Bx; B.y = By; B.z = Bz;
    BLor.x = BLorx; BLor.y = BLory; BLor.z = BLorz;
#endif
    /*  1 bohr_magneton/(1angstrom^3) = 9274009.5(amperes ∕ meter) */
    /*   mu_0 = 0.0000012566371((meter tesla) ∕ ampere) */
    /*   BLor = (mu_0/3)*M_Lor */
    /*   Note that befor this line BLor is just the sum of the magnetic moments! */
    /*      magnetic_constant * 1 bohr_magneton = 11.654064 T⋅Å^3 */

    vec3_muls(3./(4.*M_PI*pow(radius,3)),BLor);
    vec3_muls(0.33333333333*11.654064, BLor);
    /*printf("The Lorents field contribution is: %e %e %e Tesla!!\n",BLor.x,BLor.y,BLor.z); */

    vec3_get(BLor, &(out_field_lor[0]), &(out_field_lor[1]), &(out_field_lor[2]));

    /* Contact Field */
    BCont = new_vec3_zero();
    NofM = 0; /* Number of moments considered */
    SumOfWeights = 0;
    for (i=0; i < nnn_for_cont; i++) {
        if (MCont.ranks[i] > 0.0) {
            vec3_add(BCont, MCont.elements[i]);
            SumOfWeights += (1./MCont.ranks[i]);             /* Add the contribution weighted as norm^3 to the total */
            NofM++;
        }
    }

    pile_free(&MCont);

    /* (2 magnetic_constant/3)⋅1bohr_magneton   = ((2 ⋅ magnetic_constant) ∕ 3) ⋅ (1 ⋅ bohr_magneton) */
    /*   ≈ 7.769376E-27((g⋅m^3) ∕ (A⋅s^2)) */
    /*   ≈ 7.769376 T⋅Å^3 */

    if (NofM >0) {
        vec3_muls((1./SumOfWeights) * 7.769376, BCont);
    }     /* otherwise is zero anyway! */

    vec3_get(BCont, &(out_field_cont[0]), &(out_field_cont[1]), &(out_field_cont[2]));

    /* Dipolar Field */
    /* mu_0/4pi = 0.1E-6(newton ∕ ampere^2) = 0.1E-6((meter tesla) ∕ ampere) */
    vec3_muls(0.92740098, B); /* to tesla units */


    vec3_get(B, &(out_field_dip[0]), &(out_field_dip[1]), &(out_field_dip[2]));

    /* now free everything */
    vec3_free(atmpos);
    vec3_free(m);
    vec3_free(sk);
    vec3_free(isk);
    vec3_free(R);
    vec3_free(K);
    vec3_free(muonpos);
    vec3_free(B);
    vec3_free(BLor);
    vec3_free(BCont);
    mat3_free(sc_lat);



}


