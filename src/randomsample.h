/**
 * @file randomsample.h
 * @author Pietro Bonfa
 * @date 18 March 2018
 * @brief Dipolar field calculator
 *     
 */


 

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


void  RandomSample(const double *in_positions, 
          const double *in_fc, const double *in_K, const double *in_phi,
          const int * in_supercell, const double *in_cell, 
          const double radius, const unsigned int nnn_for_cont, const double cont_radius,
          const unsigned int in_natoms, const unsigned int in_nmounpos,
          const double *in_constraints, const int * in_constraint_active, 
          const unsigned int in_nconstraints,
          double* out_muonpos, double *out_field_cont, double *out_field_dip, double *out_field_lor);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

