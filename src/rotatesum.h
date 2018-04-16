#ifndef ROTATE_SUM_H
#define ROTATE_SUM_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


//Arbitrary size sum with rotations
void RotataSum(const double *in_positions, 
          const double *in_fc, const double *in_K, const double *in_phi,
          const double *in_muonpos, const int * in_supercell, const double *in_cell, 
          const double radius, const unsigned int nnn_for_cont, const double cont_radius, 
          const unsigned int in_natoms, const unsigned int in_nmuonpos,
          const double *in_axis, unsigned int in_nangles,
          double *out_field_cont, double *out_field_dip, double *out_field_lor);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif
