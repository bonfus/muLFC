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


void  RandomSample(const T *in_positions, 
          const T *in_fc, const T *in_K, const T *in_phi,
          const int * in_supercell, const T *in_cell, 
          const T radius, const unsigned int nnn_for_cont, const T cont_radius,
          const T min_radius_from_atoms, 
          const unsigned int in_natoms, unsigned int in_nmounpos,
          T* out_muonpos, T *out_field_cont, T *out_field_dip, T *out_field_lor);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

