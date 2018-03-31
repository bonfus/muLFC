#ifndef DIPOLAR_TENSOR_H
#define DIPOLAR_TENSOR_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/** @brief Dipolar Tensor
 *         
 * Evaluation of the dipolar tensor. This function constructs the dipolar tensor with the positions given in in_positions
 * 
 */void DipolarTensor(const double *in_positions, 
          const double *in_muonpos, const int * in_supercell, const double *in_cell, 
          const double radius, unsigned int size,
          double *out_field);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif
