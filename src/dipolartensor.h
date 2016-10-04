#ifndef DIPOLAR_TENSOR_H
#define DIPOLAR_TENSOR_H
// Dipolar tensor with the positions given in in_positions
void DipolarTensor(const double *in_positions, 
          const double *in_muonpos, const int * in_supercell, const double *in_cell, 
          const double radius, unsigned int size,
          double *out_field);
#endif
