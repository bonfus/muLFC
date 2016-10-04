//Simple Sum mechanism
#ifndef SIMPLE_SUM_H
#define SIMPLE_SUM_H
void  SimpleSum(const double *in_positions, 
          const double *in_fc, const double *in_K, const double *in_phi,
          const double *in_muonpos, const int * in_supercell, const double *in_cell, 
          const double radius, const unsigned int nnn_for_cont, const double cont_radius, 
          unsigned int size,
          double *out_field_cont, double *out_field_dip, double *out_field_lor);
#endif
