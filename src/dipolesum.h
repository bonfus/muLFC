/**
 * @file dipolesum.h
 * @author Pietro Bonfa
 * @date 18 March 2018
 * @brief Dipolar field calculator
 *     
 */

#include <vector>
#include "types.h"

void  DipoleSum(const MatX& atomicPositions, 
          const CMatX& FC, const Vec3 K, const VecX& phi,
          const Vec3& muonpos, unsigned int scx, unsigned int scy, unsigned int scz,
          const Mat3& lattice, 
          const double radius, const unsigned int nnn_for_cont, const double cont_radius,
          RefVec3 BC, RefVec3 BD, RefVec3 BL);

void  TransformAndSum(const MatX& atomicPositions, 
          const CMatX& FC, const Vec3 K, const VecX& phi,
          const Vec3& muonpos, unsigned int scx, unsigned int scy, unsigned int scz,
          const Mat3& lattice, 
          const double radius, const unsigned int nnn_for_cont, const double cont_radius,
          const MatX& tr_mat,
          RefMatX BC, RefMatX BD, RefMatX BL);
