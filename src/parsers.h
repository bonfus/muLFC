

#include "types.h"

int ParseAndFilterMagneticAtoms(const T *in_positions, 
          const T *in_fc, const T *in_phi,
          const unsigned int in_natoms,
          RefMatX AtomPos, RefMatX MagAtomPos, RefCMatX FC, RefVecX phi);
