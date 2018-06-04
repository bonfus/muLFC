/**
 * @file randomness.h
 * @author Pietro Bonfa
 * @date 18 March 2018
 * @brief Dipolar field calculator
 *     
 */

#include "types.h"
#include "config.h"

/* Consider using:
 https://www.numbercrunch.de/trng/
*/


class UniformRandomInsideUnitCell {
  public:
    UniformRandomInsideUnitCell(const Mat3& lat, int n);                         // constructor; initialize the list to be empty
    void GetRandomPos(RefVec3 pos); // print the list to output

  private:
    void AllocateRandomPoints();              // add k to the end of the list
    Mat3 lattice;
    T _boxSize;
    Vec3 boxOrigShift;
    MatX _buffer;
    int _buff_idx;
    int _buff_size;
};


// void GenBoundingBox(const Mat3& lattice, T& boxSize, RefVec3 boxOrigShift);

// void UniformRandomInsideUnitCell(Vec3& boxOrigShift,const T boxSize, const Mat3& lattice, RefVec3 muonPos);
