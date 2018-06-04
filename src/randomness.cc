/**
 * @file randomness.cc
 * @author Pietro Bonfa
 * @date 18 March 2018
 * @brief Dipolar field calculator
 *     
 */

#define _USE_MATH_DEFINES

#include "randomness.h"
#include "lattice.h"

#ifdef CXX11RANDOM
#include <random>
#endif

#ifdef WIN32
#include <Windows.h>
#include <stdint.h> // portable: uint64_t   MSVC: __int64 
#else
#include <sys/time.h>
#endif

/* Possibly use:
 https://www.numbercrunch.de/trng/
*/



UniformRandomInsideUnitCell::UniformRandomInsideUnitCell(const Mat3& lat, int n)
{
    srand(time(NULL));
    
    T minx, maxx, miny, maxy, minz, maxz;
    Vec3 boxMin;
    Vec3 boxMax;
    MatX Crys(3,8);
    MatX Cart(3,8);

    lattice = lat;

    /* Find box size evaluating all lattice vectors in cartesian
      coordinates */
    boxOrigShift.setZero();
    int l=0;
    for (int i=0; i<=1; i++) {
        for (int j=0; j<=1; j++) {
            for (int k=0; k<=1; k++) {
                Crys.col(l) << (T)i,(T)j,(T)k;
                l++;
            }
        }
    }

    Crys2Cart(lattice, Crys, Cart, false);

    boxMin = Cart.rowwise().minCoeff();
    boxMax = Cart.rowwise().maxCoeff();

    _boxSize = (boxMax - boxMin).maxCoeff();
    boxOrigShift = boxMin;

    // Allocate and initialize variables
    _buffer.resize(3, n);
    _buff_size = n;
    _buff_idx  = 0; // trigger allocation
    AllocateRandomPoints();
    
}

void UniformRandomInsideUnitCell::GetRandomPos(RefVec3 muonPos)
{
    if (_buff_idx < _buff_size) {
        muonPos = _buffer.col(_buff_idx);
        _buff_idx++;
    } else {
        // regenerate
        AllocateRandomPoints();
        muonPos = _buffer.col(0);
        _buff_idx=1;
    }
}
void UniformRandomInsideUnitCell::AllocateRandomPoints(){
    /* A trivial rejection samplig is used to generate a uniform distribution inside the unit cell*/
    MatX m(3, _buff_size), one(3, _buff_size);
    MatX mFrac(3, _buff_size);
    Vec3 aux;
    int filled = 0;
    int i, iters;
    
    one.setOnes();

#ifndef CXX11RANDOM
#ifndef WIN32
    timeval curTime;
    gettimeofday(&curTime, NULL);
    std::srand((unsigned int) curTime.tv_usec);
#else
    SYSTEMTIME  system_time;
    GetSystemTime(&system_time);
    std::srand((unsigned int) system_time.wMilliseconds );
#endif
#else
    const T range_from  = 0;
    const T range_to    = _boxSize;
    std::random_device                  rand_dev;
    std::mt19937                        generator(rand_dev());
    std::uniform_real_distribution<T>  distr(range_from, range_to);
#endif

    iters = 0;
    while(filled < _buff_size) { // Add safe limit to avoid infite loops
#ifndef CXX11RANDOM
        m = MatX::Random(3, _buff_size);
        m = (m + one)*(_boxSize/2.); // shift from [-1,1] to [0,2], then multiply by box_size/2
#else
        for (int i = 0; i < _buff_size; ++i) {
            m.col(i) << distr(generator), distr(generator), distr(generator);
            m.col(i) += boxOrigShift;
        }
#endif
        //std::cout << "Random value is: " << m.transpose() << std::endl;
        /* This is very stupid, every time we do a matrix inversion here! */
        Crys2Cart(lattice, m, mFrac, true);
        //std::cout << "Random value is frac: " << mFrac.transpose() << std::endl;
        for (int i = 0; (i < _buff_size) && (filled < _buff_size); ++i) {
            aux = mFrac.col(i);
            if (aux(0) < 0.0 || aux(0) >= 1.0) continue;
            if (aux(1) < 0.0 || aux(1) >= 1.0) continue;
            if (aux(2) < 0.0 || aux(2) >= 1.0) continue;
            _buffer.col(filled) = aux;
            filled++;
        }
        iters++;
    }
    _buff_idx=0;
    if (iters > 5) std::cout << "WARNING: Iterations from uniform random: "<< iters << std::endl;
}
