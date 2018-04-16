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

/* Possibly use:
 https://www.numbercrunch.de/trng/
*/



UniformRandomInsideUnitCell::UniformRandomInsideUnitCell(const Mat3& lat, int n)
{
    srand(time(NULL));
    
    T minx, maxx, miny, maxy, minz, maxz;
    Vec3 boxMin;
    Vec3 boxMax;

    lattice = lat;

    boxOrigShift.setZero();
    
    boxMin = lat.rowwise().minCoeff();
    boxMax = lat.rowwise().maxCoeff();
    
    boxSize = boxMax.maxCoeff() - boxMin.minCoeff();
    boxOrigShift = boxMin;
    
    
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
    Vec3 posFracMin, posFracMax, aux;
    int filled = 0;
    int i, iters;

    one.setOnes();
    
    
    posFracMin.setZero();
    posFracMax.setOnes();
    
    iters = 0;
    while(filled < _buff_size) { // Add safe limit to avoid infite loops
        m = MatX::Random(3, _buff_size);

        m = (m + one)*(boxSize/2.); // shift from [-1,1] to [0,2], then multiply by box_size/2

        for (int i = 0; i < _buff_size; ++i) {
            m.col(i) += boxOrigShift;
        }
        //std::cout << "Random value is: " << m.transpose() << std::endl;
        /* This is very stupid, every time we do a matrix inversion here! */
        Crys2Cart(lattice, m, mFrac, true);
        //std::cout << "Random value is frac: " << mFrac.transpose() << std::endl;
        for (int i = 0; (i < _buff_size) && (filled < _buff_size); ++i) {
            aux = mFrac.col(i);
            if ((posFracMin.array()<=aux.array()).all() && (posFracMax.array()>aux.array()).all()) {
                _buffer.col(filled) = aux;
                filled++;
            }
        }
        iters++;
    }
    _buff_idx=0;
    if (iters > 5) std::cout << "WARNING: Iterations from uniform random: "<< iters << std::endl;
}

#ifdef __TEST_RND
int main() {
    int nloops = 1000000;
    int i;
    Mat3 lat;
    Vec3 v;
    int half_x=0;
    int half_y=0;
    int half_z=0;

    lat.setZero();
    lat(0,0) = 3.0;
    lat(1,1) = 4.0;
    lat(2,2) = 5.0;

    UniformRandomInsideUnitCell a(lat, 20000);
    
    for (i=0;i<nloops; i++)
    {
        a.GetRandomPos(v);
        if (v(0) < 0.5) half_x++;
        if (v(1) < 0.5) half_y++;
        if (v(2) < 0.5) half_z++;
    }
    std::cout << ((T)half_x)/((T)nloops) << std::endl;
    std::cout << ((T)half_y)/((T)nloops) << std::endl;
    std::cout << ((T)half_z)/((T)nloops) << std::endl;
}
#endif
