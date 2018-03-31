#include <iostream>
#include "types.h"

void Crys2Cart(const Mat3& trmat, const MatX& pos, RefMatX out, bool reverseDirection)
{
  Vec3 aux;
  if (reverseDirection) {
    out = trmat.inverse() * pos;
  } else {
    out = trmat * pos;
  }
}


void recips (const Mat3& lat, RefMatX rec)
{
  int ii,lattype;
  T bzvol, ucvol;
  Vec3 a1,a2,a3,b1,b2,b3;
  T pi = 3.1415926535897932384626433;
  
  a1 = lat.col(0);
  a2 = lat.col(1);
  a3 = lat.col(2);
  
  
  ucvol=a1(0)*(a2(1)*a3(2)-a2(2)*a3(1))+a1(1)*(a2(2)*a3(0)-a2(0)*a3(2))+a1(2)*(a2(0)*a3(1)-a2(1)*a3(0));
  bzvol=8*pi*pi*pi/ucvol;

  //calculate reciprocal-space lattice vectors b1-b3
  b1(0)=2.0*pi*(a2(1)*a3(2)-a2(2)*a3(1))/ucvol;
  b1(1)=2.0*pi*(a2(2)*a3(0)-a2(0)*a3(2))/ucvol;
  b1(2)=2.0*pi*(a2(0)*a3(1)-a2(1)*a3(0))/ucvol;
  b2(0)=2.0*pi*(a3(1)*a1(2)-a3(2)*a1(1))/ucvol;
  b2(1)=2.0*pi*(a3(2)*a1(0)-a3(0)*a1(2))/ucvol;
  b2(2)=2.0*pi*(a3(0)*a1(1)-a3(1)*a1(0))/ucvol;
  b3(0)=2.0*pi*(a1(1)*a2(2)-a1(2)*a2(1))/ucvol;
  b3(1)=2.0*pi*(a1(2)*a2(0)-a1(0)*a2(2))/ucvol;
  b3(2)=2.0*pi*(a1(0)*a2(1)-a1(1)*a2(0))/ucvol;
  rec.col(0) = b1;
  rec.col(1) = b2;
  rec.col(2) = b3;
}

T GetMinDistanceFromAtoms(const Mat3& lattice, const MatX& atomicPositions, const Vec3& intPosition)
{
    Mat3 scLattice;
    Mat3 aux;
    Vec3 intPositionScaled;
    Vec3 intPositionCart;
    Vec3 vaux;
    Vec3 vauxCart;
    Vec3 apos_aux;
    Mat3 sctmp;
    sctmp << (T) 3,     0,     0, 
                 0, (T) 3,     0,
                 0,     0, (T) 3;
    
    double nrm = std::numeric_limits<double>::max();
    
    scLattice = lattice*sctmp;

    // lattice scaled by 3 in all directions
    intPositionScaled = intPosition/3. + 0.3333333333*Vec3::Ones();
    
    // pos in cartesian coordinates
    Crys2Cart(scLattice, intPositionScaled, intPositionCart, false);
    
    // scale lattice coordinates by 3

    for (int i=0; i < 3; i ++) {
        for (int j=0; j < 3; j ++) {
            for (int k=0; k < 3; k ++) {
                vaux << ((T) i)/3. ,((T) j)/3. ,((T) k)/3.;
                for (int l=0; l < atomicPositions.cols(); l++) {
                    
                    apos_aux = atomicPositions.col(l)/3. + vaux;
                    
                    // apos in Cart coordinates
                    Crys2Cart(scLattice, apos_aux, vauxCart, false);
                    std::cout << vauxCart.transpose() << "-" <<  intPositionCart.transpose() << std::endl;
                    vauxCart = vauxCart - intPositionCart;
                    std::cout << "norm: " << vauxCart.norm() << std::endl;
                    nrm = std::min(vauxCart.norm(), nrm);
                    
                }
            }
        }
    }
    return nrm;
}

#ifdef __TEST
int main()
{
  Mat3 m(3,3);
  Mat3 r(3,3);
  
  MatX a(3,24);
  MatX aux(3,24);
  Vec3 p;
  
  // Troilite
  m = Mat3::Zero();
  m(0,0) =  5.9990000; //        0.0000000000         0.0000000000
  m(0,1) = -2.9995000;  m(1,1) = 5.1952864601  ;  //   0.0000000000
  m(2,2) = 11.7100000381;
  
  
  std::cout << "V1: " << m.col(0).transpose() << std::endl;
  std::cout << "V2: " << m.col(1).transpose() << std::endl;
  std::cout << "V3: " << m.col(2).transpose() << std::endl;
  
  recips(m, r);

  std::cout << "RV1: " << r.col(0).transpose() << std::endl;
  std::cout << "RV2: " << r.col(1).transpose() << std::endl;
  std::cout << "RV3: " << r.col(2).transpose() << std::endl;

  
  // Atoms
  a.col(0)  << 0.373899978 , 0.050199998 , 0.123000002; // Fe 
  a.col(1)  << 0.949799981 , 0.323699990 , 0.123000002; // Fe 
  a.col(2)  << 0.676299987 , 0.626099981 , 0.123000002; // Fe 
  a.col(3)  << 0.373899978 , 0.050199998 , 0.376999998; // Fe 
  a.col(4)  << 0.949799981 , 0.323699990 , 0.376999998; // Fe 
  a.col(5)  << 0.676299987 , 0.626099981 , 0.376999998; // Fe 
  a.col(6)  << 0.050200007 , 0.373899983 , 0.876999957; // Fe 
  a.col(7)  << 0.323700015 , 0.949799995 , 0.876999957; // Fe 
  a.col(8)  << 0.626100001 , 0.676299951 , 0.876999957; // Fe 
  a.col(9)  << 0.050200007 , 0.373899983 , 0.623000043; // Fe 
  a.col(10) << 0.323700015 , 0.949799995 , 0.623000043; // Fe 
  a.col(11) << 0.626100001 , 0.676299951 , 0.623000043; // Fe 
  a.col(12) << 0.000000000 , 0.000000000 , 0.000000000; // S1 
  a.col(13) << 0.000000000 ,-0.000000000 , 0.500000000; // S1 
  a.col(14) << 0.333333337 , 0.666666673 , 0.018600000; // S2 
  a.col(15) << 0.333333337 , 0.666666673 , 0.481400012; // S2 
  a.col(16) << 0.666666668 , 0.333333337 , 0.981400012; // S2 
  a.col(17) << 0.666666668 , 0.333333337 , 0.518599988; // S2 
  a.col(18) << 0.665799990 , 0.996999956 , 0.250000000; // S3 
  a.col(19) << 0.002999996 , 0.668799978 , 0.250000000; // S3 
  a.col(20) << 0.331200009 , 0.334200017 , 0.250000000; // S3 
  a.col(21) << 0.996999966 , 0.665799970 , 0.750000020; // S3 
  a.col(22) << 0.668800012 , 0.003000000 , 0.750000020; // S3 
  a.col(23) << 0.334200001 , 0.331199987 , 0.750000020; // S3 

  p << 0.5, 0.5, 0.5;
  
  Crys2Cart(m, a, aux, false);
  for (int i = 0; i<aux.cols(); i++) {
    std::cout << aux.col(i).transpose() << std::endl;
  }
  
  Crys2Cart(m, aux, a, true);
  for (int i = 0; i<aux.cols(); i++) {
    std::cout << a.col(i).transpose() << std::endl;
  }
  std::cout << "Min distance: " << GetMinDistanceFromAtoms(m,a,p) << std::endl;
}
#endif
