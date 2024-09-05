# -*- coding: utf-8 -*-
import unittest
try:
    from mulfc import locfield, find_largest_sphere
except ImportError:
    from LFC import locfield, find_largest_sphere
import numpy as np

        
class TestLFCWrappers(unittest.TestCase):
    def test_one_over_r_cube(self):
        
        p  = np.array([[0.,0.,0.]])
        fc = np.array([[0.,0.,1.]],dtype=np.complex128)
        k  = np.array([0.,0.,0.])
        
        phi= np.array([0.,])
        
        mu = np.array([0.5,0.5,0.5])
        
        sc = [1,1,1]
        latpar = np.diag([2.,2.,2.])
        
        r = 10.
        
        res = locfield(latpar, p, fc, k, phi, [mu], 's', sc, r, nnn=0)[0]
                
        # zero hyperfine field since zero nnn
        np.testing.assert_array_equal(res.C,np.zeros(3))
        
        # this is 0.3333333333⋅magnetic_constant⋅(1 bohr_magneton/(4/3⋅pi⋅(10 angstrom)^3))=9.2740095E-4 tesla
        np.testing.assert_array_almost_equal(res.L,np.array([0,0,9.2740095E-4]))
        
        mu = np.array([0.5,0.,0.])
        
        res = locfield(latpar, p, fc, k, phi, [mu], 's', sc, r, nnn=0)[0]
        
        # (1/(4pi))magnetic_constant⋅(1 bohr_magneton/(1 angstrom^3)) = 0.92740095 tesla
        np.testing.assert_array_almost_equal(res.D, np.array([0,0,-0.92740095]) )
        
        mu2 = np.array([0.654,0.,0.])

        res2 = locfield(latpar, p, fc, k, phi, [mu2], 's', sc, r, nnn=0)[0]
        # ratios must be like 1/r^3
        np.testing.assert_array_almost_equal(res2.D, np.array([0,0,-0.92740095])*(1./(np.linalg.norm(mu2*2.))**3) )
    
    def test_find_largest_sphere(self):
        
        latpar = np.diag( [ 5. , 5. , 5.])
        sc = [10,10,10]
        mu = np.zeros(3)
        r = find_largest_sphere(latpar, [mu], sc)
        np.testing.assert_almost_equal(r,25.)

        mu = np.array([0.1,0.0,0.0])
        r = find_largest_sphere(latpar, [mu], sc)
        np.testing.assert_almost_equal(r,24.5)

        mu = np.array([0.1,0.1,0.1])
        r = find_largest_sphere(latpar, [mu], sc)
        np.testing.assert_almost_equal(r,24.5)

        mu = np.array([1.0,1.0,1.0])
        r = find_largest_sphere(latpar, [mu], sc)
        np.testing.assert_almost_equal(r,20.)
        
        mu = np.array([1.0,1.0,1.0])
        r = find_largest_sphere(latpar, [mu], sc)
        np.testing.assert_almost_equal(r,20.)
        
        latpar = np.array([[1.0 ,0., 0.],
                           [0., 2., 0.],
                           [0.5209445330, 1.0260604300, 2.7705264459]])
        
        mu = np.zeros(3)
        sc = [1,1,1]
        r = find_largest_sphere(latpar, [mu], sc)
        np.testing.assert_almost_equal(r,0.)
        
        sc = [2,2,2]
        r = find_largest_sphere(latpar, [mu], sc)

        assert(r < np.cos(10*np.pi/180.)) # must be less than cos(10 deg)
                                          # which is the angle btw a and c
                                          # it also depend on the angle 
                                          # btw c and b, but this I can't
                                          # easily calculate now.
        assert(r > 0.9)                   # Clearly it must also be close to 1.
        
    
    def test_mnge(self):
        # from MnGe.cif
        #1 Mn  Mn          0.13800    0.13800    0.13800    1.000    0.008    1a         1
        #2 Mn  Mn          0.36200    0.86200    0.63800    1.000    0.008    1a         1
        #3 Mn  Mn          0.86200    0.63800    0.36200    1.000    0.008    1a         1
        #4 Mn  Mn          0.63800    0.36200    0.86200    1.000    0.008    1a         1
        
        #lattice parameters
        latpar = np.diag( [ 4.76900 , 4.76900 , 4.76900])
        
        # atomic positions
        p = [[0.13800  ,  0.13800  ,  0.13800 ],
        [0.36200  ,  0.86200  ,  0.63800 ],
        [0.86200  ,  0.63800  ,  0.36200 ],
        [0.63800  ,  0.36200  ,  0.86200 ]]
        
        # propagation vector
        k = np.array([0.0,0.0,0.16710801])
        
        # phase
        phi = np.zeros(4)
        
        # calculate Fourier Components inside the unit cell
        gen_fc = lambda x: np.array([ 1.+0.0*1j , 0.0-1.0*1j,0.0])*np.exp(-2*np.pi*1j*np.dot(k,x))
        
        # scale by 1.85 i.e. the magnetic moment on Mn
        fcs = 1.85*np.array([ gen_fc([0,0,0]),
                gen_fc([0.362-0.138,  0.862-0.138 , 0.638-0.138]),
                gen_fc([0.862-0.138,  0.638-0.138 , 0.362-0.138]),
                gen_fc([0.638-0.138,  0.362-0.138 , 0.862-0.138])])
                
                
        mu = np.array([0.543 , 0.543,  0.543])
        
        sc = [10, 10, 10]
        
        r = find_largest_sphere(latpar, [mu], sc)
        
        res = locfield(latpar, p, fcs, k, phi, [mu], 'i', sc, r, nnn=3, rcont = 10.0, nangles = 300)[0]
        
        
        # Set contact field coupling (ACont) to 0.591 T/mu_b
        res.ACont = -0.0779
        
        np.testing.assert_almost_equal(np.min(np.apply_along_axis(np.linalg.norm,1, res.C / 1.85)), 
                                        np.max(np.apply_along_axis(np.linalg.norm,1, res.C / 1.85)))
                                
        np.testing.assert_almost_equal(np.min(np.apply_along_axis(np.linalg.norm,1, res.C / 1.85)), 0.59136379)
        
        # see PRB 93 174405 for values used in comparison (Fig 3 and Tab 2).
        # the simulation is not fully converged so numbers are slighlty different
        #  according to PRB, max should be 1.148
        #  according to PRB, min should be 0.468
        np.testing.assert_almost_equal(np.min(np.apply_along_axis(np.linalg.norm,1,res.T)), 0.493427, decimal=5)
        np.testing.assert_almost_equal(np.max(np.apply_along_axis(np.linalg.norm,1,res.T)), 1.1970672, decimal=5)
        
        #### Same as above, but with the slow method
        res = locfield(latpar, p, fcs, k, phi, [mu], 'r', sc, r, nnn=3, rcont = 10.0, nangles = 300, axis=[0,0,1])[0]
        
        # Set contact field coupling (ACont) to 0.591 T/mu_b
        res.ACont = -0.0779        
        
        np.testing.assert_almost_equal(np.min(np.apply_along_axis(np.linalg.norm,1, res.C / 1.85)), 
                                        np.max(np.apply_along_axis(np.linalg.norm,1, res.C / 1.85)))
                                
        np.testing.assert_almost_equal(np.min(np.apply_along_axis(np.linalg.norm,1, res.C / 1.85)), 0.59136379)
        
        # see PRB 93 174405 for values used in comparison (Fig 3 and Tab 2).
        # the simulation is not fully converged so numbers are slighlty different
        #  according to PRB, max should be 1.148
        #  according to PRB, min should be 0.468
        np.testing.assert_almost_equal(np.min(np.apply_along_axis(np.linalg.norm,1,res.T)), 0.493475, decimal=4)
        np.testing.assert_almost_equal(np.max(np.apply_along_axis(np.linalg.norm,1,res.T)), 1.1970672, decimal=4)        
        
        
if __name__ == '__main__':
    unittest.main()
