import unittest
import numpy as np
from mr_distortion_phantom.core import rigid,evaluate,sphere_coverage
from mr_distortion_phantom.__main__ import reproduce
from tempfile import TemporaryDirectory

class GeometryTests(unittest.TestCase):
    def test_rigid_and_heldout(self):
        x=np.random.default_rng(603).normal(size=(120,3))*20
        angle=.3;R=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1.]])
        y=x@R+[3,-5,8];fit=np.arange(120)<40;s,v,_,_=evaluate(y,x,fit)
        self.assertLess(s['max_mm'],1e-10)
    def test_known_local_deformation_is_not_absorbed(self):
        x=np.random.default_rng(604).normal(size=(120,3))*20;y=x.copy();y[40:,0]+=.7
        s,v,_,_=evaluate(x,y,np.arange(120)<40)
        self.assertAlmostEqual(s['mean_mm'],.7,places=10)
    def test_coverage_limits(self):
        self.assertEqual(float(sphere_coverage(10,0)),1.)
        self.assertEqual(float(sphere_coverage(10,10)),0.)
        self.assertAlmostEqual(float(sphere_coverage(10,1)),.8505)
        with self.assertRaises(ValueError):sphere_coverage(0,1)
    def test_complete_numeric_regression(self):
        with TemporaryDirectory() as d:self.assertEqual(len(reproduce(d)),16)
if __name__=='__main__':unittest.main()
