#!/usr/bin/env python3

import sys
import argparse
import os
import os.path
import copy

import scipy as sp
import scipy.linalg
from scipy import sin, cos, sqrt, pi, dot, cross, array

from os.path import join as pjoin

import functools as ft
from functools import partial as par


# ??ds shape functions should really be a property of nodes I think...

# ??ds use factory to construct nodes with the expected properties?



class Node(object):

    def __init__(self, x):
        self.x = sp.array(x)
        self.m = sp.array([0.0, 0.0, 0.0])
        self.mprev = sp.array([0.0, 0.0, 0.0])


class Element(object):

    def __init__(self):
        self.nodes = []
        self.el_dim = 0

    # helpers
    def nodei(self):
        return range(0, len(self.nodes))

    def dim():
        return len(self.nodes[0].x)

    # functions to implement
    def shape(self, s, j):
        raise NotImplementedError

    def dshapeds(self, s, j):
        raise NotImplementedError


    # derived functions
    def dshapedx(self, s, j):
        return sp.dot(self.dshapeds(s, j), self.dsdx(s))

    def test(self, s, j):
        return self.shape(s, j)

    def dtestdx(self, s, j):
        return self.dshapedx(s, j)

    # derived: jacobian of transformation and inverse
    def dxds(self, s):
        x_mat = sp.transpose(sp.array([n.x for n in self.nodes]))
        dshapeds_mat = sp.array(
            [self.dshapeds(s, j) for j in self.nodei()])

        return sp.dot(x_mat, dshapeds_mat)

    def dsdx(self, s):
        return sp.linalg.inv(self.dxds(s))

    def interpolate(self, t, s, nodal_function):
        return sum([nodal_function(n) * self.shape(s, i)
                    for i, n in enumerate(self.nodes)])

    def interpolate_dx(self, t, s, nodal_function):
        return sum([sp.outer(nodal_function(n), self.dshapedx(s, i))
                    for i, n in enumerate(self.nodes)])


class TriElement(Element):

    def __init__(self, x0, x1, x2):
        self.nodes = [Node(x0), Node(x1), Node(x2)]
        self.el_dim = 2


    def shape(self, s, j):
        if j == 0:
            return s[0]
        elif j == 1:
            return s[1]
        elif j == 2:
            return 1 - s[0] - s[1]
        else:
            raise IndexError

    def dshapeds(self, s, j):
        if j == 0:
            return array([1, 0])
        elif j == 1:
            return array([0, 1])
        elif j ==2:
            return array([-1, -1])
        else:
            raise IndexError



class Quadrature(object):

    def __init__(self, knots, weights):
        self.weights = weights
        self.knots = knots

        assert(len(self.weights) == len(self.knots))


    def integrate(self, function, t, *arg_functions):
        # eval each arg at each knot
        arg_values = [[argf(t, k) for argf in arg_functions] for k in self.knots]

        # eval function at each knot
        knot_values = sp.array([function(t, *kargs)
                                for k, kargs in zip(self.knots, arg_values)])

        # Add up and multiply by weights
        return sp.dot(knot_values, self.weights)

# nodal_quadrature = quadrature(nodes, betas, integrand)

def m_initial(x):
    l = 50

    m = array([sin(x[0]*2*pi/l)/2 + sin(x[1]*2*pi/l)/2,
               cos(x[0]*2*pi/l)/2 + cos(x[1]*2*pi/l)/2,
               0.0])
    m[2] = sqrt(1.0 - m[0]**2 - m[1]**2) # normalise

    # Check it
    assert abs(sp.dot(m, m) - 1.0) < 1e-14

    return m


def nodal_quadrature_fac(ele):
    gaussquad = Quadrature([[1/6, 1/6], [2/3, 1/6], [1/6, 2/3]],
                           [1/3, 1/3, 1/3])

    betas = []
    for i, n in enumerate(ele.nodes):
        beta = gaussquad.integrate(lambda t,s: ele.shape(s, i),
                                   0.0,
                                   lambda t,s:s)
        betas.append(beta)

    return Quadrature([n.x for n in ele.nodes], betas)




def main():

    # create element
    ele = TriElement([1,0], [0, 1], [0, 0])

    # Set initial conditions
    for n in ele.nodes:
        n.m_prev = m_initial(n.x)
        n.m = m_initial(n.x)


    # # testing
    # s = sp.array([1, 0.5])
    # print([n.x for n in ele.nodes])
    # print([ele.shape(s, j) for j in ele.nodei()])
    # print([ele.dshapedx(s, j) for j in ele.nodei()])
    # print(ele.interpolate(s, lambda n:n.x))
    # print(ele.interpolate(s, lambda n:n.m))

    def integrand(t, x, m, dmdx, test, dtestdx):
        return


    def integrand(t, x, mnext, mprev, dmdx, test, dtestdx):

        dampc = 0.5
        dt = 0.1

        # imr:
        m = (mnext + mprev)/2
        dmdt = (mnext - mprev)/dt

        print(t,x,m, dmdt,dmdx)
        print("*", test, dtestdx)

        return (sp.dot(dmdt, test)
                + dampc * sp.dot(sp.cross(m, dmdt), test)
                - sp.dot(sp.cross(m, dmdx), dtestdx))


    # Create nodal quadrature
    # ============================================================


    quad = nodal_quadrature_fac(ele)
    i_node = 0
    quad.integrate(integrand, 0.0,
                   par(ele.interpolate, nodal_function=lambda n:n.x),
                   par(ele.interpolate, nodal_function=lambda n:n.m),
                   par(ele.interpolate, nodal_function=lambda n:n.m_prev),
                   par(ele.interpolate_dx, nodal_function=lambda n:n.m),
                   lambda t,s: ele.test(s, i_node),
                   lambda t,s: ele.dtestdx(s, i_node)
                   )




    return 0


if __name__ == "__main__":
    sys.exit(main())
