#!/usr/bin/env python3

import sys
import argparse
import os
import os.path
import copy

import scipy as sp
import scipy.linalg
import scipy.integrate
import scipy.optimize
from scipy import sin, cos, sqrt, pi, dot, cross, array
from scipy.optimize import newton

from os.path import join as pjoin

import functools as ft
from functools import partial as par

import itertools as it

from pprint import pprint


# ??ds shape functions should really be a property of nodes I think...

# ??ds use factory to construct nodes with the expected properties?


# # testing for element
# s = sp.array([1, 0.5])
# print([n.x for n in ele.nodes])
# print([ele.shape(s, j) for j in ele.nodei()])
# print([ele.dshapedx(s, j) for j in ele.nodei()])
# print(ele.interpolate(s, lambda n:n.x))
# print(ele.interpolate(s, lambda n:n.m))


class Node(object):

    def __init__(self, x):
        self.x = sp.array(x)
        self.m = sp.array([0.0, 0.0, 0.0])
        self.mprev = sp.array([0.0, 0.0, 0.0])


class Element(object):

    def __init__(self):
        self.nodes = []
        self.el_dim = 0
        self.quad = None

    # helpers
    def nnode(self):
        return len(self.nodes)

    def nodei(self):
        return range(0, self.nnode())

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

    def J(self, s):
        return sp.linalg.det(self.dxds(s))

    def interpolate(self, t, s, nodal_function):
        return sum([nodal_function(n) * self.shape(s, i)
                    for i, n in enumerate(self.nodes)])

    def interpolate_dx(self, t, s, nodal_function):
        return sum([sp.outer(nodal_function(n), self.dshapedx(s, i))
                    for i, n in enumerate(self.nodes)])

    def interpolate_data(self, t, s, nodal_values):
        return sum([nodal_values[i] * self.shape(s, i)
                    for i, n in enumerate(self.nodes)])

    def interpolate_data_dx(self, t, s, nodal_values):
        return sum([sp.outer(nodal_values[i], self.dshapedx(s, i))
                    for i, n in enumerate(self.nodes)])


class TriElement(Element):

    def __init__(self, x0, x1, x2):
        self.nodes = [Node(x0), Node(x1), Node(x2)]
        self.el_dim = 2

    def local_coordinate_of_node(self, j):
        if j == 0:
            return sp.array([1, 0])
        elif j == 1:
            return sp.array([0, 1])
        elif j == 2:
            return sp.array([0, 0])
        else:
            raise IndexError


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


    def integrate(self, function, p=False):
        # eval function at each knot
        knot_values = sp.array([function(s) for s in self.knots])

        # Add up and multiply by weights
        thesum = 0.0
        for v, w in zip(knot_values, self.weights):
            thesum += v*w
            if p:
                print("integartion", v,w,thesum)


        return thesum


def m_initial(x):
    l = 50

    m = array([sin(x[0]*2*pi/l)/2 + sin(x[1]*2*pi/l)/2,
               cos(x[0]*2*pi/l)/2 + cos(x[1]*2*pi/l)/2,
               0.0])
    m[2] = sqrt(1.0 - m[0]**2 - m[1]**2) # normalise

    # Check it
    assert abs(sp.dot(m, m) - 1.0) < 1e-14

    return m


def m_initial_z(x):
    m = array([0.5, 0.1, 1.0])
    ml = sqrt(sp.dot(m,m))
    m = m/ml
    return m


def nodal_quadrature_fac(ele):
    betas = []
    for i, n in enumerate(ele.nodes):
        beta, _ = sp.integrate.dblquad(lambda x,y: ele.shape([x,y], i),
                                        0, 1,
                                        lambda x: 0,
                                        lambda x: 1-x)
        betas.append(beta)

    return Quadrature([n.x for n in ele.nodes], betas)


def skew(a):
    result = sp.zeros((3,3))
    result[0, 1] = -a[2];
    result[0, 2] = a[1];
    result[1, 0] = a[2];
    result[1, 2] = -a[0];
    result[2, 0] = -a[1];
    result[2, 1] = a[0];

    return result


def integrand(t, x, mnext, mprev, dmdx, test, dtestdx, happ, p=False):

    dampc = 0.5
    dt = 0.1

    # imr approximations
    m = (mnext + mprev)/2
    dmdt = (mnext - mprev)/dt

     # Componentwise dot product of mxdmdx and dtestdx (pretending we
    # have vector test functions)
    mxdmdx = sp.dot(skew(m), dmdx)
    exch = sp.array([sp.dot(r_mxdmdx, dtestdx) for r_mxdmdx in mxdmdx])

    if p:
        print(mnext, mprev)
        a = ((mnext - mprev) / 0.1)*test + sp.dot(sp.cross((mnext + mprev) / 2, happ), test)
        b = -0.5 * sp.cross((mnext + mprev) / 2, (mnext - mprev)/0.1) * test
        print("*", a)
        print("**", b)
        print("dmdt =", dmdt)

    return (sp.dot(dmdt, test)
            + sp.dot(sp.cross(m, happ), test)
            - dampc * sp.dot(sp.cross(m, dmdt), test)
            - exch)


def residual(ele, nodal_m_vals, time, happ, p=False):

    # Calculate residuals for each test function
    rs = []
    for i_node in ele.nodei():

        # Integrand as a function only of local coord
        def integrand2(s):
            x = ele.interpolate(time, s, lambda n:n.x)
            mnext = ele.interpolate_data(time, s, nodal_m_vals)
            mprev = ele.interpolate(time, s, lambda n:n.mprev)
            dmdx = ele.interpolate_data_dx(time, s, nodal_m_vals)
            test = ele.test(s, i_node)
            dtestdx = ele.dtestdx(s, i_node)

            return integrand(time, x, mnext, mprev, dmdx, test, dtestdx, happ, p)

        rs.append(ele.quad.integrate(integrand2, p))

    return sp.array(rs).flatten()


def mlen(m):
    return abs(1 - sqrt(sp.dot(m,m)))


def output_solution(t, ele):
    print("t =", t)
    for n in ele.nodes:
        print(n.x, n.m, n.mprev, mlen(n.m))
    print()

    return


def initialise_solution(ele, func):
    for n in ele.nodes:
        n.mprev = copy.deepcopy(func(n.x))
        n.m = copy.deepcopy(func(n.x))


def main():

    # create element and quadrature
    ele = TriElement([1,0], [0, 1], [0, 0])
    ele.quad = nodal_quadrature_fac(ele)

    # Set initial conditions

    # initialise_solution(ele, m_initial)
    # happ = array([0,0,0])

    initialise_solution(ele, m_initial_z)
    happ = array([0, 0, -1.1])

    dt = 0.1
    t = 0.0

    output_solution(t, ele)

    while t < 1.0:

        t = t + dt

        temp = sp.array([n.m for n in ele.nodes])
        print(residual(ele, temp, happ=happ, time=t, p=True))

        # minimise residual
        curr_m = sp.array([[1,1,1] for n in ele.nodes])
        out = sp.optimize.newton_krylov(lambda m: residual(ele, m,
                                                     happ=happ,
                                                     time=t),
                                        copy.deepcopy(curr_m),
                                        f_tol=1e-12,
                                        rdiff=1e-14,
                                        )

        # Update m values
        for new_m, n in zip(out, ele.nodes):
            n.mprev = copy.deepcopy(n.m)
            n.m = copy.deepcopy(new_m)

            print(n.m - n.mprev)

        # Print
        output_solution(t, ele)

    return 0


if __name__ == "__main__":
    sys.exit(main())
