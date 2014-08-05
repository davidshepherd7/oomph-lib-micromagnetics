#!/usr/bin/env python3

import sys
import argparse
import copy

import scipy as sp
import scipy.linalg
import scipy.integrate
import scipy.optimize

from scipy import sin, cos, exp, sqrt, pi, dot, cross, array

# my library
import utils

# ??ds shape functions should really be a property of nodes I think...
# ??ds use factory to construct nodes with the expected properties?


class Node(object):

    def __init__(self, x):
        self.x = array(x)
        self.m = array([0.0, 0.0, 0.0])
        self.mprev = array([0.0, 0.0, 0.0])


class Element(object):

    def __init__(self):
        self.nodes = []
        self.el_dim = None
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
        x_mat = sp.transpose(array([n.x for n in self.nodes]))
        dshapeds_mat = array([self.dshapeds(s, j) for j in self.nodei()])

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


class TriangleElement(Element):

    def __init__(self, x0, x1, x2):
        self.nodes = [Node(x0), Node(x1), Node(x2)]
        self.el_dim = 2

    def local_coordinate_of_node(self, j):
        if j == 0:
            return array([1, 0])
        elif j == 1:
            return array([0, 1])
        elif j == 2:
            return array([0, 0])
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
        knot_values = array([function(s) for s in self.knots])

        # Add up and multiply by weights, do it "manually" to make sure we do
        # the right thing with vector fuctions.
        thesum = 0.0
        for v, w in zip(knot_values, self.weights):
            thesum += v*w

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
    m = array([0.2, 0.0, 1.0])
    ml = sqrt(sp.dot(m,m))
    m = m/ml
    return m


def nodal_quadrature_factory(ele):
    """Create the local nodal quadrature scheme for an element.
    """
    betas = []
    for i, n in enumerate(ele.nodes):
        beta, _ = sp.integrate.dblquad(lambda x,y: ele.shape([x,y], i),
                                        0, 1,
                                        lambda x: 0,
                                        lambda x: 1-x)
        betas.append(beta)

    return Quadrature([n.x for n in ele.nodes], betas)


def skew(a):
    """Skew matrix of a vector a.
    """
    return array([[0.0, -a[2], a[1]],
                  [a[2], 0.0, -a[0]],
                  [-a[1], a[0], 0.0]])


def llg_residual_integrand(t, x, m, dmdt, dmdx, test, dtestdx, happ, dampc):

    # Componentwise dot product of mxdmdx and dtestdx (pretending we
    # have vector test functions)
    mxdmdx = sp.dot(skew(m), dmdx)
    exch = array([sp.dot(r_mxdmdx, dtestdx) for r_mxdmdx in mxdmdx])

    return (sp.dot(dmdt, test)
            + sp.dot(sp.cross(m, happ), test)
            - dampc * sp.dot(sp.cross(m, dmdt), test)
            - exch)


def m_length_error(m):
    return abs(1 - sqrt(sp.dot(m,m)))


def output_solution(t, ele):
    print("t =", t)
    for n in ele.nodes:
        print("x =", n.x)
        print("m =",n.m)

    m_len_errors = [m_length_error(n.m) for n in ele.nodes]
    print("max length error =", max(m_len_errors))
    print()

    return


def initialise_solution(ele, func):
    for n in ele.nodes:
        n.mprev = copy.deepcopy(func(n.x))
        n.m = copy.deepcopy(func(n.x))


def time_integrate(residual, y0, dt, tmax, actions_after_time_step=None):
    """Integrate non-linear residual function over time.
    """

    t = 0
    ts = [t]
    ys = [y0]

    while t < tmax:
        t += dt

        # Construct residual dependant only on next y value
        yn = ys[-1]
        def time_residual(ynp1):
            y = (ynp1 + yn)/2
            dydt = (ynp1 - yn)/dt
            return residual(t, y, dydt)

        # Solve system
        ynp1 = sp.optimize.newton_krylov(time_residual,
                                         yn,
                                         f_tol=1e-12,
                                         rdiff=1e-14)

        # Store values or whatever
        if actions_after_time_step is not None:
            actions_after_time_step(t, ynp1)

        # Store
        ts.append(t)
        ys.append(ynp1)


    return (ts, ys)


def main():

    # create element and quadrature
    ele = TriangleElement([1,0], [0, 0.5], [0, 0])
    ele.quad = nodal_quadrature_factory(ele)

    # Set initial conditions:

    # Varying in space, no applied field
    initialise_solution(ele, m_initial)
    happ = array([0,0,0])
    damping = 0.5
    dt = 0.1
    tmax = 1.0

    # # Constant in space, reversal under applied field
    # initialise_solution(ele, m_initial_z)
    # happ = array([0, 0, -1.1])
    # damping = 0.5
    # dt = 0.1
    # tmax = 1.0


    # Doc initial solution
    output_solution(0.0, ele)

    # Define the residual function in terms of m and dmdt.
    def residual(time, m, dmdt):

        # Calculate residuals for each test function
        rs = []
        for i_node in ele.nodei():

            # Integrand as a function only of local coord
            def integrand(s):
                x_s = ele.interpolate(time, s, lambda n:n.x)
                m_s = ele.interpolate_data(time, s, m)
                dmdx_s = ele.interpolate_data_dx(time, s, m)
                dmdt_s = ele.interpolate_data(time, s, dmdt)
                test_s = ele.test(s, i_node)
                dtestdx_s = ele.dtestdx(s, i_node)

                return llg_residual_integrand(time, x_s, m_s,
                                 dmdt_s, dmdx_s,
                                 test_s, dtestdx_s,
                                 happ, damping)

            rs.append(ele.quad.integrate(integrand))

        return array(rs)


    def my_actions_after_timestep(t, new_m):
        """Actions after time step: store new values and print.
        """
        # Update nodal values
        for new_m_node, n in zip(new_m, ele.nodes):
            n.mprev = copy.deepcopy(n.m)
            n.m = copy.deepcopy(new_m_node)

        # Print
        output_solution(t, ele)

    # Get initial condition
    m0 = array([n.m for n in ele.nodes])

    # Integrate
    ts, ms = time_integrate(residual, m0, dt, tmax,
                            actions_after_time_step=my_actions_after_timestep)

    return 0




# Self tests
# ============================================================

def test_time_integration():

    dt = 0.01
    tmax = 1

    def check_time_integration(res, y0, exact):
        ts, ys = time_integrate(res, y0, dt, tmax)

        exacts = [exact(t) for t in ts]

        print(ts, ys, exacts)

        utils.assert_list_almost_equal(ys, exacts, 1e-3)

    tests = [(lambda t,y,dy: y - dy, 1.0, lambda t: array(exp(t))),

             (lambda t,y,dy: y + dy, 1.0, lambda t: array(exp(-t))),

             (lambda t,y,dy: array([-0.1*sin(t), y[1]]) - dy,
              array([0.1*cos(0.0), exp(0.0)]),
              lambda t: array([0.1*cos(t), exp(t)])),

            ]

    for r, y0, exact in tests:
        yield check_time_integration, r, y0, exact



if __name__ == "__main__":
    sys.exit(main())
