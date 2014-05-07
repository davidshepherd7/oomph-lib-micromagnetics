#!/usr/bin/env python3

# Python 2/3 compatability
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import sys
import argparse
import os
import os.path
import scipy as sp
import itertools as it

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import proj3d

from os.path import join as pjoin

from pprint import pprint as pprint

# Struct to store data for an input point
class point():
    def __init__(self, x, y, z, boundary):
        self.x = x
        self.y = y
        self.z = z

        self.bound = boundary
        self.i = 0


    def __repr__(self):
        return "("+str(self.x) + " " + str(self.y) + " " + str(self.z)+")"


def cylinder_points(radius, length, ncircle):

    # Generate points around a circle
    t_points = sp.linspace(0, 2*sp.pi, ncircle+1)
    circle = [[radius*sp.sin(t), radius*sp.cos(t)] for t in t_points[:-1]]

    # Convert into top/bottom of cylinder
    top = [point(c[0], c[1], length, 0) for c in circle]
    bottom = [point(c[0], c[1], 0.0, 2) for c in circle]

    # Combine and number
    cylinder = top + bottom

    return cylinder, top, bottom


def number_points(points):
    """Give each point a number, starting from 1."""
    for i, p in enumerate(points):
        p.i = i + 1 # TetgenMesh in oomph is stupid and needs 1-indexed
                    # nodes



def generate_facets(l1, l2):

    # This is the hard bit: we need a list of facets around the sides, so
    # we need sets of 4 points: 2 at the top, 2 at the bottom.

    facets = []
    for a, b, c, d in zip(l1, l2,
                          l1[1:] + [l1[0]],
                          l2[1:] + [l2[0]]):
        facets.append([b, a, c, d])

    return facets


def print_points(points, boundaries, outfile):

    point_zs = [p.z for p in points]
    length = max(point_zs) - min(point_zs)

    def oprt(*args, **kwargs):
        print(*args, file=outfile, **kwargs)

    oprt("# A cylinder with aspect ratio", length)
    oprt("\n# Part 1 - node list")
    oprt("# First line: <# of points> <dimension (3)> <# of attributes> <boundary markers (0 or 1)>")
    oprt(len(points), "3 0 1")

    oprt("# <point #> <x> <y> <z> [attributes] [boundary marker]")
    for p in points:
        oprt(p.i, p.x, p.y, p.z, p.bound)


    oprt("\n# Part 2 - facet list")
    oprt("# First line: <# of facets> <boundary markers (0 or 1)>")
    oprt(len(list(it.chain(*boundaries))), "1")

    oprt("# List of facets:")
    oprt("# One line: <# of polygons> [# of holes] [boundary marker]")
    oprt("# Then list of polygons: <# of corners> <corner 1> <corner 2> ... <corner #>")
    oprt("# Then list of holes: <hole #> <x> <y> <z>")

    # Print list of list of facets
    for i_boundary, facets in enumerate(boundaries):
        for facet in facets:

            # Always one polygon, no holes, ith boundary
            oprt("1 0", i_boundary)

            # List of points making up the polygon
            oprt(len(facet), end=" ")
            for p in facet:
                oprt(p.i, end=" ")

            oprt("\n")

    oprt("\n# Part 3 - hole list")
    oprt("0            # no holes")

    oprt("\n# Part 4 - region list")
    oprt("0            # no regions")


def plot_points(points):
    """Plot the points on a labeled 3d scatter plot.
    """


    # Create 3d plot
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot it
    xs = [p.x for p in points]
    ys = [p.y for p in points]
    zs = [p.z for p in points]
    i_s = [p.i for p in points]
    ax.scatter(xs, ys, zs)

    labels = []

    for i, x, y, z in zip(i_s, xs, ys, zs):
        x2, y2, _ = proj3d.proj_transform(x, y, z, ax.get_proj())

        label = plt.annotate(i,
            xy = (x2, y2), xytext = (-20, 20),
            textcoords = 'offset points', ha = 'right', va = 'bottom',
            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

        labels.append(label)

    def update_position(e):
        for i, x, y, z, l in zip(i_s, xs, ys, zs, labels):
            x2, y2, _ = proj3d.proj_transform(x, y, z, ax.get_proj())
            l.xy = x2,y2
            l.update_positions(fig.canvas.renderer)

        fig.canvas.draw()

    fig.canvas.mpl_connect('button_release_event', update_position)


def main():

    # Parse arguments
    # ============================================================

    parser = argparse.ArgumentParser(description=main.__doc__,
                                     # Don't mess up my formating in the help message
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--npoints-circle', '-n', action = "store", default=5,
                        help="number of points around the circle ends")
    parser.add_argument('--length', '-l', action = "store", default=3,
                        help="length of cylinder (=aspect ratio)")
    parser.add_argument('--plot', '-p', action = "store_true",
                        help="Plot points generated")
    parser.add_argument('--tube', '-t', action = "store_true",
                        help="Generate a hollow tu")
    args = parser.parse_args()


    if args.tube:
        shape_name = "tube"

        outer_points, outer_top, outer_bottom \
           = cylinder_points(1, float(args.length),
                             int(args.npoints_circle))

        inner_points, inner_top, inner_bottom \
           = cylinder_points(0.3, float(args.length),
                             int(args.npoints_circle))

        points = inner_points + outer_points

        number_points(points)

        top_facets = generate_facets(outer_top, inner_top)
        bottom_facets = generate_facets(outer_bottom, inner_bottom)
        inner_facets = generate_facets(inner_top, inner_bottom)
        outer_facets = generate_facets(outer_top, outer_bottom)

        boundaries = [top_facets,
                      outer_facets,
                      inner_facets,
                      bottom_facets]

    else:
        shape_name = "cylinder"

        points, top, bottom = cylinder_points(1, float(args.length),
                                              int(args.npoints_circle))
        number_points(points)

        boundaries = [[top], generate_facets(top, bottom), [bottom]]


    # Plot the points
    if args.plot:
        plot_points(points)
        plt.show()


    # Write to file
    root = ""
    basename = shape_name + str(args.length) \
      + "_" + str(args.npoints_circle)
    name = pjoin(root, basename + ".poly")

    with open(name, 'w') as outfile:
        print_points(points, boundaries, outfile)


    return 0


if __name__ == "__main__":
    sys.exit(main())
