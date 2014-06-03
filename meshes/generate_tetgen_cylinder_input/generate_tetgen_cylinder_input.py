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

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import proj3d

from os.path import join as pjoin

# Struct to store data for an input point
class point():
    def __init__(self, x, y, z, boundary):
        self.x = x
        self.y = y
        self.z = z

        self.bound = boundary
        self.i = 0


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
    args = parser.parse_args()


    # Generate points around a circle
    t_points = sp.linspace(0, 2*sp.pi, int(args.npoints_circle)+1)
    circle = [[sp.sin(t), sp.cos(t)] for t in t_points[:-1]]

    # Convert into top/bottom of cylinder
    top = [point(c[0], c[1], float(args.length), 0) for c in circle]
    bottom = [point(c[0], c[1], 0.0, 2) for c in circle]

    # Combine and number
    cylinder = top + bottom

    # Number the points
    for i, p in enumerate(cylinder):
        p.i = i + 1 # TetgenMesh in oomph is stupid and needs 1-indexed nodes

    # Write to file
    root = ""
    name = pjoin(root, "cylinder"+str(args.length)+"_"+str(args.npoints_circle)+".poly")
    with open(name, 'w') as outfile:

        def oprt(*args, **kwargs):
            print(*args, file=outfile, **kwargs)


        oprt("# A cylinder with aspect ratio", args.length)
        oprt("\n# Part 1 - node list")
        oprt("# First line: <# of points> <dimension (3)> <# of attributes> <boundary markers (0 or 1)>")
        oprt(len(cylinder), "3 0 1")

        oprt("# <point #> <x> <y> <z> [attributes] [boundary marker]")
        for p in cylinder:
            oprt(p.i, p.x, p.y, p.z, p.bound)


        oprt("\n# Part 2 - facet list")
        oprt("# First line: <# of facets> <boundary markers (0 or 1)>")
        oprt(len(top) + 2 ,"1")

        oprt("# List of facets:")
        oprt("# One line: <# of polygons> [# of holes] [boundary marker]")
        oprt("# Then list of polygons: <# of corners> <corner 1> <corner 2> ... <corner #>")
        oprt("# Then list of holes: <hole #> <x> <y> <z>")

        oprt("1 0 1         # 1 polygon, no hole, boundary 0")
        oprt(len(bottom), end=" ")
        for p in bottom: # print list of points around the bottom
            oprt(p.i, end=" ")
        oprt("# bottom\n")

        # This is the hard bit: we need to print a list of facets around
        # the sides, so we need sets of 4 points: 2 at the top, 2 at the
        # bottom.
        for a, b, c, d in zip(top, bottom, top[1:] + [top[0]], bottom[1:] + [bottom[0]]):
            oprt("1 0 2")
            oprt("4", b.i, a.i, c.i, d.i)

        oprt("\n1 0 3")
        oprt(len(top), end=" ")
        for p in top: # print list of points around the top
            oprt(p.i, end=" ")
        oprt("# top")

        oprt("\n# Part 3 - hole list")
        oprt("0            # no holes")

        oprt("\n# Part 4 - region list")
        oprt("0            # no regions")

    if args.plot:
        # Create 3d plot
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Plot it
        xs = [p.x for p in cylinder]
        ys = [p.y for p in cylinder]
        zs = [p.z for p in cylinder]
        i_s = [p.i for p in cylinder]
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

        plt.show()

    return 0


if __name__ == "__main__":
    sys.exit(main())
