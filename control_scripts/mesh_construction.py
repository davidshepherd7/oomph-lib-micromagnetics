#!/usr/bin/env python3

# Future proofing
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

# Imports from main libraries
import subprocess as subp
from multiprocessing import Pool
import multiprocessing
import sys
import argparse
import os
import os.path

from os.path import join as pjoin


def generate_triangle_meshes(initial_poly_file):
    """Generate a set of mesh refinements using triangle.
    """

    # Chop up the filename into parts so we can construct the refined
    # filenames later.
    dirname, fname = os.path.split(initial_poly_file)
    basename, ext = os.path.splitext(fname)

    # Initial mesh generation
    initial_area = 0.05
    subp.check_call(['triangle', '-q30', '-a' + str(initial_area),
                     initial_poly_file])

    # Generate refinements
    for refine in [1, 2, 3, 4, 5]:

        # File to refine is something like square.1.poly (if we started
        # with square.poly).
        file_to_refine = (os.path.join(dirname, basename)
                          + '.' + str(refine) + ext)

        # Refine triangle areas by a factor of 4 each time (i.e. 2^dim).
        subp.check_call(['triangle', '-q30',
                         '-a' + str(initial_area / (4 ** refine)),
                         file_to_refine])

    return


def generate_tetgen_meshes(initial_poly_file):
    """Generate a set of mesh refinements using tegen.
    """

    # Chop up the filename into parts so we can construct the refined
    # filenames later.
    dirname, fname = os.path.split(initial_poly_file)
    basename, ext = os.path.splitext(fname)

    # Initial mesh generation
    initial_vol = 0.02
    subp.check_call(['tetgen', '-V', '-q', '-a' + str(initial_vol), '-p',
                     initial_poly_file])


    # Generate refinements
    for refine in [1, 2, 3, 4, 5]:

        # File to refine is something like square.1.poly (if we started with
        # square.poly).
        file_to_refine = (os.path.join(dirname, basename)
                          + '.' + str(refine))

        # Refine tet volumes by a factor of 8 each time (i.e. 2^dim).
        subp.check_call(['tetgen', '-V', '-q',
                         '-a' + str(initial_vol / (4 ** refine)),
                         '-r', file_to_refine])

    return

def generate_mumag4_meshes(initial_poly_file):
    """Generate a set of mesh refinements using tegen.
    """

    # Chop up the filename into parts so we can construct the refined
    # filenames later.
    dirname, fname = os.path.split(initial_poly_file)
    basename, ext = os.path.splitext(fname)

    # Initial mesh generation
    initial_vol = 10
    subp.check_call(['tetgen', '-V', '-q5', '-a' + str(initial_vol), '-p',
                     initial_poly_file])


    # Generate refinements
    for refine in [1, 2]:

        # File to refine is something like square.1.poly (if we started with
        # square.poly).
        file_to_refine = (os.path.join(dirname, basename)
                          + '.' + str(refine))

        # Refine tet volumes by a factor of 8 each time (i.e. 2^dim).
        subp.check_call(['tetgen', '-V', '-q',
                         '-a' + str(initial_vol / (4 ** refine)),
                         '-r', file_to_refine])

    return


def generate_sphere_mesh(mesh_dir, radius, refinement):

    # Generate the input file (a polyhedron)
    command = pjoin(mesh_dir, "generate_tetgen_sphere_input",
                    "generate_tetgen_sphere_input")

    # Name the file with "refinement-1" because tetgen will increment the
    # label (I guess it assumes we are refining something in the normal
    # way).
    mesh_input_filename = pjoin(mesh_dir, "sphere." + str(refinement - 1) + ".poly")

    with open(mesh_input_filename, 'w') as mesh_input_file:
        subp.check_call([command, str(radius), str(refinement)],
                        stdout = mesh_input_file)

    initial_vol = 0.1

    # Now use the input file to create a mesh using tetgen
    subp.check_call(['tetgen',
                     '-V', # dump info on quality
                     '-q', # "Good quality" mesh
                     '-a' + str(initial_vol / (4 ** refinement)), # specify an area
                     '-Y', # Don't split boundary faces (so that surface refinement
                           # is effectively controlled by the input file).
                     '-p', mesh_input_filename],
                    cwd = mesh_dir)

    return


def generate_sphere_meshes(mesh_dir, radius):
    """
    Note: poly files for a refine may or may not correspond to the other
    files...

    Get ~ 4x number of nodes with each refine.
    """

    # Make sure the input generation code is built
    subp.check_call(['make', '-k', '-s'],
                    cwd = pjoin(mesh_dir, "generate_tetgen_sphere_input"))

    # Make some meshes with different refinement levels
    for refine in [1, 2, 3, 4, 5, 6]:
         generate_sphere_mesh(mesh_dir, radius, refine)

    return


def main():
    """Build unstructured meshes using the poly files in "../meshes".


    ??ds -- circle input mesh is too big
    ??ds -- 2d meshes, last two refines are the same...
    """


    # Parse inputs
    # ============================================================

    parser = argparse.ArgumentParser(description=main.__doc__,

    # Don't mess up my formating in the help message
    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-C', dest='mesh_dir', default='../meshes',
                       help = 'Set the mesh directory relative to pwd.')

    args = parser.parse_args()


    # Main function
    # ============================================================

    triangle_poly_files = ['square.poly', 'circle.poly']
    tetgen_poly_files = ['cubeoid.poly']

    for poly_file in triangle_poly_files:
        generate_triangle_meshes(os.path.join(args.mesh_dir, poly_file))

    for poly_file in tetgen_poly_files:
        generate_tetgen_meshes(os.path.join(args.mesh_dir, poly_file))

    generate_mumag4_meshes(os.path.join(args.mesh_dir, 'mumag4.poly'))

    # Sphere meshes are special
    generate_sphere_meshes(args.mesh_dir, 1.0)

    return 0


# If this script is run from a shell then run main() and return the result.
if __name__ == "__main__":
    sys.exit(main())
