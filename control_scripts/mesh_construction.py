#!/usr/bin/env python

# Imports from main libraries
import subprocess as subp
from multiprocessing import Pool
import multiprocessing
import sys
import argparse
import os
import os.path



def generate_triangle_meshes(initial_poly_file):
    """Generate a set of mesh refinements using triangle.
    """

    # Chop up the filename into parts so we can construct the refined
    # filenames later.
    dirname, fname = os.path.split(initial_poly_file)
    basename, ext = os.path.splitext(fname)

    # Initial mesh generation
    initial_area = 0.1
    subp.check_call(['triangle', '-q30', '-a' + str(initial_area),
                     initial_poly_file])

    # Generate refinements
    for refine in [1, 2, 3, 4, 5]:

        # File to refine is someting like square.1.poly (if we started with
        # square.poly).
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

    raise NotYetImplementedError


def main():
    """Build unstructured meshes using the poly files in "../meshes".
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
    tetgen_poly_files = []

    for poly_file in triangle_poly_files:
        generate_triangle_meshes(os.path.join(args.mesh_dir, poly_file))

    for poly_file in tetgen_poly_files:
        generate_tetgen_meshes(os.path.join(args.mesh_dir, poly_file))

    return 0


# If this script is run from a shell then run main() and return the result.
if __name__ == "__main__":
    sys.exit(main())
