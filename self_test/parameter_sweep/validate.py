#!/usr/bin/env python3

# Python 2/3 compatability
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import sys
import argparse

import os
import shutil
import os.path
import scipy as sp
import itertools as it
import subprocess as subp

from os.path import join as pjoin

# Make sure *this* versions oomphpy is in the path (before any other
# versions in other places)
sys.path.insert(1, pjoin(os.path.dirname(__file__), "../../etc"))
import oomphpy
import oomphpy.micromagnetics as mm



def main():

    parameter_sweep_path = pjoin(mm.rootdir(), 'control_scripts/parameter-sweep.py')

    print(parameter_sweep_path)
    # Check parameter sweep exists
    assert(os.path.isfile(parameter_sweep_path))

    # Check we get output when runs correctly
    subp.check_call([parameter_sweep_path, '-p', 'self_test', '--clean',
                     '--no-build'])
    root = pjoin(mm.rootdir(), 'experiments/parameter_sweeps/self_test')
    for _, odirs, _ in os.walk(root):
        for odir in odirs:
            odir = pjoin(root, odir)
            print("checking files in", odir)
            assert os.path.isfile(pjoin(odir, "trace")), pjoin(odir, "trace")
            assert(os.path.isfile(pjoin(odir, "info")))
            assert(os.path.isfile(pjoin(odir, "run_script")))
            assert(not os.path.isfile(pjoin(odir, "FAILED")))

    assert(os.path.isfile(pjoin(root, "parameter_file")))


    # and when it fails
    print("Note that the following run is supposed to fail!")
    subp.check_call([parameter_sweep_path, '-p', 'fail_test', '--clean',
                     '--no-build'])
    root = pjoin(mm.rootdir(), 'experiments/parameter_sweeps/fail_test')
    for _, odirs, _ in os.walk(root):
        for odir in odirs:
            odir = pjoin(root, odir)
            assert(os.path.isfile(pjoin(odir, "trace")))
            assert(os.path.isfile(pjoin(odir, "info")))
            assert(os.path.isfile(pjoin(odir, "run_script")))
            assert(os.path.isfile(pjoin(odir, "FAILED")))

    assert(os.path.isfile(pjoin(root, "parameter_file")))


    return 0


if __name__ == "__main__":
    sys.exit(main())
