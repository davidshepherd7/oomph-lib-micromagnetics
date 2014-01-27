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
    for odir in os.listdir(root):
        odir = pjoin(root, odir)
        assert(os.path.isfile(pjoin(odir, "trace")))
        assert(os.path.isfile(pjoin(odir, "info")))
        assert(os.path.isfile(pjoin(odir, "run_script")))
        assert(not os.path.isfile(pjoin(odir, "FAILED")))

    # and when it fails
    subp.check_call([parameter_sweep_path, '-p', 'fail_test', '--clean',
                     '--no-build'])
    root = pjoin(mm.rootdir(), 'experiments/parameter_sweeps/fail_test')
    for odir in os.listdir(root):
        odir = pjoin(root, odir)
        assert(os.path.isfile(pjoin(odir, "trace")))
        assert(os.path.isfile(pjoin(odir, "info")))
        assert(os.path.isfile(pjoin(odir, "run_script")))
        assert(os.path.isfile(pjoin(odir, "FAILED")))

    return 0


if __name__ == "__main__":
    sys.exit(main())
