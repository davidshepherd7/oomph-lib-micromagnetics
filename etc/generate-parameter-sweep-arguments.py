#!/usr/bin/env python3

# Python 2/3 compatability
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import sys
import argparse
import os
import os.path

import ast

import oomphpy
import oomphpy.micromagnetics as mm
import oomphpy.utils as utils


from os.path import join as pjoin

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('parameter_file')
    args = parser.parse_args()

    with open(args.parameter_file, 'r') as pfile:
        args_file = ast.literal_eval(pfile.read())

    try:
        args_dict = args_file[0]
        extra = list(args_file[1:])
    except KeyError:
        args_dict = args_file
        extra = None

    argdicts = sum([mm.generate_argdicts(args_dict, extra_args_dicts=e) for e in extra], [])

    arglists, _, _ = utils.unzip([mm.argdict2list(a) for a in argdicts])

    for a in arglists:
        print("./driver", ' '.join(a), "\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
