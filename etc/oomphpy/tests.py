#!/usr/bin/env python3

# Python 2/3 compatability
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import sys
import argparse
import os
import os.path

from os.path import join as pjoin

# Make sure *this* versions oomphpy is in the path (before any other
# versions in other places)
sys.path.insert(1, pjoin(os.path.dirname(__file__), "../"))
import oomphpy
import oomphpy.micromagnetics as mm


def _default_label(data):
    """Get last two parts of the path"""
    b, c = os.path.split(data['-outdir'])
    a, b = os.path.split(b)
    return pjoin(b, c)


def check_error_norm(data, tol=1e-5, identifier=None):
    """Check that max error norm < tol."""

    if identifier is None:
        identifier = _default_label(data)

    max_err_norm = max(map(abs, data['error_norms']))
    norm_test = max_err_norm < 1e-5

    if not norm_test:
        mm.badprint("FAILED in", identifier)
        mm.badprint("with max error norm of", max_err_norm)

    else:
        mm.okprint("error norm ok in", identifier)

    return norm_test


def check_m_length(data, tol=1e-8, identifier=None):
    """Check that abs(|m| - 1) < tol."""

    if identifier is None:
        identifier = _default_label(data)

    max_length_err = max(map(abs, data['m_length_error_means']))
    length_test = max_length_err < 1e-14

    if not length_test:
        mm.badprint("FAILED in ", identifier)
        mm.badprint("|m| error of", max_length_err)

    else:
        mm.okprint("|m| ok in", identifier)

    return length_test
