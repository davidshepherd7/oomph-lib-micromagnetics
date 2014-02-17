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

from fpdiff import fpdiff

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
    norm_test = max_err_norm < tol

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


def check_ndt_less_than(data, max_ndt, identifier=None):

    if identifier is None:
        identifier = _default_label(data)

    ndt = len(data['dts'])

    passed = ndt < max_ndt

    if not passed:
        mm.badprint("FAILED in ",  identifier)
        mm.badprint("n steps is", ndt, "which is more than", max_ndt)

    else:
        mm.okprint("n steps ok in", identifier)

    return passed



def check_restarted_in_middle(restart_outdir, restart_point=20):

    last_no_soln = "soln%i.dat" % (restart_point - 1)
    first_soln = "soln%i.dat" % restart_point

    ok = (not os.path.isfile(pjoin(restart_outdir, last_no_soln))) \
      and (os.path.isfile(pjoin(restart_outdir, first_soln)))

    if ok:
        mm.okprint("Restarted at correct point in", restart_outdir)
    else:
        mm.badprint("Failed or restarted at wrong point in", restart_outdir)

    return ok


def check_solns_match(main_dir="Validation", compare_dir="validata", **kwargs):

    matches = []
    for fname in os.listdir(main_dir):
        if ("soln" in fname) and (".dat" in fname):
            match = fpdiff(pjoin(main_dir, fname),
                            pjoin(compare_dir, fname), **kwargs)
            matches.append(match)

    if len(matches) == 0:
        mm.badprint("No files in", main_dir)
        return False

    ok = all(matches)
    if ok:
        mm.okprint("Files match in", main_dir, compare_dir)
    else:
        mm.badprint("Files don't match in", main_dir, compare_dir)

    return ok
